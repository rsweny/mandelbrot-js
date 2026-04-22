/*
 * WebGPU compute-shader renderer for the Mandelbulb.
 *
 * Replaces the web-worker pool with a GPU progressive renderer that mirrors
 * the CPU pipeline: per-pixel depth march, ray-traced ambient + primary light,
 * frost neighbour samples, ray-hit goodPoints, fog/glow trace projection, and
 * shared occlusion buffer. Each compute dispatch is one full progressive pass;
 * the CPU keeps owning histogram / tone-mapping driven by max_alpha.
 *
 * Reads the same global parameters from mandelbulb.js (zoom, xcen, ycen,
 * CameraMatrix, IrotX, IrotZ, iterations, formula, etc.) so a config change
 * just rebinds uniforms and restarts the pass loop.
 */

// ---------- WGSL ----------

var GPU_RENDER_WGSL = `
const MAX_DEPTH_STEPS: u32 = 1024u;
const MAX_GOODPOINTS: u32  = 5u;

struct Params {
    M: array<vec4<f32>, 3>,
    IrotX: array<vec4<f32>, 3>,
    IrotZ: array<vec4<f32>, 3>,
    LightVector: vec4<f32>,
    fogColor: vec4<f32>,

    ximlen: u32, yimlen: u32, iterations: u32, formula: u32,
    RAY_STEPS: u32, passSeed: u32, tileOriginX: u32, tileOriginY: u32,

    power: f32, azimuth: f32, zoom: f32, xcen: f32,
    ycen: f32, cameraPersp: f32, half_ximlen: f32, half_yimlen: f32,
    AMBIENT_LIGHT: f32, primary_light: f32, shadow_darkness: f32, HORIZON: f32,
    ray_step: f32, stepDetail: f32, root_zoom: f32, opacity: f32,
    frost: f32, fog_factor: f32, cameraDOF: f32, factorDOF: f32,
    focus: f32, focus_depth: f32, min_y: f32, max_y: f32,
    useVolumetricFog: u32, _pad0: u32, _pad1: u32, _pad2: u32,
};

@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read_write> accum:    array<atomic<u32>>;
@group(0) @binding(2) var<storage, read_write> occLive:  array<atomic<u32>>;
@group(0) @binding(3) var<storage, read>       occSnap:  array<f32>;
@group(0) @binding(4) var<storage, read>       palette:  array<vec4<f32>, 256>;
@group(0) @binding(5) var<storage, read_write> maxAlpha: array<atomic<u32>>;

// f32 atomic add via compare-exchange loop on bit-cast bits. Returns the new value.
fn accAddF32(idx: u32, v: f32) -> f32 {
    var oldBits = atomicLoad(&accum[idx]);
    loop {
        let newVal  = bitcast<f32>(oldBits) + v;
        let newBits = bitcast<u32>(newVal);
        let res = atomicCompareExchangeWeak(&accum[idx], oldBits, newBits);
        if (res.exchanged) { return newVal; }
        oldBits = res.old_value;
    }
}

// PCG-ish hash; returns f32 in [0, 1).
fn nextRand(s: ptr<function, u32>) -> f32 {
    var x = *s;
    x = x * 747796405u + 2891336453u;
    let w = ((x >> ((x >> 28u) + 4u)) ^ x) * 277803737u;
    let r = (w >> 22u) ^ w;
    *s = x;
    return f32(r) * (1.0 / 4294967296.0);
}

fn rotate3(M0: vec3<f32>, M1: vec3<f32>, M2: vec3<f32>, x: f32, y: f32, z: f32) -> vec3<f32> {
    return vec3<f32>(
        M0.x*x + M0.y*y + M0.z*z,
        M1.x*x + M1.y*y + M1.z*z,
        M2.x*x + M2.y*y + M2.z*z);
}

struct InsideResult { iter: u32, color: f32 };

fn insideFractal(c: vec3<f32>) -> InsideResult {
    var p = c;
    var iter: u32 = 0u;
    var color: f32 = 0.0;
    var r: f32 = 0.0;
    let pwr = P.power; let az = P.azimuth;
    let formula = P.formula; let maxIter = P.iterations;
    for (var i: u32 = 0u; i < maxIter; i = i + 1u) {
        r = length(p);
        let theta_p = atan2(p.y, p.x) * pwr;
        let r_p = pow(r, pwr);
        var phi: f32;
        if (formula == 0u) {
            phi = asin(p.z / r);
            let pc = cos(phi * pwr);
            let ct = cos(theta_p); let st = sin(theta_p);
            let sp = sin(phi * pwr);
            p = vec3<f32>(r_p*ct*pc, r_p*st*pc, r_p*sp*az) + c;
        } else {
            phi = atan2(length(p.xy), p.z);
            let ps = sin(phi * pwr);
            let ct = cos(theta_p); let st = sin(theta_p);
            let cp = cos(phi * pwr);
            p = vec3<f32>(r_p*ct*ps, r_p*st*ps, r_p*cp*az) + c;
        }
        color = phi / 3.0;
        iter = iter + 1u;
        if (r >= 8.0) { break; }
    }
    return InsideResult(iter, color);
}

struct GoodPoints { pts: array<vec4<f32>, MAX_GOODPOINTS>, cnt: u32 };

fn pushGP(gp: ptr<function, GoodPoints>, p: vec3<f32>, color: f32) {
    let n = (*gp).cnt;
    if (n < MAX_GOODPOINTS) {
        (*gp).pts[n] = vec4<f32>(p, color);
        (*gp).cnt = n + 1u;
    }
}

// March a ray; on solid hit returns 0 and (optionally) records a goodPoint.
fn calcRayGP(start: vec3<f32>, steps: u32, step: vec3<f32>, bright: f32, fuzzy: f32, gp: ptr<function, GoodPoints>) -> f32 {
    let s = step * fuzzy;
    var p = start;
    for (var i: u32 = 1u; i < steps; i = i + 1u) {
        p = p + s * f32(i);
        let r = insideFractal(p);
        if (r.iter == P.iterations) { pushGP(gp, p, r.color); return 0.0; }

        // ray has left the general area of the solid, stop tracing.
        if (r.iter < 4u) { break; }
    }
    return bright;
}

fn calcRay(start: vec3<f32>, steps: u32, step: vec3<f32>, bright: f32, fuzzy: f32) -> f32 {
    let s = step * fuzzy;
    var p = start;
    for (var i: u32 = 1u; i < steps; i = i + 1u) {
        p = p + s * f32(i);
        let r = insideFractal(p);
        if (r.iter == P.iterations) { return 0.0; }

        // ray has left the general area of the solid, stop tracing
        // (go a bit further out here to make volume lighting look good)
        if (r.iter < 2u) { break; }
    }
    return bright;
}

fn calcRaysGP(p: vec3<f32>, fuzzy: f32, rng: ptr<function, u32>, gp: ptr<function, GoodPoints>) -> f32 {
    var lf: f32 = 1.0;
    let rs = P.ray_step;
    let lv = P.LightVector.xyz;
    let n0 = rs / (15.0 * nextRand(rng) + 1.0);
    let n1 = rs / (15.0 * nextRand(rng) + 1.0);
    let n2 = rs / (15.0 * nextRand(rng) + 1.0);
    lf = lf + calcRayGP(p, P.RAY_STEPS,      vec3<f32>(-n0, -n1, -n2),                                  P.AMBIENT_LIGHT, fuzzy, gp);
    lf = lf + calcRayGP(p, P.RAY_STEPS,      vec3<f32>( n0,  n1,  n2),                                  P.AMBIENT_LIGHT, fuzzy, gp);
    lf = lf + calcRayGP(p, P.RAY_STEPS,      vec3<f32>( n0, -n1, -rs),                                  P.AMBIENT_LIGHT, fuzzy, gp);
    lf = lf + calcRayGP(p, P.RAY_STEPS,      vec3<f32>(n0*lv.x*9.0, n1*lv.y*9.0, n2*lv.z*9.0),          P.AMBIENT_LIGHT, fuzzy, gp);
    lf = lf + calcRayGP(p, P.RAY_STEPS*4u,   vec3<f32>(rs*lv.x, rs*lv.y, rs*lv.z),                      P.primary_light, fuzzy, gp);
    return lf;
}

fn calcRaysSimple(p: vec3<f32>, fuzzy: f32, rng: ptr<function, u32>) -> f32 {
    var lf: f32 = 1.0;
    let rs = P.ray_step;
    let lv = P.LightVector.xyz;
    let n0 = rs / (15.0 * nextRand(rng) + 1.0);
    let n1 = rs / (15.0 * nextRand(rng) + 1.0);
    let n2 = rs / (15.0 * nextRand(rng) + 1.0);
    lf = lf + calcRay(p, P.RAY_STEPS,      vec3<f32>(-n0, -n1, -n2),                                    P.AMBIENT_LIGHT, fuzzy);
    lf = lf + calcRay(p, P.RAY_STEPS,      vec3<f32>( n0,  n1,  n2),                                    P.AMBIENT_LIGHT, fuzzy);
    lf = lf + calcRay(p, P.RAY_STEPS,      vec3<f32>( n0, -n1, -rs),                                    P.AMBIENT_LIGHT, fuzzy);
    lf = lf + calcRay(p, P.RAY_STEPS,      vec3<f32>(n0*lv.x*9.0, n1*lv.y*9.0, n2*lv.z*9.0),            P.AMBIENT_LIGHT, fuzzy);
    lf = lf + calcRay(p, P.RAY_STEPS*4u,   vec3<f32>(rs*lv.x, rs*lv.y, rs*lv.z),                        P.primary_light, fuzzy);
    return lf;
}

// Mirrors calcVolumetricRays in the worker: primary-light ray only, no
// ambient neighbours. Used to tint fog traces for volumetric lighting
fn calcVolumetricRays(p: vec3<f32>, fuzzy: f32) -> f32 {
    let rs = P.ray_step*5;
    let lv = P.LightVector.xyz;
    var lf: f32 = 0.05 + calcRay(p, P.RAY_STEPS*3u, vec3<f32>(rs*lv.x, rs*lv.y, rs*lv.z), P.primary_light, fuzzy)*0.035;
    return lf;
}

// Reverse-project a fractal-space point back to screen-space (pixel x, depth y, pixel z).
fn reversePoint(fp: vec3<f32>) -> vec3<f32> {
    var q = rotate3(P.IrotZ[0].xyz, P.IrotZ[1].xyz, P.IrotZ[2].xyz, fp.x, fp.y, fp.z);
    q = rotate3(P.IrotX[0].xyz, P.IrotX[1].xyz, P.IrotX[2].xyz, q.x, q.y, q.z);
    let persp = 1.0 + q.y * P.cameraPersp;
    let sx = ((q.x * persp + P.xcen) / P.zoom) * f32(P.ximlen) + P.half_ximlen;
    let sz = ((q.z * persp + P.ycen) / P.zoom) * f32(P.yimlen) + P.half_yimlen;
    return vec3<f32>(sx, q.y, sz);
}

fn pixelIdx(x: u32, y: u32) -> u32 { return (x * P.yimlen + y) * 4u; }
fn occIdx(x: u32, y: u32)   -> u32 { return x * P.yimlen + y; }

fn writeSurface(tx: i32, ty: i32, depth: f32, color: f32, lightFactor: f32) {
    if (tx < 0 || ty < 0 || tx >= i32(P.ximlen) || ty >= i32(P.yimlen)) { return; }
    let ux = u32(tx); let uy = u32(ty);
    let ci = u32(min(255.0, abs(color) * 255.0));
    let pal = palette[ci];
    let idx = pixelIdx(ux, uy);
    accAddF32(idx + 0u, (pal.r * lightFactor) / P.shadow_darkness);
    accAddF32(idx + 1u, (pal.g * lightFactor) / P.shadow_darkness);
    accAddF32(idx + 2u, (pal.b * lightFactor) / P.shadow_darkness);
    let newAlpha = accAddF32(idx + 3u, 1.0);
    // Occlusion: atomicMin on f32 bits after offsetting depth so it is positive.
    atomicMin(&occLive[occIdx(ux, uy)], bitcast<u32>(depth + P.HORIZON));
    atomicMax(&maxAlpha[0], bitcast<u32>(newAlpha));
}

fn blurOffset(x: f32, y: f32, blurFactor: f32, rng: ptr<function, u32>) -> vec2<f32> {
    var bf = blurFactor;
    if (bf > 0.0) {
        bf = bf - P.focus_depth;
        if (bf < 0.0) { return vec2<f32>(x, y); }
    }
    let r2 = nextRand(rng) * 6.2831853;
    let dr = nextRand(rng) * P.factorDOF * bf;
    return vec2<f32>(x + dr * cos(r2), y + dr * sin(r2));
}

fn plotSurface(sx: u32, sy: u32, depth: f32, color: f32, lf: f32, rng: ptr<function, u32>) {
    if (P.cameraDOF > 0.0) {
        let bp = blurOffset(f32(sx), f32(sy), P.focus - depth, rng);
        let tx = i32(floor(bp.x)); let ty = i32(floor(bp.y));
        if (tx < 0 || ty < 0 || tx >= i32(P.ximlen) || ty >= i32(P.yimlen)) { return; }
        let snap = occSnap[u32(tx) * P.yimlen + u32(ty)];
        if (depth > snap + P.stepDetail) { return; }
        writeSurface(tx, ty, depth, color, lf);
    } else {
        writeSurface(i32(sx), i32(sy), depth, color, lf);
    }
}

// Reverse-project a goodPoint and plot as a surface pixel.
fn plotGoodPoint(fp: vec3<f32>, color: f32, lf: f32, rng: ptr<function, u32>) {
    let sp = reversePoint(fp);
    let tx = i32(round(sp.x)); let ty = i32(round(sp.z));
    if (tx < 0 || ty < 0 || tx >= i32(P.ximlen) || ty >= i32(P.yimlen)) { return; }
    let snap = occSnap[u32(tx) * P.yimlen + u32(ty)];
    if (sp.y >= snap) { return; }
    plotSurface(u32(tx), u32(ty), sp.y, color, lf, rng);
}

fn plotFog(fp: vec3<f32>, factor: f32, rng: ptr<function, u32>) {
    var sp = reversePoint(fp);
    var tx = i32(round(sp.x)); var ty = i32(round(sp.z));
    if (tx < 0 || ty < 0 || tx >= i32(P.ximlen) || ty >= i32(P.yimlen)) { return; }
    let snap = occSnap[u32(tx) * P.yimlen + u32(ty)];
    if (sp.y > snap) { return; }
    if (P.cameraDOF > 0.0) {
        let bp = blurOffset(sp.x, sp.z, P.focus - sp.y, rng);
        tx = i32(round(bp.x)); ty = i32(round(bp.y));
        if (tx < 0 || ty < 0 || tx >= i32(P.ximlen) || ty >= i32(P.yimlen)) { return; }
    }
    let fc = P.fogColor.xyz;
    let idx = pixelIdx(u32(tx), u32(ty));
    accAddF32(idx + 0u, fc.r * factor);
    accAddF32(idx + 1u, fc.g * factor);
    accAddF32(idx + 2u, fc.b * factor);
    accAddF32(idx + 3u, factor);
}

// Streaming fog: re-run the iteration and plot each trace point as we go,
// avoiding the need to store a per-pixel trace history in private memory.
// When volumetric fog is enabled, each sample gets its own primary-light
// shadow factor (mirroring the CPU for (c=0; c<iter; c++) loop).
// unlike CPU this traces all points, rather than just those that escape the set
fn fogTrace(c: vec3<f32>, factor: f32, fuzzy: f32, rng: ptr<function, u32>) {
    var p = c;
    var r: f32 = 0.0;
    let pwr = P.power; let az = P.azimuth;
    let formula = P.formula; let maxIter = P.iterations;
    for (var i: u32 = 0u; i < 20u; i = i + 1u) {
        r = length(p);
        let theta_p = atan2(p.y, p.x) * pwr;
        let r_p = pow(r, pwr);
        if (formula == 0u) {
            let phi = asin(p.z / r);
            let pc = cos(phi * pwr);
            let ct = cos(theta_p); let st = sin(theta_p);
            let sp = sin(phi * pwr);
            p = vec3<f32>(r_p*ct*pc, r_p*st*pc, r_p*sp*az) + c;
        } else {
            let phi = atan2(length(p.xy), p.z);
            let ps = sin(phi * pwr);
            let ct = cos(theta_p); let st = sin(theta_p);
            let cp = cos(phi * pwr);
            p = vec3<f32>(r_p*ct*ps, r_p*st*ps, r_p*cp*az) + c;
        }

        if (r >= 8.0) { break; }

        // glow spikes fog for even iterations
        if (i > 0u || maxIter % 2 == 0) {
            var volumetricLightFactor: f32 = 1.0;
            if (P.useVolumetricFog != 0u) {
                
                // fade out the fog close to the camera unless user prefers intense fog
                if (factor < 1) {
                    let sp = reversePoint(p);
                    let fogFadeout = 0.23 + fuzzy*0.1; // 0.0 = near, 1.0 = far
                    let minFogDepth = mix(P.min_y, P.max_y, fogFadeout); 
                    if (sp.y < minFogDepth) {
                        continue;
                    }
                }

                // calc volume lighting
                volumetricLightFactor = calcVolumetricRays(p, fuzzy);
            }

            let fogGlowFactor = factor * volumetricLightFactor;
            plotFog(p, fogGlowFactor, rng); 
        }
    }
}

@compute @workgroup_size(8, 8)
fn renderMain(@builtin(global_invocation_id) gid: vec3<u32>) {
    let x = gid.x + P.tileOriginX;
    let yRow = gid.y + P.tileOriginY;
    if (x >= P.ximlen || yRow >= P.yimlen) { return; }

    var rng: u32 = (yRow * P.ximlen + x) * 1664525u + P.passSeed * 2654435761u + 1u;
    // Warm-up a few rolls so nearby seeds diverge quickly.
    let _w0 = nextRand(&rng); let _w1 = nextRand(&rng);

    let M0 = P.M[0].xyz; let M1 = P.M[1].xyz; let M2 = P.M[2].xyz;

    let jitter1 = ((0.5 - nextRand(&rng)) / f32(P.ximlen)) * 0.4;
    let jitter2 = ((0.5 - nextRand(&rng)) / f32(P.yimlen)) * 0.4;
    var stepAmount = (P.stepDetail + nextRand(&rng) * P.stepDetail) * P.root_zoom;

    let z2 = (f32(yRow) / f32(P.yimlen) - 0.5 + jitter1) * P.zoom - P.ycen;
    let x2 = (f32(x)    / f32(P.ximlen) - 0.5 + jitter2) * P.zoom - P.xcen;

    var found_limit: u32 = u32(ceil(P.frost));
    if (P.min_y == -2.0) { found_limit = 0u; }

    var y = P.min_y;
    var steps: u32 = 0u;
    loop {
        if (y >= P.max_y) { break; }
        if (steps >= MAX_DEPTH_STEPS) { break; }
        steps = steps + 1u;

        let persp = 1.0 + y * P.cameraPersp;
        let p3 = rotate3(M0, M1, M2, x2 / persp, y, z2 / persp);
        let r = insideFractal(p3);

        if (r.iter == P.iterations) {

            // point is in set so plot solid pixel
            let fuzzy = max(P.opacity * nextRand(&rng), 0.4);
            var gp: GoodPoints;
            gp.cnt = 0u;
            let lf = 1.0 + calcRaysGP(p3, fuzzy, &rng, &gp);
            plotSurface(x, yRow, y, r.color, lf, &rng);

            // Frost: sample a few extra y-values near the surface.
            var found: u32 = 0u;
            var cnt: u32 = 0u;
            let maxCnt: u32 = min(100u, u32(100.0 * P.frost));
            loop {
                if (found >= found_limit || cnt >= maxCnt) { break; }
                cnt = cnt + 1u;
                let newy = y - stepAmount * nextRand(&rng) * 2.0;
                let persp2 = 1.0 + newy * P.cameraPersp;
                let p2 = rotate3(M0, M1, M2, x2 / persp2, newy, z2 / persp2);
                let r2 = insideFractal(p2);
                if (r2.iter == P.iterations) {
                    let lf2 = 1.0 + calcRaysSimple(p2, fuzzy, &rng);
                    plotSurface(x, yRow, newy, r2.color, lf2, &rng);
                    found = found + 1u;
                }
            }

            // goodPoints scatter (only when frost > 0.9, matching CPU).
            if (P.frost > 0.9) {
                for (var i: u32 = 0u; i < gp.cnt; i = i + 1u) {
                    let tp = gp.pts[i];
                    let lfg = 1.0 + calcRaysSimple(tp.xyz, fuzzy, &rng);
                    plotGoodPoint(tp.xyz, tp.w, lfg, &rng);
                }
            }
            break;
        } else if (P.min_y != -2.0) {
            // point is not in set (and we have done a rough pass to generate the rough occlusions)
            // use much smaller and more accurate stepAmount
            let rnd = nextRand(&rng);
            stepAmount = (P.stepDetail + rnd * P.stepDetail) *
                         (f32(P.iterations) / f32(max(r.iter, 1u)) / f32(P.iterations)) *
                         P.root_zoom * 0.5;

            // and for points not in the set, optionally plot traces that act as a fog / glow
            if (P.fog_factor > 0.0 && rnd > 0.9 && r.iter > 1u) {
                let fogFuzzy = max(P.opacity * nextRand(&rng), 0.4);
                fogTrace(p3, P.fog_factor, fogFuzzy, &rng);
            }
        }
        y = y + stepAmount;
    }
}

`;

// ---- auxiliary compute modules ----

var GPU_CLEAR_WGSL = `
struct ClearParams { count: u32, fillValue: u32, _a: u32, _b: u32 };
@group(0) @binding(0) var<uniform> CP: ClearParams;
@group(0) @binding(1) var<storage, read_write> buf: array<atomic<u32>>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= CP.count) { return; }
    atomicStore(&buf[i], CP.fillValue);
}
`;

var GPU_CLEARF_WGSL = `
struct ClearParams { count: u32, fillValue: u32, _a: u32, _b: u32 };
@group(0) @binding(0) var<uniform> CP: ClearParams;
@group(0) @binding(1) var<storage, read_write> buf: array<f32>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= CP.count) { return; }
    buf[i] = bitcast<f32>(CP.fillValue);
}
`;

// Copy occLive (atomic u32 of depth+HORIZON bits) into occSnap (plain f32 depth).
var GPU_SNAP_WGSL = `
struct CopyParams { count: u32, _a: u32, horizon: f32, _b: u32 };
@group(0) @binding(0) var<uniform> YP: CopyParams;
@group(0) @binding(1) var<storage, read>       occSrc: array<u32>;
@group(0) @binding(2) var<storage, read_write> occDst: array<f32>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= YP.count) { return; }
    occDst[i] = bitcast<f32>(occSrc[i]) - YP.horizon;
}
`;

// Tone-map fragment shader: reads the accumulation buffer, applies the same
// histogram stretch as updateHistogram() in mandelbulb.js, and writes the
// displayed pixel to the WebGPU canvas.
var GPU_TONEMAP_WGSL = `
struct ToneParams {
    ximlen: u32, yimlen: u32, drawFocus: u32, _pad: u32,
    gradient: f32, brightness: f32, focus: f32, focus_depth: f32,
};
@group(0) @binding(0) var<uniform> TP: ToneParams;
@group(0) @binding(1) var<storage, read> accum:    array<u32>;
@group(0) @binding(2) var<storage, read> maxAlpha: array<u32>;
@group(0) @binding(3) var<storage, read> occSnap:  array<f32>;

struct VSOut { @builtin(position) pos: vec4<f32>, @location(0) uv: vec2<f32> };

@vertex
fn vs_main(@builtin(vertex_index) vi: u32) -> VSOut {
    // Fullscreen triangle
    var positions = array<vec2<f32>, 3>(
        vec2<f32>(-1.0, -1.0), vec2<f32>( 3.0, -1.0), vec2<f32>(-1.0,  3.0));
    var uvs = array<vec2<f32>, 3>(
        vec2<f32>(0.0, 1.0), vec2<f32>(2.0, 1.0), vec2<f32>(0.0, -1.0));
    var out: VSOut;
    out.pos = vec4<f32>(positions[vi], 0.0, 1.0);
    out.uv  = uvs[vi];
    return out;
}

@fragment
fn fs_main(in: VSOut) -> @location(0) vec4<f32> {
    let u = clamp(in.uv.x, 0.0, 0.9999);
    let v = clamp(in.uv.y, 0.0, 0.9999);
    // The CPU renderer indexes pixels as [x][y] with x=horizontal, y=vertical.
    // Our accum layout is (x * yimlen + y), so keep that mapping here.
    let x = u32(u * f32(TP.ximlen));
    let y = u32(v * f32(TP.yimlen));
    let idx = (x * TP.yimlen + y) * 4u;

    let alpha = bitcast<f32>(accum[idx + 3u]);
    if (alpha <= 0.0) { return vec4<f32>(0.0, 0.0, 0.0, 1.0); }

    let peak = bitcast<f32>(maxAlpha[0]);
    let maxA = pow(max(peak, 1.0), TP.gradient);
    let z = pow(alpha, TP.gradient) * TP.brightness / maxA;

    var r = bitcast<f32>(accum[idx + 0u]) * z / alpha;
    var g = bitcast<f32>(accum[idx + 1u]) * z / alpha;
    var b = bitcast<f32>(accum[idx + 2u]) * z / alpha;

    if (TP.drawFocus != 0u) {
        let d = occSnap[x * TP.yimlen + y];
        if (d > TP.focus)                        { b = b + 50.0; r = 20.0; g = 20.0; }
        else if (d < TP.focus - TP.focus_depth)  { r = r + 50.0; g = 20.0; b = 20.0; }
    }

    // The CPU path clamps to 0..255, so divide by 255 for the 0..1 output.
    r = clamp(r, 0.0, 255.0) / 255.0;
    g = clamp(g, 0.0, 255.0) / 255.0;
    b = clamp(b, 0.0, 255.0) / 255.0;
    return vec4<f32>(r, g, b, 1.0);
}
`;

// ---------- JS runtime ----------

var _gpu = null;          // persistent GPU state (device, pipelines, buffers)
var _gpuRunning = false;
var _gpuGen = 0;          // bumped on reset so queued frames can bail
var _gpuPassIndex = 0;
var _gpuStartMs = 0;
var _gpuLastUpdateMs = 0;
var GPU_INTERPASS_SLEEP_MS = 100;  // idle gap between passes to ease GPU load

// Tile the per-pass render dispatch into N×N pixel chunks, each submitted as
// its own command buffer. Keeps any single GPU command well under Windows'
// ~2 s GPU TDR (DXGI_ERROR_DEVICE_HUNG) limit at heavy iteration counts and
// large canvases. TILES_PER_AWAIT bounds the in-flight queue depth so the
// browser stays responsive.
const GPU_TILE_PIXELS    = 256;
const GPU_TILES_PER_AWAIT = 8;

const GPU_PARAMS_BYTES = 320;
const GPU_TONE_BYTES   = 32;
const GPU_CLEAR_BYTES  = 16;
const GPU_SNAP_BYTES   = 16;
const GPU_HORIZON_OFFSET = 20.0;  // matches HORIZON in mandelbulb.js
const GPU_TILE_ORIGIN_BYTE_OFFSET = 50 * 4;  // tileOriginX/Y in Params (u32 slots 50,51)

function isWebGPUSupported() {
    return typeof navigator !== 'undefined' && !!navigator.gpu;
}

async function initGPURenderer(destCanvas, width, height) {
    if (!isWebGPUSupported()) return null;
    const adapter = await navigator.gpu.requestAdapter();
    if (!adapter) return null;
    const device = await adapter.requestDevice();

    const gpuCanvas = document.createElement('canvas');
    gpuCanvas.width  = width;
    gpuCanvas.height = height;
    const gpuCtx = gpuCanvas.getContext('webgpu');
    const format = navigator.gpu.getPreferredCanvasFormat();
    gpuCtx.configure({ device: device, format: format, alphaMode: 'opaque' });

    const pixelCount = width * height;

    const paramsBuf = device.createBuffer({
        size: GPU_PARAMS_BYTES,
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    const toneBuf = device.createBuffer({
        size: GPU_TONE_BYTES,
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    const clearParamsBuf = device.createBuffer({
        size: GPU_CLEAR_BYTES,
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    const snapParamsBuf = device.createBuffer({
        size: GPU_SNAP_BYTES,
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    // 4 channels (R,G,B,A) of f32 stored as u32 bit-pattern for CAS atomic adds.
    // COPY_DST lets resetGPUAccumulators use queue.clearBuffer / queue.writeBuffer.
    const accumBuf = device.createBuffer({
        size: pixelCount * 4 * 4,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    const occLiveBuf = device.createBuffer({
        size: pixelCount * 4,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    const occSnapBuf = device.createBuffer({
        size: pixelCount * 4,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    const paletteBuf = device.createBuffer({
        size: 256 * 4 * 4,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    const maxAlphaBuf = device.createBuffer({
        size: 16,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC,
    });
    // CPU-mappable staging buffer for reading the per-pass peak alpha back so
    // gpuStatsCallback can update max_alpha (used by the PNG download name).
    const maxAlphaReadBuf = device.createBuffer({
        size: 16,
        usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
    });

    const renderModule  = device.createShaderModule({ code: GPU_RENDER_WGSL });
    const clearModule   = device.createShaderModule({ code: GPU_CLEAR_WGSL });
    const clearFModule  = device.createShaderModule({ code: GPU_CLEARF_WGSL });
    const snapModule    = device.createShaderModule({ code: GPU_SNAP_WGSL });
    const toneModule    = device.createShaderModule({ code: GPU_TONEMAP_WGSL });

    const renderPipeline = device.createComputePipeline({
        layout: 'auto',
        compute: { module: renderModule, entryPoint: 'renderMain' },
    });
    const clearPipeline  = device.createComputePipeline({
        layout: 'auto',
        compute: { module: clearModule, entryPoint: 'main' },
    });
    const clearFPipeline = device.createComputePipeline({
        layout: 'auto',
        compute: { module: clearFModule, entryPoint: 'main' },
    });
    const snapPipeline   = device.createComputePipeline({
        layout: 'auto',
        compute: { module: snapModule, entryPoint: 'main' },
    });
    const tonePipeline   = device.createRenderPipeline({
        layout: 'auto',
        vertex:   { module: toneModule, entryPoint: 'vs_main' },
        fragment: { module: toneModule, entryPoint: 'fs_main',
                    targets: [{ format: format }] },
        primitive: { topology: 'triangle-list' },
    });

    const renderBG = device.createBindGroup({
        layout: renderPipeline.getBindGroupLayout(0),
        entries: [
            { binding: 0, resource: { buffer: paramsBuf } },
            { binding: 1, resource: { buffer: accumBuf } },
            { binding: 2, resource: { buffer: occLiveBuf } },
            { binding: 3, resource: { buffer: occSnapBuf } },
            { binding: 4, resource: { buffer: paletteBuf } },
            { binding: 5, resource: { buffer: maxAlphaBuf } },
        ],
    });
    const toneBG = device.createBindGroup({
        layout: tonePipeline.getBindGroupLayout(0),
        entries: [
            { binding: 0, resource: { buffer: toneBuf } },
            { binding: 1, resource: { buffer: accumBuf } },
            { binding: 2, resource: { buffer: maxAlphaBuf } },
            { binding: 3, resource: { buffer: occSnapBuf } },
        ],
    });

    _gpu = {
        device, format, gpuCanvas, gpuCtx, destCanvas,
        width, height, pixelCount,
        paramsBuf, toneBuf, clearParamsBuf, snapParamsBuf,
        accumBuf, occLiveBuf, occSnapBuf, paletteBuf, maxAlphaBuf, maxAlphaReadBuf,
        renderPipeline, clearPipeline, clearFPipeline, snapPipeline, tonePipeline,
        renderBG, toneBG,
        paramsBytes: new ArrayBuffer(GPU_PARAMS_BYTES),
        toneBytes:   new ArrayBuffer(GPU_TONE_BYTES),
        clearBytes:  new ArrayBuffer(GPU_CLEAR_BYTES),
        snapBytes:   new ArrayBuffer(GPU_SNAP_BYTES),
        destCtx: destCanvas.getContext('2d'),
        paletteUploaded: false,
    };
    return _gpu;
}

// ---- uniform encoding ----

function _writeMat3(f32v, offset, m) {
    // mat3 stored as 3 × vec4; only xyz used.
    f32v[offset +  0] = m[0][0]; f32v[offset +  1] = m[0][1]; f32v[offset +  2] = m[0][2]; f32v[offset +  3] = 0;
    f32v[offset +  4] = m[1][0]; f32v[offset +  5] = m[1][1]; f32v[offset +  6] = m[1][2]; f32v[offset +  7] = 0;
    f32v[offset +  8] = m[2][0]; f32v[offset +  9] = m[2][1]; f32v[offset + 10] = m[2][2]; f32v[offset + 11] = 0;
}

function writeGPUParams(s, passSeed) {
    const f = new Float32Array(_gpu.paramsBytes);
    const u = new Uint32Array(_gpu.paramsBytes);
    _writeMat3(f,  0, s.CameraMatrix);
    _writeMat3(f, 12, s.IrotX);
    _writeMat3(f, 24, s.IrotZ);
    // LightVector (vec4)
    f[36] = s.LightVector[0]; f[37] = s.LightVector[1]; f[38] = s.LightVector[2]; f[39] = 0;
    // fogColor
    f[40] = s.fog_color.r; f[41] = s.fog_color.g; f[42] = s.fog_color.b; f[43] = 0;
    // u32 block 1
    u[44] = s.ximlen; u[45] = s.yimlen; u[46] = s.iterations; u[47] = s.formula;
    // u32 block 2 (slots 50/51 are tileOriginX/Y, overwritten per tile by writeGPUTileOrigin)
    u[48] = s.RAY_STEPS; u[49] = passSeed >>> 0; u[50] = 0; u[51] = 0;
    // f32 tail
    f[52] = s.power;         f[53] = s.azimuth;        f[54] = s.zoom;            f[55] = s.xcen;
    f[56] = s.ycen;          f[57] = s.cameraPersp;    f[58] = s.half_ximlen;     f[59] = s.half_yimlen;
    f[60] = s.AMBIENT_LIGHT; f[61] = s.primary_light;  f[62] = s.shadow_darkness; f[63] = s.HORIZON;
    f[64] = s.ray_step;      f[65] = s.stepDetail;     f[66] = s.root_zoom;       f[67] = s.opacity;
    f[68] = s.frost;         f[69] = s.fog_factor;     f[70] = s.cameraDOF;       f[71] = s.factorDOF;
    f[72] = s.focus;         f[73] = s.focus_depth;    f[74] = s.min_y;           f[75] = s.max_y;
    // useVolumetricFog + 3 padding u32s (slots 77..79) keep the struct 16-byte aligned.
    u[76] = s.useVolumetricFog ? 1 : 0; u[77] = 0; u[78] = 0; u[79] = 0;
    _gpu.device.queue.writeBuffer(_gpu.paramsBuf, 0, _gpu.paramsBytes);
}

// Patch just the two tileOrigin u32s in the params buffer between tile submits.
// queue.writeBuffer is queue-ordered with subsequent submits, so each tile
// dispatch reads the origin written immediately before it.
function writeGPUTileOrigin(originX, originY) {
    const buf = new Uint32Array(2);
    buf[0] = originX >>> 0;
    buf[1] = originY >>> 0;
    _gpu.device.queue.writeBuffer(_gpu.paramsBuf, GPU_TILE_ORIGIN_BYTE_OFFSET, buf.buffer);
}

function writeGPUTone(ximlen, yimlen, gradientVal, brightnessVal, drawFocusVal, focusVal, focusDepthVal) {
    const u = new Uint32Array(_gpu.toneBytes);
    const f = new Float32Array(_gpu.toneBytes);
    u[0] = ximlen; u[1] = yimlen; u[2] = drawFocusVal ? 1 : 0; u[3] = 0;
    f[4] = gradientVal; f[5] = brightnessVal; f[6] = focusVal; f[7] = focusDepthVal;
    _gpu.device.queue.writeBuffer(_gpu.toneBuf, 0, _gpu.toneBytes);
}

function uploadGPUPalette(pal) {
    // Palette is [[r,g,b], …] integers 0..255; shader multiplies by light_factor / shadow_darkness
    // so keep the same 0..255 range here.
    const data = new Float32Array(256 * 4);
    for (let i = 0; i < 256; i++) {
        const src = pal[i] || [0, 0, 0];
        data[i*4 + 0] = src[0];
        data[i*4 + 1] = src[1];
        data[i*4 + 2] = src[2];
        data[i*4 + 3] = 0;
    }
    _gpu.device.queue.writeBuffer(_gpu.paletteBuf, 0, data.buffer);
    _gpu.paletteUploaded = true;
}

function _writeSnapParams(count, horizon) {
    const u = new Uint32Array(_gpu.snapBytes);
    const f = new Float32Array(_gpu.snapBytes);
    u[0] = count; u[1] = 0; f[2] = horizon; u[3] = 0;
    _gpu.device.queue.writeBuffer(_gpu.snapParamsBuf, 0, _gpu.snapBytes);
}

// Build and cache CPU-side fill arrays keyed by HORIZON. These get re-used
// for both the full reset and the per-frame occLive clear.
function _ensureFillData(HORIZON) {
    if (_gpu.fillHorizon === HORIZON) return;
    const pc = _gpu.pixelCount;
    _gpu.fillHorizon = HORIZON;

    const occLiveBits = new Uint32Array(new Float32Array([HORIZON + HORIZON]).buffer)[0];
    const occLive = new Uint32Array(pc);
    occLive.fill(occLiveBits);
    _gpu.occLiveFillData = occLive.buffer;

    const occSnap = new Float32Array(pc);
    occSnap.fill(HORIZON);
    _gpu.occSnapFillData = occSnap.buffer;

    const maxAlphaFill = new Uint32Array(1);
    maxAlphaFill[0] = new Uint32Array(new Float32Array([1.0]).buffer)[0];
    _gpu.maxAlphaFillData = maxAlphaFill.buffer;
}

// Full reset: zero accum, seed occLive / occSnap to HORIZON sentinels, reset
// maxAlpha to bitcast(1.0). Uses encoder.clearBuffer + queue.writeBuffer
// instead of a compute dispatch because accumBuf's element count
// (pixelCount*4) can exceed WebGPU's maxComputeWorkgroupsPerDimension limit of
// 65535 at common canvas sizes (>1024x1024), which would silently fail
// validation and leave the buffer uncleared.
function resetGPUAccumulators(HORIZON) {
    _ensureFillData(HORIZON);
    const d = _gpu.device;
    const encoder = d.createCommandEncoder();
    encoder.clearBuffer(_gpu.accumBuf);
    d.queue.submit([encoder.finish()]);
    d.queue.writeBuffer(_gpu.occLiveBuf, 0, _gpu.occLiveFillData);
    d.queue.writeBuffer(_gpu.occSnapBuf, 0, _gpu.occSnapFillData);
    d.queue.writeBuffer(_gpu.maxAlphaBuf, 0, _gpu.maxAlphaFillData);
}

// Clear occLive at the start of every pass (fresh surface depths each pass).
// writeBuffer is queue-ordered with the subsequent submit, so the render
// dispatch sees the seeded values before any atomicMin runs.
function _clearOccLive(HORIZON) {
    _ensureFillData(HORIZON);
    _gpu.device.queue.writeBuffer(_gpu.occLiveBuf, 0, _gpu.occLiveFillData);
}

function _snapOcclusion(encoder, HORIZON) {
    _writeSnapParams(_gpu.pixelCount, HORIZON);
    const bg = _gpu.device.createBindGroup({
        layout: _gpu.snapPipeline.getBindGroupLayout(0),
        entries: [
            { binding: 0, resource: { buffer: _gpu.snapParamsBuf } },
            { binding: 1, resource: { buffer: _gpu.occLiveBuf } },
            { binding: 2, resource: { buffer: _gpu.occSnapBuf } },
        ],
    });
    const pass = encoder.beginComputePass();
    pass.setPipeline(_gpu.snapPipeline);
    pass.setBindGroup(0, bg);
    pass.dispatchWorkgroups(Math.ceil(_gpu.pixelCount / 64));
    pass.end();
}

// ---- per-frame dispatch ----

async function _gpuFrame(stateFn, onStats) {
    if (!_gpu || !_gpuRunning) return;
    // stateFn may call resetGPURender() which bumps _gpuGen; capture the
    // generation after it runs so a reset triggered by this frame still lets
    // the frame complete and reschedule the loop.
    const s = stateFn();
    if (!s) return;
    const genAtStart = _gpuGen;

    if (!_gpu.paletteUploaded || s.paletteDirty) {
        uploadGPUPalette(s.pallet);
        s.paletteDirty = false;
    }

    const passSeed = (Math.random() * 0xffffffff) >>> 0;
    writeGPUParams(s, passSeed);
    writeGPUTone(s.ximlen, s.yimlen, s.gradient, s.brightness,
                 s.drawFocus ? 1 : 0, s.focus, s.focus_depth);

    const d = _gpu.device;

    // Reset per-pass occlusion before encoding the render. queue.writeBuffer
    // is ordered with the subsequent queue.submit() so the dispatch reads the
    // seeded occLive values.
    _clearOccLive(s.HORIZON);

    // Render in tiles, each in its own command buffer, so no single GPU
    // command list runs long enough to trip Windows' GPU TDR. Atomic writes
    // to accum / occLive remain queue-ordered across tiles, so the visual
    // result matches a single full-image dispatch.
    const tilePixels = GPU_TILE_PIXELS;
    const wgPerTile  = tilePixels / 8;
    const tilesX     = Math.ceil(_gpu.width  / tilePixels);
    const tilesY     = Math.ceil(_gpu.height / tilePixels);
    let tileCount = 0;
    for (let ty = 0; ty < tilesY; ty++) {
        for (let tx = 0; tx < tilesX; tx++) {
            const originX = tx * tilePixels;
            const originY = ty * tilePixels;
            writeGPUTileOrigin(originX, originY);

            const tileWgX = Math.min(wgPerTile, Math.ceil((_gpu.width  - originX) / 8));
            const tileWgY = Math.min(wgPerTile, Math.ceil((_gpu.height - originY) / 8));

            const tileEnc  = d.createCommandEncoder();
            const tilePass = tileEnc.beginComputePass();
            tilePass.setPipeline(_gpu.renderPipeline);
            tilePass.setBindGroup(0, _gpu.renderBG);
            tilePass.dispatchWorkgroups(tileWgX, tileWgY);
            tilePass.end();
            d.queue.submit([tileEnc.finish()]);

            tileCount++;
            if (tileCount % GPU_TILES_PER_AWAIT === 0) {
                await d.queue.onSubmittedWorkDone();
                if (_gpuGen !== genAtStart || !_gpuRunning) return;
            }
        }
    }

    // Snap-occlusion + tonemap as the final small command buffer. Snap runs
    // every few passes so fog / goodPoints can filter against a populated
    // depth field (mirrors the CPU occSnap update every 50 passes).
    const finalEnc = d.createCommandEncoder();
    if (_gpuPassIndex % 8 === 0) {
        _snapOcclusion(finalEnc, s.HORIZON);
    }
    const view = _gpu.gpuCtx.getCurrentTexture().createView();
    const tonePass = finalEnc.beginRenderPass({
        colorAttachments: [{
            view: view,
            clearValue: { r: 0, g: 0, b: 0, a: 1 },
            loadOp: 'clear',
            storeOp: 'store',
        }],
    });
    tonePass.setPipeline(_gpu.tonePipeline);
    tonePass.setBindGroup(0, _gpu.toneBG);
    tonePass.draw(3);
    tonePass.end();

    // Stage maxAlpha for CPU readback. Cheap (4 bytes) and only when stats
    // are wired up.
    if (onStats) {
        finalEnc.copyBufferToBuffer(_gpu.maxAlphaBuf, 0, _gpu.maxAlphaReadBuf, 0, 4);
    }

    d.queue.submit([finalEnc.finish()]);

    await d.queue.onSubmittedWorkDone();
    if (_gpuGen !== genAtStart || !_gpuRunning) return;

    let peakAlpha = null;
    if (onStats) {
        await _gpu.maxAlphaReadBuf.mapAsync(GPUMapMode.READ, 0, 4);
        if (_gpuGen !== genAtStart || !_gpuRunning) {
            _gpu.maxAlphaReadBuf.unmap();
            return;
        }
        const mapped = _gpu.maxAlphaReadBuf.getMappedRange(0, 4);
        peakAlpha = new Float32Array(mapped, 0, 1)[0];
        _gpu.maxAlphaReadBuf.unmap();
    }

    // Blit the WebGPU canvas onto the main 2D canvas so the existing UI,
    // Save PNG, and control overlay keep working unchanged.
    _gpu.destCtx.drawImage(_gpu.gpuCanvas, 0, 0);
    _gpuPassIndex++;

    if (onStats) {
        const now = (new Date()).getTime();
        if (now - _gpuLastUpdateMs >= 1000) {
            const elapsed = (now - _gpuStartMs) / 1000.0;
            onStats({
                pass: _gpuPassIndex,
                elapsedSec: elapsed,
                pixelsPerSec: Math.floor((_gpu.pixelCount * _gpuPassIndex) / Math.max(elapsed, 0.001)),
                max_alpha: peakAlpha,
            });
            _gpuLastUpdateMs = now;
        }
    }

    // Sleep briefly, then schedule the next pass. The setTimeout yields back
    // to the browser so the GPU isn't hammered at full-tilt rAF cadence.
    setTimeout(function() {
        requestAnimationFrame(function() { _gpuFrame(stateFn, onStats); });
    }, GPU_INTERPASS_SLEEP_MS);
}

function startGPURender(stateFn, onStats) {
    if (!_gpu) return false;
    _gpuRunning = true;
    _gpuGen++;
    _gpuPassIndex = 0;
    _gpuStartMs = (new Date()).getTime();
    _gpuLastUpdateMs = _gpuStartMs;
    const s0 = stateFn();
    resetGPUAccumulators(s0.HORIZON);
    requestAnimationFrame(function() { _gpuFrame(stateFn, onStats); });
    return true;
}

function stopGPURender() { _gpuRunning = false; _gpuGen++; }

function resetGPURender(HORIZON) {
    if (!_gpu) return;
    _gpuGen++;
    _gpuPassIndex = 0;
    _gpuStartMs = (new Date()).getTime();
    _gpuLastUpdateMs = _gpuStartMs;
    resetGPUAccumulators(HORIZON);
}

function markGPUPaletteDirty() { if (_gpu) _gpu.paletteUploaded = false; }

function getGPUPassIndex() { return _gpuPassIndex; }
