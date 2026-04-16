/*
 * WebGL2 fragment-shader preview renderer for the Mandelbulb.
 *
 * Replaces the original CPU "first quick pass" with a single GPU pass that
 * produces a fully-lit image in well under a second. The CPU progressive
 * renderer still owns "final quality" mode after the preview timeout.
 *
 * Reads the same global parameters from mandelbulb.js that the CPU pipeline
 * uses (zoom, xcen, ycen, CameraMatrix, iterations, formula, power, etc.) so
 * the preview faithfully tracks whatever the user has configured.
 */

// ---------- GLSL ----------

var GPU_PREVIEW_VERT = `#version 300 es
void main() {
    // Fullscreen triangle via gl_VertexID, no VBO needed.
    vec2 pos = vec2((gl_VertexID == 1) ? 3.0 : -1.0,
                    (gl_VertexID == 2) ? 3.0 : -1.0);
    gl_Position = vec4(pos, 0.0, 1.0);
}
`;

var GPU_PREVIEW_FRAG = `#version 300 es
precision highp float;
precision highp int;

uniform vec2  u_resolution;
uniform float u_zoom, u_xcen, u_ycen, u_cameraPersp;
uniform float u_azimuth, u_power;
uniform int   u_iterations, u_formula;
uniform float u_M[9];                    // CameraMatrix, row-major
uniform vec3  u_lightVector;
uniform float u_min_y, u_max_y, u_stepDetail, u_root_zoom;
uniform float u_ray_step;
uniform int   u_RAY_STEPS;
uniform float u_ambientLight, u_primaryLight, u_shadowDarkness;
uniform float u_brightness, u_opacity;
uniform sampler2D u_palette;

out vec4 fragColor;

const int MAX_ITER        = 200;
const int MAX_DEPTH_STEPS = 1024;
const int MAX_RAY_STEPS   = 256;

// Cheap per-pixel hash in place of Math.random(); deterministic, adequate for
// the light stochastic jitter the CPU renderer uses (rndFuzzy and ray jitter).
float hash21(vec2 p) {
    p = fract(p * vec2(123.34, 456.21));
    p += dot(p, p + 45.32);
    return fract(p.x * p.y);
}

vec3 rotateVec(float x, float y, float z) {
    return vec3(
        u_M[0]*x + u_M[1]*y + u_M[2]*z,
        u_M[3]*x + u_M[4]*y + u_M[5]*z,
        u_M[6]*x + u_M[7]*y + u_M[8]*z);
}

// Returns vec2(iter, pixelColor). iter == u_iterations means "inside".
vec2 insideFractal(vec3 p) {
    vec3 c = p;
    int iter = 0;
    float pixelColor = 0.0;
    float r = 0.0;
    for (int i = 0; i < MAX_ITER; i++) {
        if (i >= u_iterations) break;
        r = length(p);
        float theta_power = atan(p.y, p.x) * u_power;
        float r_power = pow(r, u_power);
        float phi;
        if (u_formula == 0) {
            phi = asin(p.z / r);
            float pc = cos(phi * u_power);
            float ct = cos(theta_power), st = sin(theta_power);
            float sp = sin(phi * u_power);
            p = vec3(r_power*ct*pc, r_power*st*pc, r_power*sp*u_azimuth) + c;
        } else {
            phi = atan(length(p.xy), p.z);
            float ps = sin(phi * u_power);
            float ct = cos(theta_power), st = sin(theta_power);
            float cp = cos(phi * u_power);
            p = vec3(r_power*ct*ps, r_power*st*ps, r_power*cp*u_azimuth) + c;
        }
        pixelColor = phi / 3.0;
        iter++;
        if (r >= 8.0) break;
    }
    return vec2(float(iter), pixelColor);
}

// Mirrors CPU calculateRay: the step accumulates quadratically (origPoint += step*i each iter).
float calculateRay(vec3 origPoint, int steps, vec3 step, float bright, float rndFuzzy) {
    step *= rndFuzzy;
    vec3 p = origPoint;
    for (int i = 1; i < MAX_RAY_STEPS; i++) {
        if (i >= steps) break;
        p += step * float(i);
        vec2 res = insideFractal(p);
        int it = int(res.x);
        if (it == u_iterations) return 0.0;   // shadow
        if (it < 3) break;                    // ray escaped
    }
    return bright;
}

float calculateRays(vec3 p, float rndFuzzy) {
    float lf = 1.0;
    float rs = u_ray_step;
    vec2 fc = gl_FragCoord.xy;
    float h1 = hash21(fc + vec2(1.7, 2.3));
    float h2 = hash21(fc + vec2(7.1, 3.5));
    float h3 = hash21(fc + vec2(5.9, 8.2));
    float n0 = rs / (15.0*h1 + 1.0);
    float n1 = rs / (15.0*h2 + 1.0);
    float n2 = rs / (15.0*h3 + 1.0);
    vec3 lv = u_lightVector;

    lf += calculateRay(p, u_RAY_STEPS,   vec3(-n0,-n1,-n2),                                u_ambientLight, rndFuzzy);
    lf += calculateRay(p, u_RAY_STEPS,   vec3( n0, n1, n2),                                u_ambientLight, rndFuzzy);
    lf += calculateRay(p, u_RAY_STEPS,   vec3( n0,-n1,-rs),                                u_ambientLight, rndFuzzy);
    lf += calculateRay(p, u_RAY_STEPS,   vec3(n0*lv.x*9.0, n1*lv.y*9.0, n2*lv.z*9.0),      u_ambientLight, rndFuzzy);
    lf += calculateRay(p, u_RAY_STEPS*4, vec3(rs*lv.x, rs*lv.y, rs*lv.z),                  u_primaryLight, rndFuzzy);
    return lf;
}

void main() {
    // Match CPU renderer's top-left screen origin (WebGL is bottom-up).
    vec2 pix = vec2(gl_FragCoord.x, u_resolution.y - gl_FragCoord.y);
    float x2 = (pix.x / u_resolution.x - 0.5) * u_zoom - u_xcen;
    float z2 = (pix.y / u_resolution.y - 0.5) * u_zoom - u_ycen;

    float stepAmount = u_stepDetail * u_root_zoom;
    vec3 outColor = vec3(0.0);

    for (int s = 0; s < MAX_DEPTH_STEPS; s++) {
        float y = u_min_y + float(s) * stepAmount;
        if (y >= u_max_y) break;

        float persp = 1.0 + y * u_cameraPersp;
        vec3 p = rotateVec(x2 / persp, y, z2 / persp);
        vec2 res = insideFractal(p);

        if (int(res.x) == u_iterations) {
            float rndFuzzy = max(u_opacity * hash21(gl_FragCoord.xy), 0.4);
            float lightFactor = 1.0 + calculateRays(p, rndFuzzy);
            int ci = min(int(abs(res.y) * 255.0), 255);
            vec3 cn = texelFetch(u_palette, ivec2(ci, 0), 0).rgb;   // 0..1
            outColor = clamp(cn * lightFactor * u_brightness / u_shadowDarkness, 0.0, 1.0);
            break;
        }
    }
    fragColor = vec4(outColor, 1.0);
}
`;

// ---------- WebGL setup ----------

var _gpuPreviewState = null;  // cached {canvas, gl, prog, vao, paletteTex, paletteRef, uloc}

function _gpuCompile(gl, type, src, label) {
    var s = gl.createShader(type);
    gl.shaderSource(s, src);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
        console.error('gpu-preview ' + label + ' compile error:\n' + gl.getShaderInfoLog(s));
        gl.deleteShader(s);
        return null;
    }
    return s;
}

function _gpuPreviewInit(width, height) {
    var c = document.createElement('canvas');
    c.width = width; c.height = height;
    var gl = c.getContext('webgl2', { preserveDrawingBuffer: true, antialias: false, alpha: false });
    if (!gl) return null;

    var vs = _gpuCompile(gl, gl.VERTEX_SHADER,   GPU_PREVIEW_VERT, 'vertex');
    var fs = _gpuCompile(gl, gl.FRAGMENT_SHADER, GPU_PREVIEW_FRAG, 'fragment');
    if (!vs || !fs) return null;

    var prog = gl.createProgram();
    gl.attachShader(prog, vs); gl.attachShader(prog, fs);
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {
        console.error('gpu-preview link error:\n' + gl.getProgramInfoLog(prog));
        return null;
    }
    var vao = gl.createVertexArray();  // fullscreen triangle uses gl_VertexID only

    // Pre-resolve uniform locations once
    var names = ['u_resolution','u_zoom','u_xcen','u_ycen','u_cameraPersp','u_azimuth',
        'u_power','u_iterations','u_formula','u_M[0]','u_lightVector','u_min_y','u_max_y',
        'u_stepDetail','u_root_zoom','u_ray_step','u_RAY_STEPS','u_ambientLight',
        'u_primaryLight','u_shadowDarkness','u_brightness','u_opacity','u_palette'];
    var uloc = {};
    for (var i = 0; i < names.length; i++) uloc[names[i]] = gl.getUniformLocation(prog, names[i]);

    return { canvas: c, gl: gl, prog: prog, vao: vao, paletteTex: null, paletteRef: null, uloc: uloc };
}

function _gpuUpdatePaletteTexture(st, palette) {
    if (st.paletteRef === palette && st.paletteTex) return;
    var gl = st.gl;
    var bytes = new Uint8Array(256 * 4);
    for (var i = 0; i < 256; i++) {
        bytes[i*4  ] = palette[i][0];
        bytes[i*4+1] = palette[i][1];
        bytes[i*4+2] = palette[i][2];
        bytes[i*4+3] = 255;
    }
    if (!st.paletteTex) st.paletteTex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, st.paletteTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 256, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, bytes);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S,     gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T,     gl.CLAMP_TO_EDGE);
    st.paletteRef = palette;
}

// ---------- public entry point ----------

/*
 * Render a fully-lit Mandelbulb preview onto destCtx using the current global
 * render parameters from mandelbulb.js. Returns true on success, false on
 * any failure (no WebGL2, shader error, etc.) so the caller can fall back
 * to the CPU preview path.
 */
function renderPreviewGPU(destCanvas, destCtx) {
    try {
        if (!_gpuPreviewState
                || _gpuPreviewState.canvas.width  !== destCanvas.width
                || _gpuPreviewState.canvas.height !== destCanvas.height) {
            _gpuPreviewState = _gpuPreviewInit(destCanvas.width, destCanvas.height);
            if (!_gpuPreviewState) return false;
        }
        var st = _gpuPreviewState, gl = st.gl, u = st.uloc;

        _gpuUpdatePaletteTexture(st, pallet);

        // Flatten CameraMatrix (jagged 3x3) into row-major 9-float array.
        var M = CameraMatrix;
        var M9 = new Float32Array([
            M[0][0], M[0][1], M[0][2],
            M[1][0], M[1][1], M[1][2],
            M[2][0], M[2][1], M[2][2]]);

        gl.viewport(0, 0, st.canvas.width, st.canvas.height);
        gl.useProgram(st.prog);
        gl.bindVertexArray(st.vao);
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, st.paletteTex);
        gl.uniform1i(u['u_palette'], 0);

        gl.uniform2f(u['u_resolution'], st.canvas.width, st.canvas.height);
        gl.uniform1f(u['u_zoom'], zoom);
        gl.uniform1f(u['u_xcen'], xcen);
        gl.uniform1f(u['u_ycen'], ycen);
        gl.uniform1f(u['u_cameraPersp'], cameraPersp);
        gl.uniform1f(u['u_azimuth'], azimuth);
        gl.uniform1f(u['u_power'], power);
        gl.uniform1i(u['u_iterations'], iterations);
        gl.uniform1i(u['u_formula'], formula);
        gl.uniform1fv(u['u_M[0]'], M9);
        gl.uniform3f(u['u_lightVector'], LightVector[0], LightVector[1], LightVector[2]);
        gl.uniform1f(u['u_min_y'], min_y);
        gl.uniform1f(u['u_max_y'], max_y);
        gl.uniform1f(u['u_stepDetail'], stepDetail);
        gl.uniform1f(u['u_root_zoom'], root_zoom);
        gl.uniform1f(u['u_ray_step'], ray_step);
        gl.uniform1i(u['u_RAY_STEPS'], RAY_STEPS);
        gl.uniform1f(u['u_ambientLight'], AMBIENT_LIGHT);
        gl.uniform1f(u['u_primaryLight'], primary_light);
        gl.uniform1f(u['u_shadowDarkness'], shadow_darkness);
        gl.uniform1f(u['u_brightness'], brightness);
        gl.uniform1f(u['u_opacity'], opacity);

        gl.drawArrays(gl.TRIANGLES, 0, 3);

        var err = gl.getError();
        if (err !== gl.NO_ERROR) { console.warn('gpu-preview gl error:', err); return false; }

        destCtx.drawImage(st.canvas, 0, 0);
        return true;
    } catch (e) {
        console.warn('gpu-preview threw, falling back to CPU:', e);
        return false;
    }
}
