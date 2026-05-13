/*
 * The Mandelbrot Set, in HTML5 canvas and javascript.
 * https://github.com/rsweny/mandelbrot-js
 *
 * Copyright (C) 2012 Christian Stigen Larsen
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may
 * not use this file except in compliance with the License.  You may obtain
 * a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
 * License for the specific language governing permissions and limitations
 * under the License.
 *
 */

var zoomStart = 3.4;
var zoom = zoomStart;
var lookAtDefault = [-0.6, 0.01];
var lookAt = lookAtDefault;
// Parallel 256-bit fixed-point view of `lookAt`.  Built lazily from the f64
// pair the first time it is needed (see ensureLookAtFix).  Drag handlers
// update this directly so sub-LSB pan deltas survive at deep zooms.
var lookAtFix = null;
var xRange = [0, 0];
var yRange = [0, 0];
var escapeRadius = 80.0;
var interiorColor = [0, 0, 0, 255];
var reInitCanvas = true; // Whether to reload canvas size, etc
var dragToZoom = true;
var colors = [[0, 0, 0, 0]];
var renderId = 0; // To zoom before current render is finished

let doOrbit = false;
let orbitType = 3;
const orbit_trap = 0.5;
const orbitRealPoint = 0.5;
const orbitImgPoint = -0.25;
const doOrbitAverage = false;
let colourSeed = 1 + Math.random() * 0.5;

// WebGPU accelerates the per-pixel Mandelbrot iterations.  Colour mapping
// stays on the CPU so the existing palette controls and orbit-trap shaping
// stay in sync with the original renderer.
var mandelGPU = null;
var mandelGPUDisabled = false;
var mandelGPUStatus = 'CPU';
var mandelGPUInFlight = false;

var MANDEL_PIXEL_WGSL = `
struct Params {
  width: u32,
  height: u32,
  iterations: u32,
  samples: u32,
  doOrbit: u32,
  orbitType: u32,
  doOrbitAverage: u32,
  _pad: u32,
  xStart: f32,
  yStart: f32,
  dx: f32,
  dy: f32,
  escapeRadius: f32,
  orbitTrap: f32,
  orbitRealPoint: f32,
  orbitImgPoint: f32,
};

@group(0) @binding(0) var<uniform> P: Params;
// x=n, y=Tr, z=Ti, w=distance
@group(0) @binding(1) var<storage, read_write> outPixels: array<vec4<f32>>;

fn hash01(v: u32) -> f32 {
  var x = v;
  x = ((x >> 16u) ^ x) * 0x45d9f3bu;
  x = ((x >> 16u) ^ x) * 0x45d9f3bu;
  x = (x >> 16u) ^ x;
  return f32(x & 0x00ffffffu) / 16777216.0;
}

fn calcDistGpu(x1: f32, y1: f32, x2: f32, y2: f32) -> f32 {
  let dx = x1 - x2;
  let dy = y1 - y2;
  return sqrt(dx * dx + dy * dy);
}

fn iterateMandelbrotGpu(Cr: f32, Ci: f32) -> vec4<f32> {
  var Zr: f32 = 0.0;
  var Zi: f32 = 0.0;
  var Tr: f32 = 0.0;
  var Ti: f32 = 0.0;
  var n: u32 = 0u;
  var distance: f32 = 0.0;

  loop {
    if (n >= P.iterations) { break; }
    if ((Tr + Ti) > P.escapeRadius) { break; }

    Zi = 2.0 * Zr * Zi + Ci;
    Zr = Tr - Ti + Cr;
    Tr = Zr * Zr;
    Ti = Zi * Zi;

    if (P.doOrbit != 0u) {
      var dist: f32 = 0.0;
      switch (P.orbitType) {
        case 0u: {
          dist = abs((Zi + P.orbitImgPoint) * (Zr + P.orbitRealPoint));
        }
        case 1u: {
          dist = abs(Zi - P.orbitImgPoint);
          let dist2 = abs((Zr - P.orbitRealPoint) * (Zi - P.orbitImgPoint));
          if (dist2 > dist) { dist = dist2; }
        }
        case 2u: {
          dist = abs((Zi - P.orbitImgPoint) * (Zr - P.orbitRealPoint));
        }
        case 3u: {
          dist = calcDistGpu(P.orbitRealPoint, P.orbitImgPoint, Zr, Zi);
        }
        default: {
          dist = atan(abs((Zi - P.orbitImgPoint) / (Zr - P.orbitRealPoint)));
        }
      }

      if (dist < P.orbitTrap) {
        if (P.doOrbitAverage != 0u) {
          if (distance == 0.0) { distance = dist; }
          else { distance = dist * 0.25 + distance * 0.75; }
        } else {
          if (distance == 0.0 || dist < distance) { distance = dist; }
        }
      }
    }
    n = n + 1u;
  }

  return vec4<f32>(f32(n), Tr, Ti, distance);
}

@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  if (gid.x >= P.width || gid.y >= P.height * P.samples) { return; }

  let x = gid.x;
  let sampleY = gid.y;
  let y = sampleY / P.samples;
  let sampleIdx = sampleY - y * P.samples;
  let pixelIdx = y * P.width + x;
  let idx = pixelIdx * P.samples + sampleIdx;

  var Cr = P.xStart + f32(x) * P.dx;
  var Ci = P.yStart + f32(y) * P.dy;
  if (P.samples > 1u) {
    let rx = hash01(pixelIdx * 1664525u + sampleIdx * 1013904223u + 17u);
    let ry = hash01(pixelIdx * 22695477u + sampleIdx * 1103515245u + 29u);
    Cr = Cr - rx * P.dx * 0.5;
    Ci = Ci - ry * P.dy * 0.5;
  }

  outPixels[idx] = iterateMandelbrotGpu(Cr, Ci);
}
`;

// Perturbation-theory shader: each pixel iterates a small delta from a
// high-precision reference orbit computed on the CPU with BigInt fixed
// point.  This sidesteps the f32 precision wall of WGSL at deep zooms.
var MANDEL_PERTURB_WGSL = `
struct Params {
  width: u32,
  height: u32,
  iterations: u32,
  samples: u32,
  doOrbit: u32,
  orbitType: u32,
  doOrbitAverage: u32,
  refOrbitLen: u32,
  dcrStart: f32,
  dciStart: f32,
  dx: f32,
  dy: f32,
  escapeRadius: f32,
  orbitTrap: f32,
  orbitRealPoint: f32,
  orbitImgPoint: f32,
  glitchEps: f32,
  _pad0: f32,
  _pad1: f32,
  _pad2: f32,
};

@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> refOrbit: array<vec2<f32>>;
@group(0) @binding(2) var<storage, read_write> outPixels: array<vec4<f32>>;

fn hash01(v: u32) -> f32 {
  var x = v;
  x = ((x >> 16u) ^ x) * 0x45d9f3bu;
  x = ((x >> 16u) ^ x) * 0x45d9f3bu;
  x = (x >> 16u) ^ x;
  return f32(x & 0x00ffffffu) / 16777216.0;
}

fn calcDistGpu(x1: f32, y1: f32, x2: f32, y2: f32) -> f32 {
  let dx = x1 - x2;
  let dy = y1 - y2;
  return sqrt(dx * dx + dy * dy);
}

fn iteratePerturbGpu(dcr: f32, dci: f32) -> vec4<f32> {
  var dzr: f32 = 0.0;
  var dzi: f32 = 0.0;
  var Tr: f32 = 0.0;
  var Ti: f32 = 0.0;
  var n: u32 = 0u;   // total iterations, used for colouring
  var k: u32 = 0u;   // index into reference orbit; resets on rebase
  var distance: f32 = 0.0;
  let lastK = P.refOrbitLen - 1u;

  loop {
    if (n >= P.iterations) { break; }
    if ((Tr + Ti) > P.escapeRadius) { break; }
    if (k >= lastK) { break; }

    let Zr = refOrbit[k].x;
    let Zi = refOrbit[k].y;

    // delta_z_{n+1} = 2*Z_k*delta_z_n + delta_z_n^2 + delta_c
    let new_dzr = 2.0 * (Zr * dzr - Zi * dzi) + (dzr * dzr - dzi * dzi) + dcr;
    let new_dzi = 2.0 * (Zr * dzi + Zi * dzr + dzr * dzi) + dci;
    dzr = new_dzr;
    dzi = new_dzi;

    let Zr1 = refOrbit[k + 1u].x;
    let Zi1 = refOrbit[k + 1u].y;
    let zr = Zr1 + dzr;
    let zi = Zi1 + dzi;
    Tr = zr * zr;
    Ti = zi * zi;

    if (P.doOrbit != 0u) {
      var dist: f32 = 0.0;
      switch (P.orbitType) {
        case 0u: {
          dist = abs((zi + P.orbitImgPoint) * (zr + P.orbitRealPoint));
        }
        case 1u: {
          dist = abs(zi - P.orbitImgPoint);
          let dist2 = abs((zr - P.orbitRealPoint) * (zi - P.orbitImgPoint));
          if (dist2 > dist) { dist = dist2; }
        }
        case 2u: {
          dist = abs((zi - P.orbitImgPoint) * (zr - P.orbitRealPoint));
        }
        case 3u: {
          dist = calcDistGpu(P.orbitRealPoint, P.orbitImgPoint, zr, zi);
        }
        default: {
          dist = atan(abs((zi - P.orbitImgPoint) / (zr - P.orbitRealPoint)));
        }
      }

      if (dist < P.orbitTrap) {
        if (P.doOrbitAverage != 0u) {
          if (distance == 0.0) { distance = dist; }
          else { distance = dist * 0.25 + distance * 0.75; }
        } else {
          if (distance == 0.0 || dist < distance) { distance = dist; }
        }
      }
    }

    // Pauldelbrot rebasing: when |z|^2 drops well below |Z_ref|^2 the
    // delta has lost significance against the reference; adopt the full
    // position as the new delta and restart from refOrbit[0] (which is 0,
    // so the perturbation formula reduces to dz' = dz^2 + dc and the same
    // reference orbit can be reused).  Also rebase if we are about to run
    // out of reference orbit so the iteration can continue.
    let refMag2 = Zr1 * Zr1 + Zi1 * Zi1;
    let nextK = k + 1u;
    if ((Tr + Ti) < P.glitchEps * refMag2 || nextK >= lastK) {
      dzr = zr;
      dzi = zi;
      k = 0u;
    } else {
      k = nextK;
    }
    n = n + 1u;
  }

  return vec4<f32>(f32(n), Tr, Ti, distance);
}

@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  if (gid.x >= P.width || gid.y >= P.height * P.samples) { return; }

  let x = gid.x;
  let sampleY = gid.y;
  let y = sampleY / P.samples;
  let sampleIdx = sampleY - y * P.samples;
  let pixelIdx = y * P.width + x;
  let idx = pixelIdx * P.samples + sampleIdx;

  var dcr = P.dcrStart + f32(x) * P.dx;
  var dci = P.dciStart + f32(y) * P.dy;
  if (P.samples > 1u) {
    let rx = hash01(pixelIdx * 1664525u + sampleIdx * 1013904223u + 17u);
    let ry = hash01(pixelIdx * 22695477u + sampleIdx * 1103515245u + 29u);
    dcr = dcr - rx * P.dx * 0.5;
    dci = dci - ry * P.dy * 0.5;
  }

  outPixels[idx] = iteratePerturbGpu(dcr, dci);
}
`;

// BigInt fixed-point parameters for the reference-orbit computation.
// 256 fractional bits gives ~76 decimal digits, plenty of headroom for
// any zoom level the f64-encoded lookAt coordinates can express.
var MANDEL_BIGINT_SCALE_BITS = 256n;
var MANDEL_BIGINT_SCALE_BITS_NUM = 256;
var MANDEL_BIGINT_DIV_SCALE = Math.pow(2, -MANDEL_BIGINT_SCALE_BITS_NUM);
var MANDEL_PERTURB_ZOOM_THRESHOLD = 1e-3;

// Pauldelbrot glitch criterion: rebase when |z|^2 drops below
// mandelGlitchTolerance * |Z_ref|^2.  1e-3 is roughly matched to f32
// mantissa precision; lower values rebase more aggressively.
var mandelGlitchTolerance = 1e-3;

var mandelGPUPerturb = null;
var mandelRefOrbitCache = null;

var MANDEL_GPU_PARAMS_BYTES = 64;
var MANDEL_PERTURB_PARAMS_BYTES = 80;
var MANDEL_GPU_WORKGROUP_X = 16;
var MANDEL_GPU_WORKGROUP_Y = 16;

function isMandelbrotWebGPUSupported() {
  return typeof navigator !== 'undefined' && !!navigator.gpu;
}

function destroyMandelbrotGPU() {
  if (!mandelGPU) return;
  var bufs = ['paramsBuf', 'outputBuf', 'readBuf'];
  for (var i = 0; i < bufs.length; i++) {
    try { mandelGPU[bufs[i]].destroy(); } catch (e) { }
  }
  mandelGPU = null;
}

async function initMandelbrotGPU(resultCount) {
  if (mandelGPUDisabled || !isMandelbrotWebGPUSupported()) return null;
  if (mandelGPU && mandelGPU.resultCapacity >= resultCount) return mandelGPU;

  destroyMandelbrotGPU();
  try {
    var adapter = await navigator.gpu.requestAdapter();
    if (!adapter) throw new Error('No WebGPU adapter');
    var device = await adapter.requestDevice();
    var outputBytes = resultCount * 4 * 4;

    var paramsBuf = device.createBuffer({ size: MANDEL_GPU_PARAMS_BYTES, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    var outputBuf = device.createBuffer({ size: outputBytes, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC });
    var readBuf = device.createBuffer({ size: outputBytes, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
    var module = device.createShaderModule({ code: MANDEL_PIXEL_WGSL });
    var pipeline = await device.createComputePipelineAsync({ layout: 'auto', compute: { module: module, entryPoint: 'main' } });
    var bindGroup = device.createBindGroup({
      layout: pipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: paramsBuf } },
        { binding: 1, resource: { buffer: outputBuf } },
      ],
    });

    mandelGPU = {
      device: device,
      resultCapacity: resultCount,
      outputBytes: outputBytes,
      paramsBuf: paramsBuf,
      outputBuf: outputBuf,
      readBuf: readBuf,
      pipeline: pipeline,
      bindGroup: bindGroup,
      paramsBytes: new ArrayBuffer(MANDEL_GPU_PARAMS_BYTES),
    };
    mandelGPUStatus = 'GPU';
    return mandelGPU;
  } catch (e) {
    console.warn('Mandelbrot WebGPU path disabled; falling back to CPU.', e);
    mandelGPUDisabled = true;
    mandelGPUStatus = 'CPU fallback';
    destroyMandelbrotGPU();
    return null;
  }
}

function writeMandelbrotGPUParams(gpu, width, height, dx, dy, sampleCount, steps) {
  var u = new Uint32Array(gpu.paramsBytes);
  var f = new Float32Array(gpu.paramsBytes);
  u[0] = width >>> 0;
  u[1] = height >>> 0;
  u[2] = (parseInt(steps, 10) || 0) >>> 0;
  u[3] = sampleCount >>> 0;
  u[4] = doOrbit ? 1 : 0;
  u[5] = orbitType >>> 0;
  u[6] = doOrbitAverage ? 1 : 0;
  u[7] = 0;
  f[8] = xRange[0];
  f[9] = yRange[0];
  f[10] = dx;
  f[11] = dy;
  f[12] = escapeRadius;
  f[13] = orbit_trap;
  f[14] = orbitRealPoint;
  f[15] = orbitImgPoint;
  gpu.device.queue.writeBuffer(gpu.paramsBuf, 0, gpu.paramsBytes);
}

async function calcMandelbrotPixelsGPU(width, height, dx, dy, samples, steps) {
  var sampleCount = Math.max(1, Math.min(4, parseInt(samples, 10) || 1));
  var resultCount = width * height * sampleCount;
  var gpu = await initMandelbrotGPU(resultCount);
  if (!gpu) return null;

  writeMandelbrotGPUParams(gpu, width, height, dx, dy, sampleCount, steps);
  var outputBytes = resultCount * 4 * 4;
  var d = gpu.device;
  var encoder = d.createCommandEncoder();
  var pass = encoder.beginComputePass();
  pass.setPipeline(gpu.pipeline);
  pass.setBindGroup(0, gpu.bindGroup);
  pass.dispatchWorkgroups(
    Math.ceil(width / MANDEL_GPU_WORKGROUP_X),
    Math.ceil((height * sampleCount) / MANDEL_GPU_WORKGROUP_Y)
  );
  pass.end();
  encoder.copyBufferToBuffer(gpu.outputBuf, 0, gpu.readBuf, 0, outputBytes);
  d.queue.submit([encoder.finish()]);
  await d.queue.onSubmittedWorkDone();
  await gpu.readBuf.mapAsync(GPUMapMode.READ, 0, outputBytes);

  var range = gpu.readBuf.getMappedRange(0, outputBytes);
  return {
    pixels: new Float32Array(range),
    samples: sampleCount,
    release: function () {
      try { gpu.readBuf.unmap(); } catch (e) { }
    },
  };
}

// Convert a JS Number to a fixed-point BigInt with MANDEL_BIGINT_SCALE_BITS
// fractional bits.  Decomposes the IEEE 754 representation directly so the
// f64 mantissa is preserved without going through Math.round, which would
// be lossy for values whose magnitude exceeds 2^53.
function mandelNumberToFix(x) {
  if (x === 0 || !isFinite(x)) return 0n;
  var negative = x < 0;
  var absX = negative ? -x : x;
  var buf = new ArrayBuffer(8);
  new Float64Array(buf)[0] = absX;
  var u = new Uint32Array(buf);
  var lo = u[0] >>> 0;
  var hi = u[1] >>> 0;
  var rawExp = (hi >>> 20) & 0x7FF;
  var mantissa = (BigInt(hi & 0xFFFFF) << 32n) | BigInt(lo);
  var exponent;
  if (rawExp === 0) {
    exponent = -1022;
  } else {
    mantissa |= 1n << 52n;
    exponent = rawExp - 1023;
  }
  var shift = BigInt(exponent - 52) + MANDEL_BIGINT_SCALE_BITS;
  var result;
  if (shift >= 0n) result = mantissa << shift;
  else result = mantissa >> -shift;
  return negative ? -result : result;
}

// Inverse of mandelNumberToFix.  BigInt -> Number conversion already takes
// the top 53 bits, and any Mandelbrot coordinate fits comfortably below f64's
// exponent range, so this is a single multiply.
function mandelFixToNumber(fix) {
  if (fix === 0n) return 0.0;
  return Number(fix) * MANDEL_BIGINT_DIV_SCALE;
}

// Hex serialisation of a fix-point BigInt for the URL hash.  Exact and
// compact (one hex digit per 4 bits, ~64 chars per coord at 256-bit scale).
function mandelFixToHex(fix) {
  if (fix < 0n) return '-' + (-fix).toString(16);
  return fix.toString(16);
}
function mandelHexToFix(s) {
  if (!s) return 0n;
  var neg = s.charAt(0) === '-';
  var body = neg ? s.slice(1) : s;
  if (body.length === 0) return 0n;
  var v = BigInt('0x' + body);
  return neg ? -v : v;
}

// Lazily build the high-precision view of `lookAt` from its f64 mirror.
function ensureLookAtFix() {
  if (lookAtFix === null) {
    lookAtFix = [mandelNumberToFix(lookAt[0]), mandelNumberToFix(lookAt[1])];
  }
  return lookAtFix;
}

// Set both the high-precision and f64 views of the centre coordinate.
function setLookAtFix(fx, fy) {
  lookAtFix = [fx, fy];
  lookAt = [mandelFixToNumber(fx), mandelFixToNumber(fy)];
}

// Iterate Z_{n+1} = Z_n^2 + C in BigInt fixed point at the reference centre
// and capture the orbit as f32 pairs for the GPU.  cxBig/cyBig are taken
// straight in fixed-point so the centre keeps any sub-LSB precision the f64
// `lookAt` mirror would have lost.  Stops when the orbit escapes or when we
// reach maxIter
function mandelComputeReferenceOrbit(cxBig, cyBig, maxIter, escRadius) {
  var scaleBits = MANDEL_BIGINT_SCALE_BITS;
  var divScale = MANDEL_BIGINT_DIV_SCALE;
  var escBig = mandelNumberToFix(escRadius);
  var orbit = new Float32Array(2 * (maxIter + 1));
  var Zx = 0n;
  var Zy = 0n;
  orbit[0] = 0;
  orbit[1] = 0;
  var n = 0;
  for (var i = 1; i <= maxIter; i++) {
    var Zx2 = (Zx * Zx) >> scaleBits;
    var Zy2 = (Zy * Zy) >> scaleBits;
    var newZx = Zx2 - Zy2 + cxBig;
    var newZy = ((Zx * Zy) >> (scaleBits - 1n)) + cyBig;
    Zx = newZx;
    Zy = newZy;
    orbit[2 * i] = Number(Zx) * divScale;
    orbit[2 * i + 1] = Number(Zy) * divScale;
    n = i;
    if (Zx2 + Zy2 > escBig) break;
  }
  return { orbit: orbit, length: n + 1, escaped: n < maxIter };
}

// cxBig/cyBig are BigInt fixed-point coordinates; `===` between BigInts is
// value-based so the cache key works without any extra normalisation.
function getReferenceOrbit(cxBig, cyBig, steps, esc) {
  var c = mandelRefOrbitCache;
  if (c && c.cx === cxBig && c.cy === cyBig && c.steps === steps && c.esc === esc) {
    return c;
  }
  var res = mandelComputeReferenceOrbit(cxBig, cyBig, steps, esc);
  mandelRefOrbitCache = {
    cx: cxBig, cy: cyBig, steps: steps, esc: esc,
    orbit: res.orbit, length: res.length, escaped: res.escaped,
  };
  return mandelRefOrbitCache;
}

function destroyMandelbrotPerturbGPU() {
  if (!mandelGPUPerturb) return;
  var bufs = ['paramsBuf', 'orbitBuf', 'outputBuf', 'readBuf'];
  for (var i = 0; i < bufs.length; i++) {
    try { mandelGPUPerturb[bufs[i]].destroy(); } catch (e) { }
  }
  mandelGPUPerturb = null;
}

async function initMandelbrotPerturbGPU(resultCount, orbitCount) {
  if (mandelGPUDisabled || !isMandelbrotWebGPUSupported()) return null;
  if (mandelGPUPerturb &&
      mandelGPUPerturb.resultCapacity >= resultCount &&
      mandelGPUPerturb.orbitCapacity >= orbitCount) return mandelGPUPerturb;

  destroyMandelbrotPerturbGPU();
  try {
    var adapter = await navigator.gpu.requestAdapter();
    if (!adapter) throw new Error('No WebGPU adapter');
    var device = await adapter.requestDevice();
    var orbitCap = Math.max(1024, orbitCount);
    var outputBytes = resultCount * 4 * 4;
    var orbitBytes = orbitCap * 8;

    var paramsBuf = device.createBuffer({ size: MANDEL_PERTURB_PARAMS_BYTES, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    var orbitBuf = device.createBuffer({ size: orbitBytes, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });
    var outputBuf = device.createBuffer({ size: outputBytes, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC });
    var readBuf = device.createBuffer({ size: outputBytes, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
    var module = device.createShaderModule({ code: MANDEL_PERTURB_WGSL });
    var pipeline = await device.createComputePipelineAsync({ layout: 'auto', compute: { module: module, entryPoint: 'main' } });
    var bindGroup = device.createBindGroup({
      layout: pipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: paramsBuf } },
        { binding: 1, resource: { buffer: orbitBuf } },
        { binding: 2, resource: { buffer: outputBuf } },
      ],
    });

    mandelGPUPerturb = {
      device: device,
      resultCapacity: resultCount,
      orbitCapacity: orbitCap,
      outputBytes: outputBytes,
      paramsBuf: paramsBuf,
      orbitBuf: orbitBuf,
      outputBuf: outputBuf,
      readBuf: readBuf,
      pipeline: pipeline,
      bindGroup: bindGroup,
      paramsBytes: new ArrayBuffer(MANDEL_PERTURB_PARAMS_BYTES),
    };
    return mandelGPUPerturb;
  } catch (e) {
    console.warn('Mandelbrot WebGPU perturbation path failed; falling back.', e);
    destroyMandelbrotPerturbGPU();
    return null;
  }
}

function writeMandelbrotPerturbGPUParams(gpu, width, height, dx, dy, sampleCount, steps, refLen) {
  var u = new Uint32Array(gpu.paramsBytes);
  var f = new Float32Array(gpu.paramsBytes);
  u[0] = width >>> 0;
  u[1] = height >>> 0;
  u[2] = (parseInt(steps, 10) || 0) >>> 0;
  u[3] = sampleCount >>> 0;
  u[4] = doOrbit ? 1 : 0;
  u[5] = orbitType >>> 0;
  u[6] = doOrbitAverage ? 1 : 0;
  u[7] = refLen >>> 0;
  // The dcrStart/dciStart and per-pixel step values are perturbation
  // offsets from the reference centre.  Deriving them from xRange/dx
  // (which round-trip through lookAt +/- zoom/2) suffers catastrophic
  // cancellation once zoom drops below the f64 LSB of lookAt: the
  // subtraction collapses to zero or one-ULP noise.  Compute them
  // straight from zoom instead, where every operation stays in the
  // representable range of f64.
  var aspect = width / height;
  var zoomX = aspect >= 1 ? zoom * aspect : zoom;
  var zoomY = aspect >= 1 ? zoom : zoom / aspect;
  f[8] = -zoomX / 2;
  f[9] = -zoomY / 2;
  f[10] = zoomX / (0.5 + (width - 1));
  f[11] = zoomY / (0.5 + (height - 1));
  f[12] = escapeRadius;
  f[13] = orbit_trap;
  f[14] = orbitRealPoint;
  f[15] = orbitImgPoint;
  f[16] = mandelGlitchTolerance;
  f[17] = 0.0;
  f[18] = 0.0;
  f[19] = 0.0;
  gpu.device.queue.writeBuffer(gpu.paramsBuf, 0, gpu.paramsBytes);
}

async function calcMandelbrotPerturbPixelsGPU(width, height, dx, dy, samples, steps) {
  var resultCount = width * height * samples;
  var fix = ensureLookAtFix();
  var ref = getReferenceOrbit(fix[0], fix[1], steps ?? 20, escapeRadius);
  var orbitCount = ref.length;
  if (orbitCount < 2) return null;

  var gpu = await initMandelbrotPerturbGPU(resultCount, orbitCount);
  if (!gpu) return null;

  gpu.device.queue.writeBuffer(gpu.orbitBuf, 0, ref.orbit.buffer, ref.orbit.byteOffset, orbitCount * 8);
  writeMandelbrotPerturbGPUParams(gpu, width, height, dx, dy, samples, steps, orbitCount);

  var outputBytes = resultCount * 4 * 4;
  var d = gpu.device;
  var encoder = d.createCommandEncoder();
  var pass = encoder.beginComputePass();
  pass.setPipeline(gpu.pipeline);
  pass.setBindGroup(0, gpu.bindGroup);
  pass.dispatchWorkgroups(
    Math.ceil(width / MANDEL_GPU_WORKGROUP_X),
    Math.ceil((height * samples) / MANDEL_GPU_WORKGROUP_Y)
  );
  pass.end();
  encoder.copyBufferToBuffer(gpu.outputBuf, 0, gpu.readBuf, 0, outputBytes);
  d.queue.submit([encoder.finish()]);
  await d.queue.onSubmittedWorkDone();
  await gpu.readBuf.mapAsync(GPUMapMode.READ, 0, outputBytes);

  var range = gpu.readBuf.getMappedRange(0, outputBytes);
  return {
    pixels: new Float32Array(range),
    samples: samples,
    refEscaped: ref.escaped,
    refLength: ref.length,
    release: function () {
      try { gpu.readBuf.unmap(); } catch (e) { }
    },
  };
}




// Initialize canvas
var canvas = $('canvasMandelbrot');
canvas.width = window.innerWidth;
canvas.height = window.innerHeight;
var ccanvas = $('canvasControls');
ccanvas.width = window.innerWidth;
ccanvas.height = window.innerHeight;
var ctx = canvas.getContext('2d');
var img = ctx.createImageData(canvas.width, 1);

// Fetch given element, jQuery-style
function $(id) {
  return document.getElementById(id);
}

function focusOnSubmit() {
  var e = $('submitButton');
  if (e) e.focus();
}

function getSamples() {
  // more than 8x anti-alias samples can overload the GPU
  var i = parseInt($('superSamples').value, 10);
  if (i <= 0) i = 1;
  return Math.min(8, i);
}

/*
 * Main renderer equation.
 *
 * Returns number of iterations and values of Z_{n}^2 = Tr + Ti at the time
 * we either converged (n == iterations) or diverged.  We use these to
 * determined the color at the current pixel.
 *
 * The Mandelbrot set is rendered taking
 *
 *     Z_{n+1} = Z_{n} + C
 *
 * with C = x + iy, based on the "look at" coordinates.
 *
 * The Julia set can be rendered by taking
 *
 *     Z_{0} = C = x + iy
 *     Z_{n+1} = Z_{n} + K
 *
 * for some arbitrary constant K.  The point C for Z_{0} must be the
 * current pixel we're rendering, but K could be based on the "look at"
 * coordinate, or by letting the user select a point on the screen.
 */
function iterateEquation(Cr, Ci, iterations) {
  var Zr = 0;
  var Zi = 0;
  var Tr = 0;
  var Ti = 0;
  var n = 0;
  var distance = 0;

  for (; n < iterations && (Tr + Ti) <= escapeRadius; ++n) {
    Zi = 2 * Zr * Zi + Ci;
    Zr = Tr - Ti + Cr;
    Tr = Zr * Zr;
    Ti = Zi * Zi;

    if (doOrbit) {
      var dist;
      if (orbitType == 0) {
        dist = Math.abs((Zi + orbitImgPoint) * (Zr + orbitRealPoint));
      }
      else if (orbitType == 1) {
        //pointy rect
        dist = Math.abs((Zi - orbitImgPoint));
        var dist2 = Math.abs((Zr - orbitRealPoint) * (Zi - orbitImgPoint));
        if (dist2 > dist)
          dist = dist2;
      }
      else if (orbitType == 2) {
        //hyperbola
        dist = Math.abs((Zi - orbitImgPoint) * (Zr - orbitRealPoint));
      }
      else if (orbitType == 3) {
        //circle
        dist = calcDistance(orbitRealPoint, orbitImgPoint, Zr, Zi);
      }
      else {
        //spiral
        dist = Math.atan(Math.abs((Zi - orbitImgPoint) / (Zr - orbitRealPoint)));
      }

      if (dist < orbit_trap) {
        if (doOrbitAverage) {
          if (distance == 0) distance = dist;
          else distance = dist * 0.25 + distance * 0.75;
        }
        else {
          if (distance == 0 || dist < distance)
            distance = dist;
        }
      }
    }
  }
  return [n, Tr, Ti, distance];
}

function calcDistance(x1, y1, x2, y2) {
  return Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

function getColorPicker() {
  var p = $("colorScheme").value;
  if (p == "pickGradientColorBands") return pickGradientColorBands;
  return pickSharpColorBands;
}


/*
 * Render the Mandelbrot set
 */
function draw(pickColor, superSamples) {
  if (lookAt === null) lookAt = [-0.6, 0];
  if (zoom === null) zoom = zoomStart;

  smoothColor = $("smooth").checked;

  initialColor = $("colorSlider").value / 100.0;

  if (smoothColor)
    contrast = $("contrastSlider").value / 200.0;
  else
    contrast = 1.0 - $("contrastSlider").value / 200.0;

  if (reInitCanvas) {
    reInitCanvas = false;

    canvas = $('canvasMandelbrot');
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;

    ccanvas = $('canvasControls');
    ccanvas.width = window.innerWidth;
    ccanvas.height = window.innerHeight;

    ctx = canvas.getContext('2d');
    img = ctx.createImageData(canvas.width, 1);
  }

  var aspect = canvas.width / canvas.height;
  var zoomX = aspect >= 1 ? zoom * aspect : zoom;
  var zoomY = aspect >= 1 ? zoom : zoom / aspect;
  xRange = [lookAt[0] - zoomX / 2, lookAt[0] + zoomX / 2];
  yRange = [lookAt[1] - zoomY / 2, lookAt[1] + zoomY / 2];

  var steps = Math.min(100000, parseInt($('steps').value));
  if ($('autoIterations').checked) {
    var f = Math.sqrt(
      0.001 + 4.0 * Math.min(
        Math.abs(xRange[0] - xRange[1]),
        Math.abs(yRange[0] - yRange[1])));

    steps = Math.floor(250.0 / f);
    $('steps').value = String(steps);
  }

  // `spectrum` controls how many colour cycles span the iteration range
  // (hue is roughly n/steps before this multiplier).  Scale by log2(steps)
  // so banding density stays comparable whether steps is ~50 or ~50000
  const sliderSpectrum = parseFloat($("spectrumSlider").value) || 1;
  const factor = $("colorScheme").value == "pickGradientColorBands" && spectrum > 11 ? 1 : 100;
  spectrum = (sliderSpectrum/factor) * Math.log2(Math.min(steps,10000)) / Math.log2(250);
  console.log("initialColor: " + initialColor + " spectrum: " + spectrum);

  var dx = (xRange[1] - xRange[0]) / (0.5 + (canvas.width - 1));
  var dy = (yRange[1] - yRange[0]) / (0.5 + (canvas.height - 1));
  var Ci_step = (yRange[1] - yRange[0]) / (0.5 + (canvas.height - 1));

  updateHashTag(superSamples, steps);
  updateInfoBox();

  // Only enable one render at a time
  renderId += 1;

  function drawLineSuperSampled(Ci, off, Cr_init, Cr_step) {
    var Cr = Cr_init;

    for (var x = 0; x < canvas.width; ++x, Cr += Cr_step) {
      var color = [0, 0, 0, 255];

      for (var s = 0; s < superSamples; ++s) {
        var rx = Math.random() * Cr_step;
        var ry = Math.random() * Ci_step;
        var p = iterateEquation(Cr - rx / 2, Ci - ry / 2, steps);
        color = addRGB(color, pickColor(steps, p[0], p[1], p[2], p[3]));
      }

      color = divRGB(color, superSamples);

      img.data[off++] = color[0];
      img.data[off++] = color[1];
      img.data[off++] = color[2];
      img.data[off++] = 255;
    }
  }

  function drawLine(Ci, off, Cr_init, Cr_step) {
    var Cr = Cr_init;

    for (var x = 0; x < canvas.width; ++x, Cr += Cr_step) {
      var p = iterateEquation(Cr, Ci, steps);
      var color = pickColor(steps, p[0], p[1], p[2], p[3]);
      img.data[off++] = color[0];
      img.data[off++] = color[1];
      img.data[off++] = color[2];
      img.data[off++] = 255;
    }
  }

  function drawSolidLine(y, color) {
    var off = y * canvas.width;

    for (var x = 0; x < canvas.width; ++x) {
      img.data[off++] = color[0];
      img.data[off++] = color[1];
      img.data[off++] = color[2];
      img.data[off++] = color[3];
    }
  }

  function render() {
    var start = (new Date).getTime();
    var startHeight = canvas.height;
    var startWidth = canvas.width;
    var lastUpdate = start;
    var updateTimeout = 200;
    var pixels = 0;
    var Ci = yRange[0];
    var sy = 0;
    var drawLineFunc = superSamples > 1 ? drawLineSuperSampled : drawLine;
    var ourRenderId = renderId;

    var scanline = function () {
      if (renderId != ourRenderId ||
        startHeight != canvas.height ||
        startWidth != canvas.width) {
        // Stop drawing
        return;
      }

      drawLineFunc(Ci, 0, xRange[0], dx);
      Ci += Ci_step;
      pixels += canvas.width;
      ctx.putImageData(img, 0, sy);

      var now = (new Date).getTime();

      /*
       * Javascript is inherently single-threaded, and the way
       * you yield thread control back to the browser is MYSTERIOUS.
       *
       * People seem to use setTimeout() to yield, which lets us
       * make sure the canvas is updated, so that we can do animations.
       *
       * But if we do that for every scanline, it will take 100x longer
       * to render everything, because of overhead.  So therefore, we'll
       * do something in between.
       */
      if (sy++ < canvas.height) {
        if ((now - lastUpdate) >= updateTimeout) {
          // show the user where we're rendering
          drawSolidLine(0, [255, 59, 3, 255]);
          ctx.putImageData(img, 0, sy);

          // Update speed and time taken
          var elapsedMS = now - start;
          $('renderTime').innerHTML = (elapsedMS / 1000.0).toFixed(1); // 1 comma

          var speed = Math.floor(pixels / elapsedMS);

          if (metric_units(speed).substring(0, 3) == "NaN") {
            speed = Math.floor(60.0 * pixels / elapsedMS);
            $('renderSpeedUnit').innerHTML = 'minute';
          } else
            $('renderSpeedUnit').innerHTML = 'second';

          $('renderSpeed').innerHTML = metric_units(speed);

          // yield control back to browser, so that canvas is updated
          lastUpdate = now;
          setTimeout(scanline, 0);
        } else
          scanline();
      }
    };

    // Disallow redrawing while rendering
    scanline();
  }

  async function renderGPU() {
    if (mandelGPUInFlight) return false;
    mandelGPUInFlight = true;

    var start = (new Date).getTime();
    var startHeight = canvas.height;
    var startWidth = canvas.width;
    var ourRenderId = renderId;

    function renderIsCurrent() {
      return renderId == ourRenderId && startHeight == canvas.height && startWidth == canvas.width;
    }

    function paintGPUBatch(batch) {
      if (!renderIsCurrent()) return true;

      var gpuImg = ctx.createImageData(canvas.width, canvas.height);
      var data = gpuImg.data;
      var pixels = batch.pixels;
      var sampleCount = batch.samples;

      for (var idx = 0; idx < canvas.width * canvas.height; idx++) {
        var colorSum = [0, 0, 0, 255];

        for (var s = 0; s < sampleCount; s++) {
          var poff = (idx * sampleCount + s) * 4;
          var n = Math.round(pixels[poff]);
          var Tr = pixels[poff + 1];
          var Ti = pixels[poff + 2];
          var dist = pixels[poff + 3];

          if (!isFinite(Tr) || !isFinite(Ti) || !isFinite(dist)) {
            // GPU overflow can produce NaN/Infinity; treat those samples as interior.
            colorSum = addRGB(colorSum, interiorColor);
            continue;
          }

          colorSum = addRGB(colorSum, pickColor(steps, n, Tr, Ti, dist));
        }

        var color = divRGB(colorSum, sampleCount);
        var outOff = idx * 4;
        data[outOff] = color[0];
        data[outOff + 1] = color[1];
        data[outOff + 2] = color[2];
        data[outOff + 3] = 255;
      }

      if (!renderIsCurrent()) return true;
      ctx.putImageData(gpuImg, 0, 0);
      return true;
    }

    try {
      var batch = null;
      var usePerturb = zoom < MANDEL_PERTURB_ZOOM_THRESHOLD;
      try {
        if (usePerturb) {
          batch = await calcMandelbrotPerturbPixelsGPU(canvas.width, canvas.height, dx, dy, superSamples, steps);
        }
        if (!batch) {
          batch = await calcMandelbrotPixelsGPU(canvas.width, canvas.height, dx, dy, superSamples, steps);
          usePerturb = false;
        }
      } catch (e) {
        console.warn('Mandelbrot WebGPU render failed; falling back to CPU.', e);
        mandelGPUDisabled = true;
        mandelGPUStatus = 'CPU fallback';
        destroyMandelbrotGPU();
        destroyMandelbrotPerturbGPU();
        return false;
      }

      if (!batch) return false;
      try {
        if (!paintGPUBatch(batch)) return false;
      } finally {
        batch.release();
      }

      if (!renderIsCurrent()) return true;

      var elapsedMS = Math.max(1, (new Date).getTime() - start);
      mandelGPUStatus = usePerturb ? 'GPU perturb' : 'GPU';
      $('renderTime').innerHTML = (elapsedMS / 1000.0).toFixed(1);
      $('renderSpeed').innerHTML = metric_units(Math.floor((canvas.width * canvas.height) / elapsedMS));
      $('renderSpeedUnit').innerHTML = 'second (' + mandelGPUStatus + ')';
      return true;
    } finally {
      mandelGPUInFlight = false;
    }
  }

  if (!mandelGPUDisabled && isMandelbrotWebGPUSupported()) {
    var fallbackRenderId = renderId;
    renderGPU().then(function (usedGPU) {
      if (!usedGPU && renderId == fallbackRenderId) {
        render();
      }
    });
  } else {
    render();
  }
}

// Some constants used with smoothColor
var logBase = 1.0 / Math.log(2.0);
var smoothColor = true;
var band = 1.0;
var contrast = 0.45;
var glow = 1.7;
var spectrum = 0.3;
var initialColor = 0.95;


function HSVtoRGB(h, s, v) {
  var r, g, b, i, f, p, q, t;
  if (h && s === undefined && v === undefined) {
    s = h.s, v = h.v, h = h.h;
  }
  i = Math.floor(h * 6);
  f = h * 6 - i;
  p = v * (1 - s);
  q = v * (1 - f * s);
  t = v * (1 - (1 - f) * s);
  switch (i % 6) {
    case 0: r = v, g = t, b = p; break;
    case 1: r = q, g = v, b = p; break;
    case 2: r = p, g = v, b = t; break;
    case 3: r = p, g = q, b = v; break;
    case 4: r = t, g = p, b = v; break;
    case 5: r = v, g = p, b = q; break;
  }

  var rgb = [0, 0, 0];
  rgb[0] = Math.floor(r * 255);
  rgb[1] = Math.floor(g * 255);
  rgb[2] = Math.floor(b * 255);

  return rgb;
}

function getHue(steps, n, Tr, Ti, distance) {
  var hue = (n - logBase * Math.log(Math.log(Tr + Ti))) / steps;

  //orbit trap
  if (doOrbit) {
    // The circle trap is a Euclidean distance: its values spread fairly
    // uniformly across [0, orbit_trap], so `dx` regularly dips to its
    // 0.03 floor and amplifies hue ~33x.  The product-style traps' values
    // cluster near 0, so they rarely hit that floor.  Squaring (and
    // renormalising back into [0, orbit_trap]) skews the circle trap's
    // distribution toward 0 to match the other traps, without altering
    // the trap region or the no-hit (distance==0) case.
    var d = distance;
    if (orbitType == 3 && orbit_trap > 0) d = (d * d) / orbit_trap;
    var dx = (orbit_trap - d + 0.03);
    hue = hue / dx;
  }

  return hue;
}


function pickGradientColorBands(steps, n, Tr, Ti, distance) {
  if (n == steps) return interiorColor;

  var hue = getHue(steps, n, Tr, Ti, distance);
  var huesat = hue;
  hue = hue * 16000 * spectrum;
  hue = (hue % 16000) / 16000;

  huesat = huesat * 16000 * band * colourSeed;
  huesat = (huesat % 16000) / 16000;
  var new_sat = 1.0 - huesat;

  if (smoothColor) {
    if (huesat > 0.5) huesat = 1.0 - huesat;
    huesat *= 2.0;
  }

  let gama_hue = huesat;
  if (steps > 10000) {
    gama_hue = gama_hue * 16000 * colourSeed * 2;
    gama_hue = (gama_hue % 16000) / 16000;
  }

  var gamma = Math.min(Math.pow(gama_hue, contrast), 1.0);
  if (smoothColor) {
    if (gamma > 0.5) gamma = 1.0 - gamma;
    gamma *= 2.0;
  }
  gamma = Math.min(1.0, gamma+0.1);

  var c = HSVtoRGB(hue + initialColor, new_sat, gamma);
  c.push(255); // alpha
  return c;
}

function pickSharpColorBands(steps, n, Tr, Ti, distance) {
  if (n == steps) return interiorColor;

  let hue = getHue(steps, n, Tr, Ti, distance);
  hue = hue * 16000 * spectrum;
  hue = (hue % 16000) / 16000;

  let hue2 = hue * 16000 * spectrum * colourSeed;
  hue2 = (hue2 % 16000) / 16000;


  let gamma = Math.min(glow * Math.pow(hue2, contrast), 0.99);
  if (smoothColor) {
    if (gamma > 0.5) gamma = 1 - gamma;
    gamma *= 2;
    gamma = Math.min(1.0, gamma+0.1);
  }
  
  new_sat = 1.0 - gamma;
  new_sat = Math.pow(new_sat / band, 0.75);
  if (new_sat > 1.0) new_sat = 1.0;

  if (smoothColor) {
    if (new_sat > 0.5) new_sat = 1.0 - new_sat;
    new_sat *= 2.0;
  }

  const c = HSVtoRGB(hue + initialColor, new_sat, gamma);
  c.push(255); // alpha
  return c;
}

/*
 * Update URL's hash with render parameters so we can pass it around.
 */
var lastWrittenHash = null;

// Map a canvas pixel to the complex point under it, returned as a fix-point
// BigInt pair.  The pan delta `zoomX * nx` is computed in f64 (which is fine:
// |nx|<=0.5 and zoomX is small, so the product is well within f64 precision),
// converted exactly to fix-point, then added to lookAtFix.  This bypasses
// both the xRange[1]-xRange[0] cancellation and the f64 LSB on `lookAt`
// itself, so drag-to-zoom keeps tracking the cursor past zoom ~ 1e-16.
function pixelToLookAtFix(px, py) {
  var aspect = canvas.width / canvas.height;
  var zoomX = aspect >= 1 ? zoom * aspect : zoom;
  var zoomY = aspect >= 1 ? zoom : zoom / aspect;
  var nx = px / (canvas.width - 0.5) - 0.5;
  var ny = py / (canvas.height - 0.5) - 0.5;
  var fix = ensureLookAtFix();
  return [fix[0] + mandelNumberToFix(zoomX * nx),
          fix[1] + mandelNumberToFix(zoomY * ny)];
}

function updateHashTag(samples, iterations) {
  var scheme = $('colorScheme').value;

  // Number.prototype.toString already produces the shortest decimal that
  // round-trips through parseFloat at full f64 precision, so zoom and the
  // f64 `lookAt` restore exactly at any zoom f64 can express.  `lookAtHi`
  // carries the BigInt fixed-point centre hex-encoded so sub-LSB pan offsets
  // accumulated by deep-zoom drags survive a refresh too.
  var fix = ensureLookAtFix();
  var h = 'zoom=' + zoom + '&' +
    'lookAt=' + lookAt + '&' +
    'lookAtHi=' + mandelFixToHex(fix[0]) + ',' + mandelFixToHex(fix[1]) + '&' +
    'iterations=' + iterations + '&' +
    'superSamples=' + samples + '&' +
    'escapeRadius=' + $('escapeRadius').value + '&' +
    'color=' + $('colorSlider').value + '&' +
    'spectrum=' + $('spectrumSlider').value + '&' +
    'contrast=' + $('contrastSlider').value + '&' +
    'smooth=' + ($('smooth').checked ? 1 : 0) + '&' +
    'colorScheme=' + scheme + '&' +
    'orbitType=' + $('orbitType').value;
  lastWrittenHash = h;
  location.hash = h;
}

/*
 * Update small info box in lower right hand side
 */
function updateInfoBox() {
  // Update infobox
  $('infoBox').innerHTML =
    'x<sub>0</sub>=' + xRange[0] + ' y<sub>0</sub>=' + yRange[0] + ' ' +
    'x<sub>1</sub>=' + xRange[1] + ' y<sub>1</sub>=' + yRange[1] + ' ' +
    'wxh=' + canvas.width + 'x' + canvas.height + ' '
    + (canvas.width * canvas.height / 1000000.0).toFixed(1) + 'MP';
}

/*
 * Parse URL hash tag, returns whether we should redraw.
 */
function readHashTag() {
  var redraw = false;
  var hash = location.hash.replace(/^#/, '');
  var tags = hash.split('&');
  // `lookAtHi` (exact BigInt) wins over the f64 `lookAt` regardless of which
  // tag appears first in the URL.
  var lookAtHiSeen = false;

  for (var i = 0; i < tags.length; ++i) {
    var tag = tags[i].split('=');
    var key = tag[0];
    var val = tag[1];

    switch (key) {
      case 'zoom': {
        zoom = parseFloat(val);
        redraw = true;
      } break;

      case 'lookAt': {
        if (!lookAtHiSeen) {
          var l = val.split(',');
          lookAt = [parseFloat(l[0]), parseFloat(l[1])];
          lookAtFix = null; // re-derive lazily from the new f64 mirror
        }
        redraw = true;
      } break;

      case 'lookAtHi': {
        var l = val.split(',');
        setLookAtFix(mandelHexToFix(l[0]), mandelHexToFix(l[1]));
        lookAtHiSeen = true;
        redraw = true;
      } break;

      case 'iterations': {
        $('steps').value = val;
        $('autoIterations').checked = false;
        redraw = true;
      } break;

      case 'superSamples': {
        $('superSamples').value = String(parseInt(val, 10));
        redraw = true;
      } break;

      case 'colorScheme': {
        $('colorScheme').value = String(val);
        redraw = true;
      } break;

      case 'escapeRadius': {
        $('escapeRadius').value = val;
        redraw = true;
      } break;

      case 'color': {
        $('colorSlider').value = String(parseInt(val, 10));
        redraw = true;
      } break;

      case 'spectrum': {
        $('spectrumSlider').value = String(parseInt(val, 10));
        redraw = true;
      } break;

      case 'contrast': {
        $('contrastSlider').value = String(parseInt(val, 10));
        redraw = true;
      } break;

      case 'smooth': {
        $('smooth').checked = (val === '1' || val === 'true');
        redraw = true;
      } break;

      case 'orbitType': {
        var ot = parseInt(val, 10);
        if (isFinite(ot)) {
          $('orbitType').value = String(ot);
          doOrbit = (ot >= 0);
          if (doOrbit) orbitType = ot;
          redraw = true;
        }
      } break;
    }
  }

  if (redraw)
    reInitCanvas = true;

  return redraw;
}

/*
 * Return number with metric units
 */
function metric_units(number) {
  var unit = ["", "k", "M", "G", "T", "P", "E"];
  var mag = Math.ceil((1 + Math.log(number) / Math.log(10)) / 3);
  return "" + (number / Math.pow(10, 3 * (mag - 1))).toFixed(2) + unit[mag];
}


function addRGB(v, w) {
  v[0] += w[0];
  v[1] += w[1];
  v[2] += w[2];
  v[3] += w[3];
  return v;
}

function divRGB(v, div) {
  v[0] /= div;
  v[1] /= div;
  v[2] /= div;
  v[3] /= div;
  return v;
}

function main() {
  $('viewPNG').onclick = function (event) {
    var link = document.createElement('a');
    link.download = `mandelbrot-${lookAt[0]}-${lookAt[1]}.png`;
    link.href = canvas.toDataURL('image/png');
    link.click();
  };

  $('steps').onblur = function (event) {
    // disable auto-iterations when user edits it manually
    $('autoIterations').checked = false;
    draw(getColorPicker(), getSamples());
  };

  $("contrastSlider").onchange = function () {
    draw(getColorPicker(), getSamples());
  };

  $("colorSlider").onchange = function () {
    draw(getColorPicker(), getSamples());
  };

  $("spectrumSlider").onchange = function () {
    draw(getColorPicker(), getSamples());
  };

  $("colorScheme").onchange = function () {
    draw(getColorPicker(), getSamples());
  };

  $("orbitType").onchange = function () {
    var v = parseInt($("orbitType").value, 10);
    doOrbit = (v >= 0);
    if (doOrbit) orbitType = v;
    draw(getColorPicker(), getSamples());
  };

  $("smooth").onchange = function () {
    if (steps > 100) colourSeed = 1 + Math.random() * 3.5;
    draw(getColorPicker(), getSamples());
  };

  $('resetButton').onclick = function (even) {
    $('settingsForm').reset();
    setTimeout(function () { location.hash = ''; }, 1);
    zoom = zoomStart;
    lookAt = lookAtDefault;
    lookAtFix = null;
    reInitCanvas = true;
    draw(getColorPicker(), getSamples());
  };

  if (dragToZoom == true) {
    var box = null;

    $('canvasControls').onmousedown = function (e) {
      if (e.button == 2) return;
      if (box == null)
        box = [e.clientX, e.clientY, 0, 0];
    };

    $('canvasControls').oncontextmenu = function (e) {
      e.preventDefault();
      var d = $('description');
      d.style.display = (d.style.display == 'none') ? '' : 'none';
      return false;
    };

    $('canvasControls').onmousemove = function (e) {
      if (box != null) {
        var c = ccanvas.getContext('2d');
        c.lineWidth = 1;

        // clear out old box first
        c.clearRect(0, 0, ccanvas.width, ccanvas.height);

        // draw new box
        c.strokeStyle = '#FF3B03';
        box[2] = e.clientX;
        box[3] = e.clientY;
        c.strokeRect(box[0], box[1], box[2] - box[0], box[3] - box[1]);
      }
    };

    var zoomOut = function (event) {
      var nf = pixelToLookAtFix(event.clientX, event.clientY);
      setLookAtFix(nf[0], nf[1]);

      if (event.shiftKey) {
        zoom /= 0.5;
      }

      draw(getColorPicker(), getSamples());
    };

    $('canvasControls').onmouseup = function (e) {
      if (box != null) {
        // Zoom out
        if (e.shiftKey) {
          box = null;
          zoomOut(e);
          return;
        }

        // Abort if the zoom box is too small
        if (Math.abs(box[0] - box[2]) < 30) {
          var c = ccanvas.getContext('2d');
          c.clearRect(0, 0, ccanvas.width, ccanvas.height);
          box = null;
          return;
        }

        var c = ccanvas.getContext('2d');
        c.clearRect(0, 0, ccanvas.width, ccanvas.height);

        // Calculate new rectangle to render.  Compute the box centre in
        // complex space directly from `zoom`/`lookAtFix` (see
        // pixelToLookAtFix) so deep-zoom drags don't lose precision via
        // xRange cancellation or via the f64 LSB of `lookAt` itself.
        var bx = Math.min(box[0], box[2]) + Math.abs(box[0] - box[2]) / 2.0;
        var by = Math.min(box[1], box[3]) + Math.abs(box[1] - box[3]) / 2.0;

        var nf = pixelToLookAtFix(bx, by);
        setLookAtFix(nf[0], nf[1]);

        var xf = Math.abs(Math.abs(box[0] - box[2]) / canvas.width);
        var yf = Math.abs(Math.abs(box[1] - box[3]) / canvas.height);

        zoom *= Math.max(xf, yf); // retain aspect ratio

        box = null;
        draw(getColorPicker(), getSamples());
      }
    };
  }

  /*
   * Enable zooming (currently, the zooming is inexact!) Click to zoom;
   * perfect to mobile phones, etc.
   */
  if (dragToZoom == false) {
    $('canvasMandelbrot').onclick = function (event) {
      var nf = pixelToLookAtFix(event.clientX, event.clientY);
      setLookAtFix(nf[0], nf[1]);

      if (event.shiftKey) {
        zoom /= 0.5;
      } else {
        zoom *= 0.5;
      }

      draw(getColorPicker(), getSamples());
    };
  }

  /*
   * When resizing the window, be sure to update all the canvas stuff.
   */
  window.onresize = function (event) {
    reInitCanvas = true;
  };

  /*
   * Re-render when the URL hash is edited manually.
   */
  window.onhashchange = function () {
    if (location.hash.replace(/^#/, '') === lastWrittenHash) return;
    if (readHashTag())
      draw(getColorPicker(), getSamples());
  };

  /*
   * Read hash tag and render away at page load.
   */
  readHashTag();
  draw(getColorPicker(), getSamples());
}

main();
