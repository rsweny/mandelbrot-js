/*
 * Taylor Series Fractal, in HTML5 canvas and javascript.
 * Ported from Taylor.java by Ryan Sweny
 * https://github.com/rsweny/mandelbrot-js
 *
 * Algorithm reference: http://fractorama.com/doc/taylor.html
 *
 * Rendering approach: accumulation / histogram over multiple passes.
 */

// ---------------------------------------------------------------------------
// Complex number arithmetic
// ---------------------------------------------------------------------------

function Complex(re, im) {
  this.re = (re !== undefined) ? re : 0;
  this.im = (im !== undefined) ? im : 0;
}

Complex.prototype.add = function(w) {
  return new Complex(this.re + w.re, this.im + w.im);
};
Complex.prototype.sub = function(w) {
  return new Complex(this.re - w.re, this.im - w.im);
};
Complex.prototype.mul = function(w) {
  return new Complex(
    this.re * w.re - this.im * w.im,
    this.re * w.im + this.im * w.re
  );
};
Complex.prototype.div = function(w) {
  var d = w.re * w.re + w.im * w.im;
  if (d === 0) return new Complex(NaN, NaN);
  return new Complex(
    (this.re * w.re + this.im * w.im) / d,
    (this.im * w.re - this.re * w.im) / d
  );
};
Complex.prototype.abs = function() {
  return Math.sqrt(this.re * this.re + this.im * this.im);
};
Complex.prototype.arg = function() {
  return Math.atan2(this.im, this.re);
};
Complex.prototype.log = function() {
  return new Complex(Math.log(this.abs()), this.arg());
};
Complex.prototype.sinh = function() {
  return new Complex(
    Math.sinh(this.re) * Math.cos(this.im),
    Math.cosh(this.re) * Math.sin(this.im)
  );
};
Complex.prototype.cot = function() {
  var a = this.re, b = this.im;
  var sinZ = new Complex(Math.sin(a) * Math.cosh(b),  Math.cos(a) * Math.sinh(b));
  var cosZ = new Complex(Math.cos(a) * Math.cosh(b), -Math.sin(a) * Math.sinh(b));
  return cosZ.div(sinZ);
};

// z^n where n is real (polar form)
Complex.pow = function(z, n) {
  if (n === 0) return new Complex(1, 0);
  var r = z.abs();
  if (r === 0) return new Complex(0, 0);
  var rn    = Math.pow(r, n);
  var theta = z.arg() * n;
  return new Complex(rn * Math.cos(theta), rn * Math.sin(theta));
};

// ---------------------------------------------------------------------------
// Globals
// ---------------------------------------------------------------------------

function $(id) { return document.getElementById(id); }

var canvas, ctx, ccanvas, cctx;
var ximlen, yimlen;
var img_red, img_green, img_blue, img_alpha;
var imageData;

// Render-loop state
var quad, lasti, points;
var running      = false;
var renderTimer  = null;
var lastUpdateTime = 0;
var lastPassTime   = 0;
var passCount      = 0;

// Fractal parameters (matching Taylor.java defaults)
var zoom         = 9;
var xcen         = 0.0;
var ycen         = 0.01;
var depth        = 20;
var glow         = 0.015;
var brightness   = 1.28;
var complexScale = 0;
var cx = 0, cy = 0;
var scalePoint    = false;
var compoundPoint = false;
var recurse       = false;
var mandelbrot    = false;
var alg           = 0;

// Drag-to-zoom
var dragBox = null;

// ---------------------------------------------------------------------------
// UI helpers
// ---------------------------------------------------------------------------

function readControls() {
  try { depth        = parseInt($('depth').value, 10)     || 20;   } catch(e) {}
  try { zoom         = parseFloat($('zoom').value)         || 9;    } catch(e) {}
  try { xcen         = parseFloat($('xcen').value)         || 0;    } catch(e) {}
  try { ycen         = parseFloat($('ycen').value)         || 0.01; } catch(e) {}
  glow         = $('glowSlider').value / 1000;
  brightness   = $('brightnessSlider').value / 100;
  try { complexScale = parseFloat($('complexScale').value) || 0;    } catch(e) {}
  scalePoint    = $('scalePoint').checked;
  compoundPoint = $('compoundPoint').checked;
  recurse       = $('recurse').checked;
  mandelbrot    = $('mandelbrotMode').checked;
  alg           = parseInt($('algorithm').value, 10);
}

function writeControls() {
  $('depth').value        = depth;
  $('zoom').value         = zoom;
  $('xcen').value         = xcen;
  $('ycen').value         = ycen;
  $('glowSlider').value        = Math.round(glow * 1000);
  $('brightnessSlider').value  = Math.round(brightness * 100);
  $('complexScale').value = complexScale;
  $('scalePoint').checked    = scalePoint;
  $('compoundPoint').checked = compoundPoint;
  $('recurse').checked       = recurse;
  $('mandelbrotMode').checked = mandelbrot;
  $('algorithm').value = alg;
}

function updateStatus(msg) {
  $('status').textContent = msg || '';
}

// ---------------------------------------------------------------------------
// URL hash — persist fractal parameters so views can be bookmarked / shared.
// ---------------------------------------------------------------------------

var lastWrittenHash = null;

function updateHashTag() {
  var h = 'zoom='          + zoom +
          '&xcen='         + xcen +
          '&ycen='         + ycen +
          '&depth='        + depth +
          '&glow='         + Math.round(glow * 1000) +
          '&brightness='   + Math.round(brightness * 100) +
          '&complexScale=' + complexScale +
          '&alg='          + alg +
          '&scalePoint='   + (scalePoint    ? '1' : '0') +
          '&compoundPoint='+ (compoundPoint ? '1' : '0') +
          '&recurse='      + (recurse       ? '1' : '0') +
          '&mandelbrot='   + (mandelbrot    ? '1' : '0') +
          '&cx='           + cx +
          '&cy='           + cy;
  lastWrittenHash = h;
  location.hash = h;
}

function readHashTag() {
  var tags = location.hash.replace(/^#/, '').split('&');
  var redraw = false;
  for (var i = 0; i < tags.length; i++) {
    var kv  = tags[i].split('=');
    var key = kv[0];
    var val = decodeURIComponent(kv[1] !== undefined ? kv[1] : '');
    var fv  = parseFloat(val);
    var iv  = parseInt(val, 10);
    switch (key) {
      case 'zoom':          if (!isNaN(fv)) zoom         = fv;        redraw = true; break;
      case 'xcen':          if (!isNaN(fv)) xcen         = fv;        redraw = true; break;
      case 'ycen':          if (!isNaN(fv)) ycen         = fv;        redraw = true; break;
      case 'depth':         if (!isNaN(iv)) depth        = iv;        redraw = true; break;
      case 'glow':          if (!isNaN(iv)) glow         = iv / 1000; redraw = true; break;
      case 'brightness':    if (!isNaN(iv)) brightness   = iv / 100;  redraw = true; break;
      case 'complexScale':  if (!isNaN(fv)) complexScale = fv;        redraw = true; break;
      case 'alg':           if (!isNaN(iv)) alg          = iv;        redraw = true; break;
      case 'scalePoint':    scalePoint    = val === '1';               redraw = true; break;
      case 'compoundPoint': compoundPoint = val === '1';               redraw = true; break;
      case 'recurse':       recurse       = val === '1';               redraw = true; break;
      case 'mandelbrot':    mandelbrot    = val === '1';               redraw = true; break;
      case 'cx':            if (!isNaN(fv)) cx = fv;                   redraw = true; break;
      case 'cy':            if (!isNaN(fv)) cy = fv;                   redraw = true; break;
    }
  }
  return redraw;
}

// ---------------------------------------------------------------------------
// Canvas / reset
// ---------------------------------------------------------------------------

function initCanvas() {
  canvas  = $('taylorCanvas');
  ccanvas = $('canvasControls');
  canvas.width  = ccanvas.width  = window.innerWidth;
  canvas.height = ccanvas.height = window.innerHeight;
  ctx  = canvas.getContext('2d');
  cctx = ccanvas.getContext('2d');
}

function clearAndReset() {
  readControls();

  ximlen = canvas.width;
  yimlen = canvas.height;
  points = ximlen;
  // Ensure quad is smaller than points so the loop body executes
  quad   = Math.min(10, Math.max(1, Math.floor(points / 4)));
  lasti  = 0;
  passCount = 0;

  var total = ximlen * yimlen;
  img_red   = new Float64Array(total);
  img_green = new Float64Array(total);
  img_blue  = new Float64Array(total);
  img_alpha = new Float64Array(total);

  imageData = ctx.createImageData(ximlen, yimlen);
  // Black background (alpha=255)
  for (var i = 3; i < imageData.data.length; i += 4) imageData.data[i] = 255;
  ctx.putImageData(imageData, 0, 0);

  lastUpdateTime = 0;
  lastPassTime   = Date.now();
  updateStatus('Rendering...');
}

// ---------------------------------------------------------------------------
// Core fractal computation (direct port of Java nextpoints logic)
// ---------------------------------------------------------------------------

function computePixel(xi, yi) {
  var xf = xi / ximlen;
  var yf = yi / yimlen;

  var aspect = ximlen / yimlen;
  var a = 2 * (xf - 0.5) * zoom * aspect + xcen;
  var b = 2 * (yf - 0.5) * zoom + ycen;

  var current = new Complex(a, b);
  var one = new Complex(1, 0);
  var two = new Complex(2, 0);

  // Compute the "real value" for this pixel based on the chosen algorithm
  var zvalue;
  if (alg === 0) {
    zvalue = Complex.pow(current, Math.E);
  } else if (alg === 1) {
    zvalue = one.add(current).div(one.sub(current)).log();
  } else {
    zvalue = current.sinh();
  }

  var xavg = 0, yavg = 0, mavg = 0;
  var xmax = 0, ymax = 0, mmax = 0;
  var xtot = 0, ytot = 0, mtot = 0;

  var den  = 1;
  var den2 = 1;

  // Expansion point: f(current)^e + (cx,cy)
  var point = Complex.pow(current, Math.E).add(new Complex(cx, cy));
  var zSum  = new Complex(0, 0);

  for (var count = 0; count < depth; count++) {
    var scale      = new Complex(1, complexScale);
    var oldcurrent = current;

    if (scalePoint) {
      if (alg === 0) {
        scale = Complex.pow(point, Math.E);
      } else if (alg === 1) {
        scale = one.add(point).div(one.sub(point)).log();
      } else {
        scale = point.sinh();
      }
    }

    if (compoundPoint) {
      current = current.sub(point);
    }

    // Factorial denominators
    if (count > 0) {
      den  *= count;
      den2 *= (count * 2 + 1);
    }

    var denominator  = new Complex(den);
    var denominator2 = new Complex(den2);

    if (alg === 0) {
      // f(x) = x^n / n!   (Taylor series for e^x)
      var num = scale.mul(Complex.pow(current, count));
      zSum = zSum.add(num.div(denominator));
    } else if (alg === 1) {
      // f(x) = 2x^(2n-1) / (2n-1)   (Taylor series for log((1+x)/(1-x)))
      var twoN1 = new Complex(2 * count - 1);
      if (twoN1.re !== 0) {
        var num = scale.mul(two.mul(Complex.pow(current, twoN1)));
        zSum = zSum.add(num.div(twoN1));
      }
    } else {
      // f(x) = x^(2n+1) / (2n+1)!   (Taylor series for sinh(x))
      var num = scale.mul(Complex.pow(current, den2));
      zSum = zSum.add(num.div(denominator2));
    }

    if (mandelbrot) {
      zSum = zSum.mul(zvalue.sub(zSum));
    }

    if (!recurse) {
      current = oldcurrent;
    }

    var value = zSum.cot();

    if (isFinite(value.re) && isFinite(value.im) && !isNaN(value.re)   && !isNaN(value.im)) {

      var x = Complex.pow(value, 0.1).abs();
      var y = Complex.pow(value, 0.5).abs();
      var m = Complex.pow(value, 0.9).abs();

      // m = 1 + m - Math.log(m)/Math.log(2.0);

      if (xavg === 0) { xavg = x; yavg = y; mavg = m; }

      var xdev = Math.abs(x - xavg);
      var ydev = Math.abs(y - yavg);
      var mdev = Math.abs(m - mavg);

      xavg = x * 0.5 + xavg * (0.5 + 0.01 * glow);
      yavg = y * 0.5 + yavg * (0.5 + 0.1  * glow);
      mavg = m * 0.5 + mavg * (0.5 + 1.0  * glow);

      if (xdev > xmax) xmax = xdev;
      if (ydev > ymax) ymax = ydev;
      if (mdev > mmax) mmax = mdev;

      xtot += xdev;
      ytot += ydev;
      mtot += mdev;
    }
  }

  const brightness_factor = Math.pow(depth, 2.1 - brightness);
  let red = 0, green = 0, blue = 0;
  if (brightness_factor > 0) {
    if (xmax > 0) red   = (xtot * 255 / brightness_factor) / xmax;
    if (ymax > 0) green = (ytot * 255 / brightness_factor) / ymax;
    if (mmax > 0) blue  = (mtot * 255 / brightness_factor) / mmax;
  }

  return [red, green, blue];
}

// ---------------------------------------------------------------------------
// Render batch (direct port of Java nextpoints structure)
// ---------------------------------------------------------------------------

function nextPoints() {
  var totalPixels = ximlen * yimlen;
  var i;

  for (i = lasti; i < lasti + points - quad; i += quad) {
    const xi = i % ximlen;
    const yi = Math.floor(i / ximlen);
    const rgb = computePixel(xi, yi);
    const idx = yi * ximlen + xi;

    img_red[idx]   += rgb[0];
    img_green[idx] += rgb[1];
    img_blue[idx]  += rgb[2];
    img_alpha[idx]++;

    const alpha = img_alpha[idx];
    const r = Math.min(255, img_red[idx]   / alpha) | 0;
    const g = Math.min(255, img_green[idx] / alpha) | 0;
    const b = Math.min(255, img_blue[idx]  / alpha) | 0;

    const di = idx * 4;
    imageData.data[di]     = r;
    imageData.data[di + 1] = g;
    imageData.data[di + 2] = b;
    // alpha channel already 255 from init
  }

  // Advance lasti (port of Java end-of-loop logic)
  if (i <= totalPixels - points * quad) {
    lasti += points * quad;
  } else if (quad > 1) {
    lasti = 0;
    quad  = 1;
  } else {
    lasti += points;
    if (lasti % totalPixels === 0) {
      passCount++;
      var diff = ((Date.now() - lastPassTime) / 1000).toFixed(1);
      lastPassTime = Date.now();
      if (passCount >= 1) {
        running = false;
        if (renderTimer) { clearTimeout(renderTimer); renderTimer = null; }
        ctx.putImageData(imageData, 0, 0);
        updateStatus('Done — ' + diff + 's');
      } else {
        updateStatus('Pass: ' + passCount + ' — ' + diff + 's');
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Render loop
// ---------------------------------------------------------------------------

function renderFrame() {
  var batchStart = Date.now();
  do {
    nextPoints();
  } while (running && Date.now() - batchStart < 50);

  var now = Date.now();
  if (now - lastUpdateTime > 500) {
    ctx.putImageData(imageData, 0, 0);
    lastUpdateTime = now;
  }

  if (running) renderTimer = setTimeout(renderFrame, 0);
}

function startRendering() {
  running = true;
  if (renderTimer) clearTimeout(renderTimer);
  renderTimer = setTimeout(renderFrame, 0);
}

function stopRendering() {
  running = false;
  if (renderTimer) { clearTimeout(renderTimer); renderTimer = null; }
}

// ---------------------------------------------------------------------------
// WebGPU GPU path — mirrors the structure used by newton.js.
// TAYLOR_PIXEL_WGSL is a direct WGSL port of computePixel().
// ---------------------------------------------------------------------------

var TAYLOR_PIXEL_WGSL = `
struct Params {
  width:        u32,
  height:       u32,
  depth:        u32,
  alg:          u32,
  scalePoint:   u32,
  compoundPoint:u32,
  recurse:      u32,
  mandelbrot:   u32,
  zoom:         f32,
  xcen:         f32,
  ycen:         f32,
  complexScale: f32,
  cx:           f32,
  cy:           f32,
  glow:         f32,
  brightness:   f32,
};

@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read_write> outPixels: array<vec4<f32>>;

const MATH_E: f32 = 2.718281828459045;
const ONE = vec2<f32>(1.0, 0.0);

fn cmul(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {
  return vec2<f32>(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}
fn cdiv(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {
  let d = dot(b, b);
  if (d == 0.0) { return vec2<f32>(0.0, 0.0); }
  return vec2<f32>((a.x*b.x + a.y*b.y)/d, (a.y*b.x - a.x*b.y)/d);
}
fn clog(v: vec2<f32>) -> vec2<f32> {
  let r = length(v);
  if (r == 0.0) { return vec2<f32>(-1e30, 0.0); }
  return vec2<f32>(log(r), atan2(v.y, v.x));
}
fn csinh(v: vec2<f32>) -> vec2<f32> {
  return vec2<f32>(sinh(v.x)*cos(v.y), cosh(v.x)*sin(v.y));
}
fn ccot(v: vec2<f32>) -> vec2<f32> {

  // Guard against f32 sinh/cosh overflow (happens at |y| > ~88, much sooner
  // than the f64 CPU path at ~710).  For large |Im(v)|, cot(x+iy) → -i·sign(y),
  //  if (abs(v.y) > 80) {
  //    return vec2<f32>(0.0, select(-1.0, 1.0, v.y > 0.0));
  //  }

  let sinZ = vec2<f32>(sin(v.x)*cosh(v.y),  cos(v.x)*sinh(v.y));
  let cosZ = vec2<f32>(cos(v.x)*cosh(v.y), -sin(v.x)*sinh(v.y));
  return cdiv(cosZ, sinZ);
}
fn cpow_real(z: vec2<f32>, n: f32) -> vec2<f32> {
  if (n == 0.0) { return ONE; }
  let r = length(z);
  if (r == 0.0) { return vec2<f32>(0.0, 0.0); }
  let rn = pow(r, n);
  let theta = atan2(z.y, z.x) * n;
  return vec2<f32>(rn*cos(theta), rn*sin(theta));
}
fn isFiniteVec(v: vec2<f32>) -> bool {
  return all(abs(v) < vec2<f32>(1e20, 1e20));
}
fn zfunc(z: vec2<f32>) -> vec2<f32> {
  if (P.alg == 0u) { return cpow_real(z, MATH_E); }
  if (P.alg == 1u) { return clog(cdiv(ONE + z, ONE - z)); }
  return csinh(z);
}

@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  if (gid.x >= P.width || gid.y >= P.height) { return; }
  let xi = gid.x; let yi = gid.y;
  let idx = yi * P.width + xi;
  let xf = f32(xi) / f32(P.width);
  let yf = f32(yi) / f32(P.height);
  let aspect = f32(P.width) / f32(P.height);
  let a = 2.0*(xf - 0.5)*P.zoom*aspect + P.xcen;
  let b = 2.0*(yf - 0.5)*P.zoom + P.ycen;
  let current_init = vec2<f32>(a, b);
  let two  = vec2<f32>(2.0, 0.0);
  let zvalue  = zfunc(current_init);
  let point   = cpow_real(current_init, MATH_E) + vec2<f32>(P.cx, P.cy);
  var zSum    = vec2<f32>(0.0, 0.0);
  var current = current_init;
  var xavg = 0.0; var yavg = 0.0; var mavg = 0.0;
  var xmax = 0.0; var ymax = 0.0; var mmax = 0.0;
  var xtot = 0.0; var ytot = 0.0; var mtot = 0.0;
  var den  = 1.0;
  var den2 = 1.0;
  for (var count = 0u; count < P.depth; count++) {
    var scale    = vec2<f32>(1.0, P.complexScale);
    let oldcurrent = current;
    if (P.scalePoint    != 0u) { scale = zfunc(point); }
    if (P.compoundPoint != 0u) { current = current - point; }
    if (count > 0u) {
      den  *= f32(count);
      den2 *= f32(count * 2u + 1u);
    }
    if (P.alg == 0u) {
      zSum += cdiv(cmul(scale, cpow_real(current, f32(count))), vec2<f32>(den, 0.0));
    } else if (P.alg == 1u) {
      let twoN1 = f32(count) * 2.0 - 1.0;
      if (twoN1 != 0.0) {
        zSum += cdiv(cmul(scale, cmul(two, cpow_real(current, twoN1))), vec2<f32>(twoN1, 0.0));
      }
    } else {
      zSum += cdiv(cmul(scale, cpow_real(current, den2)), vec2<f32>(den2, 0.0));
    }
    if (P.mandelbrot != 0u) { zSum = cmul(zSum, zvalue - zSum); }
    if (P.recurse    == 0u) { current = oldcurrent; }
    let value = ccot(zSum);
    if (isFiniteVec(value)) {
      let x = length(cpow_real(value, 0.1));
      let y = length(cpow_real(value, 0.5));
      let m = length(cpow_real(value, 0.9));
      if (xavg == 0.0) { xavg = x; yavg = y; mavg = m; }
      let xdev = abs(x - xavg); let ydev = abs(y - yavg); let mdev = abs(m - mavg);
      xavg = x * 0.5 + xavg * (0.5 + 0.01 * P.glow);
      yavg = y * 0.5 + yavg * (0.5 + 0.1  * P.glow);
      mavg = m * 0.5 + mavg * (0.5 + 1.0  * P.glow);
      if (xdev > xmax) { xmax = xdev; }
      if (ydev > ymax) { ymax = ydev; }
      if (mdev > mmax) { mmax = mdev; }
      xtot += xdev; ytot += ydev; mtot += mdev;
    }
  }
  let bf = pow(f32(P.depth), 2.01 - P.brightness);
  var red = 0.0; var green = 0.0; var blue = 0.0;
  if (bf > 0.0) {
    if (xmax > 0.0) { red   = (xtot * 255.0 / bf) / xmax; }
    if (ymax > 0.0) { green = (ytot * 255.0 / bf) / ymax; }
    if (mmax > 0.0) { blue  = (mtot * 255.0 / bf) / mmax; }
  }
  outPixels[idx] = vec4<f32>(red, green, blue, 0.0);
}
`;

var TAYLOR_GPU_PARAMS_BYTES = 64;
var TAYLOR_GPU_WORKGROUP_X  = 16;
var TAYLOR_GPU_WORKGROUP_Y  = 16;

var taylorGPU         = null;
let taylorGPUDisabled = false;

// ---------------------------------------------------------------------------
// Public draw / reset
// ---------------------------------------------------------------------------

function isTaylorWebGPUSupported() {
  return typeof navigator !== 'undefined' && !!navigator.gpu;
}

function destroyTaylorGPU() {
  if (!taylorGPU) return;
  var bufs = ['paramsBuf', 'outputBuf', 'readBuf'];
  for (var i = 0; i < bufs.length; i++) {
    try { taylorGPU[bufs[i]].destroy(); } catch(e) {}
  }
  taylorGPU = null;
}

async function initTaylorGPU(resultCount) {
  if (taylorGPUDisabled || !isTaylorWebGPUSupported()) return null;
  if (taylorGPU && taylorGPU.resultCapacity >= resultCount) return taylorGPU;
  destroyTaylorGPU();
  try {
    var adapter = await navigator.gpu.requestAdapter();
    if (!adapter) throw new Error('No WebGPU adapter');
    var device = await adapter.requestDevice();
    var outputBytes = resultCount * 4 * 4;
    var paramsBuf = device.createBuffer({ size: TAYLOR_GPU_PARAMS_BYTES, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    var outputBuf = device.createBuffer({ size: outputBytes, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC });
    var readBuf   = device.createBuffer({ size: outputBytes, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
    var module    = device.createShaderModule({ code: TAYLOR_PIXEL_WGSL });
    var pipeline  = await device.createComputePipelineAsync({ layout: 'auto', compute: { module: module, entryPoint: 'main' } });
    var bindGroup = device.createBindGroup({
      layout: pipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: paramsBuf } },
        { binding: 1, resource: { buffer: outputBuf } },
      ],
    });
    taylorGPU = {
      device:         device,
      resultCapacity: resultCount,
      paramsBuf:      paramsBuf,
      outputBuf:      outputBuf,
      readBuf:        readBuf,
      pipeline:       pipeline,
      bindGroup:      bindGroup,
      paramsBytes:    new ArrayBuffer(TAYLOR_GPU_PARAMS_BYTES),
    };
    return taylorGPU;
  } catch(e) {
    console.warn('Taylor WebGPU path disabled; falling back to CPU.', e);
    taylorGPUDisabled = true;
    destroyTaylorGPU();
    return null;
  }
}

function writeTaylorGPUParams(gpu, width, height) {
  var u = new Uint32Array(gpu.paramsBytes);
  var f = new Float32Array(gpu.paramsBytes);
  u[0] = width  >>> 0;
  u[1] = height >>> 0;
  u[2] = (depth | 0) >>> 0;
  u[3] = (alg   | 0) >>> 0;
  u[4] = scalePoint    ? 1 : 0;
  u[5] = compoundPoint ? 1 : 0;
  u[6] = recurse       ? 1 : 0;
  u[7] = mandelbrot    ? 1 : 0;
  f[8]  = zoom;
  f[9]  = xcen;
  f[10] = ycen;
  f[11] = complexScale;
  f[12] = cx;
  f[13] = cy;
  f[14] = glow;
  f[15] = brightness;
  gpu.device.queue.writeBuffer(gpu.paramsBuf, 0, gpu.paramsBytes);
}

async function calcTaylorPixelsGPU(gpu, width, height) {
  writeTaylorGPUParams(gpu, width, height);
  var outputBytes = width * height * 4 * 4;
  var d = gpu.device;
  var encoder = d.createCommandEncoder();
  var pass = encoder.beginComputePass();
  pass.setPipeline(gpu.pipeline);
  pass.setBindGroup(0, gpu.bindGroup);
  pass.dispatchWorkgroups(
    Math.ceil(width  / TAYLOR_GPU_WORKGROUP_X),
    Math.ceil(height / TAYLOR_GPU_WORKGROUP_Y)
  );
  pass.end();
  encoder.copyBufferToBuffer(gpu.outputBuf, 0, gpu.readBuf, 0, outputBytes);
  d.queue.submit([encoder.finish()]);
  await d.queue.onSubmittedWorkDone();
  await gpu.readBuf.mapAsync(GPUMapMode.READ, 0, outputBytes);
  var range = gpu.readBuf.getMappedRange(0, outputBytes);
  return {
    pixels:  new Float32Array(range),
    release: function() { try { gpu.readBuf.unmap(); } catch(e) {} },
  };
}

async function renderGPU() {
  try {
    var resultCount = ximlen * yimlen;
    var gpu = await initTaylorGPU(resultCount);
    if (!gpu) return false;
    var batch = await calcTaylorPixelsGPU(gpu, ximlen, yimlen);
    if (!batch) return false;
    var gpuImg = ctx.createImageData(ximlen, yimlen);
    var data   = gpuImg.data;
    var pixels = batch.pixels;
    for (var idx = 0; idx < ximlen * yimlen; idx++) {
      var off = idx * 4;
      data[off]     = Math.min(255, Math.max(0, pixels[off]    )) | 0;
      data[off + 1] = Math.min(255, Math.max(0, pixels[off + 1])) | 0;
      data[off + 2] = Math.min(255, Math.max(0, pixels[off + 2])) | 0;
      data[off + 3] = 255;
    }
    batch.release();
    ctx.putImageData(gpuImg, 0, 0);
    updateStatus('Done (GPU)');
    return true;
  } catch(e) {
    console.warn('Taylor WebGPU render failed; falling back to CPU.', e);
    taylorGPUDisabled = true;
    destroyTaylorGPU();
    return false;
  }
}

function draw() {
  stopRendering();
  initCanvas();
  clearAndReset();
  updateHashTag();
  if (!taylorGPUDisabled && isTaylorWebGPUSupported()) {
    updateStatus('Rendering (GPU)…');
    renderGPU().then(function(ok) {
      if (!ok) startRendering();
    });
    return;
  }
  startRendering();
}

function resetAndDraw() {
  zoom         = 9;
  xcen         = 0.0;
  ycen         = 0.01;
  depth        = 20;
  glow         = 0.01;
  brightness   = 1.28;
  complexScale = 0;
  cx = 0; cy = 0;
  scalePoint    = false;
  compoundPoint = false;
  recurse       = false;
  mandelbrot    = false;
  alg           = 0;
  writeControls();
  draw();
}

// ---------------------------------------------------------------------------
// Mouse / zoom interaction (port of Java mouse handlers)
// ---------------------------------------------------------------------------

function setupMouse() {
  ccanvas.oncontextmenu = function(e) { e.preventDefault(); };

  ccanvas.onmousedown = function(e) {

    if (e.button === 2) {
      // right click hide control panel
      var desc = $('description');
      desc.style.display = (desc.style.display === 'none') ? 'block' : 'none';
      return;
    }

    if (e.button === 0 && e.shiftKey) {
      // Shift+left: set compound expansion point
      cx =  ((e.clientX / ximlen) - 0.5) * 2 * zoom * (ximlen / yimlen) + xcen;
      cy = -((e.clientY / yimlen) - 0.5) * 2 * zoom + ycen;  // y flipped (imaginary axis)
      compoundPoint = true;
      mandelbrot = true;
      $('compoundPoint').checked = true;
      $('mandelbrotMode').checked = mandelbrot;
      writeControls();
      draw();
      return;
    }

    if (e.button === 0) {
      dragBox = [e.clientX, e.clientY, e.clientX, e.clientY];
    }
  };

  ccanvas.onmousemove = function(e) {
    if (dragBox) {
      dragBox[2] = e.clientX;
      dragBox[3] = e.clientY;
      cctx.clearRect(0, 0, ccanvas.width, ccanvas.height);
      cctx.strokeStyle = '#FF3B03';
      cctx.lineWidth = 1;
      cctx.strokeRect(dragBox[0], dragBox[1], dragBox[2] - dragBox[0], dragBox[3] - dragBox[1]);
    }
  };

  ccanvas.onmouseup = function(e) {
    if (dragBox && e.button === 0) {
      cctx.clearRect(0, 0, ccanvas.width, ccanvas.height);

      var dx   = Math.abs(dragBox[2] - dragBox[0]);
      var dy   = Math.abs(dragBox[3] - dragBox[1]);
      var dmax = Math.max(dx, dy);

      if (dmax >= 10) {
        var midX = (dragBox[0] + dragBox[2]) / 2;
        var midY = (dragBox[1] + dragBox[3]) / 2;

        // Port of Java mouseReleased coordinate mapping
        var aspect  = ximlen / yimlen;
        var newXcen = ((midX / ximlen) - 0.5) * 2 * zoom * aspect + xcen;
        var newYcen = ((midY / yimlen) - 0.5) * 2 * zoom + ycen;
        // zoom = y half-range; derive from whichever drag dimension is larger
        var newZoom = Math.max(dx / ximlen, dy / yimlen) * zoom;

        xcen = newXcen;
        ycen = newYcen;
        zoom = newZoom;

        writeControls();
        draw();
      }
      dragBox = null;
    }
  };

  window.onresize = function() {
    initCanvas();
    draw();
  };
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

function main() {
  initCanvas();
  setupMouse();

  // CPU rendering sometimes contains more detail and less rendering artifacts than GPU
  // likely due to float 32 limits
  $('drawButton').onclick = () => { taylorGPUDisabled = true; draw(); };
  $('resetButton').onclick = () => { taylorGPUDisabled = false; resetAndDraw(); };

  // Re-render automatically when any control changes
  var textIds     = ['depth', 'zoom', 'xcen', 'ycen', 'complexScale'];
  var toggleIds   = ['scalePoint', 'compoundPoint', 'recurse', 'mandelbrotMode'];
  var selectIds   = ['algorithm'];

  textIds.forEach(function(id) { $(id).onkeyup = draw; });
  $('brightnessSlider').onchange = function() { draw(); };
  $('glowSlider').onchange = function() { draw(); };
  toggleIds.forEach(function(id) { $(id).onchange = draw; });
  selectIds.forEach(function(id) { $(id).onchange = draw; });

  $('viewPNG').onclick = function() {
    var link = document.createElement('a');
    link.download = 'taylor.png';
    link.href = canvas.toDataURL('image/png');
    link.click();
  };

  window.onhashchange = function() {
    if (location.hash.replace(/^#/, '') === lastWrittenHash) return;
    if (readHashTag()) {
      writeControls();
      draw();
    }
  };

  readHashTag();
  writeControls();
  draw();
}

main();
