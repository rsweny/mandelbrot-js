/*
 * Newton's Method fractal, in HTML5 canvas and javascript.
 * https://github.com/rsweny/mandelbrot-js
 *
 * Copyright (C) 2018 Ryan Sweny
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

//constants
var mode = 0;
var iterations = 20;
var zero = Complex["0"];
var two = Complex(2, 0);
var one = Complex(1, 0);
var pointone = Complex(0.1, 0);

var complex_error = 0;
var oneError = Complex(1,complex_error);

var order = 3;
var img_order = 0;

var zoomStart = 3.4;
var zoom = [zoomStart, zoomStart];
var lookAtDefault = [-0.6, 0];
var lookAt = lookAtDefault;
var xRange = [0, 0];
var yRange = [0, 0];
var interiorColor = [0, 0, 0, 255];
var reInitCanvas = true; // Whether to reload canvas size, etc
var dragToZoom = true;
var colors = [[0,0,0,0]];
var renderId = 0; // To zoom before current render is finished

// WebGPU accelerates the per-pixel Newton iterations.  Root discovery and
// colour mapping stay on the CPU so the existing scanline ordering and palette
// side effects are preserved.
var newtonGPU = null;
var newtonGPUDisabled = false;
var newtonGPUStatus = 'CPU';
var newtonGPUInFlight = false;

var NEWTON_PIXEL_WGSL = `
struct Params {
  width: u32,
  height: u32,
  iterations: u32,
  mode: u32,
  mandelbrotAddition: u32,
  samples: u32,
  _pad1: u32,
  _pad2: u32,
  xStart: f32,
  yStart: f32,
  dx: f32,
  dy: f32,
  order: f32,
  imgOrder: f32,
  complexError: f32,
  rootBoundary: f32,
};

@group(0) @binding(0) var<uniform> P: Params;
// x=hue, y=iteration count, z=final real, w=final imaginary
@group(0) @binding(1) var<storage, read_write> outPixels: array<vec4<f32>>;

fn cadd(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> { return a + b; }
fn csub(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> { return a - b; }
fn cmul(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {
  return vec2<f32>(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}
fn cdiv(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {
  let d = b.x * b.x + b.y * b.y;
  if (d == 0.0) { return vec2<f32>(0.0 / d, 0.0 / d); }
  return vec2<f32>((a.x * b.x + a.y * b.y) / d, (a.y * b.x - a.x * b.y) / d);
}
fn clog(z: vec2<f32>) -> vec2<f32> {
  return vec2<f32>(log(length(z)), atan2(z.y, z.x));
}
fn cpow(z: vec2<f32>, w: vec2<f32>) -> vec2<f32> {
  if (z.x == 0.0 && z.y == 0.0) { return vec2<f32>(0.0, 0.0); }
  let logr = log(length(z));
  let theta = atan2(z.y, z.x);
  let newR = exp(w.x * logr - w.y * theta);
  let newTheta = w.y * logr + w.x * theta;
  return vec2<f32>(newR * cos(newTheta), newR * sin(newTheta));
}

fn hash01(v: u32) -> f32 {
	var x = v;
	x = ((x >> 16u) ^ x) * 0x45d9f3bu;
	x = ((x >> 16u) ^ x) * 0x45d9f3bu;
	x = (x >> 16u) ^ x;
	return f32(x & 0x00ffffffu) / 16777216.0;
}

fn newtonStepGpu(z: vec2<f32>, seed: vec2<f32>) -> vec2<f32> {
  let expo = vec2<f32>(P.order, P.imgOrder);
  let one = vec2<f32>(1.0, 0.0);
  let two = vec2<f32>(2.0, 0.0);
  let pointOne = vec2<f32>(0.1, 0.0);
  let oneError = vec2<f32>(1.0, P.complexError);
  let expM1 = csub(expo, oneError);

  switch (P.mode) {
    case 0u: {
      let f = csub(cpow(z, expo), one);
      let df = cmul(expo, cpow(z, expM1));
      return csub(z, cdiv(f, df));
    }
    case 1u: {
      let f = csub(cpow(z, expo), cdiv(one, z));
      let df = cadd(cmul(expo, cpow(z, expM1)), cpow(z, vec2<f32>(-2.0, 0.0)));
      return csub(z, cdiv(f, df));
    }
    case 2u: {
      let ten = vec2<f32>(10.0, 0.0);
      let nine = vec2<f32>(9.0, 0.0);
      let p2i = vec2<f32>(0.0, 0.2);
      let f = csub(cadd(cpow(z, ten), cmul(p2i, cpow(z, expo))), one);
      let df = cadd(cmul(ten, cpow(z, nine)), cmul(cmul(p2i, expo), cpow(z, expM1)));
      return csub(z, cdiv(f, df));
    }
    case 3u: {
      let f = cadd(csub(cmul(two, cpow(z, expo)), seed), one);
      let df = csub(cmul(cmul(two, expo), cpow(z, expM1)), one);
      return csub(z, cdiv(f, df));
    }
    case 4u: {
      let c = vec2<f32>(1.0 + P.order, P.imgOrder);
      let cMinus1Im = select(P.imgOrder - 1.0, 0.0, P.imgOrder == 0.0);
      let cMinus1 = vec2<f32>(P.order, cMinus1Im);
      let f = cadd(csub(cpow(z, c), z), pointOne);
      let df = csub(cmul(c, cpow(z, cMinus1)), one);
      return csub(z, cdiv(f, df));
    }
    case 5u: {
      let three = vec2<f32>(3.0, 0.0);
      let six = vec2<f32>(6.0, 0.0);
      let five = vec2<f32>(5.0, 0.0);
      let four = vec2<f32>(4.0, 0.0);
      let fifteen = vec2<f32>(15.0, 0.0);
      let eighteen = vec2<f32>(18.0, 0.0);
      let f = cadd(csub(cadd(csub(cpow(z, expo), cmul(three, cpow(z, five))), cmul(six, cpow(z, three))), cmul(three, z)), three);
      let df = csub(cadd(csub(cmul(expo, cpow(z, expM1)), cmul(fifteen, cpow(z, four))), cmul(eighteen, cpow(z, two))), three);
      return csub(z, cdiv(f, df));
    }
    default: {
      let con = vec2<f32>(P.order, P.imgOrder);
      let zToZ = cpow(z, z);
      let f = csub(zToZ, cmul(con, z));
      let df = csub(cmul(zToZ, cadd(one, clog(z))), con);
      return csub(z, cdiv(f, df));
    }
  }
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
	var seed = vec2<f32>(P.xStart + f32(x) * P.dx, P.yStart + f32(y) * P.dy);
	if (P.samples > 1u) {
	let sx = hash01(pixelIdx * 1664525u + sampleIdx * 1013904223u + 17u);
	let sy = hash01(pixelIdx * 22695477u + sampleIdx * 1103515245u + 29u);
	seed = seed - vec2<f32>(sx * P.dx * 0.5, sy * P.dy * 0.5);
	}

  var n = 0u;
  var old = seed;
  var z = newtonStepGpu(seed, seed);
  var hue = 0.0;

  loop {
    if (n >= P.iterations) { break; }
    if (!(distance(old, z) > P.rootBoundary)) { break; }
    old = z;
    z = newtonStepGpu(z, seed);
    if (P.mandelbrotAddition != 0u) {
      z = z + seed * 0.5;
    }
    let delta = length(z - old);
    let w = 1.0 / delta;
    hue = hue + pow(1.05, -w);
    n = n + 1u;
  }

  outPixels[idx] = vec4<f32>(hue, f32(n), z.x, z.y);
}
`;

var NEWTON_GPU_PARAMS_BYTES = 64;
var NEWTON_GPU_WORKGROUP_X = 16;
var NEWTON_GPU_WORKGROUP_Y = 16;

function isNewtonWebGPUSupported()
{
  return typeof navigator !== 'undefined' && !!navigator.gpu;
}

function destroyNewtonGPU()
{
  if (!newtonGPU) return;
  var bufs = ['paramsBuf', 'outputBuf', 'readBuf'];
  for (var i = 0; i < bufs.length; i++) {
    try { newtonGPU[bufs[i]].destroy(); } catch(e) {}
  }
  newtonGPU = null;
}

async function initNewtonGPU(resultCount)
{
  if (newtonGPUDisabled || !isNewtonWebGPUSupported()) return null;
	  if (newtonGPU && newtonGPU.resultCapacity >= resultCount) return newtonGPU;

  destroyNewtonGPU();
  try {
    var adapter = await navigator.gpu.requestAdapter();
    if (!adapter) throw new Error('No WebGPU adapter');
    var device = await adapter.requestDevice();
	var outputBytes = resultCount * 4 * 4;

    var paramsBuf = device.createBuffer({ size: NEWTON_GPU_PARAMS_BYTES, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    var outputBuf = device.createBuffer({ size: outputBytes, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC });
    var readBuf = device.createBuffer({ size: outputBytes, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
    var module = device.createShaderModule({ code: NEWTON_PIXEL_WGSL });
    var pipeline = await device.createComputePipelineAsync({ layout: 'auto', compute: { module: module, entryPoint: 'main' } });
    var bindGroup = device.createBindGroup({
      layout: pipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: paramsBuf } },
        { binding: 1, resource: { buffer: outputBuf } },
      ],
    });

    newtonGPU = {
      device: device,
	  resultCapacity: resultCount,
      outputBytes: outputBytes,
      paramsBuf: paramsBuf,
      outputBuf: outputBuf,
      readBuf: readBuf,
      pipeline: pipeline,
      bindGroup: bindGroup,
      paramsBytes: new ArrayBuffer(NEWTON_GPU_PARAMS_BYTES),
    };
    newtonGPUStatus = 'GPU';
    return newtonGPU;
  } catch (e) {
    console.warn('Newton WebGPU path disabled; falling back to CPU.', e);
    newtonGPUDisabled = true;
    newtonGPUStatus = 'CPU fallback';
    destroyNewtonGPU();
    return null;
  }
}

function writeNewtonGPUParams(gpu, width, height, dx, dy, sampleCount)
{
  var u = new Uint32Array(gpu.paramsBytes);
  var f = new Float32Array(gpu.paramsBytes);
  u[0] = width >>> 0;
  u[1] = height >>> 0;
  u[2] = (parseInt(iterations, 10) || 0) >>> 0;
  u[3] = mode >>> 0;
  u[4] = mandelbrotAddition ? 1 : 0;
  u[5] = sampleCount >>> 0;
  u[6] = u[7] = 0;
  f[8] = xRange[0];
  f[9] = yRange[0];
  f[10] = dx;
  f[11] = dy;
  f[12] = parseFloat(order) || 0;
  f[13] = parseFloat(img_order) || 0;
  f[14] = complex_error;
  f[15] = rootBoundry;
  gpu.device.queue.writeBuffer(gpu.paramsBuf, 0, gpu.paramsBytes);
}

async function calcNewtonPixelsGPU(width, height, dx, dy, samples)
{
	var sampleCount = Math.max(1, parseInt(samples, 10) || 1);
	var resultCount = width * height * sampleCount;
	var gpu = await initNewtonGPU(resultCount);
 	if (!gpu) return null;

	writeNewtonGPUParams(gpu, width, height, dx, dy, sampleCount);
	var outputBytes = resultCount * 4 * 4;
  var d = gpu.device;
  var encoder = d.createCommandEncoder();
  var pass = encoder.beginComputePass();
  pass.setPipeline(gpu.pipeline);
  pass.setBindGroup(0, gpu.bindGroup);
	  pass.dispatchWorkgroups(
	    Math.ceil(width / NEWTON_GPU_WORKGROUP_X),
	    Math.ceil((height * sampleCount) / NEWTON_GPU_WORKGROUP_Y)
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
    release: function() {
      try { gpu.readBuf.unmap(); } catch(e) {}
    },
  };
}

/*
 * Initialize canvas
 */
var canvas = $('canvasMandelbrot');
canvas.width  = window.innerWidth;
canvas.height = window.innerHeight;
//
var ccanvas = $('canvasControls');
ccanvas.width  = window.innerWidth;
ccanvas.height = window.innerHeight;
//
var ctx = canvas.getContext('2d');
var img = ctx.createImageData(canvas.width, 1);

/*
 * Just a shorthand function: Fetch given element, jQuery-style
 */
function $(id)
{
	return document.getElementById(id);
}

function focusOnSubmit()
{
	var e = $('submitButton');
	if ( e ) e.focus();
}

function getSamples()
{
	// more than 4x anti-alias samples can overload the GPU
	var i = Math.min(4, parseInt($('superSamples').value));
	return !isFinite(i) || i <= 0 ? 2 : i;
}

function iterateEquation(i, j, equation, derivative) 
{
	var n = 0;
	var z = new Complex(i,j);
	var old = new Complex(i,j);
	
	z = z.sub( equation(z,i,j).div(derivative(z,i,j)) );
	
	var hue = 0.0;
	var w = 0.0;
	while(n < iterations && distance(old,z) > rootBoundry) 
	{
		old = z;
		z = z.sub( equation(z,i,j).div(derivative(z,i,j)) );
		if (mandelbrotAddition) z = z.add(Complex(i/2.0,j/2.0));	
		n++;
		
		//normal smoothing
		w = 1.0 / distance(z.sub(old), zero);
		hue += Math.pow(1.05, -w);
	}
	
	var vals = new Array();
	if (n != iterations) 
	{
		vals[0] = hue;
		//console.log(iterations + " " + n + " new root at " + i + " " + j + " " + distance(old,z));
		vals[1] = addRoot(z);
	}
	else
	{
		vals[0] = n;
		vals[1] = 0;
	}
	
	vals[2] = z.re;
	vals[3] = z.i;
	
	return vals;
}

var mandelbrotAddition = false;
var rootBoundry = 0.0001;
var roots = [];

/* Newton functions */
function addRoot(root)
{
	var colorFactor = 10/rootBoundry;
	for (var i = 0; i < roots.length; i++)
	{
		var c = roots[i];
		if ( distance(root, c) < rootBoundry )
		{
			//console.log("use existing root: " + distance(root, c));
			var retVal = i + colorFactor*distance(root, c);
			return Math.min(roots.length - 1, retVal);
		}
	}

	if (roots.length < 20)
	{
		roots.push(root);
		console.log("add root:" + roots.length);
	}
	return roots.length - 1;
}

function distance(a, b) 
{
	return Math.sqrt((a.re-b.re)*(a.re-b.re) + (a.i-b.i)*(a.i-b.i));
}


//Mode0: z^n - 1 = 0
function UnityFunction(z, i, j) 
{
	var exponent = Complex(order, img_order);
	return z.cPow(exponent).sub( one );
}
//Mode0: n*z^(n-1)
function UnityDerivative(z, i, j) 
{
	var exponent = Complex(order, img_order);
	var exponentLessOne = exponent.sub(oneError);
	return exponent.mult( z.cPow(exponentLessOne) );
}
	

//Mode1: z^n - 1 / z
function DivZFunction(z, i, j)
{
	var exponent = Complex(order, img_order);
	return (z.cPow(exponent)).sub( one.div(z) );
}
//Mode1: n*z^(n-1) + z^(-2)
function DivZDerivative(z, i, j)
{
	var minusTwo = new Complex(-2,0);
	var exponent = new Complex(order, img_order);
	var exponentLessOne = exponent.sub(oneError);
	return (exponent.mult( z.cPow(exponentLessOne) )).add( z.cPow(minusTwo) );
}


//Mode2: z^10 + 0.2i * z^5 - 1
function Poly2Function(z, i, j)
{
	var exponent = Complex(order, img_order);
	var point2i  = Complex(0, 0.2);
	var ten      = Complex(10, 0);
	return z.cPow(ten).add( point2i.mult(z.cPow(exponent)) ).sub(one);
}
//Mode2: 10z^9 + 0.2i*5*z^4
function Poly2Derivative(z, i, j)
{
	var exponent        = Complex(order, img_order);
	var exponentLessOne = exponent.sub(oneError);
	var ten             = Complex(10, 0);
	var nine            = Complex(9, 0);
	var point2i         = Complex(0, 0.2);
	return ten.mult(z.cPow(nine)).add( point2i.mult(exponent).mult(z.cPow(exponentLessOne)) );
}

//Mode3: 2z^n - c + 1  (c = pixel coordinate)
function PolyMFunction(z, i, j)
{
	var c = Complex(i, j);
	var exponent = Complex(order, img_order);
	return two.mult(z.cPow(exponent)).sub(c).add(one);
}
//Mode3: 2*n*z^(n-1) - 1
function PolyMDerivative(z, i, j)
{
	var exponent = Complex(order, img_order);
	var exponentLessOne = exponent.sub(oneError);
	return two.mult(exponent).mult(z.cPow(exponentLessOne)).sub(one);
}

//Mode4: z^c - z + 0.1  (c = 1+order)
function Poly3Function(z, i, j)
{
	var c = Complex(1 + parseFloat(order), img_order);
	return z.cPow(c).sub(z).add(pointone);
}

//Mode4: c*z^(c-1) - 1
function Poly3Derivative(z, i, j)
{
	var n = parseFloat(order);
	var c = Complex(1 + n, img_order);
	var cminus1 = Complex(n, (img_order === 0) ? 0 : img_order - 1);
	return c.mult(z.cPow(cminus1)).sub(one);
}


//Mode5: z^n - 3z^5 + 6z^3 - 3z + 3
function PolyFunction(z, i, j)
{
	var exponent = Complex(order, img_order);
	var three    = Complex(3, 0);
	var six      = Complex(6, 0);
	var five     = Complex(5, 0);
	return z.cPow(exponent)
		.sub( three.mult(z.cPow(five)) )
		.add( six.mult(z.cPow(three)) )
		.sub( three.mult(z) )
		.add( three );
}
//Mode5: n*z^(n-1) - 15z^4 + 18z^2 - 3
function PolyDerivative(z, i, j)
{
	var exponent        = Complex(order, img_order);
	var exponentLessOne = exponent.sub(oneError);
	var three    = Complex(3, 0);
	var four     = Complex(4, 0);
	var fifteen  = Complex(15, 0);
	var eighteen = Complex(18, 0);
	return exponent.mult(z.cPow(exponentLessOne))
		.sub( fifteen.mult(z.cPow(four)) )
		.add( eighteen.mult(z.cPow(two)) )
		.sub( three );
}


//Mode6: z^z - c*z  (c = order + img_order*i)
function ZZFunction(z, i, j)
{
	var con = Complex(order, img_order);
	return z.cPow(z).sub( con.mult(z) );
}
//Mode6: z^z * (1 + ln z) - c
function ZZDerivative(z, i, j)
{
	var con = Complex(order, img_order);
	var lnz = Complex.log(z).add(one);
	return z.cPow(z).mult(lnz).sub(con);
}

/*
 * Render the Mandelbrot set
 */
function draw(superSamples)
{
	mode = $("mode").selectedIndex;

	var equation, derivative;
	if (mode == 0) {
		equation = UnityFunction;
		derivative = UnityDerivative;
	} else if (mode == 1) {
		equation = DivZFunction;
		derivative = DivZDerivative;
	} else if (mode == 2) {
		equation = Poly2Function;
		derivative = Poly2Derivative;
	} else if (mode == 3) {
		equation = PolyMFunction;
		derivative = PolyMDerivative;
	} else if (mode == 4) {
		equation = Poly3Function;
		derivative = Poly3Derivative;
	} else if (mode == 5) {
		equation = PolyFunction;
		derivative = PolyDerivative;
	} else if (mode == 6) {
		equation = ZZFunction;
		derivative = ZZDerivative;
	}

	console.log("mode: " + mode);
	order = $("real_exponent").value;
	img_order = parseFloat($("complex_exponent").value);
	iterations = $("txtIterations").value;
	roots = [];


	if ( lookAt === null ) lookAt = [-0.6, 0];
	if ( zoom === null ) zoom = [zoomStart, zoomStart];

	xRange = [lookAt[0]-zoom[0]/2, lookAt[0]+zoom[0]/2];
	yRange = [lookAt[1]-zoom[1]/2, lookAt[1]+zoom[1]/2];

	mandelbrotAddition = $("mandelbrot").checked;

	initialColor = $("colorSlider").value / 100.0;
	contrast = $("contrastSlider").value / 100.0;
	brightness = Math.max(1.0, Math.min(10.0, parseFloat($("brightnessSlider").value) || 2.0));
	console.log(order + " " + img_order + " iterations: " + iterations + " initialColor: "  + initialColor + " mode: " + mode);

	if ( reInitCanvas ) {
		reInitCanvas = false;

		canvas = $('canvasMandelbrot');
		canvas.width  = window.innerWidth;
		canvas.height = window.innerHeight;

		ccanvas = $('canvasControls');
		ccanvas.width  = window.innerWidth;
		ccanvas.height = window.innerHeight;

		ctx = canvas.getContext('2d');
		img = ctx.createImageData(canvas.width, 1);

		adjustAspectRatio(xRange, yRange, canvas);
	}

	var dx = (xRange[1] - xRange[0]) / (0.5 + (canvas.width-1));
	var dy = (yRange[1] - yRange[0]) / (0.5 + (canvas.height-1));
	var Ci_step = (yRange[1] - yRange[0]) / (0.5 + (canvas.height-1));

	updateHashTag(superSamples);

	// Only enable one render at a time
	renderId += 1;

	function drawLineSuperSampled(Ci, off, Cr_init, Cr_step)
	{
		var Cr = Cr_init;

		for ( var x=0; x < canvas.width; ++x, Cr += Cr_step ) {
			var color = [255, 255, 255, 255];

			for ( var s=0; s<superSamples; ++s ) {
				var rx = Math.random()*Cr_step;
				var ry = Math.random()*Ci_step;
				var p = iterateEquation(Cr - rx/2, Ci - ry/2, equation, derivative);
				color = addRGB(color, pickColor(p));
			}

			color = divRGB(color, superSamples);

			img.data[off++] = color[0];
			img.data[off++] = color[1];
			img.data[off++] = color[2];
			img.data[off++] = 255;
		}
	}

	function drawLine(Ci, off, Cr_init, Cr_step) {
		let Cr = Cr_init;
		for (let x = 0; x < canvas.width; ++x, Cr += Cr_step) {
			let p = iterateEquation(Cr, Ci, equation, derivative);
			let color = pickColor(p);
			img.data[off++] = color[0];
			img.data[off++] = color[1];
			img.data[off++] = color[2];
			img.data[off++] = 255;
		}
	}

	function pickColor(iter)
	{
		return pickColorValues(iter[0], iter[1], iter[2], iter[3]);
	}

	function pickColorValues(iter0, rootIndex, finalRe, finalIm)
	{
		if (iter0 == iterations) {
			if (mandelbrotAddition || mode == 3) {
				return [0,0,0,255];
			} else {
				return [255,255,255,255];
			}
		}
		var hue;
		if (mandelbrotAddition) {
			//get Color based on angle of z
			var angle = Math.atan(finalIm/finalRe);
			angle = angle + (Math.PI / 2);
			angle = angle / Math.PI;
			hue = angle;
		} else {
			//limit the colors to 30% of the spectrum, unless a complex root
			var limitfactor = 0.3;
			if (img_order != 0) limitfactor = 1.0;
			hue = initialColor + (rootIndex*limitfactor) / roots.length;
		}

		//normal shading
		huesat = iter0/(iterations/brightness);
		if (huesat >= 1.0) huesat = 0.9999;

		huesat = huesat * 16000000;
		if (huesat > 16000000)  {
			huesat = Math.abs(16000000 + (16000000-huesat));
		}
		huesat = (huesat % 16000000)/16000000.0;

		new_gamma = Math.pow(huesat, contrast);

		var c = HSVtoRGB(hue, 1.0 - huesat, new_gamma);
		c.push(255); // alpha
		return c;
	}

	function drawSolidLine(y, color)
	{
		var off = y*canvas.width;

		for ( var x=0; x<canvas.width; ++x ) {
			img.data[off++] = color[0];
			img.data[off++] = color[1];
			img.data[off++] = color[2];
			img.data[off++] = color[3];
		}
	}

	function render()
	{
		var start  = (new Date).getTime();
		var startHeight = canvas.height;
		var startWidth = canvas.width;
		var lastUpdate = start;
		var updateTimeout = 200;
		var pixels = 0;
		var Ci = yRange[0];
		var sy = 0;
		var drawLineFunc = superSamples>1? drawLineSuperSampled : drawLine;
		var ourRenderId = renderId;

		var scanline = function()
		{
			if (renderId != ourRenderId || startHeight != canvas.height || startWidth != canvas.width) {
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
			if ( sy++ < canvas.height ) {
				if ( (now - lastUpdate) >= updateTimeout ) {
					// show the user where we're rendering
					drawSolidLine(0, [255,59,3,255]);
					ctx.putImageData(img, 0, sy);

					// Update speed and time taken
					var elapsedMS = now - start;
					$('renderTime').innerHTML = (elapsedMS/1000.0).toFixed(1); // 1 comma

					var speed = Math.floor(pixels / elapsedMS);

					if ( metric_units(speed).substr(0,3)=="NaN" ) {
						speed = Math.floor(60.0*pixels / elapsedMS);
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

		async function renderGPU()
		{
			if (newtonGPUInFlight) return false;
			newtonGPUInFlight = true;

			var start = (new Date).getTime();
			var startHeight = canvas.height;
			var startWidth = canvas.width;
			var ourRenderId = renderId;
			var requestedSamples = superSamples;

			function renderIsCurrent()
			{
				return renderId == ourRenderId && startHeight == canvas.height && startWidth == canvas.width;
			}

			function paintGPUBatch(batch)
			{
				if (!renderIsCurrent()) return true;

				var gpuImg = ctx.createImageData(canvas.width, canvas.height);
				var data = gpuImg.data;
				var pixels = batch.pixels;
				var iterationLimit = parseInt(iterations, 10) || 0;
				var sampleCount = batch.samples;

				for (var idx = 0; idx < canvas.width * canvas.height; idx++) {
					var colorSum = [0, 0, 0, 255];

					for (var s = 0; s < sampleCount; s++) {
						var poff = (idx * sampleCount + s) * 4;
						var hue = pixels[poff];
						var n = Math.round(pixels[poff + 1]);
						var zre = pixels[poff + 2];
						var zim = pixels[poff + 3];
						var rootIndex = 0;
						var iter0 = n;

						if (!isFinite(hue) || !isFinite(n) || !isFinite(zre) || !isFinite(zim)) {
							// GPU overflow/poles can produce NaN/Infinity; ImageData turns NaN colour channels into black, so treat these as non-converged white samples.
							colorSum = addRGB(colorSum, [255,255,255,255]);
							continue;
						}

						if (n == iterationLimit) {
							// Match the CPU renderer: samples that hit the iteration limit are interior/non-converged and should be white.
							if (mandelbrotAddition || mode == 3) {
								colorSum = addRGB(colorSum, [0,0,0,255]);
							} else {
								colorSum = addRGB(colorSum, [255,255,255,255]);
							}
							continue;
						}

						if (n != iterationLimit) {
							iter0 = hue;
							rootIndex = addRoot({re: zre, i: zim});
						}

						colorSum = addRGB(colorSum, pickColorValues(iter0, rootIndex, zre, zim));
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

			async function runGPUPass(samplesForPass)
			{
				var batch = null;
				try {
					batch = await calcNewtonPixelsGPU(canvas.width, canvas.height, dx, dy, samplesForPass);
				} catch (e) {
					console.warn('Newton WebGPU render failed; falling back to CPU.', e);
					newtonGPUDisabled = true;
					newtonGPUStatus = 'CPU fallback';
					destroyNewtonGPU();
					return false;
				}

				if (!batch) return false;
				try {
					return paintGPUBatch(batch);
				} finally {
					batch.release();
				}
			}

			try {
				// First pass deliberately uses one sample to discover roots quickly.
				// The second pass keeps those roots and recolours the full image.
				if (!await runGPUPass(1)) return false;
				if (!renderIsCurrent()) return true;
				if (!await runGPUPass(requestedSamples)) return false;

				var elapsedMS = Math.max(1, (new Date).getTime() - start);
				$('renderTime').innerHTML = (elapsedMS/1000.0).toFixed(1);
				$('renderSpeed').innerHTML = metric_units(Math.floor((canvas.width * canvas.height) / elapsedMS));
				$('renderSpeedUnit').innerHTML = 'second (GPU)';
				newtonGPUStatus = 'GPU';
				return true;
			} finally {
				newtonGPUInFlight = false;
			}
		}

			if (!newtonGPUDisabled && isNewtonWebGPUSupported()) {
			var fallbackRenderId = renderId;
			renderGPU().then(function(usedGPU) {
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
var logHalfBase = Math.log(0.5)*logBase;
var smoothColor = true;
var band = 1.0;
var contrast = 0.35;
var brightness = 2.0;
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

	var rgb = [0,0,0];
	rgb[0] = Math.floor(r * 255)
	rgb[1] = Math.floor(g * 255)
	rgb[2] = Math.floor(b * 255)

	return rgb;
}

/*
 * Update URL's hash with render parameters so we can pass it around.
 */
function updateHashTag(samples) {
	var scheme = $('mode').value;
	console.log("updateHashTag() " + iterations);
	location.hash = 'zoom=' + zoom + '&lookAt=' + lookAt + '&iterations=' + iterations + '&superSamples=' + samples +
		'&realExponent=' + encodeURIComponent($('real_exponent').value) +
		'&complexExponent=' + encodeURIComponent($('complex_exponent').value) +
		'&mandelbrot=' + ($('mandelbrot').checked ? '1' : '0') +
		'&mode=' + scheme;
}

/*
 * Parse URL hash tag, returns whether we should redraw.
 */
function readHashTag()
{
	var redraw = false;
	var tags = location.hash.split('&');

	for ( var i=0; i<tags.length; ++i ) {
		var tag = tags[i].split('=');
		var key = tag[0];
		var val = tag[1];

		switch ( key ) {
			case '#zoom': {
				var z = val.split(',');
				zoom = [parseFloat(z[0]), parseFloat(z[1])];
				redraw = true;
			} break;

			case 'lookAt': {
				var l = val.split(',');
				lookAt = [parseFloat(l[0]), parseFloat(l[1])];
				redraw = true;
			} break;

			case 'iterations': {
				$('txtIterations').value = val;
				console.log("readHashTag() " + $('txtIterations').value);
				redraw = true;
			} break;

			case 'superSamples': {
				$('superSamples').value = String(parseInt(val));
				redraw = true;
			} break;

			case 'realExponent': {
				$('real_exponent').value = decodeURIComponent(val);
				redraw = true;
			} break;

			case 'complexExponent': {
				$('complex_exponent').value = decodeURIComponent(val);
				redraw = true;
			} break;

			case 'mandelbrot': {
				$('mandelbrot').checked = val === '1' || val === 'true';
				redraw = true;
			} break;

			case 'mode': {
				$('mode').value = String(val);
				redraw = true;
			} break;
		}
	}

	if ( redraw )
		reInitCanvas = true;

	return redraw;
}

/*
 * Return number with metric units
 */
function metric_units(number)
{
	var unit = ["", "k", "M", "G", "T", "P", "E"];
	var mag = Math.ceil((1+Math.log(number)/Math.log(10))/3);
	return "" + (number/Math.pow(10, 3*(mag-1))).toFixed(2) + unit[mag];
}


/*
 * Adjust aspect ratio based on plot ranges and canvas dimensions.
 */
function adjustAspectRatio(xRange, yRange, canvas)
{
	var ratio = Math.abs(xRange[1]-xRange[0]) / Math.abs(yRange[1]-yRange[0]);
	var sratio = canvas.width/canvas.height;
	if ( sratio>ratio ) {
		var xf = sratio/ratio;
		xRange[0] *= xf;
		xRange[1] *= xf;
			zoom[0] *= xf;
	} else {
		var yf = ratio/sratio;
		yRange[0] *= yf;
		yRange[1] *= yf;
			zoom[1] *= yf;
	}
}

function addRGB(v, w)
{
	v[0] += w[0];
	v[1] += w[1];
	v[2] += w[2];
	v[3] += w[3];
	return v;
}

function divRGB(v, div)
{
	v[0] /= div;
	v[1] /= div;
	v[2] /= div;
	v[3] /= div;
	return v;
}


function main()
{
	$('viewPNG').onclick = function(event)
	{
		var link = document.createElement('a');
		link.download = `newton-${mode}.png`;
		link.href = canvas.toDataURL('image/png');
		link.click();
	};

	$('resetButton').onclick = function(even)
	{
		$('settingsForm').reset();
		setTimeout(function() { location.hash = ''; }, 1);
		zoom = [zoomStart, zoomStart];
		lookAt = lookAtDefault;
		reInitCanvas = true;
		draw(getSamples());
	};

	$("txtIterations").onchange = function() {
		draw(getSamples());
	}

	$("superSamples").onchange = function() {
		draw(getSamples());
	}

	$("mode").onchange = function() {
		draw(getSamples());
	}

	$("real_exponent").onchange = function() {
		draw(getSamples());
	}

	$("complex_exponent").onchange = function() {
		draw(getSamples());
	}

	$("mandelbrot").onchange = function() {
		draw(getSamples());
	}

	$("contrastSlider").onchange = function() {
		contrast = $("contrastSlider").value / 100.0;
		draw(getSamples());
	}

	$("brightnessSlider").onchange = function() {
		brightness = Math.max(1.0, Math.min(5.0, parseFloat($("brightnessSlider").value) || 2.0));
		draw(getSamples());
	}

	$("colorSlider").onchange = function() {
		initialColor = $("colorSlider").value / 100.0;
		draw(getSamples());
	}

	if ( dragToZoom == true ) {
		var box = null;

		$('canvasControls').oncontextmenu = function(e) {
			e.preventDefault();
		};

		$('canvasControls').onmousedown = function(e)
		{
			if ( e.button === 2 ) {
				var desc = $('description');
				desc.style.display = (desc.style.display === 'none') ? 'block' : 'none';
				return;
			}
			if ( box == null ) box = [e.clientX, e.clientY, 0, 0];
		}

		$('canvasControls').onmousemove = function(e)
		{
			if ( box != null ) {
				var c = ccanvas.getContext('2d');
				c.lineWidth = 1;

				// clear out old box first
				c.clearRect(0, 0, ccanvas.width, ccanvas.height);

				// draw new box
				c.strokeStyle = '#FF3B03';
				box[2] = e.clientX;
				box[3] = e.clientY;
				c.strokeRect(box[0], box[1], box[2]-box[0], box[3]-box[1]);
			}
		}

		var zoomOut = function(event) {
			var x = event.clientX;
			var y = event.clientY;

			var w = window.innerWidth;
			var h = window.innerHeight;

			var dx = (xRange[1] - xRange[0]) / (0.5 + (canvas.width-1));
			var dy = (yRange[1] - yRange[0]) / (0.5 + (canvas.height-1));

			x = xRange[0] + x*dx;
			y = yRange[0] + y*dy;

			lookAt = [x, y];

			if ( event.shiftKey ) {
				zoom[0] /= 0.5;
				zoom[1] /= 0.5;
			}

			draw(getSamples());
		};

		$('canvasControls').onmouseup = function(e)
		{
			if ( box != null ) {
				// Zoom out?
				if ( e.shiftKey ) {
					box = null;
					zoomOut(e);
					return;
				}

				/*
				 * Cleaer entire canvas
				 */
				var c = ccanvas.getContext('2d');
				c.clearRect(0, 0, ccanvas.width, ccanvas.height);

				/*
				 * Calculate new rectangle to render
				 */
				var x = Math.min(box[0], box[2]) + Math.abs(box[0] - box[2]) / 2.0;
				var y = Math.min(box[1], box[3]) + Math.abs(box[1] - box[3]) / 2.0;

				var dx = (xRange[1] - xRange[0]) / (0.5 + (canvas.width-1));
				var dy = (yRange[1] - yRange[0]) / (0.5 + (canvas.height-1));

				x = xRange[0] + x*dx;
				y = yRange[0] + y*dy;

				lookAt = [x, y];

				var xf = Math.abs(Math.abs(box[0]-box[2])/canvas.width);
				var yf = Math.abs(Math.abs(box[1]-box[3])/canvas.height);

				zoom[0] *= Math.max(xf, yf); // retain aspect ratio
				zoom[1] *= Math.max(xf, yf);

				box = null;
				draw(getSamples());
			}
		}
	}

	/*
	 * Enable zooming (currently, the zooming is inexact!) Click to zoom;
	 * perfect to mobile phones, etc.
	 */
	if ( dragToZoom == false ) {
		$('canvasMandelbrot').onclick = function(event)
		{
			var x = event.clientX;
			var y = event.clientY;
			var w = window.innerWidth;
			var h = window.innerHeight;

			var dx = (xRange[1] - xRange[0]) / (0.5 + (canvas.width-1));
			var dy = (yRange[1] - yRange[0]) / (0.5 + (canvas.height-1));

			x = xRange[0] + x*dx;
			y = yRange[0] + y*dy;

			lookAt = [x, y];

			if ( event.shiftKey ) {
				zoom[0] /= 0.5;
				zoom[1] /= 0.5;
			} else {
				zoom[0] *= 0.5;
				zoom[1] *= 0.5;
			}

			draw(getSamples());
		};
	}

	/*
	 * When resizing the window, be sure to update all the canvas stuff.
	 */
	window.onresize = function(event)
	{
		reInitCanvas = true;
	};

	/*
	 * Read hash tag and render away at page load.
	 */
	readHashTag();
	draw(getSamples());
}

main();
