/*
 * The Mandelbulb, in HTML5 canvas and javascript.
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

var pallet = palettes["Geyser"];
var mode = 0;

// detail level
var rayDetail = 0.003;
var stepDetail = 0.028;

// higher values here can reduce cloudiness at higher zooms since it biases points 
// near the already-found surface
var frost = 2.0;

// ray traced lighting
var LightVector = [ 0.12, 0.15, -0.19 ];
var AMBIENT_LIGHT = 8.0;
var RAY_STEPS = 32;
var MAX_DEPTH_STEPS = 4096;
var ray_step;
var primary_light = 28.0;
var shadow_darkness = 30.0;
var HORIZON = 20;

// fog based on path traces
var fog_factor = 0.0;
var fog_color = {r: 250, g: 250, b: 250};
var useVolumetricFog = false;
var min_y, max_y;

// formula variation
var formula = 0;
var azimuth = 1;

// settings that need saving
let iterations = 20;
var ximlen = 0;
var yimlen = 0;
var pal = 1;
var power = 8;
var gradient = 0.5;
var brightness = 1.8;
var zoom = 3.0;
var xcen = 0.0;
var ycen = 0.0;
var cameraPersp = 0.0;
var cameraYaw = 0.1;
var cameraPitch = 0.75;

var cameraDOF = 0.0;
var focus_depth = 1.0;
var focus = -0.15; // higher values for further distance in focus.
const FOCUS_WINDOW = 0.5 // 50% of depth of field is in focus

var opacity = 1.4;


// data
var root_zoom;
var occlusionPositions = null;
var img_alpha = null;
var img_red = null;
var img_green = null;
var img_blue = null;
var lastDir = "";
var reset = 0;

var CameraMatrix;
var IrotX;
var IrotZ;
var rotX;
var rotZ;

// mouse
var xcurr,ycurr,xanchor,yanchor;
var m_down;
var drawFocus = false;

// stats
var t1 = 0;
var visiblePixels, allPixels, rayPoints;
var renderpass = 0;
var max_alpha = 1;

// worker pool (web workers, used as fallback when WebGPU is unavailable)
var workerPool = [];
var workerGen   = 0;   // incremented on each restart; stale results are discarded
var renderY     = 0;   // next row to dispatch

// GPU renderer state. Populated by initGPURenderer() during startup.
var useGPU = false;

// Build a Blob URL from the inert <script type="text/mandelbulb-worker"> block so that
// workers load correctly on file:// origins (no HTTP server required).
var WORKER_BLOB_URL = (function() {
    var src = document.getElementById('workerScript').textContent;
    var blob = new Blob([src], { type: 'application/javascript' });
    return URL.createObjectURL(blob);
}());

// initialize canvas
var canvas = $('canvasMandelbrot');
canvas.width  = 1800;
canvas.height = 1800;
var ccanvas = $('canvasControls');
ccanvas.width  = 1800;
ccanvas.height = 1800;
var ctx = canvas.getContext('2d');
var ctx_img = ctx.createImageData(canvas.width, canvas.height);

// JQuery shorthand
function $(id) { return document.getElementById(id) }

function updateMinMaxY()
{
	console.log("--------------------------updateMinMaxY()");
	min_y = HORIZON;
	max_y = -HORIZON;
	for (var i=0; i<ximlen; i++)
	{
		for (var j=0; j<yimlen; j++)
		{
			var _occ = occlusionPositions[i * yimlen + j];
			if (_occ < min_y)
			{
				min_y = _occ - stepDetail;
			}
			else if (_occ != HORIZON && _occ > max_y)
			{
				max_y = _occ + stepDetail;
			}
		}
	}
	focus_depth = (max_y - min_y) * FOCUS_WINDOW;
	console.log(focus_depth + " Y Bounds: " + min_y + " to " + max_y);
}
	
function findPeak(arr)
{
	var max = 0;
	for (var i = 1; i < ximlen; i++)
	{
		for (var j = 1; j < yimlen; j++)
		{
			if (arr[i][j] > max)
			{
				max = arr[i][j];
			}
		}
	}
	return max;
}

function clearScreenAndReset()
{
	renderpass = 0;
	max_alpha = 1;
	min_y = -2.0;
	max_y = 2.0;
	
	occlusionPositions = new Float32Array(ximlen * yimlen); // initialised to 0
	img_alpha = matrix(canvas.width, canvas.height, 0.0);
	img_red = matrix(canvas.width, canvas.height, 0.0);
	img_green = matrix(canvas.width, canvas.height, 0.0);
	img_blue = matrix(canvas.width, canvas.height, 0.0);

	m_down = false;
}
	
function setCamera()
{
	// 3D camera precalc
	console.log("-----------------------setCamera");
	rotX = RotateX(cameraPitch);
	rotZ = RotateZ(cameraYaw);
	CameraMatrix = matrixMult(rotZ, rotX);
	IrotX = RotateX(-cameraPitch);
	IrotZ = RotateZ(-cameraYaw);
}
	
function round7(v)
{
	return Math.round(v * 1e7) / 1e7;
}

function setZoom(z)
{
	zoom = z;
	root_zoom = Math.pow(zoom, 0.5);

	// at zoom < 0.02 higher values result in better shadows, lower values result in sharp images
	const adjRayStep = zoom < 0.02 ? 0.0015 : 0.008;
	ray_step = rayDetail*Math.max(adjRayStep, zoom);
	console.log("ray_step " + ray_step);
}
	
function reDraw()
{
	readHashTag();
	updateHashTag();
	draw(false);
}

/*
 * Render the Mandelbrot set
 */
function draw(startScanning)
{
	console.log("--------------------------draw()");
	updateHashTag();
	render(startScanning);
}

function generateConfig() {
	const half_ximlen = ximlen / 2;
	const half_yimlen = yimlen / 2;
	return {
		ximlen,
		yimlen,
		half_ximlen,
		half_yimlen,
		zoom: zoom, xcen: xcen, ycen: ycen,
		cameraPersp: cameraPersp,
		iterations: iterations, formula: formula, azimuth: azimuth, power: power,
		stepDetail: stepDetail, frost: frost,
		root_zoom: root_zoom, opacity: opacity,
		focus: focus, focus_depth: focus_depth,
		cameraDOF: cameraDOF, factorDOF: cameraDOF * (ximlen / 3),
		LightVector: LightVector, RAY_STEPS: RAY_STEPS,
		AMBIENT_LIGHT: AMBIENT_LIGHT, primary_light: primary_light,
		shadow_darkness: shadow_darkness, HORIZON: HORIZON,
		ray_step: ray_step,
		CameraMatrix: CameraMatrix, IrotX: IrotX, IrotZ: IrotZ,
		fog_factor: fog_factor, fog_color: fog_color,
		useVolumetricFog: useVolumetricFog,
		pallet: pallet
	};
}

// Push the latest config to every worker so non-reset controls (fog, focus,
// DOF, primary light, etc.) take effect on the next dispatched row without
// terminating and re-initialising the pool.
function broadcastConfig() {
	if (workerPool.length === 0) return;
	var config = generateConfig();
	workerPool.forEach(function(w) {
		w.postMessage({ type: 'configUpdate', config: config });
	});
}

// Snapshot of all globals the WebGPU renderer reads each frame. Any global
// update (camera, palette, iterations, etc.) is picked up here; setting
// reset = 1 drains the accumulators before the next frame is encoded.
function gpuStateFn() {
	if (reset === 1) {
		reset = 0;
		// Stop the progressive render, paint the WebGL preview onto the main
		// canvas, then restart the GPU loop. Returning null aborts this frame before any blit
		// so the preview isn't overwritten.
		stopGPURender();
		renderPreview(function() {
			setTimeout(function() {
				clearScreenAndReset();
				startGPURender(gpuStateFn, gpuStatsCallback);
			}, 4000);
		});
		return null;
	} else if (reset === 2) {
		reset = 0;
	}
	var passIdx = getGPUPassIndex();
	// Every 8 passes, pad the depth-march range outwards so the GPU briefly
	// explores a wider slab. Paired with _clearYBounds on the GPU side, this
	// gives the atomic min/max reduction a chance to discover surface points
	// just outside the previously-tight bounds, then re-converge over the
	// next several passes.
	var passMinY = min_y;
	var passMaxY = max_y;
	if (passIdx > 0 && passIdx % 8 === 0) {
		var pad = 10 * stepDetail;
		passMinY = min_y - pad;
		passMaxY = max_y + pad;
	}
	return {
		ximlen: ximlen, yimlen: yimlen,
		half_ximlen: ximlen / 2, half_yimlen: yimlen / 2,
		zoom: zoom, xcen: xcen, ycen: ycen,
		cameraPersp: cameraPersp,
		iterations: iterations, formula: formula, azimuth: azimuth, power: power,
		stepDetail: stepDetail, frost: frost,
		root_zoom: root_zoom, opacity: opacity,
		focus: focus, focus_depth: focus_depth,
		cameraDOF: cameraDOF, factorDOF: cameraDOF * (ximlen / 3),
		LightVector: LightVector, RAY_STEPS: RAY_STEPS,
		MAX_DEPTH_STEPS: zoom < 0.01 ? MAX_DEPTH_STEPS : MAX_DEPTH_STEPS*2,
		AMBIENT_LIGHT: AMBIENT_LIGHT, primary_light: primary_light,
		shadow_darkness: shadow_darkness, HORIZON: HORIZON,
		ray_step: ray_step,
		CameraMatrix: CameraMatrix, IrotX: IrotX, IrotZ: IrotZ,
		fog_factor: fog_factor, fog_color: fog_color,
		useVolumetricFog: useVolumetricFog,
		pallet: pallet,
		gradient: gradient, brightness: brightness, drawFocus: drawFocus,
		// Pass 0 is the rough depth-finder (frost/fog disabled via min_y == -2.0).
		// Subsequent passes use the live min_y/max_y refined by gpuStatsCallback
		// from the GPU-side atomicMin/atomicMax reductions over surface depths,
		// with a periodic 10*stepDetail widening (see passMinY/passMaxY above).
		min_y: (passIdx === 0 ? -2.0 : passMinY),
		max_y: (passIdx === 0 ?  2.0 : passMaxY)
	};
}

function gpuStatsCallback(s) {
	renderpass = s.pass
	if (typeof s.max_alpha === 'number') {
		max_alpha = Math.pow(Math.max(s.max_alpha, 1.0), gradient);
	}
	// Tighten the depth-march range using bounds reduced on the GPU. The
	// shader fills minY/maxY only when a surface is actually plotted, so we
	// ignore the seed sentinels (HORIZON / -HORIZON) until the first hit.
	if (typeof s.min_y === 'number' && typeof s.max_y === 'number' &&
		s.min_y < HORIZON && s.max_y > -HORIZON && s.max_y > s.min_y) {
		min_y = s.min_y - stepDetail;
		max_y = s.max_y + stepDetail;
		focus_depth = (max_y - min_y) * FOCUS_WINDOW;
		console.log("GPU pass " + s.pass + " Y bounds: " + min_y.toFixed(4) + " to " + max_y.toFixed(4) + " (focus_depth " + focus_depth.toFixed(4) + ")");
	}
	$('renderSpeed').innerHTML = Math.round(s.pixelsPerSec/1000) + "k px/sec pass: " + s.pass;
}

function render(startScanning)
{
	if (!startScanning) return;  // workers handle reset signals themselves

	if (useGPU) {
		startGPURender(gpuStateFn, gpuStatsCallback);
		return;
	}

	var start      = (new Date).getTime();
	var lastUpdate = start;
	var pixels     = 0;

	// Terminate any previously running workers
	workerPool.forEach(function(w) { w.terminate(); });
	workerPool = [];
	workerGen++;
	renderY = 0;

	var gen = workerGen;  // captured for stale-result detection

	function applyResult(data) {
		renderpass++;
		pixels += ximlen;
		visiblePixels += data.stats.visiblePixels;
		allPixels     += data.stats.allPixels;
		rayPoints     += data.stats.rayPoints;

		// Replicate the per-pixel occlusion reset that the original gridpoints() did
		// for every x in the processed row before searching for a surface.
		var yRow = data.yRow;
		for (var xd = 0; xd < ximlen; xd++) {
			occlusionPositions[xd * yimlen + yRow] = HORIZON;
		}

		// Apply pixel writes: flat Float64Array [x, y, depth, r, g, b, a, setOcc, ...]
		// setOcc=1 → surface pixel, update occlusionPositions; setOcc=0 → fog, colour only.
		var pw = data.pixelWrites;
		for (var i = 0; i < pw.length; i += 8) {
			var px = pw[i], py = pw[i+1], depth = pw[i+2];
			img_red[px][py]   += pw[i+3];
			img_green[px][py] += pw[i+4];
			img_blue[px][py]  += pw[i+5];
			img_alpha[px][py] += pw[i+6];
			if (pw[i+7]) occlusionPositions[px * yimlen + py] = depth;
		}

		if (renderpass % ximlen == 0) updateMinMaxY();

		if (renderpass % 50 == 0) {
			const t2 = (new Date()).getTime();
			const completeness = Math.round(max_alpha*100)/100.0;
			const strStatus = (new Date()).toISOString().substring(0,19) + " " + data.yRow + " " + visiblePixels + " " + rayPoints + " " + allPixels + " Pass: " + renderpass + " max value: " + completeness + " in " + (t2-t1);
			t1 = (new Date()).getTime();
			console.log(strStatus);
			updateHistogram();
			visiblePixels = 0;
			allPixels = 0;
			rayPoints = 0;
			// Send a fresh occlusion snapshot to all workers so they can filter
			// fog and ray-traced goodPoints against up-to-date surface depths.
			var snap = occlusionPositions.slice();
			workerPool.forEach(function(w) {
				w.postMessage({ type: 'occlusion', data: snap });
			});
		}

		var now = (new Date).getTime();
		if ((now - lastUpdate) >= 10000) {
			var elapsed = (now - start) / 1000.0;
			$('renderSpeed').innerHTML = Math.floor(pixels / elapsed) + " px/sec";
			lastUpdate = now;
		}
	}

	function dispatchRow(worker) {
		if (workerGen !== gen) return;  // this worker belongs to an old generation

		// Send updated depth bounds so workers benefit from updateMinMaxY() refinements
		worker.postMessage({ type: 'compute', yRow: renderY, min_y: min_y, max_y: max_y });
		renderY++;
		if (renderY > yimlen - 1) renderY = 0;
	}

	function onResult(worker) {
		return function(e) {
			if (workerGen !== gen) return;  // stale result from a superseded render

			applyResult(e.data);

			if (reset === 1) {
				console.log("reset render");
				reset = 0;
				workerPool.forEach(function(w) { w.terminate(); });
				workerPool = [];
				renderPreview(function() {
				 	setTimeout(function() { clearScreenAndReset(); render(true); }, 9000);
				});
				return;
			}
			if (reset === 2) {
				reset = 0;
				updateHistogram();
			}

			// When mouse is held down, yield to the browser for 2 s before next row
			if (m_down) {
				setTimeout(function() { dispatchRow(worker); }, 2000);
			} else {
				dispatchRow(worker);
			}
		};
	}

	var NUM_WORKERS = 3;
	var config = generateConfig();
	for (var i = 0; i < NUM_WORKERS; i++) {
		var w = new Worker(WORKER_BLOB_URL);
		w.postMessage({ type: 'init', config: config });
		w.onmessage = onResult(w);
		workerPool.push(w);
		dispatchRow(w);
	}
}


/*
 * Render a preview pass using the WebGL2 fragment-shader renderer from
 * gpu-preview.js. Produces a fully-lit image in well under a second, then
 * hands off to the CPU progressive renderer for final quality.
 */
function renderPreview(onComplete) {
	renderPreviewGPU(canvas, ctx);
	if (onComplete) onComplete();
}

function updateHistogram()
{
	var off = 0;
	var red = 0;
	var green = 0;
	var blue = 0;
	
	max_alpha = Math.pow(findPeak(img_alpha), gradient);
	
	for (var y=0; y<canvas.height; y++)
	{
		for (var x=0; x<canvas.width; x++)
		{
			var z = Math.pow(img_alpha[x][y], gradient)*brightness/max_alpha;
			red = ( (img_red[x][y]*z)/img_alpha[x][y] );
			green = ( (img_green[x][y]*z)/img_alpha[x][y] );
			blue = ( (img_blue[x][y]*z)/img_alpha[x][y] );

			//tool to view focus bounds
			if (drawFocus) {
				var _occ2 = occlusionPositions[x * yimlen + y];
				if (_occ2 > focus) {
					green = 20; red = 20; blue += 50;
				} else if (_occ2 < focus - focus_depth) {
					green = 20; blue = 20; red += 50;
				}
			}

			if (red > 255) red = 255;
			if (green > 255) green = 255;
			if (blue > 255) blue = 255;

			ctx_img.data[off++] = red;
			ctx_img.data[off++] = green;
			ctx_img.data[off++] = blue;
			ctx_img.data[off++] = 255;
		}
	}
	ctx.putImageData(ctx_img, 0, 0);
	drawFocus = false;
}

/*
 * When resizing the window, be sure to update all the canvas stuff.
 */
window.onresize = function(_event)
{
	reInitCanvas = true;
};

var redrawLightWidget = function() {};  // set by setupLightWidget

// Small trackball-style widget for rotating LightVector. Horizontal drag
// rotates about the Y axis (mixes x and z), vertical drag rotates about the
// X axis (mixes y and z). Magnitude is preserved.
function setupLightWidget()
{
	var widget = $('lightWidget');
	if (!widget) return;
	var wctx = widget.getContext('2d');
	var dragging = false;
	var lastX = 0, lastY = 0;
	var ROT_SPEED = 0.012;  // radians per pixel

	function drawWidget()
	{
		var w = widget.width, h = widget.height;
		var cx = w / 2, cy = h / 2;
		var R = Math.min(cx, cy) - 6;

		wctx.clearRect(0, 0, w, h);

		// sphere outline
		wctx.strokeStyle = '#555';
		wctx.lineWidth = 1;
		wctx.beginPath();
		wctx.arc(cx, cy, R, 0, Math.PI * 2);
		wctx.stroke();

		// crosshair
		wctx.strokeStyle = '#3a3a3a';
		wctx.beginPath();
		wctx.moveTo(cx - R, cy); wctx.lineTo(cx + R, cy);
		wctx.moveTo(cx, cy - R); wctx.lineTo(cx, cy + R);
		wctx.stroke();

		var lx = LightVector[0], ly = LightVector[1], lz = LightVector[2];
		var mag = Math.sqrt(lx*lx + ly*ly + lz*lz) || 1;

		// Transform fractal-space light vector into the renderer's screen-aligned
		// frame using the inverse camera rotations (same convention as reversePoint
		// in the worker): v[0] = screen horizontal, v[1] = depth (smaller = closer
		// to camera), v[2] = screen vertical (positive points down on screen).
		var v = [0, 0, 0];
		if (IrotZ && IrotX) {
			var tmp = [0, 0, 0];
			rotateVector(tmp, IrotZ, lx, ly, lz);
			rotateVector(v, IrotX, tmp[0], tmp[1], tmp[2]);
		} else {
			v[0] = lx; v[1] = ly; v[2] = lz;
		}
		var nx = v[0] / mag, ny = v[1] / mag, nz = v[2] / mag;

		// project onto widget: x right, screen-z down (renderer's +z is screen-down)
		var px = cx + nx * R;
		var py = cy + nz * R;

		// dot size/brightness modulated by depth (ny): smaller y = in front of bulb = bright
		var front = (1 - ny) * 0.5;  // 0..1, brighter when ny is negative (closer to camera)
		var radius = 4 + front * 4;
		var shade = Math.round(80 + front * 175);
		wctx.fillStyle = 'rgb(' + shade + ',' + Math.round(shade * 0.85) + ',' + Math.round(shade * 0.3) + ')';
		wctx.beginPath();
		wctx.arc(px, py, radius, 0, Math.PI * 2);
		wctx.fill();
		wctx.strokeStyle = ny <= 0 ? '#FF3B03' : '#222';
		wctx.lineWidth = 1;
		wctx.stroke();

		$('lightVecLabel').textContent = lx.toFixed(3) + ', ' + ly.toFixed(3) + ', ' + lz.toFixed(3);
	}

	function rotate(dxPx, dyPx)
	{
		var ay = -dxPx * ROT_SPEED;  // yaw: rotate around Y -> mixes x and z
		var ax = dyPx * ROT_SPEED;  // pitch: rotate around X -> mixes y and z

		var x = LightVector[0], y = LightVector[1], z = LightVector[2];

		var cosY = Math.cos(ay), sinY = Math.sin(ay);
		var x1 =  cosY * x + sinY * z;
		var z1 = -sinY * x + cosY * z;

		var cosX = Math.cos(ax), sinX = Math.sin(ax);
		var y1 =  cosX * y + sinX * z1;
		var z2 = -sinX * y + cosX * z1;

		LightVector[0] = x1;
		LightVector[1] = y1;
		LightVector[2] = z2;
	}

	widget.addEventListener('mousedown', function(e) {
		dragging = true;
		lastX = e.clientX; lastY = e.clientY;
		widget.style.cursor = 'grabbing';
		e.preventDefault();
		e.stopPropagation();
	});
	window.addEventListener('mousemove', function(e) {
		if (!dragging) return;
		var dx = e.clientX - lastX;
		var dy = e.clientY - lastY;
		lastX = e.clientX; lastY = e.clientY;
		if (dx === 0 && dy === 0) return;
		rotate(dx, dy);
		drawWidget();
		broadcastConfig();
	});
	window.addEventListener('mouseup', function() {
		if (!dragging) return;
		dragging = false;
		widget.style.cursor = 'grab';
		updateHashTag();
	});

	redrawLightWidget = drawWidget;
	drawWidget();
}

function init()
{
	yimlen = canvas.height;
	ximlen = canvas.height;

	CameraMatrix = matrix(3, 3, 0.0);
	IrotX = matrix(3, 3, 0.0);
	IrotZ = matrix(3, 3, 0.0);
	rotX = matrix(3, 3, 0.0);
	rotZ = matrix(3, 3, 0.0);

	clearScreenAndReset();
	setCamera();
}

function main()
{
	$('savePNG').onclick = function(_e)
	{
		var link = document.createElement('a');
		link.download = 'mandelbulb-' + renderpass + '-' + Math.round(max_alpha) + '.png';
		link.href = canvas.toDataURL('image/png');
		link.click();
	};

	$("zoomInput").onchange = function() {
		setZoom(parseFloat($("zoomInput").value));
		updateHashTag();
		reset = 1;
		draw(false);
	}

	$("xcenInput").onchange = function() {
		xcen = round7(parseFloat($("xcenInput").value));
		updateHashTag();
		reset = 1;
		draw(false);
	}

	$("ycenInput").onchange = function() {
		ycen = round7(parseFloat($("ycenInput").value));
		updateHashTag();
		reset = 1;
		draw(false);
	}

	$("contrastSlider").onchange = function() {
		gradient =  $("contrastSlider").value / 100.0;
		updateHistogram();
		updateHashTag();
		console.log("gradient: " + gradient);
		m_down = false;
		broadcastConfig();
	}

	$("brightnessSlider").onchange = function() {
		brightness =  $("brightnessSlider").value / 100.0;
		updateHistogram();
		updateHashTag();
		console.log("brightness: " + brightness);
		m_down = false;
		broadcastConfig();
	}

	$("primary_light").onchange = function() {
		primary_light = parseFloat($("primary_light").value);
		updateHashTag();
		broadcastConfig();
	}

	setupLightWidget();

	$("fog").onchange = function() {
		fog_factor = parseFloat($("fog").value);
		updateHashTag();
		broadcastConfig();
	}

	$("fogColor").oninput = function() {
		var hex = $("fogColor").value;
		fog_color.r = parseInt(hex.slice(1,3), 16);
		fog_color.g = parseInt(hex.slice(3,5), 16);
		fog_color.b = parseInt(hex.slice(5,7), 16);
		broadcastConfig();
	}

	$("volumetricFog").onchange = function() {
		useVolumetricFog = $("volumetricFog").checked;
		broadcastConfig();
	}

	$("power").onkeyup = function() {
		power = parseFloat($("power").value);
		updateHashTag();
		reset = 1;
	}

	$("DOF").onkeyup = function() {
		cameraDOF = parseFloat($("DOF").value);
		updateHashTag();
		drawFocus = true;
		broadcastConfig();
	}

	$("focus").onkeyup = function() {
		focus = parseFloat($("focus").value);
		updateHashTag();
		drawFocus = true;
		broadcastConfig();
	}

	$("cameraYaw").onchange = function() {
		cameraYaw = parseFloat($("cameraYaw").value);
		updateHashTag();
		reset = 1;
		setCamera();
		redrawLightWidget();
		draw(false);
	}

	$("cameraPitch").onchange = function() {
		cameraPitch = parseFloat($("cameraPitch").value);
		updateHashTag();
		reset = 1;
		setCamera();
		redrawLightWidget();
		draw(false);
	}

	$("colorPalette").onchange = function() {
		pallet = palettes[$("colorPalette").value];
		updateHashTag();
		markGPUPaletteDirty();
		reset = 1;
	}

	$("formulaSelect").onchange = function() {
		formula = parseInt($("formulaSelect").value);
		reset = 1;
		zoom = 2.9;
		draw(false);
	}

	$("iterationsInput").onchange = function() {
		iterations = Math.min(100, parseInt($("iterationsInput").value));
		updateHashTag();
		reset = 1;
	}

	$("inverseAzimuth").onchange = function() {
		azimuth = $("inverseAzimuth").checked ? -1 : 1;
		reset = 1;
	}

	$('canvasControls').onmousedown = function(e)
	{
		m_down = true;
		var rect = ccanvas.getBoundingClientRect();
		xanchor = e.clientX - rect.left;
		yanchor = e.clientY - rect.top;
		console.log(xanchor + " " + yanchor + " - " + this.offsetLeft + " " + this.offsetTop);
	}

	$('canvasControls').oncontextmenu = function(e)
	{
		e.preventDefault();
		var panel = $('description');
		panel.style.display = (panel.style.display === 'none') ? '' : 'none';
	}

	$('canvasControls').onmousemove = function(e)
	{
		if (m_down)
		{
			var c = ccanvas.getContext('2d');
			c.lineWidth = 1;

			// clear out old box first
			c.clearRect(0, 0, ccanvas.width, ccanvas.height);

			var rect = ccanvas.getBoundingClientRect();
			xcurr = e.clientX - rect.left;
			ycurr = e.clientY - rect.top;

			var dx = Math.abs(xcurr - xanchor);
			var dy = Math.abs(ycurr - yanchor);

			// draw new box
			c.strokeStyle = '#FF3B03';
			c.strokeRect(xanchor, yanchor, dx, dy);
			console.log(xanchor + " " + yanchor + " " + dx + " " + dy);
		}
	}

	$('canvasControls').onmouseup = function(e)
	{
		// Ignore mouseup events that didn't originate from a canvas mousedown
		// (e.g. releasing over the canvas after dragging the light widget).
		if (!m_down) return;

		console.log("mouse up!");
		const half_ximlen = ximlen / 2;
		const half_yimlen = yimlen / 2;

		// clear entire canvas
		var c = ccanvas.getContext('2d');
		c.clearRect(0, 0, ccanvas.width, ccanvas.height);

		// do the zoom and restart render
		m_down = false;
		var rect = ccanvas.getBoundingClientRect();
		xcurr = e.clientX - rect.left;
		ycurr = e.clientY - rect.top;

		var dx = Math.abs(xcurr - xanchor);
		var dy = Math.abs(ycurr - yanchor);
		if (dy > dx)  dx = dy;
		
	    // make sure zoom isn't too small
	   	if (dx > 10)
	   	{
			var newxcen = xanchor + dx/2.0;
			newxcen = ((newxcen - half_ximlen)/ximlen)*zoom;
			xcen = round7(xcen - newxcen);

			var newycen = yanchor + dx/2.0;
			newycen = ((newycen - half_yimlen)/yimlen)*zoom;
			ycen = round7(ycen - newycen);

			console.log(dx + " " + dy + " Xcen is " + xcen + " Ycen is " + ycen);
			primary_light = primary_light*1.2;
			$("primary_light").value = primary_light;
			setZoom( (dx/ximlen)*zoom );
			reset = 1;
			setCamera();
			draw(false);
		}
		else if (e.shiftKey)
		{
			console.log("toggle focus planes")
			drawFocus = !drawFocus;
			reset = 2;
		}
	}

	// read hash tag and render away at page load.
	readHashTag();
	setZoom(zoom);
	init();
	updateHashTag();
	redrawLightWidget();

	// Fragment-shader GPU preview first, then kick off the WebGPU progressive
	// renderer once it has finished initialising. If WebGPU is unavailable we
	// fall back to the web-worker CPU renderer on a longer warm-up timer.
	// Firefox's WebGPU implementation is skipped; use the worker path instead.
	var isFirefox = navigator.userAgent.indexOf("Firefox") !== -1;
	renderPreview(function() {
		if (!isFirefox && isWebGPUSupported()) {
			initGPURenderer(canvas, canvas.width, canvas.height).then(function(g) {
				if (g) {
					useGPU = true;
					setTimeout(function() { clearScreenAndReset(); draw(true); }, 4000);
				} else {
					setTimeout(function() { clearScreenAndReset(); draw(true); }, 9000);
				}
			}).catch(function(err) {
				console.warn("WebGPU init failed, falling back to workers:", err);
				setTimeout(function() { clearScreenAndReset(); draw(true); }, 9000);
			});
		} else {
			setTimeout(function() { clearScreenAndReset(); draw(true); }, 9000);
		}
	});
}

main();


/*
 * Parse URL hash tag
 */
function readHashTag()
{
	var tags = location.hash.split('&');
	for ( var i=0; i<tags.length; ++i ) {
		var tag = tags[i].split('=');
		var key = tag[0];
		var val = tag[1];

		switch ( key ) {
			case '#zoom': {
				zoom = parseFloat(val);
				console.log("readHashTag() zoom : " + zoom);
				break;
			} 
			case 'xcen': {
				xcen = round7(parseFloat(val));
				console.log("readHashTag() xcen : " + xcen);
				break;
			}
			case 'ycen': {
				ycen = round7(parseFloat(val));
				console.log("readHashTag() ycen : " + ycen);
				break;
			}
			case 'contrast': {
				gradient = parseFloat(val);
				$("contrastSlider").value = gradient * 100;
				console.log("readHashTag() contrast : " + gradient);
				break;
			}
			case 'brightness': {
				brightness = parseFloat(val);
				$("brightnessSlider").value = brightness * 100;
				console.log("readHashTag() brightness : " + brightness);
				break;
			}
			case 'power': {
				power = parseFloat(val);
				$("power").value = power;
				console.log("readHashTag() power : " + power);
				break;
			}
			case 'fog': {
				fog_factor = parseFloat(val);
				$("fog").value = fog_factor;
				console.log("readHashTag() fog : " + fog_factor);
				break;
			} 
			case 'primary_light': {
				primary_light = parseFloat(val);
				$("primary_light").value = primary_light;
				console.log("readHashTag() primary_light : " + primary_light);
				break;
			}
			case 'light': {
				var parts = val.split(',');
				if (parts.length === 3) {
					var lx = parseFloat(parts[0]);
					var ly = parseFloat(parts[1]);
					var lz = parseFloat(parts[2]);
					if (isFinite(lx) && isFinite(ly) && isFinite(lz)) {
						LightVector[0] = lx;
						LightVector[1] = ly;
						LightVector[2] = lz;
						console.log("readHashTag() light : " + LightVector);
					}
				}
				break;
			}
			case 'dof': {
				cameraDOF = parseFloat(val);
				$("DOF").value = cameraDOF;
				console.log("readHashTag() dof : " + cameraDOF);
				break;
			}
			case 'focus': {
				focus = parseFloat(val);
				$("focus").value = focus;
				console.log("readHashTag() focus : " + focus);
				break;
			}
			case 'yaw': {
				cameraYaw = parseFloat(val);
				$("cameraYaw").value = cameraYaw;
				console.log("readHashTag() yaw : " + cameraYaw);
				break;
			}
			case 'pitch': {
				cameraPitch = parseFloat(val);
				$("cameraPitch").value = cameraPitch;
				console.log("readHashTag() pitch : " + cameraPitch);
				break;
			}
			case 'azimuth': {
				azimuth = parseFloat(val);
				$("inverseAzimuth").checked = (azimuth === -1);
				console.log("readHashTag() azimuth : " + azimuth);
				break;
			}
			case 'formula': {
				formula = parseInt(val);
				$("formulaSelect").value = formula;
				console.log("readHashTag() formula : " + formula);
				break;
			}
			case 'iterations': {
				iterations = parseInt(val);
				$("iterationsInput").value = iterations;
				console.log("readHashTag() iterations : " + iterations);
				break;
			}
			case 'palette': {
				var paletteName = decodeURIComponent(val);
				if (palettes[paletteName]) {
					pallet = palettes[paletteName];
					$("colorPalette").value = paletteName;
					console.log("readHashTag() palette : " + paletteName);
				}
				break;
			}
		}
	}
}

/*
 * Update URL's hash with render parameters so we can pass it around.
 */
var lastSetHash = "";
function updateHashTag()
{
	console.log("updateHashTag(): " + zoom);
	$("zoomInput").value = zoom;
	$("xcenInput").value = xcen;
	$("ycenInput").value = ycen;
	var lightStr = LightVector[0].toFixed(4) + ',' + LightVector[1].toFixed(4) + ',' + LightVector[2].toFixed(4);
	var newHash = 'zoom=' + zoom + '&xcen=' + xcen + '&ycen=' + ycen + '&contrast=' + gradient + '&brightness=' + brightness + "&fog=" +  fog_factor + "&primary_light=" + primary_light + "&power=" + power + "&dof=" + cameraDOF + "&focus=" + focus + "&yaw=" + cameraYaw + "&pitch=" + cameraPitch + "&azimuth=" + azimuth + "&formula=" + formula + "&iterations=" + iterations + "&light=" + lightStr + "&palette=" + encodeURIComponent($("colorPalette").value);
	lastSetHash = '#' + newHash;
	location.hash = newHash;
}

/*
 * Re-render when the user edits the URL hash and presses Enter.
 * Skip the event when it was triggered by our own updateHashTag() call.
 */
window.onhashchange = function()
{
	if (location.hash === lastSetHash) return;
	console.log("hashchange: re-reading parameters from URL");
	readHashTag();
	setZoom(zoom);
	reset = 1;
	setCamera();
	draw(false);
}