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
var frost = 1.0;

// ray traced lighting
var LightVector = [ 0.12, 0.15, -0.19 ];
var AMBIENT_LIGHT = 8.0;
var RAY_STEPS = 25;
var ray_step;
var primary_light = 28.0;
var shadow_darkness = 30.0;
var HORIZON = 20;

// fog based on path traces
var fog_factor = 0.01;
var fog_color = {r: 250, g: 250, b: 250};
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
var cameraPersp = -0.1;
var cameraYaw = 0.1;
var cameraPitch = 0.75;
var cameraDOF = 0.0;
var factorDOF = 0.0;
var focus_depth = 1.0;
var opacity = 1.4;
var focus = -0.15; //higher values for further distance in focus.

// data
var root_zoom;
var half_ximlen, half_yimlen;
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

// worker pool (web workers for parallel scanline rendering)
var workerPool = [];
var workerGen   = 0;   // incremented on each restart; stale results are discarded
var renderY     = 0;   // next row to dispatch

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
	focus_depth = (max_y - min_y) * 0.33;
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
	
/*
function gridpoints(yRow)
{
	//make first pass quick
	var found_limit = Math.ceil(frost);
	if (min_y == -2.0) found_limit = 0;
	
	var trace_history = matrix(iterations+1, 5, 0.0);

	var jitter1 = ((0.5 - Math.random()) / ximlen)*0.4;
	var jitter2 = ((0.5 - Math.random()) / yimlen)*0.4;
	
	//step amount for surface detection
	var stepAmount = (stepDetail + Math.random()*stepDetail)*root_zoom;
	
	var z2 = (yRow / yimlen - 0.5 + jitter1) * zoom - ycen;
	var red, green, blue;
	for (var x = 0; x < ximlen; x++)
	{
		//convert screen to fractal coordinates
		var x2 = (x / ximlen - 0.5 + jitter2) * zoom - xcen;
		
		//reset occlusion for this pixel
		occlusionPositions[x][yRow] = HORIZON;
		
		//main loop for y values (depth)
		for (var y = min_y; y < max_y; y += stepAmount)
		{
			allPixels++;
	
			//rotate coordinate system from screen to fractal coordinates
			var persp = 1.0 + y * cameraPersp;
			var fractal_x = x2 / persp;
			var fractal_z = z2 / persp;
			var point3D = [0,0,0,0,0];
			rotateVector(point3D, CameraMatrix, fractal_x, y, fractal_z);

			//do the calculation
			insideFractal(point3D, trace_history);
			var iter = point3D[3];
			var color = point3D[4];

			//if (allPixels%1000 == 0) console.log("checking " + point3D[0] + " " + point3D[1] + " " + point3D[2] + " " + point3D[3]);

			if (iter == iterations)
			{
				//plot pixel
				var goodPoints = [];
				var rndFuzzy = Math.max(opacity*Math.random(), 0.4);
				var light_factor = 1.0 + calculateRays(point3D, rndFuzzy, goodPoints);
				plotPixel(x, yRow, y, color, light_factor);
				visiblePixels++;
				
				//plot more pixels near the surface
				var found = 0;
				var counter = 0;
				var newy = y;
				while (found < found_limit && counter < 100*frost) {
					newy = y - stepAmount*Math.random()*2.0;
					persp = 1.0 + newy * cameraPersp;
					fractal_x = x2 / persp;
					fractal_z = z2 / persp;
					rotateVector(point3D, CameraMatrix, fractal_x, newy, fractal_z);
					insideFractal(point3D, null);
					iter = Math.floor(point3D[3]);
					color = point3D[4];
					if (iter == iterations) {
						light_factor = 1.0 + calculateRays(point3D, rndFuzzy, goodPoints);
						plotPixel(x, yRow, y, color, light_factor);
						visiblePixels++;
						found++;
					}
					counter++;
				}
				
				//plot additional points found from ray tracing
				if (frost > 0.9)
				{
					for (let i = 0; i < goodPoints.length; i++)
					{
						const traced_point = goodPoints[i]; //double[]
						color = traced_point[4];
						light_factor = 1.0 + calculateRays(traced_point, rndFuzzy, null);
						const screen_point = reversePoint(traced_point);
						
						//if the ray point is not occluded, draw it
						const tempx = Math.round(screen_point[0]);
						const tempy = Math.round(screen_point[2]);
						if ( tempx >= 0 && tempy >= 0 && tempx < ximlen && tempy < yimlen && screen_point[1] < occlusionPositions[tempx][tempy]) {
							plotPixel(tempx, tempy, screen_point[1], color, light_factor);
							rayPoints++;
						}
					}
				}
				break;
			}
			else if (min_y != -2.0) //draw glow after first pass
			{
				//increase step accuracy after first pass and near surface
				var rnd = Math.random();
				stepAmount = (stepDetail + rnd*stepDetail) * ((iterations/iter) / iterations) * root_zoom * 0.5;
				
				const plot_start = 1;
				if (fog_factor > 0 && rnd > 0.9 && iter > plot_start) {
					red = fog_color.r;
					green = fog_color.g;
					blue = fog_color.b;
					for (var c = 0; c < iter; c++) {
						plotFogPixel(trace_history[c], fog_factor, red, green, blue);
					}
				}
			}
		}
	}
}
	
function plotPixel(x, yRow, depth, color, light_factor)
{
	if (cameraDOF > 0)
	{
		var blur_factor = focus-depth;
		
		var blurred_point = [ x, 0, yRow ];
		blur(blurred_point, blur_factor);
		
		var tempx = Math.floor(blurred_point[0]);
		var tempy = Math.floor(blurred_point[2]);
		
		if (tempx >= 0 && tempy >=0 && tempx < ximlen && tempy < yimlen)
		{
			if (depth > occlusionPositions[tempx][tempy] + stepDetail) return;
			plotShadowPixel(tempx, tempy, depth, color, light_factor);
		}
	}
	else
	{
		plotShadowPixel(x, yRow, depth, color, light_factor);
	}
}
	
function blur(orig, blur_factor)
{
	if (blur_factor > 0)
	{
		//increase the range that is in focus
		blur_factor -= focus_depth;
		if (blur_factor < 0) return;
	}
	
	var r2 = Math.random()*2*Math.PI;
	var dsin = Math.sin(r2);
	var dcos = Math.cos(r2);
	
	var dr = Math.random() * factorDOF * blur_factor;
	orig[0] = orig[0] + dr*dcos;
	orig[2] = orig[2] + dr*dsin;
}
	
function insideFractal(data, trace_history)
{
	var magnitude, r, theta_power, r_power, phi, phi_cos, phi_sin;
	var iter = 0;
	var x = data[0];
	var y = data[1];
	var z = data[2];
	
	var pixelColor = 0.0;
	
	do
	{	
		magnitude = x*x + y*y + z*z;
		r = Math.sqrt(magnitude);
		theta_power = Math.atan2(y,x)*power;
		r_power = Math.pow(r,power);
		
		if (formula == 0)
		{
			//2D compatible / sin
			phi = Math.asin(z / r);
			phi_cos = Math.cos(phi*power);
			x = r_power * Math.cos(theta_power) * phi_cos + data[0];
			y = r_power * Math.sin(theta_power) * phi_cos + data[1];
			z = r_power * Math.sin(phi*power)*azimuth + data[2];
			pixelColor = phi / 3.0;
		}
		else if (formula == 1)
		{
			//wikipedia / original / cos
			phi = Math.atan2(Math.sqrt(x*x + y*y), z);
			phi_sin = Math.sin(phi*power);
			x = r_power * Math.cos(theta_power) * phi_sin + data[0];
			y = r_power * Math.sin(theta_power) * phi_sin + data[1];
			z = r_power * Math.cos(phi*power)*azimuth + data[2];
			pixelColor = phi / 3.0;
		}

		//used for nebula / fog
		if (trace_history != null)
		{
			trace_history[iter][0] = x;
			trace_history[iter][1] = y;
			trace_history[iter][2] = z;
			trace_history[iter][3] = iter;
			trace_history[iter][4] = pixelColor;
		}
		
		iter++;
	}
	while ( iter < iterations && r < 8 );

	data[3] = iter;
	data[4] = pixelColor;
}

function calculateRays(origPoint, rndFuzzy, goodPoints)
{
	var light_factor = 1.0;

	var n0 = ray_step / (15.0*Math.random() + 1);
	var n1 = ray_step / (15.0*Math.random() + 1);
	var n2 = ray_step / (15.0*Math.random() + 1);
	
	//ambient light
	light_factor += calculateRay(origPoint, RAY_STEPS, -n0, -n1, -n2, AMBIENT_LIGHT, rndFuzzy, goodPoints);
	light_factor += calculateRay(origPoint, RAY_STEPS, n0, n1, n2, AMBIENT_LIGHT, rndFuzzy, goodPoints);
	light_factor += calculateRay(origPoint, RAY_STEPS, n0, -n1, -ray_step, AMBIENT_LIGHT, rndFuzzy, goodPoints);
	light_factor += calculateRay(origPoint, RAY_STEPS, n0*LightVector[0]*9.0,  n1*LightVector[1]*9.0, n2*LightVector[2]*9.0, AMBIENT_LIGHT, rndFuzzy, goodPoints);

	//direct light
	light_factor += calculateRay(origPoint, RAY_STEPS*4, ray_step*LightVector[0],  ray_step*LightVector[1], ray_step*LightVector[2], primary_light, rndFuzzy, goodPoints);
	
	return light_factor;
}
	

function calculateRay(origPoint, steps,stepx,stepy, stepz, bright, rndFuzzy, goodPoints)
{
	stepx *= rndFuzzy;
	stepy *= rndFuzzy;
	stepz *= rndFuzzy;

	var x = origPoint[0];
	var y = origPoint[1];
	var z = origPoint[2];
	
	for (var i = 1; i < steps; i++)
	{
		origPoint[0] += stepx*i;
		origPoint[1] += stepy*i;
		origPoint[2] += stepz*i;
		
		insideFractal( origPoint, null );
		if ( origPoint[3] == iterations )
		{
			//ray hit solid, this is in shadow
			bright = 0.0;
			if (goodPoints != null)
			{
				//save points that are moderately close to the known surface

				//double[] clone_data = new double[5];
				//System.arraycopy(origPoint, 0, clone_data, 0, clone_data.length);

				var clone_data = [];
				for (var k = 0; k < origPoint.length; k++) {
					clone_data[k] = origPoint[k];
				}

				goodPoints.push(clone_data);
			}
			break;
		}
		else if (origPoint[3] < 3)
		{
			//ray has left the general area of the solid, stop tracing.
			break;
		}
	}
	origPoint[0] = x;
	origPoint[1] = y;
	origPoint[2] = z;
	
	return bright;
}
	
function plotShadowPixel(tempx, tempy, depth, colorVal, light_factor)
{
	//get color
	var colorIndex = Math.floor(Math.abs(colorVal)*255.0);
	if (colorIndex > 255) colorIndex = 255;
	
	//apply lighting
	var red = (pallet[colorIndex][0]*light_factor) / shadow_darkness;
	var green = (pallet[colorIndex][1]*light_factor) / shadow_darkness;
	var blue = (pallet[colorIndex][2]*light_factor) / shadow_darkness;
	
	img_red[tempx][tempy] += red;
	img_green[tempx][tempy] += green;
	img_blue[tempx][tempy] += blue;
	img_alpha[tempx][tempy] += 1.0;
	
	//record depth to mask fog
	occlusionPositions[tempx][tempy] = depth;
}
	
function plotFogPixel(origPoint, factor, r, g, b)
{
	var screenPoint = reversePoint(origPoint);
	var tempx = Math.round(screenPoint[0]);
	var tempy = Math.round(screenPoint[2]);
	var depth = screenPoint[1];
	if (tempx >= 0 && tempy >=0 && tempx < ximlen && tempy < yimlen) {
		//solid occludes glow
		if (depth > occlusionPositions[tempx][tempy]) return;

		if (cameraDOF > 0) {
			var blur_factor = focus-depth;
			blur(screenPoint, blur_factor);
			tempx = Math.round(screenPoint[0]);
			tempy = Math.round(screenPoint[2]);
			if (tempx < 0 || tempy < 0 || tempx > ximlen -1 || tempy > yimlen-1) return;
		}
		img_red[tempx][tempy] += (r*factor);
		img_green[tempx][tempy] += (g*factor);
		img_blue[tempx][tempy] += (b*factor);
		img_alpha[tempx][tempy] += factor;
	}
}
	
function reversePoint(fracPoint)
{
	var point3D = [];
	point3D[3] = fracPoint[3];
	point3D[4] = fracPoint[4];
	
	rotateVector(point3D, IrotZ, fracPoint[0], fracPoint[1], fracPoint[2]);
	var xx = point3D[0];
	var yy = point3D[1];
	var zz = point3D[2];
	rotateVector(point3D, IrotX, xx, yy, zz);
	
	var persp = 1.0 + point3D[1] * cameraPersp;
	point3D[0] = point3D[0]*persp;
	point3D[2] = point3D[2]*persp;

	point3D[0] = ((point3D[0] + xcen)/zoom)*ximlen + half_ximlen ;
	point3D[2] = ((point3D[2] + ycen)/zoom)*yimlen + half_yimlen ;
	
	return point3D;
}
*/
	
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
	
function setZoom(z)
{
	zoom = z;
	root_zoom = Math.pow(zoom, 0.5);
	ray_step = rayDetail*zoom;
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

function render(startScanning)
{
	if (!startScanning) return;  // workers handle reset signals themselves

	var start      = (new Date).getTime();
	var lastUpdate = start;
	var pixels     = 0;

	// Terminate any previously running workers
	workerPool.forEach(function(w) { w.terminate(); });
	workerPool = [];
	workerGen++;
	renderY = 0;

	var gen = workerGen;  // captured for stale-result detection

	function generateConfig() {
		return {
			ximlen: ximlen, yimlen: yimlen,
			zoom: zoom, xcen: xcen, ycen: ycen,
			half_ximlen: half_ximlen, half_yimlen: half_yimlen,
			cameraPersp: cameraPersp,
			iterations: iterations, formula: formula, azimuth: azimuth, power: power,
			stepDetail: stepDetail, frost: frost,
			root_zoom: root_zoom, opacity: opacity,
			focus: focus, focus_depth: focus_depth,
			cameraDOF: cameraDOF, factorDOF: factorDOF,
			LightVector: LightVector, RAY_STEPS: RAY_STEPS,
			AMBIENT_LIGHT: AMBIENT_LIGHT, primary_light: primary_light,
			shadow_darkness: shadow_darkness, HORIZON: HORIZON,
			ray_step: ray_step,
			CameraMatrix: CameraMatrix, IrotX: IrotX, IrotZ: IrotZ,
			fog_factor: fog_factor, fog_color: fog_color,
			pallet: pallet
		};
	}

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
			$('renderTime').innerHTML = elapsed.toFixed(1) + " pass: " + renderpass;
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
				clearScreenAndReset();
				render(true);
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

function init()
{
	yimlen = canvas.height;
	ximlen = canvas.height;

	half_ximlen = ximlen / 2;
	half_yimlen = yimlen / 2;

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
	$('viewPNG').onclick = function(_e)
	{
		var link = document.createElement('a');
		link.download = 'mandelbulb-' + renderpass + '-' + max_alpha + '.png';
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
		xcen = parseFloat($("xcenInput").value);
		updateHashTag();
		reset = 1;
		draw(false);
	}

	$("ycenInput").onchange = function() {
		ycen = parseFloat($("ycenInput").value);
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
	}

	$("brightnessSlider").onchange = function() {
		brightness =  $("brightnessSlider").value / 100.0; 
		updateHistogram();
		updateHashTag();
		console.log("brightness: " + brightness);
		m_down = false;
	}

	$("primary_light").onchange = function() {
		primary_light = parseFloat($("primary_light").value);
		updateHashTag();
	}

	$("fog").onchange = function() {
		fog_factor = parseFloat($("fog").value);
		updateHashTag();
	}

	$("fogColor").oninput = function() {
		var hex = $("fogColor").value;
		fog_color.r = parseInt(hex.slice(1,3), 16);
		fog_color.g = parseInt(hex.slice(3,5), 16);
		fog_color.b = parseInt(hex.slice(5,7), 16);
	}

	$("power").onchange = function() {
		power = parseFloat($("power").value);
		updateHashTag();
		reset = 1;
	}

	$("DOF").onchange = function() {
		cameraDOF = parseFloat($("DOF").value);
		factorDOF = cameraDOF * (ximlen/3);
		updateHashTag();
		drawFocus = true;
	}

	$("focus").onchange = function() {
		focus = parseFloat($("focus").value);
		updateHashTag();
		drawFocus = true;
	}

	$("cameraYaw").onchange = function() {
		cameraYaw = parseFloat($("cameraYaw").value);
		updateHashTag();
		reset = 1;
		setCamera();
		draw(false);
	}

	$("cameraPitch").onchange = function() {
		cameraPitch = parseFloat($("cameraPitch").value);
		updateHashTag();
		reset = 1;
		setCamera();
		draw(false);
	}

	$("colorPalette").onchange = function() {
		pallet = palettes[$("colorPalette").value];
		reset = 1;
	}

	$("formulaSelect").onchange = function() {
		formula = parseInt($("formulaSelect").value);
		reset = 1;
		zoom = 2.9;
		draw(false);
	}

	$("iterationsInput").onchange = function() {
		iterations = parseInt($("iterationsInput").value);
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
		console.log("mouse up!");

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
			xcen = xcen - newxcen;

			var newycen = yanchor + dx/2.0;
			newycen = ((newycen - half_yimlen)/yimlen)*zoom;
			ycen = ycen - newycen;

			console.log(dx + " " + dy + " Xcen is " + xcen + " Ycen is " + ycen);
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
	draw(true);
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
				xcen = parseFloat(val);
				console.log("readHashTag() xcen : " + xcen);
				break;
			} 
			case 'ycen': {
				ycen = parseFloat(val);
				console.log("readHashTag() ycen : " + ycen);
				reDraw = true;
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
			case 'dof': {
				cameraDOF = parseFloat(val);
				factorDOF = cameraDOF * (ximlen/3);
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
		}
	}
}


/*
 * Update URL's hash with render parameters so we can pass it around.
 */
function updateHashTag()
{
	console.log("updateHashTag(): " + zoom);
	$("zoomInput").value = zoom;
	$("xcenInput").value = xcen;
	$("ycenInput").value = ycen;
	location.hash = 'zoom=' + zoom + '&xcen=' + xcen + '&ycen=' + ycen + '&contrast=' + gradient + '&brightness=' + brightness + "&fog=" +  fog_factor + "&primary_light=" + primary_light + "&power=" + power + "&dof=" + cameraDOF + "&focus=" + focus + "&yaw=" + cameraYaw + "&pitch=" + cameraPitch + "&azimuth=" + azimuth + "&formula=" + formula;
}

// file:///C:/devel/mandelbrot-js/mandelbulb.html#zoom=2.9&xcen=1.1019999999999999&ycen=-0.17883333333333334&contrast=0.5&brightness=1.8&fog=0.01&primary_light=38&power=2&dof=0&focus=-0.15&yaw=-0.8&pitch=0.1&azimuth=-1&formula=1
// file:///C:/devel/mandelbrot-js/mandelbulb.html#zoom=0.6673544444444444&xcen=-0.02137944444444445&ycen=0.9455127777777776&contrast=1.08&brightness=2.62&fog=0.51&primary_light=28&power=3&dof=0&focus=-0.15&yaw=0.1&pitch=0.75&azimuth=-1&formula=1
