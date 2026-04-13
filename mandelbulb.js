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

function matrix(rows, cols, defaultValue) {
	var arr = [];
	for(var i = 0; i < rows; i++) {
		arr.push([]);
		arr[i].push(new Array(cols));

		for(var j=0; j < cols; j++) {
			arr[i][j] = defaultValue;
		}
	}
	return arr;
}
	
function updateMinMaxY()
{
	console.log("updateMinMaxY()");
	min_y = HORIZON;
	max_y = -HORIZON;
	for (var i=0; i<ximlen; i++)
	{
		for (var j=0; j<yimlen; j++)
		{
			if (occlusionPositions[i][j] < min_y)
			{
				min_y = occlusionPositions[i][j]-stepDetail;
			}
			else if (occlusionPositions[i][j] != HORIZON && occlusionPositions[i][j] > max_y)
			{
				max_y = occlusionPositions[i][j]+stepDetail;
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
	
	occlusionPositions = matrix(canvas.width,canvas.height,0.0);
	img_alpha = matrix(canvas.width, canvas.height, 0.0);
	img_red = matrix(canvas.width, canvas.height, 0.0);
	img_green = matrix(canvas.width, canvas.height, 0.0);
	img_blue = matrix(canvas.width, canvas.height, 0.0);

	m_down = false;
}
	
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
	
function rotateVector(point3D, rot, x, y, z)
{
	point3D[0] = rot[0][0] * x + rot[0][1] * y + rot[0][2] * z;
	point3D[1] = rot[1][0] * x + rot[1][1] * y + rot[1][2] * z;
	point3D[2] = rot[2][0] * x + rot[2][1] * y + rot[2][2] * z;
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
			if (goodPoints != null /*&& i > 1*/)
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
	
function RotateX(angle)
{
	var rot = matrix(3,3,0.0);
	var s = Math.sin(angle);
	var c = Math.cos(angle);
	rot[0][0] = 1.0;
	rot[1][1] = c;
	rot[2][2] = c;
	rot[1][2] = -s;
	rot[2][1] = s;
	return rot;
}

function RotateY(angle)
{
	var rot = matrix(3,3,0.0);
	var s = Math.sin(angle);
	var c = Math.cos(angle);
	rot[1][1] = 1.0;
	rot[2][2] = c;
	rot[0][0] = c;
	rot[2][0] = -s;
	rot[0][2] = s;
	return rot;
}
	
function RotateZ(angle)
{
	var rot = matrix(3,3,0.0);
	var s = Math.sin(angle);
	var c = Math.cos(angle);
	rot[2][2] = 1.0;
	rot[0][0] = c;
	rot[1][1] = c;
	rot[0][1] = -s;
	rot[1][0] = s;
	return rot;
}
	
function matrixMult(m, matrixArr)
{
	var result = matrix(3,3,0.0);
	result[0][0] = m[0][0] * matrixArr[0][0] + m[0][1] * matrixArr[1][0] + m[0][2] * matrixArr[2][0];
	result[0][1] = m[0][0] * matrixArr[0][1] + m[0][1] * matrixArr[1][1] + m[0][2] * matrixArr[2][1];
	result[0][2] = m[0][0] * matrixArr[0][2] + m[0][1] * matrixArr[1][2] + m[0][2] * matrixArr[2][2];
	result[1][0] = m[1][0] * matrixArr[0][0] + m[1][1] * matrixArr[1][0] + m[1][2] * matrixArr[2][0];
	result[1][1] = m[1][0] * matrixArr[0][1] + m[1][1] * matrixArr[1][1] + m[1][2] * matrixArr[2][1];
	result[1][2] = m[1][0] * matrixArr[0][2] + m[1][1] * matrixArr[1][2] + m[1][2] * matrixArr[2][2];
	result[2][0] = m[2][0] * matrixArr[0][0] + m[2][1] * matrixArr[1][0] + m[2][2] * matrixArr[2][0];
	result[2][1] = m[2][0] * matrixArr[0][1] + m[2][1] * matrixArr[1][1] + m[2][2] * matrixArr[2][1];
	result[2][2] = m[2][0] * matrixArr[0][2] + m[2][1] * matrixArr[1][2] + m[2][2] * matrixArr[2][2];
	return result;
}

/* Matrix Ops */
function determinant(mat) 
{
	var result = 0;
	
	if(mat.length == 1) {
		result = mat[0][0];
		return result;
	}

	if(mat.length == 2) {
		result = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
		return result;
	} 

	for(var i = 0; i < mat[0].length; i++) 
	{
		var temp = matrix(mat.length - 1, mat[0].length - 1, 0.0);
		for(var j = 1; j < mat.length; j++)
		{
			for(var k = 0; k < mat[0].length; k++) 
			{
				if(k < i) {
				temp[j - 1][k] = mat[j][k];
				} else if(k > i) {
				temp[j - 1][k - 1] = mat[j][k];
				}
			}
		}

		result += mat[0][i] * Math.pow(-1, i) * determinant(temp);
	}
	return result;
} 

function cofactor3x3T(m) 
{
	var temp = matrix(m[0].length, m[0].length, 0.0);
	temp[0][0] = m[1][1]*m[2][2] - m[1][2]*m[2][1];
	temp[1][0] = m[0][2]*m[2][1] - m[0][1]*m[2][2];
	temp[2][0] = m[0][1]*m[1][2] - m[0][2]*m[1][1];
	temp[0][1] = m[1][2]*m[2][0] - m[1][0]*m[2][2];
	temp[1][1] = m[0][0]*m[2][2] - m[0][2]*m[2][0];	
	temp[2][1] = m[0][2]*m[1][0] - m[0][0]*m[1][2];
	temp[0][2] = m[1][0]*m[2][1] - m[1][1]*m[2][0];
	temp[1][2] = m[0][1]*m[2][0] - m[0][0]*m[2][1];
	temp[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0];
	return temp;
} 

function inverse(m)
{
	var d = determinant(m);
	if (d == 0) d = 0.00001;  //dubious
	
	var cofactorT = cofactor3x3T(m);
	d = 1/d;
	constMul(d, cofactorT);
	return cofactorT;
}
	
function constMul(d, m)
{
	m[0][0] *= d;
	m[0][1] *= d;
	m[0][2] *= d;
	m[1][0] *= d;
	m[1][1] *= d;
	m[1][2] *= d;
	m[2][0] *= d;
	m[2][1] *= d;
	m[2][2] *= d;
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
	var start  = (new Date).getTime();
	var lastUpdate = start;
	var pixels = 0;

	var y = 0;
	var scanline = function()
	{
		if (reset == 1)
		{
			console.log("reset render at first line");
			reset = 0;
			y = 0;
			clearScreenAndReset();
		}
		else if (reset == 2)
		{
			reset = 0;
			updateHistogram();
		}
		
		renderpass++;
		
		gridpoints(y);
		if (renderpass % ximlen == 0) updateMinMaxY();
		pixels += ximlen;
		
		if (renderpass % 50 == 0)
		{
			const t2 = (new Date()).getTime();
			const completeness = Math.round(max_alpha*100)/100.0;
			const strStatus = (new Date()).toISOString().substring(0,19) + " " + y + " " + visiblePixels + " " + rayPoints + " " + allPixels + " Pass: " + renderpass + " max value: " + completeness + " in " + (t2-t1);
			t1 = (new Date()).getTime();
			console.log(strStatus);
			updateHistogram();

			//autosave png and last 5 states.
			//if (lastDir != null && lastDir.length() > 0 && (currentPass % 5000 == 0) )
			//{
			//	writePNG(lastDir + "Mandelbulb-" + currentPass + "-" + completeness + "-" + ximlen + "x" + yimlen + ".png");
			//	writeFractalData(lastDir + "Mandelbulb-" + (currentPass%30000) + "-" + ximlen + "x" + yimlen + ".fractal");
			//}
			
			visiblePixels = 0;
			allPixels = 0;
			rayPoints = 0;
		}

		var now = (new Date).getTime();
		if ( (now - lastUpdate) >= 10000) {
			// Update speed and time taken
			var elapsed = (now - start)/1000.0;
			var speed = Math.floor(pixels / elapsed);

			$('renderTime').innerHTML = elapsed.toFixed(1) + " pass: " + renderpass;
			$('renderSpeed').innerHTML = speed + " px/sec";

			lastUpdate = now;
        }

		// yield control back to browser, so that canvas is updated
		const sleepTime = m_down ? 2000 : 1;
		y++;
		if (y > yimlen-1) y = 0;
        setTimeout(scanline, sleepTime);
	}

	if (startScanning) scanline();
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

			//tool to view focus bouds
			if (drawFocus && occlusionPositions[x][y] > focus) {
				green = 20;
				red = 20;
				blue += 50;
			}
			else if (drawFocus && occlusionPositions[x][y] < focus - focus_depth) {
				green = 20;
				blue = 20;
				red += 50;
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
