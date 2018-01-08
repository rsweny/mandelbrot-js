/*
 * The Mandelbulb, in HTML5 canvas and javascript.
 * https://github.com/rsweny/mandelbrot-js
 *
 * Copyright (C) 2018 Ryan Swney
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

 var mode = 0;


//detail level
var rayDetail = 0.003;
var stepDetail = 0.028;
var frost = 1.0;

//ray traced lighting
var LightVector = [ 0.12, 0.15, -0.19 ];
var AMBIENT_LIGHT = 8.0;
var RAY_STEPS = 25;
var ray_step;
var primary_light = 28.0;
var shadow_darkness = 30.0;
var HORIZON = 20;

//fog based on path traces
var fog_factor = 0.04;
var min_y, max_y;

//formula variation
var formula = 0;
var inverse_azimuth = 1;

//settings that need saving
var depth = 20;
var ximlen = 0;
var yimlen = 0;
var pal = 1;
var power = 8;
var gradient = 0.6;
var brightness = 1.5;
var zoom = 2.4;
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

//data
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


//mouse
var xcurr,ycurr,xanchor,yanchor;
var m_down;
var drawFocus = false;


//Stats
var t1 = 0;
var visiblePixels, allPixels, rayPoints;
var renderpass = 0;
var max_alpha = 1;

//Screen
var zoomStart = 2.4;
var screenZoom = [zoomStart, zoomStart];
var lookAtDefault = [0.0, 0.0];
var lookAt = lookAtDefault;
var xRange = [0, 0];
var yRange = [0, 0];
var range_x;
var range_y;

var reInitCanvas = true; // Whether to reload canvas size, etc
var dragToZoom = true;

var renderId = 0; // To zoom before current render is finished

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
var ctx_img = ctx.createImageData(canvas.width, canvas.height);

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

function matrix(rows, cols, defaultValue)
{
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
	focus_depth = (max_y - min_y) * 0.25; //0.3
	console.log(focus_depth + " Y Bounds: " + min_y + " to " + max_y);
}
	
/*
public void paint(Graphics g)
{
	g.drawImage(img, 0, 0, this);
	if (m_down) 
    {
		g.setColor(Color.white);
		g.drawRect(xanchor,yanchor,xcurr-xanchor,ycurr-yanchor);
    }
}
*/
	
function findPeak(arr)
{
	var max = 0;
	for (var i = 1; i < canvas.width; i++)
	{
		for (var j = 1; j < canvas.height; j++)
		{
			if (arr[i][j] > max)
			{
				max = arr[i][j];
			}
		}
	}
	return max;
}

function clearScreenAndReset(fresh)
{
	renderpass = 0;
	max_alpha = 1;
	min_y = -2.0;
	max_y = 2.0;
	
	if (fresh)
	{
		occlusionPositions = matrix(ximlen,yimlen,0.0);
		img_alpha =  matrix(ximlen,yimlen,0.0);
		img_red = matrix(ximlen,yimlen,0.0);
		img_green = matrix(ximlen,yimlen,0.0);
		img_blue = matrix(ximlen,yimlen,0.0);
	}

	m_down = false;
}
	
function gridpoints(yRow)
{
	//make first pass quick
	var found_limit = Math.ceil(frost);
	if (min_y == -2.0) found_limit = 0;
	
	var trace_history = matrix(depth+1, 5, 0.0);

	var jitter1 = (0.5 - Math.random()) / ximlen;
	var jitter2 = (0.5 - Math.random()) / yimlen;
	
	//step amount for surface detection
	var stepAmount = (stepDetail + Math.random()*stepDetail)*root_zoom;
	
	var z2 = (yRow / yimlen - 0.5 + jitter1) * zoom - ycen;
	var red, green, blue;

	//console.log(root_zoom + " " + stepAmount + " lims " + ximlen + " " + yimlen);
	
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
			

			if (iter == depth)
			{
				//console.log("inside for " + fractal_x + " " + y + " " + fractal_z);
				//plot pixel
				var goodPoints = []; //java Vector
				var rnd1 = Math.random();
				var light_factor = 1.0 + calculateRays(point3D, rnd1, goodPoints);
				plotPixel(x, yRow, y, color, light_factor);
				visiblePixels++;
				
				//plot more pixels near the surface
				var found = 0;
				var counter = 0;
				var newy = y;
				while (found < found_limit && counter < 200*frost)
				{
					newy = y - stepAmount*Math.random()*2.0;
					persp = 1.0 + newy * cameraPersp;
					fractal_x = x2 / persp;
					fractal_z = z2 / persp;
					rotateVector(point3D, CameraMatrix, fractal_x, newy, fractal_z);
					insideFractal(point3D, null);
					iter = Math.floor(point3D[3]);
					color = point3D[4];
					if (iter == depth)
					{
						light_factor = 1.0 + calculateRays(point3D, rnd1, goodPoints);
						plotPixel(x, yRow, y, color, light_factor);
						visiblePixels++;
						found++;
					}
					counter++;
				}
				
				//plot additional points found from ray tracing
				if (frost > 0.9)
				{
					for (var i = 0; i < goodPoints.length; i++)
					{
						var traced_point = goodPoints[i]; //double[]
						color = traced_point[4];
						light_factor = 1.0 + calculateRays(traced_point, rnd1, null);
						var screen_point = reversePoint(traced_point);
						
						//if the ray point is not occluded, draw it
						var tempx = Math.floor(screen_point[0]);
						var tempy = Math.floor(screen_point[2]);
						
						if ( tempx >= 0 && tempy >= 0 && tempx < ximlen && tempy < yimlen && screen_point[1] < occlusionPositions[tempx][tempy])
						{
							//plot the ray pixel
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
				stepAmount = (stepDetail + rnd*stepDetail) * ((depth/iter) / depth) * root_zoom * 0.5;
				
				var plot_start = 1;
				if (cameraPersp > 0) plot_start = 4;
				if (fog_factor > 0 && rnd > 0.9 && iter > plot_start)
				{
					blue = 250;
					green = 90 + iter*2;
					red = 20 + iter;
					
					for (var c = 0; c < iter; c++)
					{
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
			z = r_power * Math.sin(phi*power)*inverse_azimuth + data[2];
			pixelColor = phi / 3.0;
		}
		else if (formula == 1)
		{
			//wikipedia / original / cos
			phi = Math.atan2(Math.sqrt(x*x + y*y), z);
			phi_sin = Math.sin(phi*power);
			x = r_power * Math.cos(theta_power) * phi_sin + data[0];
			y = r_power * Math.sin(theta_power) * phi_sin + data[1];
			z = r_power * Math.cos(phi*power)*inverse_azimuth + data[2];
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
	while ( iter < depth && r < 8 );

	data[3] = iter;
	data[4] = pixelColor;
}
	

function calculateRays(origPoint, rnd, goodPoints)
{
	var light_factor = 1.0;

	if (rnd > 0.9) rnd *= 2.0;
	var rndFuzzy = Math.max(opacity*rnd, 0.4);
	
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
		if ( origPoint[3] == depth )
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
	var screenPoint = reversePoint( origPoint );
	var tempx = Math.floor(screenPoint[0]);
	var tempy =  Math.floor(screenPoint[2]);
	var depth = screenPoint[1];
	

	if (tempx >= 0 && tempy >=0 && tempx < ximlen && tempy < yimlen)
	{
		//solid occludes glow
		if (depth > occlusionPositions[tempx][tempy]) return;
	
		if (cameraDOF > 0)
		{
			var blur_factor = focus-depth;
			blur(screenPoint, blur_factor);
			tempx = Math.floor(screenPoint[0]);
			tempy = Math.floor(screenPoint[2]);
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
	

function getConstants()
{
		/*
		power = Double.parseDouble(txtPower.getText());
		depth = Integer.parseInt(txtDepth.getText());
		brightness = Double.parseDouble(txtBrightness.getText());
		gradient = Double.parseDouble(txtGradient.getText());
		setZoom( Double.parseDouble(txtZoom.getText()) );

		cameraPersp = Double.parseDouble(txtPerspective.getText());
		stepDetail = Double.parseDouble(txtZPos.getText());
		cameraPitch = Double.parseDouble(txtPitch.getText());
		cameraYaw = Double.parseDouble(txtYaw.getText());
		cameraDOF = Double.parseDouble(txtDOF.getText());
		opacity = Double.parseDouble(txtOpacity.getText());
		focus = Double.parseDouble(txtFocus.getText());
		frost = Double.parseDouble(txtZRes.getText());
		fog_factor = Double.parseDouble(txtFog.getText());
		primary_light = Double.parseDouble(txtLight.getText());
		*/
		
		setCamera();

}
	
function setConstants()
{
		/*
		txtPower.setText(""+ power);
		txtDepth.setText(""+ depth);
		txtBrightness.setText(""+ brightness);
		txtGradient.setText(""+ gradient);
		txtZoom.setText(""+ zoom);
		
		txtPerspective.setText("" + cameraPersp);
		txtZPos.setText("" + stepDetail);
		txtPitch.setText("" + cameraPitch);
		txtYaw.setText("" + cameraYaw);
		txtDOF.setText("" + cameraDOF);
		txtOpacity.setText("" + opacity);
		txtFocus.setText("" + focus);
		txtZRes.setText("" + frost);
		txtFog.setText("" + fog_factor);
		txtLight.setText("" + primary_light);
		*/
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
	
	factorDOF = cameraDOF * (ximlen/3);

	console.log(CameraMatrix[0][0] + " ----------------------------- " + CameraMatrix[0][1] + " " + CameraMatrix[0][2] + " " + CameraMatrix[1][0] + " " + CameraMatrix[1][1] + " " + CameraMatrix[1][2]);
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

	//ray tracing
	ray_step = rayDetail*zoom;
}
	
/*
public boolean mouseDown(Event evt, int x, int y)
{
	xanchor = x;
	yanchor = y;

	if (evt.modifiers == Event.ALT_MASK)
	{
		formula  = (formula + 1)%2;
		if (formula == 0) System.out.println("sin mandelbulb");
		else if (formula == 1) System.out.println("cos mandelbulb");
		else System.out.println("symetric mandelbulb");
		reset = 1;
	}
	else if (evt.modifiers == Event.SHIFT_MASK)
	{
		pal  = (pal + 1)%(Pallet.fpalette.length);
		reset = 1;
	}
	else if (evt.modifiers == Event.CTRL_MASK)
	{
		inverse_azimuth = -inverse_azimuth;
		System.out.println("inverse azimuth is " + inverse_azimuth);
		reset = 1;
	}
	else if (evt.modifiers == Event.META_MASK)
	{
		panelMain.setVisible(!panelMain.isVisible());
	}

	return true;
}

	public boolean mouseDrag(Event evt, int x, int y)
	{
		m_down = true;
		xcurr = x;
		ycurr = y;
		repaint();
		
		int dx = Math.abs(xcurr - xanchor);
		int dy = Math.abs(ycurr - yanchor);
		double newxcen = xanchor + dx/2.0;
		double fx = ((newxcen - (double)half_ximlen)/(double)ximlen)*zoom;
		
		double newycen = yanchor + dy/2.0;
		double fy = ((newycen - (double)half_yimlen)/(double)yimlen)*zoom;
		
		String strStatus = newxcen + " " + newycen + "( " + fx + ","+fy+")";
		showStatus(strStatus);
		
		return true;
	}

	public boolean mouseUp(Event evt, int x, int y)
	{
		xcurr = x;
		ycurr = y;
		m_down = false;

		int dx = Math.abs(xcurr - xanchor);
		int dy = Math.abs(ycurr - yanchor);
		if (dy > dx)  dx = dy;
		
        //make sure zoom isn't too small
       	if (dx > 10)
       	{
			double newxcen = xanchor + dx/2.0;
			newxcen = ((newxcen - half_ximlen)/(double)ximlen)*zoom;
			xcen = xcen - newxcen;

			double newycen = yanchor + dx/2.0;
			newycen = ((newycen - half_yimlen)/(double)yimlen)*zoom;
			ycen = ycen - newycen;
			
			System.out.println("Xcen is " + xcen + " Ycen is " + ycen);
			setZoom( ((double)dx/(double)(ximlen))*zoom);
			reset = 1;
		}

		setConstants();
		getConstants();
		return true;
	}

	
	public void keyReleased(KeyEvent e)
	{
		char c = e.getKeyChar();
		int keyCode = (int)c;
		
		//22: ctrl-v
		//8: backspace
		if (keyCode == 22 || keyCode == 8 || c == ' ' || c == '-' || c == '|' || c == '.' || (c >= '0' && c <= '9') )
		{
		
			if (e.getSource() != txtGradient && e.getSource() != txtBrightness && e.getSource() != txtFocus && e.getSource() != txtDOF && e.getSource() != txtFog && e.getSource() != txtLight && e.getSource() != txtZRes && e.getSource() != txtOpacity && e.getSource() != txtZPos) 
				reset = 1;
			else
				reset = 2;
	
			getConstants();
		}
	}
	
	public void keyTyped(KeyEvent e)
	{
		char c = e.getKeyChar();
		int keyCode = (int)c;
		
		//draw focus border for DOF
		if ( c == 'f')
		{
			drawFocus = !drawFocus;
			reset = 2;
			e.consume();
		}
	}
	
	public void actionPerformed(ActionEvent e)
	{
		m_down = true;
		try
		{
			if (e.getSource() == saveFractal)
			{
				String filename = saveFile("Save Fractal", null, "mandelbrot3d_" + ximlen + "x" + yimlen + ".fractal");
				writeFractalData(filename);
			}
			else if (e.getSource() == openFractal)
			{
				String filename = loadFile("Open Fractal", null, ".fractal");
				readFractalData(filename);
			}
			else
			{
				String filename = saveFile("Export PNG", null, ".png");
				writePNG(filename);
			}
		}
		catch(Exception ex)
		{
		}
		m_down = false;
	}
	
	
	public String loadFile(String title, String defDir, String fileType) throws Exception
	{
		Frame parent = new Frame();
		FileDialog fd = new FileDialog(parent, title, FileDialog.LOAD);
		fd.setFile(fileType);
		fd.setDirectory(defDir);
		fd.setLocation(50, 50);
		fd.show();
		lastDir = fd.getDirectory();
		if (lastDir == null) throw new Exception();
		return lastDir + fd.getFile();
	}

	public String saveFile(String title, String defDir, String fileType) throws Exception
	{
		Frame parent = new Frame();
		FileDialog fd = new FileDialog(parent, title, FileDialog.SAVE);
		fd.setFile(fileType);
		fd.setDirectory(defDir);
		fd.setLocation(50, 50);
		fd.show();
		lastDir = fd.getDirectory();
		if (lastDir == null) throw new Exception();
		return lastDir + fd.getFile();
	} 
	
	public void writePNG(String filename)
	{
		try 
		{
			BufferedImage bi = new BufferedImage(ximlen, yimlen, BufferedImage.TYPE_INT_ARGB); 
			bi.setRGB(0, 0, ximlen, yimlen, pixels, 0, ximlen);
			File outputfile = new File(filename);
			ImageIO.write(bi, "png", outputfile);
			System.out.println("Export PNG success: " + outputfile.getAbsolutePath());
			showStatus("Export PNG success: " + outputfile.getAbsolutePath() );
		} 
		catch (Exception e) 
		{
			System.out.println("export: " + e.toString());
		}
	}
	
	public void writeFractalData(String filename)
	{
		try 
		{
			FileOutputStream fos = new FileOutputStream(filename);
			ObjectOutputStream oos = new ObjectOutputStream(fos);

			oos.writeObject(img_alpha);
			oos.writeObject(img_red);
			oos.writeObject(img_green);
			oos.writeObject(img_blue);
			
			//write int pref array
			int[] intPrefs = { depth, ximlen, yimlen, pal, inverse_azimuth, formula };
			oos.writeObject(intPrefs);
			
			//write double pref array
			double[] doublePrefs = { power, gradient, brightness, zoom, xcen, ycen, cameraPersp, stepDetail, cameraYaw, cameraPitch, cameraDOF, opacity, focus, frost, fog_factor, primary_light };
			oos.writeObject(doublePrefs);
			
			oos.flush();
			fos.close();
			System.out.println("write Fractal Data success: " + filename);
			showStatus("write Fractal Data success: " + filename);
		}
		catch (Throwable e) 
		{
			System.out.println("write: " + e.toString());
		} 
	}
	
	public void readFractalData(String filename)
	{
		try 
		{
			FileInputStream fis = new FileInputStream(filename);
			ObjectInputStream ois = new ObjectInputStream(fis);

			img_alpha = (float[][])ois.readObject();
			img_red = (float[][])ois.readObject();
			img_green = (float[][])ois.readObject();
			img_blue = (float[][])ois.readObject();

			//read int pref array
			int[] intPrefs = (int[])ois.readObject();
			depth = intPrefs[0];
			pal = intPrefs[3];
			
			
			//read double pref array
			double[] doublePrefs = (double[])ois.readObject();
			power = doublePrefs[0];
			gradient = doublePrefs[1];
			brightness = doublePrefs[2];
			setZoom( doublePrefs[3] );
			xcen = doublePrefs[4];
			ycen = doublePrefs[5];
			cameraPersp = doublePrefs[6];
			stepDetail = doublePrefs[7];
			cameraYaw = doublePrefs[8];
			cameraPitch = doublePrefs[9];
			cameraDOF = doublePrefs[10];
			opacity = doublePrefs[11];
			focus = doublePrefs[12];
			frost = doublePrefs[13];
			
			try
			{
				inverse_azimuth= intPrefs[4];
				fog_factor = doublePrefs[14];
				primary_light = doublePrefs[15];
				formula = intPrefs[5];
			}
			catch(Exception e)
			{
			}
			
			
			fis.close();
			
			
			if (yimlen != img_alpha[0].length)
			{
				initVars(true);
			}
			else
			{
				ximlen = intPrefs[1];
				yimlen = intPrefs[2];
				initVars(false);
			}
			System.out.println("open success! " + filename);
		}
		catch (Throwable e) 
		{
			System.out.println("read: " + e.toString());
		} 
	}

	*/
	
	

/*
 * Render the Mandelbrot set
 */
function draw()
{
	renderPass = 0;

	if ( lookAt === null ) lookAt = lookAtDefault;
	if ( screenZoom === null ) screenZoom = [zoomStart, zoomStart];

	xRange = [lookAt[0]-screenZoom[0]/2, lookAt[0]+screenZoom[0]/2];
	yRange = [lookAt[1]-screenZoom[1]/2, lookAt[1]+screenZoom[1]/2];
	range_x = xRange[1] - xRange[0];
	range_y = yRange[1] - yRange[0];

	if ( reInitCanvas ) {
		reInitCanvas = false;

		canvas = $('canvasMandelbrot');
		canvas.width  = window.innerWidth;
		canvas.height = window.innerHeight;

		ccanvas = $('canvasControls');
		ccanvas.width  = window.innerWidth;
		ccanvas.height = window.innerHeight;

		ctx = canvas.getContext('2d');
		ctx_img = ctx.createImageData(canvas.width, canvas.height);
	}

	adjustAspectRatio(xRange, yRange, canvas);
	//console.log("screenZoom: " + screenZoom[0] + " " + screenZoom[1]  + " mode: " + mode);

	img_alpha = matrix(canvas.width, canvas.height, 0.0);
	img_red = matrix(canvas.width, canvas.height, 0.0);
	img_green = matrix(canvas.width, canvas.height, 0.0);
	img_blue = matrix(canvas.width, canvas.height, 0.0);

	updateHashTag();
	updateInfoBox();

	// Only enable one render at a time
	renderId += 1;

	render();
}


function render2()
{
	var start  = (new Date).getTime();
	var startHeight = canvas.height;
	var startWidth = canvas.width;
	var lastUpdate = start;
	var pixels = 0;
	var Ci = yRange[0];
	var sy = 0;
	var ourRenderId = renderId;

	var dx = (xRange[1] - xRange[0]) / canvas.width;
	var dy = (yRange[1] - yRange[0]) / canvas.height;

	var scanline = function()
	{
		if (renderId != ourRenderId || startHeight != canvas.height || startWidth != canvas.width)
		{
			// Stop drawing
			console.log("window changed, stopping.")
			return;
		}



		 // Javascript is inherently single-threaded, and the way
		 // you yield thread control back to the browser is MYSTERIOUS.

		 // People seem to use setTimeout() to yield, which lets us
		 // make sure the canvas is updated, so that we can do animations.

		 // But if we do that for every scanline, it will take 100x longer
		 // to render everything, because of overhead.  So therefore, we'll
		 // do something in between.


		var now = (new Date).getTime();
		gridpoints(Ci);
		Ci += dy;
		if (Ci >= canvas.height) Ci = 0;
		pixels += canvas.width;
		
		if ( (now - lastUpdate) >= 10000  || renderPass == 5) {
			// show the user where we're rendering
			updateHistogram();
			updateInfoBox();

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

			lastUpdate = now;
        }



	};

	scanline();
}


	
function render()
{
	// for(var y = 0; y < yimlen; y++)
	var y = 0;

	var scanline = function()
	{
		if (reset == 1)
		{
			reset = 0;
			y = 0;
			clearScreenAndReset(true);
			return;
		}
		else if (reset == 2)
		{
			reset = 0;
			updateHistogram();
			//repaint();
		}
		
		//while (m_down)
		//{
		//	try { Thread.sleep(1000); } catch(Exception e) {}
		//}

		renderpass++;
		
		//local variable to keep consistent across threads
		var currentPass = renderpass;
		
		gridpoints(y);
		if (currentPass % ximlen == 0) updateMinMaxY();
		
		//Thread.sleep(10);
		
		if (currentPass < 20 || currentPass % 50 == 0)
		{
			var t2 = (new Date()).getTime();
			var completeness = Math.round(max_alpha*100)/100.0;
			var strStatus = y + " " + visiblePixels + " " + rayPoints + " " + allPixels + " Pass: " + currentPass + " max value: " + completeness + " in " + (t2-t1);
			t1 = (new Date()).getTime();
			console.log(strStatus);
			//showStatus(strStatus);
			updateHistogram();

			//repaint();

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

		// yield control back to browser, so that canvas is updated
		y++;
		if (y > yimlen-1) y = 0;
        setTimeout(scanline, 100);
	}

	scanline();
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
		for (var  x=0; x<canvas.width; x++)
		{
			var z = Math.pow(img_alpha[x][y], gradient)*brightness/max_alpha;
			red = ( (img_red[x][y]*z)/img_alpha[x][y] );
			green = ( (img_green[x][y]*z)/img_alpha[x][y] );
			blue = ( (img_blue[x][y]*z)/img_alpha[x][y] );

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
}





/*
 * When resizing the window, be sure to update all the canvas stuff.
 */
window.onresize = function(event)
{
	reInitCanvas = true;
};


function init()
{
	yimlen = canvas.height;
	ximlen = canvas.width;

	half_ximlen = ximlen / 2;
	half_yimlen = yimlen / 2;

	CameraMatrix = matrix(3, 3, 0.0);
	IrotX = matrix(3, 3, 0.0);
	IrotZ = matrix(3, 3, 0.0);
	rotX = matrix(3, 3, 0.0);
	rotZ = matrix(3, 3, 0.0);

	clearScreenAndReset(true);
	setConstants();
	getConstants();
}

function main()
{
	$('viewPNG').onclick = function(event)
	{
		window.open(canvas.toDataURL('image/png'));
	};


	/*
	 * Read hash tag and render away at page load.
	 */
	readHashTag();
	setZoom(zoom);
	init();
	draw();
}

main();





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
		screenZoom[0] *= xf;
	} else {
		var yf = ratio/sratio;
		yRange[0] *= yf;
		yRange[1] *= yf;
		screenZoom[1] *= yf;
	}
	xRange = [lookAt[0]-screenZoom[0]/2, lookAt[0]+screenZoom[0]/2];
	yRange = [lookAt[1]-screenZoom[1]/2, lookAt[1]+screenZoom[1]/2];
	range_x = xRange[1] - xRange[0];
	range_y = yRange[1] - yRange[0];
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
				screenZoom = [parseFloat(z[0]), parseFloat(z[1])];
				console.log("readHashTag(): " + screenZoom);
				redraw = true;
			} break;
		}
	}

	if ( redraw )
		reInitCanvas = true;

	return redraw;
}


/*
 * Update URL's hash with render parameters so we can pass it around.
 */
function updateHashTag()
{
	var alg = 0;

	console.log("updateHashTag(): " + screenZoom);

	location.hash = 'zoom=' + screenZoom + '&' +
									'lookAt=' + lookAt + '&' +
									'contrast=' + gradient + '&' +
									'brightness=' + brightness;
}

/*
 * Update small info box in lower right hand side
 */
function updateInfoBox()
{
	// Update infobox
	$('infoBox').innerHTML =
		'x<sub>0</sub>=' + xRange[0] + ' y<sub>0</sub>=' + yRange[0] + ' ' +
		'x<sub>1</sub>=' + xRange[1] + ' y<sub>1</sub>=' + yRange[1] + ' ' +
		'wxh=' + canvas.width + 'x' + canvas.height + ' ' + (canvas.width*canvas.height/1000000.0).toFixed(1) + 'MP pass: ' + renderPass;
}



// geyser27.ppm 
var pallet =
[
 [100, 0, 0], [26, 5, 5], [77, 17, 6], [98, 45, 29],
 [124, 78, 66], [114, 129, 146], [115, 168, 193], [132, 179, 210],
 [152, 186, 220], [153, 176, 202], [142, 171, 181], [157, 169, 176],
 [165, 162, 158], [154, 131, 103], [149, 104, 73], [152, 98, 60],
 [140, 85, 46], [138, 67, 23], [163, 57, 6], [173, 55, 8],
 [179, 56, 3], [190, 61, 8], [208, 86, 7], [210, 110, 25],
 [221, 127, 29], [222, 150, 44], [226, 161, 53], [222, 159, 52],
 [222, 152, 60], [213, 148, 66], [193, 134, 70], [164, 146, 113],
 [138, 141, 159], [133, 144, 165], [154, 139, 122], [148, 108, 84],
 [158, 115, 69], [189, 121, 57], [198, 124, 48], [212, 136, 45],
 [224, 146, 39], [223, 143, 35], [222, 137, 31], [222, 136, 24],
 [221, 137, 18], [205, 101, 19], [198, 76, 4], [176, 60, 14],
 [148, 86, 43], [154, 110, 61], [141, 115, 89], [125, 168, 156],
 [126, 173, 199], [109, 156, 212], [105, 144, 195], [96, 120, 161],
 [71, 74, 121], [37, 58, 125], [26, 49, 119], [16, 42, 117],
 [17, 41, 115], [18, 38, 106], [35, 42, 95], [53, 36, 45],
 [59, 35, 35], [84, 25, 14], [98, 35, 15], [116, 48, 18],
 [121, 56, 23], [135, 72, 28], [142, 82, 36], [163, 103, 51],
 [178, 108, 56], [196, 123, 61], [214, 145, 68], [202, 165, 99],
 [182, 166, 149], [177, 174, 165], [175, 176, 165], [157, 168, 171],
 [135, 153, 167], [120, 125, 129], [126, 87, 69], [125, 82, 64],
 [126, 81, 63], [127, 75, 46], [137, 83, 37], [142, 82, 36],
 [142, 80, 35], [134, 70, 20], [116, 51, 9], [107, 38, 11],
 [97, 29, 7], [86, 25, 13], [87, 26, 8], [93, 28, 8],
 [103, 31, 6], [110, 40, 6], [111, 43, 7], [141, 28, 8],
 [143, 31, 7], [145, 34, 10], [148, 36, 6], [145, 39, 3],
 [128, 54, 22], [115, 51, 25], [99, 50, 25], [99, 47, 28],
 [102, 44, 27], [95, 42, 26], [53, 35, 43], [23, 38, 100],
 [14, 34, 103], [18, 28, 89], [18, 28, 89], [18, 28, 89],
 [34, 13, 12], [75, 15, 6], [85, 15, 2], [98, 23, 4],
 [101, 24, 2], [105, 32, 2], [106, 34, 2], [104, 37, 12],
 [100, 33, 14], [91, 30, 11], [88, 29, 11], [93, 38, 18],
 [96, 39, 24], [98, 47, 30], [69, 63, 76], [36, 54, 121],
 [49, 63, 118], [100, 68, 69], [114, 71, 55], [131, 86, 57],
 [159, 109, 61], [186, 112, 45], [200, 123, 46], [217, 143, 50],
 [221, 149, 55], [221, 147, 50], [211, 139, 54], [188, 119, 49],
 [175, 103, 38], [152, 85, 32], [165, 61, 23], [160, 43, 7],
 [163, 41, 6], [153, 31, 3], [156, 35, 7], [162, 41, 7],
 [164, 43, 5], [170, 53, 9], [182, 67, 12], [195, 110, 40],
 [213, 134, 48], [224, 152, 58], [223, 161, 67], [222, 154, 64],
 [195, 134, 61], [181, 119, 60], [156, 111, 61], [145, 104, 74],
 [126, 81, 65], [91, 91, 112], [114, 123, 157], [118, 163, 183],
 [127, 182, 196], [128, 186, 211], [132, 187, 192], [154, 193, 164],
 [167, 204, 164], [170, 202, 163], [180, 203, 188], [163, 198, 225],
 [139, 193, 222], [136, 192, 219], [119, 177, 225], [115, 173, 225],
 [112, 166, 225], [109, 154, 202], [113, 118, 144], [119, 77, 72],
 [112, 67, 52], [102, 53, 32], [102, 53, 33], [111, 67, 54],
 [90, 85, 91], [44, 63, 127], [32, 52, 121], [26, 46, 119],
 [18, 38, 109], [16, 31, 95], [16, 31, 95], [16, 31, 95],
 [16, 31, 95], [23, 21, 51], [33, 11, 21], [46, 37, 52],
 [55, 66, 120], [98, 107, 144], [133, 149, 176], [137, 158, 190],
 [121, 132, 164], [76, 86, 128], [40, 58, 125], [30, 50, 118],
 [27, 47, 121], [39, 50, 109], [107, 61, 46], [104, 57, 36],
 [110, 60, 40], [130, 82, 55], [153, 108, 60], [155, 114, 66],
 [186, 132, 68], [216, 177, 98], [221, 182, 114], [220, 213, 126],
 [223, 205, 147], [194, 197, 178], [174, 198, 177], [183, 184, 159],
 [192, 145, 102], [178, 111, 64], [170, 81, 28], [173, 58, 7],
 [173, 52, 7], [173, 50, 5], [167, 48, 2], [160, 47, 3],
 [158, 48, 2], [163, 51, 6], [174, 58, 8], [189, 97, 22],
 [209, 124, 41], [219, 144, 52], [219, 148, 62], [200, 138, 62],
 [187, 134, 72], [171, 156, 133], [167, 172, 160], [138, 150, 173],
 [126, 139, 157], [142, 106, 81], [117, 71, 55], [108, 50, 28],
 [98, 42, 18], [89, 28, 10], [81, 23, 6], [78, 14, 2],
 [78, 14, 1], [88, 14, 1], [98, 14, 1], [108, 14, 1],
];

