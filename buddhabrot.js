/*
 * The Mandelbrot Set, in HTML5 canvas and javascript.
 * https://github.com/cslarsen/mandelbrot-js
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

/*
 * Global variables:
 */

var mode = 0;

var pointsPerPass = 1500;
var renderPass = 0;
var max_depth;
var img_cache;
var img_red;
var img_gree;
var img_blue;
var depth_blue = 60;
var depth_green = 100;
var depth_red = 500;
var border_buffer = 40;
var orbitWell = 0.001;

//for Julia
var juliax = 0;
var juliay = 0;

var gradient = 0.5;
var brightness = 1.1;

var normalZoom = true;
var doHolo = false;
var doInverse = false;


//adaptive algorithm for zooms
var accepted = 0;
var rejected = 0;
var smart_random_points = new matrix(pointsPerPass,2,0.0);
var smart_random_points_score = new Array(pointsPerPass);
var smart_point_factor = 0.01;

//var avg_score = new double[3];
//var replacePoint[][] = new int[3][10];


var zoomStart = 3.4;
var zoom = [zoomStart, zoomStart];
var lookAtDefault = [-0.6, 0.01];
var lookAt = lookAtDefault;
var xRange = [0, 0];
var yRange = [0, 0];
var range_x;
var range_y;
var interiorColor = [0, 0, 0, 255];
var reInitCanvas = true; // Whether to reload canvas size, etc
var dragToZoom = true;
var colors = [[0,0,0,0]];
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

function resetScore()
{
	console.log("reset score");
	accepted = 0;
	rejected = 0;

	var range_x = xRange[1] - xRange[0];
	var range_y = yRange[1] - yRange[0];
	for (var i = 0; i < pointsPerPass; i++)
	{
		if (i % 10 == 0)
		{
			smart_random_points_score[i] = 0;
			smart_random_points[i][0] = xRange[0] + Math.random()*range_x*2;
			smart_random_points[i][1] = yRange[0] + Math.random()*range_y*2;
		}
		else
		{
			smart_random_points_score[i] = 0;
			smart_random_points[i][0] = 4*(Math.random()-0.5);
			smart_random_points[i][1] = 4*(Math.random()-0.5);
		}
	}
}


// file:///C:/devel/mandelbrot/html5/buddhabrot.html#zoom=0.3144193571865245,0.16146744072182975&lookAt=-0.14251218380407338,-1.0336744869574346
// initial bench = 0.06%
// best = 80%
function updateScore(counter, i)
{
	//score how interesting this point is
	if (smart_random_points_score[i] > 0.01)
		smart_random_points_score[i] = 0.05*counter + 0.95*smart_random_points_score[i];
	else if (smart_random_points_score[i] == 0.01)
		smart_random_points_score[i] = counter*10;
	else
		smart_random_points_score[i] = counter;
		
	if (smart_random_points_score[i] > 3.1)
	{
		//this is an interesting point, try another one near it next time
		accepted++;
		smart_random_points[i][0] += (Math.random() - 0.5)*smart_point_factor*zoom[0];
		smart_random_points[i][1] += (Math.random() - 0.5)*smart_point_factor*zoom[0];
	}
	else
	{
		//try another random point
		rejected++;
		if (i > pointsPerPass*0.9)
		{
			smart_random_points_score[i] = 0.01;
			smart_random_points[i][0] = 4*(Math.random()-0.5);
			smart_random_points[i][1] = 4*(Math.random()-0.5);
	    }
	    else
	    {
	    	smart_random_points_score[i] = 0;
			smart_random_points[i][0] = xRange[0] + Math.random()*range_x;
			smart_random_points[i][1] = yRange[0] + Math.random()*range_y;
		}
	}
}

function updateHistogram()
{
	var off = 0;
	var red,green,blue,off;
	var max_red = Math.pow(findPeak(img_red), gradient);
	var max_green = Math.pow(findPeak(img_green), gradient);
	var max_blue = Math.pow(findPeak(img_blue), gradient);
	//console.log("histogram peaks: " + max_red + " " + max_green + " " + max_blue);
	
	for (var y=0; y<canvas.height; y++)
	{
		for (var  x=0; x<canvas.width; x++)
		{
			red = (( brightness*Math.pow(img_red[x][y], gradient)/max_red)*255);
			green = (( brightness*Math.pow(img_green[x][y], gradient)/max_green)*255);
			blue = (( brightness*Math.pow(img_blue[x][y], gradient)/max_blue)*255);

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





function calcReal(aa, bb, a, b)
{
	//for Julia
	if (juliax != 0)
	{
		a = juliax;
		b = juliay;
	}
	

	if (mode == 1)
	{
		return (aa*aa*aa - 3*aa*bb*bb) + a;
	}
	else if (mode == 2)
	{
		return aa*aa - Math.sin(bb)*(bb) + Math.atan(a*b);
	}
	else if (mode == 3)
	{
		return aa*aa - bb*bb + a;
	}
	else if (mode == 4)
	{
		return Math.sin(aa*aa - bb*bb) + a;
	}
	else if (mode == 5)
	{
		return -aa*aa + bb*bb + a;
	}
	else if (mode == 6)
	{
		return Math.tan(aa*aa - bb*bb) - a;
	}
	else
	{
		return aa*aa - bb*bb + a;
	}
	
}

function calcImg(aa, bb, a, b)
{
	//for Julia
	if (juliax != 0)
	{
		a = juliax;
		b = juliay;
	}
	

	if (mode == 1)
	{
		return (3*aa*aa*bb - bb*bb*bb) + b;
	}
	else if (mode == 2)
	{
		return 2*aa*bb + b - a;
	}
	else if (mode == 3)
	{
		return Math.atan(2*aa*bb) + b;
	}
	else if (mode == 4)
	{
		return 2*aa*bb + b;
	}
	else if (mode == 5)
	{
		return 2*aa*bb + b;
	}
	else
	{
		return 2*aa*bb + b;
	}
	
}

function calcAllOrbits(Cr, Ci, highlightPoints, i)
{
	//// red  ///////////////////////////////////////////////////////////////////// 
	iter = calcOrbit(Cr, Ci, depth_red, 0);
	if ( (iter == depth_red && doInverse) || (iter < depth_red && !doInverse) ) 
	{
		counter = drawOrbit(iter, img_red, highlightPoints ? '#FF0000' : null);
	}
	else
	{
		counter = 0;
	}
	if (i > -1) updateScore(counter, i);


	//// green  ///////////////////////////////////////////////////////////////////// 
	iter = calcOrbit(Cr, Ci, depth_green, 1);
	if ( (iter == depth_green && doInverse) || (iter < depth_green && !doInverse) ) 
	{ 
		counter = drawOrbit(iter, img_green,  highlightPoints ? '#00FF00' : null);
	}
	else
	{
		counter = 0;
	}
	//if (i > -1) updateScore(counter, i);

	//// blue  ///////////////////////////////////////////////////////////////////// 
	iter = calcOrbit(Cr, Ci, depth_blue, 2);
	if ( (iter == depth_blue && doInverse) || (iter < depth_blue && !doInverse) ) 
	{ 
		counter = drawOrbit(iter, img_blue,  highlightPoints ? '#0000FF' : null);
	}
	else
	{
		counter = 0;
	}
	//if (i > -1) updateScore(counter, i);


}

function drawRandomPass()
{
	console.log("Random Pass");
	for (var i = 0; i < pointsPerPass; i++)
	{
		var Cr = smart_random_points[i][0];
		var Ci = smart_random_points[i][1];

		calcAllOrbits(Cr, Ci, renderPass < 800, i);
	}

	//get a fresh random sample of points after a while
	renderPass++;
	if (renderPass % 500 == 300) {
		resetScore();
	}
}

function drawLinePass(Ci, Cr_init, Cr_step)
{
	var Cr = Cr_init;
	for ( var x=0; x<canvas.width; ++x, Cr += Cr_step ) {
		calcAllOrbits(Cr, Ci, true, -1);
	}
	renderPass++;
}

function calcOrbit(Cr, Ci, depth, color)
{
	var mag2,aanew,bbnew;
	var aa = Cr;
	var bb = Ci;
	var iter = 0;
	
	do
	{
		aanew = calcReal(aa, bb, Cr, Ci);
		bbnew = calcImg(aa, bb, Cr, Ci);
		aa = aanew;
		bb = bbnew;
		mag2 = aa*aa + bb*bb;
		
		img_cache[iter][0] = aa;
		img_cache[iter][1] = bb;
		
		iter++;
	} 
	while ( (iter < depth) && (mag2 < 4) );
	
	return iter;
}

function drawOrbit(iter, img_color, string_color)
{
	var pointx,pointy;
	var counter = 0;
	
	if (string_color != null)
		ctx.fillStyle = string_color;

	for (var j = 0; j < iter; j++)
	{
		if (normalZoom)
		{
			pointx = Math.floor( ((img_cache[j][0]-lookAt[0])*canvas.width)/zoom[0] + canvas.width/2.0 );
			pointy = Math.floor( ((img_cache[j][1]-lookAt[1])*canvas.height)/zoom[1] + canvas.height/2.0 );
			//console.log(zoom[0] + " " + img_cache[j][0] + " " + img_cache[j][1] + "   " + pointx + " " + pointy);
		}
		else
		{
			pointx = Math.floor(((img_cache[j][0]+2.0)/4.0)*canvas.width);
			pointy = Math.floor(((img_cache[j][1]+2.0)/4.0)*canvas.height);
		}
		
		if (pointx >= -border_buffer && pointy >= -border_buffer && pointx < canvas.width + border_buffer && pointy < canvas.height + border_buffer)
		{
			counter++;
			if (pointx >= 0 && pointy >= 0 && pointx < canvas.width && pointy < canvas.height)
			{
				img_color[pointx][pointy]++;

				if (string_color != null) ctx.fillRect(pointx,pointy,1,1);

				//stop accumulating if the point is just repeating itself
				if ( j > 15 && doInverse &&
					(  (Math.abs(img_cache[j][0] - img_cache[j-1][0]) < orbitWell && Math.abs(img_cache[j][1] - img_cache[j-1][1]) < orbitWell ) 
					|| (Math.abs(img_cache[j][0] - img_cache[j-2][0]) < orbitWell && Math.abs(img_cache[j][1] - img_cache[j-2][1]) < orbitWell ) 
					|| (Math.abs(img_cache[j][0] - img_cache[j-3][0]) < orbitWell && Math.abs(img_cache[j][1] - img_cache[j-3][1]) < orbitWell )
					)) 
					j = iter;
			}
		}
	}
	
	//return the number of points we actually plotted to the screen
	return counter;
}

	/*
 * Render the Mandelbrot set
 */
function draw()
{
	renderPass = 0;
	mode = $("algorithm").selectedIndex;
	console.log("alg: " + mode);
	doInverse = $("inverse").checked;

	depth_red = parseInt($("redText").value);
	depth_green = parseInt($("greenText").value);
	depth_blue = parseInt($("blueText").value);

	max_depth = Math.max(depth_green, Math.max(depth_blue, depth_red));
	img_cache = matrix(max_depth, 2, 0.0);

	resetScore();

	if ( lookAt === null ) lookAt = lookAtDefault;
	if ( zoom === null ) zoom = [zoomStart, zoomStart];

	xRange = [lookAt[0]-zoom[0]/2, lookAt[0]+zoom[0]/2];
	yRange = [lookAt[1]-zoom[1]/2, lookAt[1]+zoom[1]/2];
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
	console.log("zoom: " + zoom[0] + " " + zoom[1]  + " mode: " + mode);

	img_red = new matrix(canvas.width, canvas.height, 0.0);
	img_green = new matrix(canvas.width, canvas.height, 0.0);
	img_blue = new matrix(canvas.width, canvas.height, 0.0);

	updateHashTag();
	updateInfoBox();

	// Only enable one render at a time
	renderId += 1;

	render();
}



function render()
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
		//drawRandomPass();
		var now = (new Date).getTime();
		if (sy++ < canvas.height && zoom[0] > 3) {
			drawLinePass(Ci, xRange[0], dx);
			Ci += dy;
			pixels += canvas.width;
		}
		else {
			drawRandomPass();
			pixels += pointsPerPass;
		}

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

        // yield control back to browser, so that canvas is updated
        setTimeout(scanline, 30);

	};

	scanline();
}


/*
 * Update URL's hash with render parameters so we can pass it around.
 */
function updateHashTag()
{
	var alg = $('algorithm').value;

	console.log("updateHashTag(): " + zoom);

	location.hash = 'zoom=' + zoom + '&' +
									'lookAt=' + lookAt + '&' +
									'contrast=' + gradient + '&' +
									'brightness=' + brightness + '&' +
									'red=' + depth_red + '&' +
									'green=' + depth_green + '&' +
									'blue=' + depth_blue + '&' +
									'juliax=' + juliax + '&' +
									'juliay=' + juliay + '&' +
									'algorithm=' + alg;
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
		'wxh=' + canvas.width + 'x' + canvas.height + ' '
				+ (canvas.width*canvas.height/1000000.0).toFixed(1) + 'MP pass: ' + renderPass + ' accepted %' + ( Math.floor( (accepted*10000.0001)/(rejected+accepted+1.1) ) / 100.0);
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
				console.log("readHashTag(): " + zoom);
				redraw = true;
			} break;

			case 'lookAt': {
				var l = val.split(',');
				lookAt = [parseFloat(l[0]), parseFloat(l[1])];
				redraw = true;
			} break;

			case 'red': {
				$('redText').value = val;
				redraw = true;
			} break;

			case 'green': {
				$('greenText').value = val;
				redraw = true;
			} break;

			case 'blue': {
				$('blueText').value = val;
				redraw = true;
			} break;

			case 'juliax': {
				juliax = parseFloat(val);
				redraw = true;
			} break;

			case 'juliay': {
				juliay = parseFloat(val);
				redraw = true;
			} break;

			case 'contrast': {
				gradient = parseFloat(val);
			} break;

			case 'brightness': {
				brightness = parseFloat(val);
			} break;

			case 'algorithm': {
				$('algorithm').value = String(val);
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
	xRange = [lookAt[0]-zoom[0]/2, lookAt[0]+zoom[0]/2];
	yRange = [lookAt[1]-zoom[1]/2, lookAt[1]+zoom[1]/2];
	range_x = xRange[1] - xRange[0];
	range_y = yRange[1] - yRange[0];
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

function doJulia (event) {
	var x = event.clientX;
	var y = event.clientY;

	console.log(x + "x" + y + "   " + xRange[0] + "x" + yRange[0]);

	var dx = (xRange[1] - xRange[0]) / canvas.width;
	var dy = (yRange[1] - yRange[0]) / canvas.height;

	juliax = xRange[0] + x*dx;
	juliay = yRange[0] + y*dy;

	$("contrastSlider").value = 35;
	gradient = 0.35;

	draw();
}

function zoomOut (event) {
	var x = event.clientX;
	var y = event.clientY;

	var dx = (xRange[1] - xRange[0]) / canvas.width;
	var dy = (yRange[1] - yRange[0]) / canvas.height;

	x = xRange[0] + x*dx;
	y = yRange[0] + y*dy;

	lookAt = [x, y];

	if ( event.shiftKey ) {
		zoom[0] /= 0.5;
		zoom[1] /= 0.5;
	}
	draw();
}

/*
 * When resizing the window, be sure to update all the canvas stuff.
 */
window.onresize = function(event)
{
	reInitCanvas = true;
};

function main()
{
	$('viewPNG').onclick = function(event)
	{
		window.open(canvas.toDataURL('image/png'));
	};

	$("contrastSlider").onchange = function() {
		gradient =  $("contrastSlider").value / 100.0; 
		updateHistogram();
		updateHashTag();
		console.log("gradient: " + gradient);
	}

	$("brightnessSlider").onchange = function() {
		brightness =  $("brightnessSlider").value / 100.0; 
		updateHistogram();
		updateHashTag();
		console.log("brightness: " + brightness);
	}

	$("redText").onchange = function() {

		var tmp_red = parseInt($("redText").value);
		max_depth = Math.max(max_depth, tmp_red);
		if(max_depth > img_cache.length)
		{
			img_cache = matrix(max_depth, 2, 0.0);
		}
		depth_red = tmp_red;
		updateHashTag();
	}

	$("greenText").onchange = function() {
		var tmp_green = parseInt($("greenText").value);
		max_depth = Math.max(max_depth, tmp_green);
		if(max_depth > img_cache.length)
		{
			img_cache = matrix(max_depth, 2, 0.0);
		}
		depth_green = tmp_green;
		updateHashTag();
	}

	$("blueText").onchange = function() {
		var tmp_blue = parseInt($("blueText").value);
		max_depth = Math.max(max_depth, tmp_blue);
		if(max_depth > img_cache.length)
		{
			img_cache = matrix(max_depth, 2, 0.0);
		}
		depth_blue = tmp_blue;
		updateHashTag();
	}

	$("inverse").onchange = function() {
		doInverse = $("inverse").checked;
		$("contrastSlider").value = 35;
		gradient = 0.35;
		draw();
	}

	$("algorithm").onchange = function() {
		mode = $("algorithm").selectedIndex;
		draw();
	}

	$('resetButton').onclick = function(even)
	{
		$('settingsForm').reset();
		setTimeout(function() { location.hash = ''; }, 1);
		zoom = [zoomStart, zoomStart];
		lookAt = lookAtDefault;
		reInitCanvas = true;
		juliax = 0;
		draw();
	};

	if ( dragToZoom == true ) {
		var box = null;

		$('canvasControls').onmousedown = function(e)
		{
			if ( e.shiftKey ) {
				doJulia(e);
				return;
			}
			else if ( box == null )
				box = [e.clientX, e.clientY, 0, 0];
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

		$('canvasControls').onmouseup = function(e)
		{
			console.log("mouse up!");
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

				var dx = (xRange[1] - xRange[0]) / canvas.width;
				var dy = (yRange[1] - yRange[0]) / canvas.height;

				x = xRange[0] + x*dx;
				y = yRange[0] + y*dy;

				lookAt = [x, y];

				/*
				 * This whole code is such a mess ...
				 */

				var xf = Math.abs(Math.abs(box[0]-box[2])/canvas.width);
				var yf = Math.abs(Math.abs(box[1]-box[3])/canvas.height);

				zoom[0] *= Math.max(xf, yf); // retain aspect ratio
				zoom[1] *= Math.max(xf, yf);

				box = null;
				draw();
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

			var dx = (xRange[1] - xRange[0]) / canvas.width;
			var dy = (yRange[1] - yRange[0]) / canvas.height;

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

			draw();
		};
	}


	/*
	 * Read hash tag and render away at page load.
	 */
	readHashTag();
	draw();
}

main();


//file:///C:/devel/mandelbrot/html5/buddhabrot.html#zoom=6.613981762917935,3.4068895643363732&lookAt=-0.6,0&red=500&green=100&blue=600&colorScheme=pickGradientColorBands

//file:///C:/devel/mandelbrot/html5/buddhabrot.html#zoom=6.600606673407482,3.4&lookAt=-0.6,0&red=500&green=100&blue=60&juliax=-0.2965337124625109&juliay=-0.6234193222053617&colorScheme=pickGradientColorBands

//file:///C:/devel/mandelbrot/html5/buddhabrot.html#zoom=1.9755101873899035,1.0175935288169868&lookAt=-1.1699664267347174,-0.3964087000505816&red=239&green=100&blue=60&juliax=-0.541557128412538&juliay=-0.5930232558139534&colorScheme=pickGradientColorBands



//file:///C:/devel/mandelbrot/html5/buddhabrot.html#zoom=0.08666409798710273,0.044641037973564786&lookAt=-1.7583815408236587,0.007818101974909841&contrast=0.21&brightness=2.82&red=255&green=133&blue=60&juliax=-1.7646859921203775&juliay=0.009424879575539336
//file:///C:/devel/mandelbrot/html5/buddhabrot.html#zoom=0.2986438777325588,0.15383270576953162&lookAt=-1.7400615009158655,0.01104726600598499&contrast=0.35&brightness=1.68&red=255&green=133&blue=60&juliax=-1.7646859921203775&juliay=0.009424879575539336&colorScheme=pickGradientColorBands


//file:///C:/devel/mandelbrot/html5/buddhabrot.html#zoom=6.7,3.4&lookAt=-0.1,0&contrast=0.64&brightness=1.6&red=1500&green=250&blue=50&juliax=0&juliay=0&algorithm=5
