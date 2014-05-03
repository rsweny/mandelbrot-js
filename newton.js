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

//constants
var mode = 0;
var iterations = 20;
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
	var i = parseInt($('superSamples').value, 10);
	return i<=0? 1 : i;
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
		
		if (mandelbrotAddition)
			z = z.add(Complex(i/2.0,j/2.0));
			
		n++;
		
		//normal smoothing
		w = 1.0 / distance(z.sub(old), Complex["0"]);
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
var rootBoundry = 0.00000001;
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

	if (roots.length < 20 /*&& quad > 4*/)
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


/*
	
	//z^10 + 0.2 i * z^5 - 1.
	public Complex Poly2Function(Complex z, double i, double j) 
	{
		Complex exponent = new Complex(order, img_order);
		Complex point2i = new Complex(0, 0.2);
		Complex ten = new Complex(10, 0);
		return (z.pow(ten)).add (point2i.mul(z.pow(exponent))).sub ( one );
	}
	
	//10z^9 + 0.2i*5*z^4
	public Complex Poly2Derivative(Complex z, double i, double j) 
	{
		Complex exponent = new Complex(order, img_order);
		Complex exponentLessOne = exponent.sub(oneError);
		Complex ten = new Complex(10, 0);
		Complex nine = new Complex(9, 0);
		Complex point2i = new Complex(0, 0.2);

		return (ten.mul(z.pow(nine))).add(  point2i.mul(exponent).mul(z.pow(exponentLessOne)) );
	}
	

	//2z^3 - c + 1
	public Complex PolyMFunction(Complex z, double i, double j) 
	{
		Complex c = new Complex(i, j);
		Complex exponent = new Complex(order, img_order);
		return two.mul( z.pow(exponent) ).sub(c).add(one);
	}
	
	//6z^2 - 1
	public Complex PolyMDerivative(Complex z, double i, double j) 
	{
		Complex exponent = new Complex(order, img_order);
		Complex exponentLessOne = exponent.sub(oneError);
		return two.mul(exponent).mul( z.pow(exponentLessOne) ).sub(one);
	}
	
	
	
	//z^c - z + 0.1
	public Complex Poly3Function(Complex z, double i, double j) 
	{
		Complex c = new Complex(1+order, img_order);
		return z.pow(c).sub(z).add(pointone);
	}
	
	//c*z^(c-1) - 1
	public Complex Poly3Derivative(Complex z, double i, double j) 
	{
		Complex cminus1 = new Complex(1+order-1, (img_order == 0) ? 0 : img_order-1);
		Complex c = new Complex(1+order, img_order);
		return c.mul(z.pow(cminus1)).sub(one);
	}
	
	
	//z^n - 3z^5 + 6z^3 - 3z + 3
	public Complex PolyFunction(Complex z, double i, double j) 
	{
		Complex exponent = new Complex(order, img_order);
		Complex three = new Complex(3, 0);
		Complex six = new Complex(6, 0);
		Complex five = new Complex(5, 0);
		
		return (z.pow(exponent)).sub (three.mul(z.pow(five))).add ( six.mul(z.pow(three)) ).sub (three.mul(z)).add(three);
	}
	
	public Complex PolyDerivative(Complex z, double i, double j) 
	{
		Complex exponent = new Complex(order, img_order);
		Complex exponentLessOne = exponent.sub(oneError);
		Complex three = new Complex(3, 0);
		Complex eighteen = new Complex(18, 0);
		Complex fifteen = new Complex(15, 0);
		Complex four = new Complex(4, 0);
		return (exponent.mul(z.pow(exponentLessOne))).sub (fifteen.mul(z.pow(four))).add ( eighteen.mul(z.pow(two)) ).sub(three);
	}

	//z^z - cz
	public Complex ZZFunction(Complex z, double i, double j) 
	{
		Complex con = new Complex(order, img_order);
		return (z.pow(z)).sub( con.mul(z) );
	}
	
	//z^z * (1 + lnz) - c
	public Complex ZZDerivative(Complex z, double i, double j) 
	{
		Complex con = new Complex(order, img_order);
		Complex lnz = z.log().add(one);
		return ((z.pow(z)).mul( lnz )).sub(con);
	}

	*/






/*
 * Render the Mandelbrot set
 */
function draw(superSamples)
{
	mode = $("mode").selectedIndex;

	var equation, derivative;
	if (mode == 0)
	{
		equation = UnityFunction;
		derivative = UnityDerivative;
	}
	else if (mode == 1)
	{
		equation = DivZFunction
		derivative = DivZDerivative
	}

	console.log("mode: " + mode);
	order = $("real_exponent").value;
	img_order = parseFloat($("complex_exponent").value);
	iterations = $("txtSteps").value;
	roots = [];


	if ( lookAt === null ) lookAt = [-0.6, 0];
	if ( zoom === null ) zoom = [zoomStart, zoomStart];

	xRange = [lookAt[0]-zoom[0]/2, lookAt[0]+zoom[0]/2];
	yRange = [lookAt[1]-zoom[1]/2, lookAt[1]+zoom[1]/2];

	mandelbrotAddition = $("mandelbrot").checked;

	initialColor = $("colorSlider").value / 100.0;
	contrast = $("contrastSlider").value / 100.0;
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
	updateInfoBox();

	// Only enable one render at a time
	renderId += 1;

	function drawLineSuperSampled(Ci, off, Cr_init, Cr_step)
	{
		var Cr = Cr_init;

		for ( var x=0; x<canvas.width; ++x, Cr += Cr_step ) {
			var color = [0, 0, 0, 255];

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

	function drawLine(Ci, off, Cr_init, Cr_step)
	{
		var Cr = Cr_init;

		for ( var x=0; x<canvas.width; ++x, Cr += Cr_step ) {
			var p = iterateEquation(Cr, Ci, equation, derivative);
			var color = pickColor(p);
			img.data[off++] = color[0];
			img.data[off++] = color[1];
			img.data[off++] = color[2];
			img.data[off++] = 255;
		}
	}

	function pickColor(iter)
	{
		var hue;
		if (mandelbrotAddition)
		{
			//get Color based on angle of z	
			var angle = Math.atan(iter[3]/iter[2]);

			angle = angle + (Math.PI / 2);
			angle = angle / Math.PI;
			hue = angle;
		}
		else
		{
			//limit the colors to 30% of the spectrum, unless a complex root
			var limitfactor = 0.3;
			if (img_order != 0) limitfactor = 1.0;	
			hue = initialColor + (iter[1]*limitfactor) / roots.length;
		}

		//normal shading
		huesat = iter[0]/(iterations/brightness);
		if (huesat >= 1.0) huesat = 0.9999;

		huesat = (iter[0] == iterations) ? 1 : huesat * 16000000;
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
			if (    renderId != ourRenderId ||
					 startHeight != canvas.height ||
						startWidth != canvas.width )
			{
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

	render();
}

// Some constants used with smoothColor
var logBase = 1.0 / Math.log(2.0);
var logHalfBase = Math.log(0.5)*logBase;
var smoothColor = true;
var band = 1.0;
var contrast = 0.35;
var brightness = 2.0;
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

		var rgb = [0,0,0];
	rgb[0] = Math.floor(r * 255)
	rgb[1] = Math.floor(g * 255)
	rgb[2] = Math.floor(b * 255)

		return rgb;
}





/*
 * Update URL's hash with render parameters so we can pass it around.
 */
function updateHashTag(samples)
{
	var scheme = $('mode').value;

	console.log("updateHashTag() " + iterations);

	location.hash = 'zoom=' + zoom + '&' +
									'lookAt=' + lookAt + '&' +
									'iterations=' + iterations + '&' +
									'superSamples=' + samples + '&' +
									'mode=' + scheme;
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
				+ (canvas.width*canvas.height/1000000.0).toFixed(1) + 'MP';
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
				$('txtSteps').value = val;
				console.log("readHashTag() " + $('txtSteps').value);
				redraw = true;
			} break;

			case 'superSamples': {
				$('superSamples').value = String(parseInt(val, 10));
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
		window.location = canvas.toDataURL('image/png');
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

	if ( dragToZoom == true ) {
		var box = null;

		$('canvasControls').onmousedown = function(e)
		{
			if ( box == null )
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

				/*
				 * This whole code is such a mess ...
				 */

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
