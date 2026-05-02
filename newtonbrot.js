/*
 * Newton's Method accumulation fractal, in HTML5 canvas and javascript.
 * Ported from Newtonbrot.java by Ryan Sweny
 * https://github.com/rsweny/mandelbrot-js
 *
 * Newton's method applied to complex polynomial functions.
 * Buddhabrot-style rendering: orbits traced through the plane,
 * accumulated in a histogram, then tone-mapped to colour.
 */

// ---------------------------------------------------------------------------
// Complex arithmetic
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
// z^w  (complex base, complex exponent) via z^w = exp(w * ln z)
Complex.prototype.cpow = function(w) {
  if (this.re === 0 && this.im === 0) return new Complex(0, 0);
  var logr  = Math.log(this.abs());
  var theta = this.arg();
  var newR     = Math.exp(w.re * logr - w.im * theta);
  var newTheta = w.im * logr + w.re * theta;
  if (!isFinite(newR)) return new Complex(NaN, NaN);
  return new Complex(newR * Math.cos(newTheta), newR * Math.sin(newTheta));
};

function cdist(a, b) {
  var dr = a.re - b.re, di = a.im - b.im;
  return Math.sqrt(dr * dr + di * di);
}

// ---------------------------------------------------------------------------
// Colour utilities
// ---------------------------------------------------------------------------

function hsvToRgb(h, s, v) {
  var i = Math.floor(h * 6);
  var f = h * 6 - i;
  var p = v * (1 - s);
  var q = v * (1 - f * s);
  var t = v * (1 - (1 - f) * s);
  var r, g, b;
  switch (i % 6) {
    case 0: r=v; g=t; b=p; break;
    case 1: r=q; g=v; b=p; break;
    case 2: r=p; g=v; b=t; break;
    case 3: r=p; g=q; b=v; break;
    case 4: r=t; g=p; b=v; break;
    default: r=v; g=p; b=q; break;
  }
  return [Math.floor(r*255), Math.floor(g*255), Math.floor(b*255)];
}

// Generate 12 distinct 256-entry palettes procedurally
function buildPalettes() {
  var pals = [];
  for (var p = 0; p < 12; p++) {
    var pal = [];
    for (var i = 0; i <= 255; i++) {
      var t = i / 255.0;
      var h, s, v;
      switch (p % 6) {
        case 0: h = t;                s = 0.90; v = 0.5 + 0.5*t;  break; // rainbow
        case 1: h = 0.55 + 0.45*t;    s = 1.00; v = 0.4 + 0.6*t;  break;
        case 2: h = 0.10*t;           s = 1.00; v = 0.3 + 0.7*t;  break;
        case 3: h = 0.48 + 0.22*t;    s = 0.80; v = 0.3 + 0.7*Math.sin(t*Math.PI); break;
        case 4: h = 0.25 + 0.50*t;    s = 0.80; v = 0.70;          break;
        default:h = 0.85 + 0.20*t;    s = 1.00; v = 0.5 + 0.5*t;  break;
      }
      if (p >= 6) h = (h + 0.5) % 1.0; // complementary second set
      pal.push(hsvToRgb(h, s, v));
    }
    pals.push(pal);
  }
  return pals;
}
var palettes = buildPalettes();

// ---------------------------------------------------------------------------
// Globals
// ---------------------------------------------------------------------------

function $(id) { return document.getElementById(id); }

var canvas, ctx, ccanvas, cctx;
var ximlen, yimlen;

// Flat histogram accumulators (indexed y*ximlen+x)
var imgRed, imgGreen, imgBlue, imgAlpha;
var imageData;
var imgCacheRe, imgCacheIm; // shared orbit scratch buffer

// Fractal parameters (matching Java defaults)
var zoom          = 6;
var xcen          = 0.1;
var ycen          = 0;
var order         = 4.9242613;
var imgOrder      = 0.4841751596941611;
var complexError  = 1.0;
var gradient      = 0.11;
var brightness    = 3.0;
var depthRed      = 200;
var depthGreen    = 140;
var depthBlue     = 40;
var numPoints     = 1000;
var rootBoundary  = 1e-13;
var algMode       = 0;
var doInverse     = false;
var byStructure   = true;
var mandelbrotAdd = false;
var palIndex      = 0;

// Render-loop state
var runtotal      = 1;
var quadSize      = 36;
var LOW_RES_PASSES = 340;
var borderBuffer  = 40;
var running       = false;
var renderTimer   = null;
var lastUpdateTime = 0;
var lastStatusTime = 0;
var accepted = 0, rejected = 0;

// Smart random points  [3 channels][numPoints][xy]
var smartX, smartY, smartScore;
var replacePoint; // [3][10] circular buffer
var avgScore;     // [3]

// Root list for addRoot()
var roots = [];

// Drag-to-zoom
var dragBox = null;

// ---------------------------------------------------------------------------
// UI helpers
// ---------------------------------------------------------------------------

function readControls() {
  try { zoom         = parseFloat($('zoom').value)         || 6;    } catch(e){}
  try { xcen         = parseFloat($('xcen').value)         || 0.1;  } catch(e){}
  try { ycen         = parseFloat($('ycen').value)         || 0;    } catch(e){}
  try { order        = parseFloat($('order').value)        || 4.9;  } catch(e){}
  try { imgOrder     = parseFloat($('imgOrder').value)     || 0.48; } catch(e){}
  try { complexError = parseFloat($('complexError').value) || 1;    } catch(e){}
  gradient   = $('contrastSlider').value / 100;
  brightness = $('brightnessSlider').value / 100;
  try { depthRed     = parseInt($('depthRed').value, 10)   || 200;  } catch(e){}
  try { depthGreen   = parseInt($('depthGreen').value, 10) || 140;  } catch(e){}
  try { depthBlue    = parseInt($('depthBlue').value, 10)  || 40;   } catch(e){}
  try { numPoints    = parseInt($('numPoints').value, 10)  || 1000; } catch(e){}
  algMode       = parseInt($('algMode').value, 10);
  doInverse     = $('doInverse').checked;
  byStructure   = $('byStructure').checked;
  mandelbrotAdd = $('mandelbrotAdd').checked;
  palIndex      = parseInt($('palIndex').value, 10);
}

function writeControls() {
  $('zoom').value         = zoom;
  $('xcen').value         = xcen;
  $('ycen').value         = ycen;
  $('order').value        = order;
  $('imgOrder').value     = imgOrder;
  $('complexError').value = complexError;
  $('contrastSlider').value  = Math.round(gradient * 100);
  $('brightnessSlider').value = Math.round(brightness * 100);
  $('depthRed').value     = depthRed;
  $('depthGreen').value   = depthGreen;
  $('depthBlue').value    = depthBlue;
  $('numPoints').value    = numPoints;
  $('algMode').value      = algMode;
  $('doInverse').checked  = doInverse;
  $('byStructure').checked   = byStructure;
  $('mandelbrotAdd').checked = mandelbrotAdd;
  $('palIndex').value     = palIndex;
}

function updateStatus(msg) { $('status').textContent = msg || ''; }

// ---------------------------------------------------------------------------
// Canvas init
// ---------------------------------------------------------------------------

function initCanvas() {
  canvas  = $('newtonCanvas');
  ccanvas = $('canvasControls');
  canvas.width  = ccanvas.width  = window.innerWidth;
  canvas.height = ccanvas.height = window.innerHeight;
  ctx  = canvas.getContext('2d');
  cctx = ccanvas.getContext('2d');
}

// ---------------------------------------------------------------------------
// Reset / clear
// ---------------------------------------------------------------------------

function resetSmartPoints() {
  for (var ch = 0; ch < 3; ch++) {
    avgScore[ch] = 0;
    for (var i = 0; i < numPoints; i++) {
      smartX[ch][i]     = Math.random();
      smartY[ch][i]     = Math.random();
      smartScore[ch][i] = 0;
    }
  }
}

function clearAndReset(newRoots) {
  readControls();
  if (newRoots) roots = [];

  ximlen = canvas.width;
  yimlen = canvas.height;

  var total = ximlen * yimlen;
  imgRed   = new Float64Array(total);
  imgGreen = new Float64Array(total);
  imgBlue  = new Float64Array(total);
  imgAlpha = new Float64Array(total);

  imageData = ctx.createImageData(ximlen, yimlen);
  for (var i = 3; i < imageData.data.length; i += 4) imageData.data[i] = 255;
  ctx.putImageData(imageData, 0, 0);

  var maxDepth = Math.max(depthRed, depthGreen, depthBlue);
  imgCacheRe = new Float64Array(maxDepth);
  imgCacheIm = new Float64Array(maxDepth);

  smartX     = [new Float64Array(numPoints), new Float64Array(numPoints), new Float64Array(numPoints)];
  smartY     = [new Float64Array(numPoints), new Float64Array(numPoints), new Float64Array(numPoints)];
  smartScore = [new Float64Array(numPoints), new Float64Array(numPoints), new Float64Array(numPoints)];
  avgScore   = [0, 0, 0];
  replacePoint = [[],[],[]];
  for (var ch = 0; ch < 3; ch++) {
    for (var k = 0; k < 10; k++) replacePoint[ch][k] = 0;
  }

  resetSmartPoints();

  runtotal     = 1;
  quadSize     = 36;
  borderBuffer = 40;
  accepted     = 0;
  rejected     = 0;
  lastUpdateTime = 0;
  lastStatusTime = Date.now();
}

// ---------------------------------------------------------------------------
// Newton iteration functions (6 modes)
// Each returns z_{n+1} = z_n - f(z_n)/f'(z_n)
// ---------------------------------------------------------------------------

function newtonStep(z, aa, bb) {
  var exp  = new Complex(order, imgOrder);
  var one  = new Complex(1, 0);
  var oerr = new Complex(1, complexError);   // "oneError" — derivative distortion
  var expM1 = exp.sub(oerr);

  switch (algMode) {
    case 1: { // z^n - 1 = 0  (roots of unity)
      var f  = z.cpow(exp).sub(one);
      var df = exp.mul(z.cpow(expM1));
      return z.sub(f.div(df));
    }
    case 2: { // z^10 + 0.2i*z^n - 1
      var ten  = new Complex(10, 0);
      var nine = new Complex(9, 0);
      var p2i  = new Complex(0, 0.2);
      var f  = z.cpow(ten).add(p2i.mul(z.cpow(exp))).sub(one);
      var df = ten.mul(z.cpow(nine)).add(p2i.mul(exp).mul(z.cpow(expM1)));
      return z.sub(f.div(df));
    }
    case 3: { // (z^n - 1) / z
      var two  = new Complex(2, 0);
      var expM2 = exp.sub(two);
      var f  = z.cpow(exp).sub(one).div(z);
      var df = expM1.mul(z.cpow(expM2)).add(z.cpow(new Complex(-2,0)));
      return z.sub(f.div(df));
    }
    case 4: { // z^c - z + 0.1,  c = Complex(10+order, img_order)
      var c    = new Complex(10 + order, imgOrder);
      var cm1  = new Complex(10 + order - 1, (imgOrder === 0) ? 0 : imgOrder - 1);
      var p1   = new Complex(0.1, 0);
      var f  = z.cpow(c).sub(z).add(p1);
      var df = c.mul(z.cpow(cm1)).sub(one);
      return z.sub(f.div(df));
    }
    case 5: { // z^n - 1/z
      var f  = z.cpow(exp).sub(one.div(z));
      var df = exp.mul(z.cpow(expM1)).add(z.cpow(new Complex(-2,0)));
      return z.sub(f.div(df));
    }
    default: { // z^n - 3z^5 + 6z^3 - 3z + 3
      var three   = new Complex(3, 0);
      var six     = new Complex(6, 0);
      var five    = new Complex(5, 0);
      var four    = new Complex(4, 0);
      var eighteen = new Complex(18, 0);
      var fifteen  = new Complex(15, 0);
      var f  = z.cpow(exp)
        .sub(three.mul(z.cpow(five)))
        .add(six.mul(z.cpow(three)))
        .sub(three.mul(z))
        .add(three);
      var df = exp.mul(z.cpow(expM1))
        .sub(fifteen.mul(z.cpow(four)))
        .add(eighteen.mul(z.cpow(new Complex(2,0))))
        .sub(three);
      return z.sub(f.div(df));
    }
  }
}

// ---------------------------------------------------------------------------
// Root catalogue (for palette-based colouring)
// ---------------------------------------------------------------------------

function addRoot(z) {
  var colorFactor = 2.0 / rootBoundary;
  var minDist = 1e10;
  var closest = 0;

  for (var idx = 0; idx < roots.length; idx++) {
    var d = cdist(z, roots[idx]);
    if (d < minDist) { minDist = d; closest = idx; }
    if (d <= rootBoundary * 1.1) {
      var retVal = idx + colorFactor * d;
      return Math.min(roots.length - 1, retVal);
    }
  }

  if (roots.length < 200 && runtotal < 25) {
    roots.push({re: z.re, im: z.im});
  }

  return doInverse ? closest : (roots.length - 1);
}

// ---------------------------------------------------------------------------
// Orbit computation
// Stores orbit in imgCacheRe/Im, returns [iter, solutionColor]
// ---------------------------------------------------------------------------

function calcOrbit(pointIdx, depth, channel) {
  var xf = smartX[channel][pointIdx];
  var yf = smartY[channel][pointIdx];

  // Starting point in complex plane: [-1,1] x [-1,1]
  var aa = 2 * (xf - 0.5);
  var bb = 2 * (yf - 0.5);

  var z   = new Complex(aa, bb);
  var old = new Complex(aa, bb);
  var iter = 0;

  z = newtonStep(z, aa, bb);

  while (iter < depth && cdist(old, z) > rootBoundary) {
    old.re = z.re; old.im = z.im;
    z = newtonStep(z, aa, bb);
    if (mandelbrotAdd) {
      z.re += aa * 0.5;
      z.im += bb * 0.5;
    }
    if (iter < imgCacheRe.length) {
      imgCacheRe[iter] = z.re;
      imgCacheIm[iter] = z.im;
    }
    iter++;
  }

  var solutionColor = 0;
  if (iter !== depthRed || (doInverse && runtotal > 25)) {
    if (mandelbrotAdd) {
      z = newtonStep(z, aa, bb);
      var angle = Math.atan2(z.im, z.re) + Math.PI / 2;
      solutionColor = (angle / Math.PI) * 254.0;
    } else {
      var curroot = addRoot(z);
      solutionColor = (curroot * 254) / Math.max(1, roots.length);
    }
  }

  return [iter, solutionColor];
}

// ---------------------------------------------------------------------------
// Convert a cached orbit to pixel coordinates and accumulate
// ---------------------------------------------------------------------------

function pixelFromComplex(re, im) {
  // aspect-corrected: both axes use yimlen/zoom as scale
  var scale = yimlen / zoom;
  var px = Math.floor((re + xcen) * scale + ximlen * 0.5);
  var py = Math.floor((im + ycen) * scale + yimlen * 0.5);
  return [px, py];
}

function drawOrbitChannel(iter, imgChannel) {
  var close = 0.001;
  var counter = 0;
  for (var j = 0; j < iter; j++) {
    var pxy = pixelFromComplex(imgCacheRe[j], imgCacheIm[j]);
    var px = pxy[0], py = pxy[1];
    if (px >= -borderBuffer && py >= -borderBuffer &&
        px <  ximlen + borderBuffer && py < yimlen + borderBuffer) {
      counter++;
      if (px >= 0 && py >= 0 && px < ximlen && py < yimlen) {
        counter += 0.001; // tag: actually on screen
        imgChannel[py * ximlen + px]++;
        // Early-exit if orbit is converging (when in inverse mode)
        if (j > 5 && doInverse) {
          var close1 = Math.abs(imgCacheRe[j]-imgCacheRe[j-1]) < close && Math.abs(imgCacheIm[j]-imgCacheIm[j-1]) < close;
          var close2 = j>1 && Math.abs(imgCacheRe[j]-imgCacheRe[j-2]) < close && Math.abs(imgCacheIm[j]-imgCacheIm[j-2]) < close;
          var close3 = j>2 && Math.abs(imgCacheRe[j]-imgCacheRe[j-3]) < close && Math.abs(imgCacheIm[j]-imgCacheIm[j-3]) < close;
          if (close1 || close2 || close3) break;
        }
      }
    }
  }
  return counter;
}

function drawOrbitStructure(iter, amtR, amtG, amtB) {
  var close = 0.001;
  var s2 = quadSize >> 1;
  var counter = 0;

  for (var j = 0; j < iter; j++) {
    var pxy = pixelFromComplex(imgCacheRe[j], imgCacheIm[j]);
    var px = pxy[0], py = pxy[1];

    if (px >= -borderBuffer && py >= -borderBuffer &&
        px <  ximlen + borderBuffer && py < yimlen + borderBuffer) {
      counter++;
      if (px >= 0 && py >= 0 && px < ximlen && py < yimlen) {
        counter += 0.001;

        if (j > 5 && doInverse) {
          var c1 = Math.abs(imgCacheRe[j]-imgCacheRe[j-1]) < close && Math.abs(imgCacheIm[j]-imgCacheIm[j-1]) < close;
          var c2 = j>1 && Math.abs(imgCacheRe[j]-imgCacheRe[j-2]) < close && Math.abs(imgCacheIm[j]-imgCacheIm[j-2]) < close;
          var c3 = j>2 && Math.abs(imgCacheRe[j]-imgCacheRe[j-3]) < close && Math.abs(imgCacheIm[j]-imgCacheIm[j-3]) < close;
          if (c1 || c2 || c3) break;
        }

        var idx = py * ximlen + px;
        imgAlpha[idx] += 10;
        imgRed[idx]   += amtR * 10;
        imgGreen[idx] += amtG * 10;
        imgBlue[idx]  += amtB * 10;

        // Low-res block fill during early passes
        if (runtotal > 50 && runtotal < LOW_RES_PASSES && quadSize > 14) {
          var lo = 7, hi = quadSize - 7;
          for (var k = lo; k < hi; k++) {
            for (var m = lo; m < hi; m++) {
              var nx = px + k - s2;
              var ny = py + m - s2;
              if (nx >= 0 && ny >= 0 && nx < ximlen && ny < yimlen) {
                var nidx = ny * ximlen + nx;
                var sparseColor = 2 - imgAlpha[nidx] > 0;
                if (sparseColor) {
                  imgAlpha[nidx] ++;
                  imgRed[nidx]   += amtR;
                  imgGreen[nidx] += amtG;
                  imgBlue[nidx]  += amtB;
                }
              }
            }
          }
        }
      }
    }
  }
  return counter;
}

// ---------------------------------------------------------------------------
// Adaptive sampling score update (port of Java updateScore)
// ---------------------------------------------------------------------------

function updateScore(counter, channel, pointIdx) {
  var sc = smartScore[channel][pointIdx];
  if (sc > 0)
    smartScore[channel][pointIdx] = 0.5 * counter + 0.5 * sc;
  else
    smartScore[channel][pointIdx] = counter;

  if (borderBuffer > 40) {
    if ((counter - Math.floor(counter)) !== 0)
      smartScore[channel][pointIdx] += 1;
  }

  var score = smartScore[channel][pointIdx];
  if (score > 1.1) {
    accepted++;
    var factor = 0.01;
    if (score > avgScore[channel] * (1.2 + zoom)) {
      var repIdx = replacePoint[channel][0];
      smartX[channel][repIdx] = smartX[channel][pointIdx] + (Math.random()-0.5)*factor*zoom;
      smartY[channel][repIdx] = smartY[channel][pointIdx] + (Math.random()-0.5)*factor*zoom;
      // rotate circular buffer
      for (var k = 0; k < replacePoint[channel].length-1; k++)
        replacePoint[channel][k] = replacePoint[channel][k+1];
      replacePoint[channel][replacePoint[channel].length-1] = pointIdx;
    }
    if (avgScore[channel] > 3 && pointIdx % 3 === 0) factor *= 5;
    smartX[channel][pointIdx] += (Math.random()-0.5)*factor*zoom;
    smartY[channel][pointIdx] += (Math.random()-0.5)*factor*zoom;
  } else {
    rejected++;
    smartScore[channel][pointIdx] = 0;
    smartX[channel][pointIdx] = Math.random();
    smartY[channel][pointIdx] = Math.random();
    // push this weak index onto the replace buffer
    for (var k = replacePoint[channel].length-1; k > 0; k--) {
      replacePoint[channel][k] = replacePoint[channel][k-1];
    }
    replacePoint[channel][0] = pointIdx;
  }
}

function calcAvgScore() {
  var totalC = 0;
  for (var ch = 0; ch < 3; ch++) {
    var c = 0;
    for (var i = 0; i < numPoints; i++) c += smartScore[ch][i];
    avgScore[ch] = c / numPoints;
    totalC += c;

    // Propagate successful channel to struggling channels
    if (avgScore[ch] > 1 && borderBuffer === 40) {
      for (var other = 0; other < 3; other++) {
        if (other !== ch && avgScore[other] < 0.4) {
          var start = Math.floor(Math.random() * (numPoints/2));
          for (var j = start; j < start + numPoints/4 && j < numPoints; j++) {
            smartX[other][j]     = smartX[ch][j];
            smartY[other][j]     = smartY[ch][j];
            smartScore[other][j] = smartScore[ch][j];
          }
        }
      }
    }
  }

  if (totalC < 2) {
    borderBuffer = Math.min(80000, borderBuffer * 2);
  } else {
    if (borderBuffer === 40 || (rejected > 0 && accepted*100/(accepted+rejected) > 15 && zoom > 0.01)) {
      borderBuffer = 40;
    } else {
      borderBuffer = Math.max(40, (borderBuffer / 1.2) | 0);
    }
  }
}

// ---------------------------------------------------------------------------
// Per-frame computation (port of Java nextpoints)
// ---------------------------------------------------------------------------

function nextPoints() {
  for (var i = 0; i < numPoints; i++) {
    var counter;

    if (byStructure) {
      var data = calcOrbit(i, depthRed, 0);
      var iter = data[0];
      var solColor = data[1];

      var converged = (iter === depthRed && doInverse) || (iter < depthRed && !doInverse);
      if (converged) {
        // Interpolate palette colour
        var baseColor = Math.floor(solColor) | 0;
        var frac = solColor - baseColor;
        var pal = palettes[palIndex % palettes.length];
        var c0 = pal[baseColor] ?? pal[254];
        var c1 = pal[baseColor+1] ?? pal[0];
        var amtR = ((1-frac)*c0[0] + frac*c1[0]) | 0;
        var amtG = ((1-frac)*c0[1] + frac*c1[1]) | 0;
        var amtB = ((1-frac)*c0[2] + frac*c1[2]) | 0;
        counter = drawOrbitStructure(iter, amtR, amtG, amtB);
      } else {
        counter = 0;
      }
      updateScore(counter, 0, i);

    } else {
      // Red channel
      var dataR = calcOrbit(i, depthRed, 0);
      var iterR = dataR[0];
      var convR = (iterR === depthRed && doInverse) || (iterR < depthRed && !doInverse);
      counter = convR ? drawOrbitChannel(iterR, imgRed) : 0;
      updateScore(counter, 0, i);

      // Green channel
      var dataG = calcOrbit(i, depthGreen, 1);
      var iterG = dataG[0];
      var convG = (iterG === depthGreen && doInverse) || (iterG < depthGreen && !doInverse);
      counter = convG ? drawOrbitChannel(iterG, imgGreen) : 0;
      updateScore(counter, 1, i);

      // Blue channel
      var dataB = calcOrbit(i, depthBlue, 2);
      var iterB = dataB[0];
      var convB = (iterB === depthBlue && doInverse) || (iterB < depthBlue && !doInverse);
      counter = convB ? drawOrbitChannel(iterB, imgBlue) : 0;
      updateScore(counter, 2, i);
    }
  }

  runtotal++;

  // Fading passes during early refinement
  if (runtotal > 100 && runtotal <= LOW_RES_PASSES+100 && runtotal % 50 === 40) {
    ensureBlack();
  }

  // Periodic score recalibration
  if (runtotal % 5 === 0) {
    var denom = accepted + rejected;
    var pctAccepted = denom > 0 ? Math.round(accepted * 100 / denom) : 0;
    accepted = rejected = 0;
    calcAvgScore();

    var now = Date.now();
    if (now - lastStatusTime > 1000) {
      updateStatus('Pass: ' + runtotal + '  Accept: ' + pctAccepted + '%  Roots: ' + roots.length);
      lastStatusTime = now;
    }
  }

  if (runtotal % 5000 === 0) resetSmartPoints();
}

// ---------------------------------------------------------------------------
// Periodic fading (port of Java ensureBlack)
// ---------------------------------------------------------------------------

function ensureBlack() {
  var total = ximlen * yimlen;
  for (var idx = 0; idx < total; idx++) {
    var a = imgAlpha[idx];
    if (a > 0) {
      imgRed[idx]   -= imgRed[idx]   / a;
      imgGreen[idx] -= imgGreen[idx] / a;
      imgBlue[idx]  -= imgBlue[idx]  / a;
    }
    imgRed[idx]   /= 1.01;
    imgGreen[idx] /= 1.01;
    imgBlue[idx]  /= 1.01;
    imgAlpha[idx] /= 1.01;
  }
  quadSize = Math.max(0, quadSize - 2);
}

// ---------------------------------------------------------------------------
// Histogram → pixels (port of Java updateHistogram)
// ---------------------------------------------------------------------------

function findPeak(arr) {
  var max = 0;
  for (var i = 0; i < arr.length; i++) {
    if (arr[i] > max) max = arr[i];
  }
  if (gradient <= 0.05) return Math.pow(max, 0.3);
  if (gradient <= 0.1)  return Math.pow(max, 0.4);
  if (gradient <= 0.2)  return Math.pow(max, 0.7);
  if (gradient <= 0.3)  return Math.pow(max, 0.9);
  return max;
}

function updateHistogram() {
  var maxR = findPeak(imgRed);
  var maxG = findPeak(imgGreen);
  var maxB = findPeak(imgBlue);
  var maxA = findPeak(imgAlpha);

  var gmR = Math.pow(maxR, gradient);
  var gmG = Math.pow(maxG, gradient);
  var gmB = Math.pow(maxB, gradient);
  var gmA = Math.pow(maxA, gradient);
  if (gmA === 0) gmA = 1;

  var total = ximlen * yimlen;
  var d = imageData.data;

  for (var idx = 0; idx < total; idx++) {
    var r = 0, g = 0, b = 0;
    var a = imgAlpha[idx];

    if (byStructure) {
      if (a > 0) {
        var z = Math.pow(a, gradient) / gmA;
        r = (imgRed[idx]   * z * brightness / a) | 0;
        g = (imgGreen[idx] * z * brightness / a) | 0;
        b = (imgBlue[idx]  * z * brightness / a) | 0;
      }
    } else {
      if (gmR > 0) r = (brightness * Math.pow(imgRed[idx],   gradient) / gmR * 255) | 0;
      if (gmG > 0) g = (brightness * Math.pow(imgGreen[idx], gradient) / gmG * 255) | 0;
      if (gmB > 0) b = (brightness * Math.pow(imgBlue[idx],  gradient) / gmB * 255) | 0;
    }

    if (r > 255) r = 255;
    if (g > 255) g = 255;
    if (b > 255) b = 255;

    var di = idx * 4;
    d[di]     = r;
    d[di + 1] = g;
    d[di + 2] = b;
    d[di + 3] = 255;
  }

  ctx.putImageData(imageData, 0, 0);
}

// ---------------------------------------------------------------------------
// Render loop
// ---------------------------------------------------------------------------

function renderFrame() {
  nextPoints();

  var now = Date.now();
  if (now - lastUpdateTime > 200) {
    updateHistogram();
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
// Public draw / reset
// ---------------------------------------------------------------------------

function draw(newRoots) {
  stopRendering();
  initCanvas();
  clearAndReset(newRoots !== false);
  startRendering();
}

function resetAndDraw() {
  zoom         = 6;
  xcen         = 0.1;
  ycen         = 0;
  order        = 4.9242613;
  imgOrder     = 0.4841751596941611;
  complexError = 1.0;
  gradient     = 0.11;
  brightness   = 3.0;
  depthRed     = 200;
  depthGreen   = 140;
  depthBlue    = 40;
  numPoints    = 1000;
  algMode      = 0;
  doInverse    = false;
  byStructure  = true;
  mandelbrotAdd = false;
  palIndex     = 0;
  writeControls();
  draw(true);
}

// ---------------------------------------------------------------------------
// Mouse / zoom (same aspect-corrected coordinate mapping as newtonbrot)
// pixel ↔ complex:  re = (px - ximlen/2) * zoom/yimlen - xcen
//                   im = (py - yimlen/2) * zoom/yimlen - ycen
// ---------------------------------------------------------------------------

function setupMouse() {
  ccanvas.oncontextmenu = function(e) { e.preventDefault(); };

  ccanvas.onmousedown = function(e) {
    if (e.button === 2) {
      if (e.ctrlKey && e.shiftKey) {
        mandelbrotAdd = !mandelbrotAdd;
        $('mandelbrotAdd').checked = mandelbrotAdd;
        draw(true);
      } else if (e.ctrlKey) {
        doInverse = !doInverse;
        $('doInverse').checked = doInverse;
        draw(false);
      } else {
        var desc = $('description');
        desc.style.display = (desc.style.display === 'none') ? 'block' : 'none';
      }
      return;
    }
    if (e.button === 0 && e.altKey) {
      algMode = (algMode + 1) % 6;
      $('algMode').value = algMode;
      draw(true);
      return;
    }
    if (e.button === 0 && e.shiftKey) {
      palIndex = (palIndex + 1) % 12;
      $('palIndex').value = palIndex;
      draw(false);
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
      cctx.lineWidth   = 1;
      cctx.strokeRect(dragBox[0], dragBox[1], dragBox[2]-dragBox[0], dragBox[3]-dragBox[1]);
    }
  };

  ccanvas.onmouseup = function(e) {
    if (dragBox && e.button === 0) {
      cctx.clearRect(0, 0, ccanvas.width, ccanvas.height);
      var dx   = Math.abs(dragBox[2] - dragBox[0]);
      var dy   = Math.abs(dragBox[3] - dragBox[1]);
      if (Math.max(dx, dy) >= 10) {
        var midX = (dragBox[0] + dragBox[2]) / 2;
        var midY = (dragBox[1] + dragBox[3]) / 2;
        var scale = zoom / yimlen;
        xcen = xcen - (midX - ximlen * 0.5) * scale;
        ycen = ycen - (midY - yimlen * 0.5) * scale;
        zoom = Math.max(dx / ximlen, dy / yimlen) * zoom;
        writeControls();
        draw(false);
      }
      dragBox = null;
    }
  };

  window.onresize = function() {
    initCanvas();
    draw(false);
  };
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

function main() {
  initCanvas();
  setupMouse();

  $('drawButton').onclick  = function() { draw(true); };
  $('resetButton').onclick = resetAndDraw;
  $('viewPNG').onclick = function() {
    const link = document.createElement('a');
    link.download = `newtonbrot-${runtotal}.png`;
    link.href = canvas.toDataURL('image/png');
    link.click();
  };

  // Re-render on any control change
  const redrawWithRoots  = ['order','imgOrder','complexError','algMode'];
  const redrawNoRoots    = ['zoom','xcen','ycen','numPoints','palIndex'];
  const toggleWithRoots  = ['mandelbrotAdd'];
  const toggleNoRoots    = ['doInverse','byStructure'];
  const noReset = ['depthRed','depthGreen','depthBlue'];

  redrawWithRoots.forEach(function(id) { $(id).onchange = function() { draw(true); }; });
  redrawNoRoots.forEach(function(id)   { $(id).onchange = function() { draw(false); }; });
  toggleWithRoots.forEach(function(id) { $(id).onchange = function() { draw(true); }; });
  toggleNoRoots.forEach(function(id)   { $(id).onchange = function() { draw(false); }; });
  noReset.forEach(function(id)   { $(id).onkeydown = function() { readControls() }; });
  $('contrastSlider').onchange  = function() { readControls(); updateHistogram(); };
  $('brightnessSlider').onchange = function() { readControls(); updateHistogram(); };

  writeControls();
  draw(true);
}

main();
