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
var brightness   = 1.2;
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

    if (isFinite(value.re) && isFinite(value.im) &&
        !isNaN(value.re)   && !isNaN(value.im)) {

      var x = Complex.pow(value, 0.1).abs();
      var y = Complex.pow(value, 0.5).abs();
      var m = Complex.pow(value, 0.9).abs();

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

  const brightness_factor = Math.pow(depth, 2.01 - brightness);
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
    var xi, yi;

    if (i >= totalPixels) {
      // Random anti-aliasing pass
      xi = Math.floor(Math.random() * ximlen);
      yi = Math.floor(Math.random() * yimlen);
    } else {
      xi = i % ximlen;
      yi = Math.floor(i / ximlen);
    }

    var rgb = computePixel(xi, yi);
    var idx  = yi * ximlen + xi;

    img_red[idx]   += rgb[0];
    img_green[idx] += rgb[1];
    img_blue[idx]  += rgb[2];
    img_alpha[idx]++;

    var alpha = img_alpha[idx];
    var r = Math.min(255, img_red[idx]   / alpha) | 0;
    var g = Math.min(255, img_green[idx] / alpha) | 0;
    var b = Math.min(255, img_blue[idx]  / alpha) | 0;

    var di = idx * 4;
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
      updateStatus('Pass: ' + passCount + ' — ' + diff + 's');
      lastPassTime = Date.now();
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
  } while (Date.now() - batchStart < 50);

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
// Public draw / reset
// ---------------------------------------------------------------------------

function draw() {
  stopRendering();
  initCanvas();
  clearAndReset();
  startRendering();
}

function resetAndDraw() {
  zoom         = 9;
  xcen         = 0.0;
  ycen         = 0.01;
  depth        = 20;
  glow         = 0.01;
  brightness   = 1.0;
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
    if (e.button === 1) {
      // Middle click: reset zoom and cycle algorithm
      xcen = 0; ycen = 0; zoom = 4;
      compoundPoint = false;
      $('compoundPoint').checked = false;
      alg = (alg + 1) % 3;
      $('algorithm').value = alg;
      writeControls();
      draw();
      return;
    }

    if (e.button === 2) {
      if (e.ctrlKey) {
        recurse = !recurse;
        $('recurse').checked = recurse;
        draw();
      } else if (e.shiftKey) {
        mandelbrot = !mandelbrot;
        $('mandelbrotMode').checked = mandelbrot;
        draw();
      } else {
        var desc = $('description');
        desc.style.display = (desc.style.display === 'none') ? 'block' : 'none';
      }
      return;
    }

    if (e.button === 0 && e.ctrlKey) {
      scalePoint = !scalePoint;
      $('scalePoint').checked = scalePoint;
      draw();
      return;
    }

    if (e.button === 0 && e.shiftKey) {
      // Shift+left: set compound expansion point
      cx =  ((e.clientX / ximlen) - 0.5) * 2 * zoom * (ximlen / yimlen) + xcen;
      cy = -((e.clientY / yimlen) - 0.5) * 2 * zoom + ycen;  // y flipped (imaginary axis)
      compoundPoint = true;
      $('compoundPoint').checked = true;
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
      cctx.lineWidth   = 1;
      cctx.strokeRect(dragBox[0], dragBox[1],
                      dragBox[2] - dragBox[0],
                      dragBox[3] - dragBox[1]);
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

  $('drawButton').onclick  = draw;
  $('resetButton').onclick = resetAndDraw;

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

  writeControls();
  draw();
}

main();
