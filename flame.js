/*
 * 3D Flame Fractal in HTML5 Canvas and JavaScript
 * Based on Flame.java by Ryan Sweny (2007)
 * Ported to JavaScript
 *
 * Based on Apophysis: http://sourceforge.net/projects/apophysis/
 */

var NVARS = 27;
var FUSE = 50;
var CHOOSE_XFORM_GRAIN = 1000;
var EPS = 1e-12;

// IFS xforms array (size NVARS - extra slots for symmetry additions)
var xforms = [];

// Current iterated point: [x, y, z, color]
var pt = [0.0, 0.0, 0.0, 0.0];

// Gaussian RNG state (box-muller approximation via 4 uniforms)
var gauss_rnd = [Math.random(), Math.random(), Math.random(), Math.random()];
var gauss_N = 0;

// 3D Camera parameters
var CameraMatrix = [[0,0,0],[0,0,0],[0,0,0]];
var cameraPersp  = 0.9;
var cameraZpos   = 0.0;
var cameraYaw    = 0.3;
var cameraPitch  = 0.86;
var cameraDOF    = 0.01;
var focus        = 0.0;
var opacity      = 0.9;
var translucency = 0.1;
var horizon      = 1.0;
var zres         = 0.96;

// View parameters
var zoom      = 1.0;
var xcen      = 0.0;
var ycen      = 0.0;
var gradient  = 0.3;
var brightness = 1.5;
var ximlen    = 0;
var yimlen    = 0;

// Histogram accumulation buffers
var zPositions = null;
var img_alpha  = null;
var img_red    = null;
var img_green  = null;
var img_blue   = null;

// Render state
var renderpass    = 0;
var max_alpha     = 1;
var visiblePixels = 0;
var occludedPixels = 0;
var running       = false;
var reset         = 0;
var iterCount     = 0;

// Palette
var pallet = null;

// Mouse / UI state
var xcurr, ycurr, xanchor, yanchor;
var m_down   = false;

// Canvas
var canvas, ccanvas, ctx, ctx_img;

// Stats
var t1 = 0;

// --------------------------------------------------------------------------
// Utilities
// --------------------------------------------------------------------------

function $(id) { return document.getElementById(id); }

function createXForm() {
    return {
        vars:     new Array(NVARS).fill(0.0),
        c:        [[0,0],[0,0],[0,0]],
        density:  0.0,
        color:    0.0,
        symmetry: 0.0
    };
}

function matrix2D(rows, cols, val) {
    var arr = new Array(rows);
    for (var i = 0; i < rows; i++) {
        arr[i] = new Array(cols).fill(val);
    }
    return arr;
}

function qsrandom11() { return Math.random() * 2 - 1; }

function randomDistrib(items) {
    return items[Math.floor(Math.random() * items.length) % items.length];
}

// --------------------------------------------------------------------------
// Camera
// --------------------------------------------------------------------------

function setCamera() {
    var cy = Math.cos(-cameraYaw);
    var sy = Math.sin(-cameraYaw);
    var cp = Math.cos(cameraPitch);
    var sp = Math.sin(cameraPitch);

    CameraMatrix[0][0] = cy;        CameraMatrix[1][0] = -sy;        CameraMatrix[2][0] = 0;
    CameraMatrix[0][1] = cp * sy;   CameraMatrix[1][1] = cp * cy;    CameraMatrix[2][1] = -sp;
    CameraMatrix[0][2] = sp * sy;   CameraMatrix[1][2] = sp * cy;    CameraMatrix[2][2] = cp;
}

function projectPitchYawDOF(p) {
    var z_local = p[2] - cameraZpos;
    var x  = CameraMatrix[0][0]*p[0] + CameraMatrix[1][0]*p[1];
    var y  = CameraMatrix[0][1]*p[0] + CameraMatrix[1][1]*p[1] + CameraMatrix[2][1]*z_local;
    var z  = CameraMatrix[0][2]*p[0] + CameraMatrix[1][2]*p[1] + CameraMatrix[2][2]*z_local;
    var zr = 1 - cameraPersp * z;

    var angle = Math.random() * 2 * Math.PI;
    var dsin  = Math.sin(angle);
    var dcos  = Math.cos(angle);
    var dr    = Math.random() * cameraDOF * (z + focus);

    p[0] = (x + dr * dcos) / zr;
    p[1] = (y + dr * dsin) / zr;
    p[2] = z + dr * dsin;
}

// --------------------------------------------------------------------------
// Rendering
// --------------------------------------------------------------------------

function clearScreenAndReset() {
    renderpass     = 0;
    max_alpha      = 1;
    visiblePixels  = 0;
    occludedPixels = 0;
    iterCount      = 0;

    // img_alpha initialised to 1.0 (like Java) so unaccumulated pixels render black
    zPositions = matrix2D(ximlen, yimlen, -100.0);
    img_alpha  = matrix2D(ximlen, yimlen,   1.0);
    img_red    = matrix2D(ximlen, yimlen,   0.0);
    img_green  = matrix2D(ximlen, yimlen,   0.0);
    img_blue   = matrix2D(ximlen, yimlen,   0.0);

    pt[0] = pt[1] = pt[2] = pt[3] = 0.0;
    m_down = false;
}

function findPeak(arr) {
    var maxVal = 0;
    for (var i = 0; i < ximlen; i++) {
        for (var j = 0; j < yimlen; j++) {
            if (arr[i][j] > maxVal) maxVal = arr[i][j];
        }
    }
    return maxVal + 1;
}

function updateHistogram() {
    if (renderpass % 4 === 2) {
        max_alpha = Math.pow(findPeak(img_alpha), gradient);
    }

    var off = 0;
    for (var y = 0; y < yimlen; y++) {
        for (var x = 0; x < ximlen; x++) {
            var alpha = img_alpha[x][y];
            var z     = Math.pow(alpha, gradient) * brightness / max_alpha;
            var red   = Math.floor((img_red[x][y]   * z) / alpha);
            var green = Math.floor((img_green[x][y] * z) / alpha);
            var blue  = Math.floor((img_blue[x][y]  * z) / alpha);

            if (red   > 255) red   = 255;
            if (green > 255) green = 255;
            if (blue  > 255) blue  = 255;
            if (red   < 0)   red   = 0;
            if (green < 0)   green = 0;
            if (blue  < 0)   blue  = 0;

            ctx_img.data[off++] = red;
            ctx_img.data[off++] = green;
            ctx_img.data[off++] = blue;
            ctx_img.data[off++] = 255;
        }
    }
    ctx.putImageData(ctx_img, 0, 0);
}

// --------------------------------------------------------------------------
// Core IFS iteration
// --------------------------------------------------------------------------

function nextpoints(numIter) {
    var i, j, k, fn;
    var a, c1, c2, tx, ty, tz, v, nx, ny, r, r2;

    // Build weighted xform distribution
    var xform_distrib = new Array(CHOOSE_XFORM_GRAIN);
    var dr = 0.0;
    for (j = 0; j < NVARS; j++) dr += xforms[j].density;
    dr = dr / CHOOSE_XFORM_GRAIN;
    j = 0;
    var t = xforms[0].density;
    r = 0.0;
    for (k = 0; k < CHOOSE_XFORM_GRAIN; k++) {
        while (r >= t) {
            j++;
            if (j >= NVARS) { j = NVARS - 1; break; }
            t += xforms[j].density;
        }
        xform_distrib[k] = j;
        r += dr;
    }

    // Per-batch curl constants (Curl3D variation)
    var curlx  = Math.random();
    var curly  = Math.random();
    var curlz  = Math.random();
    var curlc2 = curlx*curlx + curly*curly + curlz*curlz + EPS;

    var palLen = pallet.length;
    var half   = ximlen / 2;

    for (i = -FUSE; i < numIter; i++) {

        // Explosion guard
        if (Math.abs(pt[0]) > 1e6 || Math.abs(pt[1]) > 1e6 || Math.abs(pt[2]) > 1e6) {
            pt[0] = pt[1] = pt[2] = 0.0;
        }

        fn = xform_distrib[Math.floor(Math.random() * CHOOSE_XFORM_GRAIN) % CHOOSE_XFORM_GRAIN];
        var xf = xforms[fn];

        // Blend color coordinate
        var s = xf.symmetry;
        pt[3] = (pt[3] + xf.color) * 0.5 * (1.0 - s) + s * pt[3];

        // Affine pre-transform
        tx = xf.c[0][0] * pt[0] + xf.c[1][0] * pt[1] + xf.c[2][0];
        ty = xf.c[0][1] * pt[0] + xf.c[1][1] * pt[1] + xf.c[2][1];
        tz = pt[2];

        pt[0] = pt[1] = pt[2] = 0.0;

        // ---- VAR 26: PRE-BLUR ----
        v = xf.vars[26];
        if (v !== 0.0) {
            var ang26  = Math.random() * 2 * Math.PI;
            var mag26  = v * (gauss_rnd[0] + gauss_rnd[1] + gauss_rnd[2] + gauss_rnd[3] - 2);
            gauss_rnd[gauss_N] = Math.random();
            gauss_N = (gauss_N + 1) % 4;
            tx += mag26 * Math.cos(ang26);
            ty += mag26 * Math.sin(ang26);
        }

        // ---- VAR 0: LINEAR ----
        v = xf.vars[0];
        if (v > 0.0) {
            pt[0] += v * tx;
            pt[1] += v * ty;
            pt[2] += v * tz;
        }

        // ---- VAR 1: SINUSOIDAL ----
        v = xf.vars[1];
        if (v > 0.0) {
            pt[0] += v * Math.sin(tx);
            pt[1] += v * Math.sin(ty);
        }

        // ---- VAR 2: SPHERICAL ----
        v = xf.vars[2];
        if (v > 0.0) {
            r2 = tx*tx + ty*ty + 1e-6;
            pt[0] += v * tx / r2;
            pt[1] += v * ty / r2;
        }

        // ---- VAR 3: SWIRL ----
        v = xf.vars[3];
        if (v > 0.0) {
            r2 = tx*tx + ty*ty;
            c1 = Math.sin(r2);
            c2 = Math.cos(r2);
            pt[0] += v * (c1*tx - c2*ty);
            pt[1] += v * (c2*tx + c1*ty);
        }

        // ---- VAR 4: HORSESHOE ----
        v = xf.vars[4];
        if (v > 0.0) {
            a  = (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS) ? Math.atan2(tx, ty) : 0.0;
            c1 = Math.sin(a);
            c2 = Math.cos(a);
            pt[0] += v * (c1*tx - c2*ty);
            pt[1] += v * (c2*tx + c1*ty);
        }

        // ---- VAR 5: POLAR ----
        v = xf.vars[5];
        if (v > 0.0) {
            nx = (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS) ? Math.atan2(tx, ty) / Math.PI : 0.0;
            ny = Math.sqrt(tx*tx + ty*ty) - 1.0;
            pt[0] += v * nx;
            pt[1] += v * ny;
        }

        // ---- VAR 6: HANDKERCHIEF ----
        v = xf.vars[6];
        if (v > 0.0) {
            a  = (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS) ? Math.atan2(tx, ty) : 0.0;
            r2 = Math.sqrt(tx*tx + ty*ty);
            pt[0] += v * Math.sin(a + r2) * r2;
            pt[1] += v * Math.cos(a - r2) * r2;
        }

        // ---- VAR 7: HEART ----
        v = xf.vars[7];
        if (v > 0.0) {
            a  = (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS) ? Math.atan2(tx, ty) : 0.0;
            r2 = Math.sqrt(tx*tx + ty*ty);
            a  = a * r2;
            pt[0] += v *  Math.sin(a) * r2;
            pt[1] += v * -Math.cos(a) * r2;
        }

        // ---- VAR 8: DISC ----
        v = xf.vars[8];
        if (v > 0.0) {
            nx = tx * Math.PI;
            ny = ty * Math.PI;
            a  = (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS) ? Math.atan2(nx, ny) : 0.0;
            r2 = Math.sqrt(nx*nx + ny*ny);
            pt[0] += v * Math.sin(r2) * a / Math.PI;
            pt[1] += v * Math.cos(r2) * a / Math.PI;
        }

        // ---- VAR 9: SPIRAL ----
        v = xf.vars[9];
        if (v > 0.0) {
            a  = (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS) ? Math.atan2(tx, ty) : 0.0;
            r2 = Math.sqrt(tx*tx + ty*ty) + 1e-6;
            pt[0] += v * (Math.cos(a) + Math.sin(r2)) / r2;
            pt[1] += v * (Math.sin(a) - Math.cos(r2)) / r2;
        }

        // ---- VAR 10: HYPERBOLIC ----
        v = xf.vars[10];
        if (v > 0.0) {
            a  = (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS) ? Math.atan2(tx, ty) : 0.0;
            r2 = Math.sqrt(tx*tx + ty*ty) + 1e-6;
            pt[0] += v * Math.sin(a) / r2;
            pt[1] += v * Math.cos(a) * r2;
        }

        // ---- VAR 11: DIAMOND ----
        v = xf.vars[11];
        if (v > 0.0) {
            a  = (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS) ? Math.atan2(tx, ty) : 0.0;
            r2 = Math.sqrt(tx*tx + ty*ty);
            pt[0] += v * Math.sin(a) * Math.cos(r2);
            pt[1] += v * Math.cos(a) * Math.sin(r2);
        }

        // ---- VAR 12: EX ----
        v = xf.vars[12];
        if (v > 0.0) {
            a  = (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS) ? Math.atan2(tx, ty) : 0.0;
            r2 = Math.sqrt(tx*tx + ty*ty);
            var n0 = Math.sin(a + r2);
            var n1 = Math.cos(a - r2);
            var m0 = n0 * n0 * n0 * r2;
            var m1 = n1 * n1 * n1 * r2;
            pt[0] += v * (m0 + m1);
            pt[1] += v * (m0 - m1);
        }

        // ---- VAR 13: JULIA ----
        v = xf.vars[13];
        if (v > 0.0) {
            a  = (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS) ? Math.atan2(tx, ty) / 2.0 : 0.0;
            if (Math.random() > 0.5) a += Math.PI;
            r2 = Math.pow(tx*tx + ty*ty, 0.25);
            pt[0] += v * r2 * Math.cos(a);
            pt[1] += v * r2 * Math.sin(a);
        }

        // ---- VAR 14: BENT ----
        v = xf.vars[14];
        if (v > 0.0) {
            nx = (tx < 0.0) ? tx * 2.0 : tx;
            ny = (ty < 0.0) ? ty / 2.0 : ty;
            pt[0] += v * nx;
            pt[1] += v * ny;
        }

        // ---- VAR 15: WAVES ----
        v = xf.vars[15];
        if (v !== 0.0) {
            var dx15 = xf.c[2][0];
            var dy15 = xf.c[2][1];
            pt[0] += v * (tx + xf.c[1][0] * Math.sin(ty / (dx15*dx15 + EPS)));
            pt[1] += v * (ty + xf.c[1][1] * Math.sin(tx / (dy15*dy15 + EPS)));
        }

        // ---- VAR 16: FISHEYE ----
        v = xf.vars[16];
        if (v !== 0.0) {
            r  = Math.sqrt(tx*tx + ty*ty);
            a  = Math.atan2(tx, ty);
            r  = 2 * r / (r + 1);
            pt[0] += v * r * Math.cos(a);
            pt[1] += v * r * Math.sin(a);
        }

        // ---- VAR 17: POPCORN ----
        v = xf.vars[17];
        if (v !== 0.0) {
            pt[0] += v * (tx + xf.c[2][0] * Math.sin(Math.tan(3 * ty)));
            pt[1] += v * (ty + xf.c[2][1] * Math.sin(Math.tan(3 * tx)));
        }

        // ---- 3D VARIATIONS ----

        // VAR 18: Z TRANSLATE
        v = xf.vars[18];
        if (v !== 0.0) {
            pt[2] += v;
        }

        // VAR 19: Z SCALE
        v = xf.vars[19];
        if (v !== 0.0) {
            pt[2] = pt[2] + pt[2] * v;
        }

        // VAR 20: BLUR 3D
        v = xf.vars[20];
        if (v !== 0.0) {
            var ang20  = Math.random() * 2 * Math.PI;
            var mag20  = v * (gauss_rnd[0] + gauss_rnd[1] + gauss_rnd[2] + gauss_rnd[3] - 2);
            gauss_rnd[gauss_N] = Math.random();
            gauss_N = (gauss_N + 1) % 4;
            var ang20b = Math.random() * Math.PI;
            var sinb20 = Math.sin(ang20b);
            pt[0] += mag20 * sinb20 * Math.cos(ang20);
            pt[1] += mag20 * sinb20 * Math.sin(ang20);
            pt[2] += mag20 * Math.cos(ang20b);
        }

        // VAR 21: Z BLUR
        v = xf.vars[21];
        if (v !== 0.0) {
            pt[2] = pt[2] + v * (gauss_rnd[0] + gauss_rnd[1] + gauss_rnd[2] + gauss_rnd[3] - 2);
            gauss_rnd[gauss_N] = Math.random();
            gauss_N = (gauss_N + 1) % 4;
        }

        // VAR 22: Z CONE
        v = xf.vars[22];
        if (v !== 0.0) {
            pt[2] += v * Math.sqrt(tx*tx + ty*ty);
        }

        // VAR 23: BUBBLE
        v = xf.vars[23];
        if (v !== 0.0) {
            var mag23 = (tx*tx + ty*ty) / 4 + 1;
            pt[2] = pt[2] + v * (2 / mag23 - 1);
            var bfac = v / mag23;
            pt[0] = pt[0] + bfac * tx;
            pt[1] = pt[1] + bfac * ty;
        }

        // VAR 24: JULIA 3D
        v = xf.vars[24];
        if (v !== 0.0) {
            var r2d24  = Math.sqrt(tx*tx + ty*ty);
            var mag24  = v * Math.sqrt(r2d24);
            pt[2] = pt[2] + mag24 * tz / (r2d24 + EPS) / 2;
            var rnd24  = Math.floor(Math.random() * 2);
            var ang24  = Math.atan2(ty, tx) / 2 + Math.PI * rnd24;
            pt[0] = pt[0] + mag24 * Math.cos(ang24);
            pt[1] = pt[1] + mag24 * Math.sin(ang24);
        }

        // VAR 25: CURL 3D
        v = xf.vars[25];
        if (v !== 0.0) {
            var cr2  = tx*tx + ty*ty + tz*tz;
            var cmag = v / (cr2*curlc2 + curlx*2*tx - curly*2*ty + curlz*2*tz + 1);
            pt[0] += cmag * (tx + curlx * cr2);
            pt[1] += cmag * (ty - curly * cr2);
            pt[2] += cmag * (tz + curlz * cr2);
        }

        // Skip fuse burn-in iterations
        if (i < 0) continue;

        // Project through camera
        var proj = [pt[0], pt[1], pt[2], pt[3]];
        projectPitchYawDOF(proj);

        var realx = proj[0] + xcen;
        var realy = proj[1] + ycen;

        if (realx > -zoom && realy > -zoom && realx < zoom && realy < zoom) {
            var tempx = Math.floor((realx / zoom) * half + half);
            var tempy = Math.floor((realy / zoom) * half + (yimlen / 2));

            if (tempx < 0 || tempx >= ximlen || tempy < 0 || tempy >= yimlen) continue;

            var horizon_dist = proj[2] + Math.random() / 2;
            if (horizon_dist < -horizon) continue;

            var colorIndex = Math.floor(Math.abs(pt[3])) % palLen;

            if (proj[2] > zPositions[tempx][tempy]) {
                // Visible – overwrite depth buffer
                zPositions[tempx][tempy] = proj[2];
                img_red[tempx][tempy]   += pallet[colorIndex][0] * opacity;
                img_green[tempx][tempy] += pallet[colorIndex][1] * opacity;
                img_blue[tempx][tempy]  += pallet[colorIndex][2] * opacity;
                img_alpha[tempx][tempy] += opacity;
                visiblePixels++;
            } else {
                // Occluded – add with reduced weight
                var totalColor = pallet[colorIndex][0] + pallet[colorIndex][1] + pallet[colorIndex][2];
                var intensity  = totalColor / 765.0;
                var occFac     = 1.0 - opacity + translucency;
                zPositions[tempx][tempy] = zPositions[tempx][tempy]*zres + proj[2]*(1.0 - zres);
                img_red[tempx][tempy]   += pallet[colorIndex][0] * occFac;
                img_green[tempx][tempy] += pallet[colorIndex][1] * occFac;
                img_blue[tempx][tempy]  += pallet[colorIndex][2] * occFac;
                img_alpha[tempx][tempy] += (1.0 - opacity + translucency/2) * intensity;
                occludedPixels++;
            }
            iterCount++;
        }
    }
}

// --------------------------------------------------------------------------
// Formula parsing
// --------------------------------------------------------------------------

function loadFormula(xform, formulaStr) {
    try {
        var tokens  = formulaStr.trim().split(/\s+/).filter(function(s) { return s.length > 0; });
        var values  = tokens.map(parseFloat);

        xform.density  = values[0];
        xform.symmetry = values[1];
        xform.color    = values[2];

        var counter = 0;
        for (var i = 0; i < 3; i++) {
            for (var j = 0; j < 2; j++) {
                xform.c[i][j] = values[3 + counter];
                counter++;
            }
        }
        for (var idx = 9; idx < values.length && (idx - 9) < NVARS; idx++) {
            xform.vars[idx - 9] = values[idx];
        }
    } catch(e) {
        console.log("Problem parsing xform: " + e);
    }
}

function parseFormula(s) {
    for (var i = 0; i < NVARS; i++) xforms[i].density = 0.0;
    var parts = s.split('|');
    var idx   = 0;
    for (var p = 0; p < parts.length; p++) {
        var part = parts[p].trim();
        if (part.length > 0) {
            if (idx >= NVARS) break;
            loadFormula(xforms[idx], part);
            idx++;
        }
    }
}

function buildFormulaText() {
    var sb = '';
    for (var i = 0; i < NVARS; i++) {
        if (xforms[i].density > 0) {
            sb += xforms[i].density + ' ' + xforms[i].symmetry + ' ' + xforms[i].color + '\n';
            for (var rr = 0; rr < 3; rr++)
                for (var cc = 0; cc < 2; cc++)
                    sb += xforms[i].c[rr][cc] + ' ';
            sb += '\n';
            for (var j = 0; j < NVARS; j++) {
                sb += xforms[i].vars[j] + ' ';
                if (j === 17) sb += '\n';
            }
            sb += '\n|\n';
        }
    }
    return sb;
}

// --------------------------------------------------------------------------
// Random generation
// --------------------------------------------------------------------------

function randomXform(ivar) {
    var xform_dist   = [1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6];
    var var_dist     = [24, 24, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17];
    var mixed_dist   = [24, 24, 0, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 14, 14, 15, 16, 17];

    var active  = randomDistrib(xform_dist);
    var varType = (ivar < 0) ? randomDistrib(var_dist) : ivar;

    for (var i = 0; i < NVARS; i++) xforms[i].density = 0.0;

    for (i = 0; i < active; i++) {
        xforms[i].symmetry = 0.0;
        xforms[i].density  = 1.0 / active;
        xforms[i].color    = (i / active) * 255;
        for (var j = 0; j < 3; j++)
            for (var k = 0; k < 2; k++)
                xforms[i].c[j][k] = qsrandom11();
        for (j = 0; j < NVARS; j++) xforms[i].vars[j] = 0.0;
        if (varType >= 0) xforms[i].vars[varType] = 1.0;
        else              xforms[i].vars[randomDistrib(mixed_dist)] = 1.0;
    }
}

function randomizeC() {
    for (var i = 0; i < 2; i++)
        for (var j = 0; j < 3; j++)
            for (var k = 0; k < 2; k++)
                xforms[i].c[j][k] = qsrandom11();
}

// --------------------------------------------------------------------------
// UI constants sync
// --------------------------------------------------------------------------

function getConstants() {
    zoom         = parseFloat($('zoomInput').value);
    xcen         = parseFloat($('xcenInput').value);
    ycen         = parseFloat($('ycenInput').value);
    cameraPersp  = parseFloat($('perspInput').value);
    cameraZpos   = parseFloat($('zposInput').value);
    cameraYaw    = parseFloat($('yawInput').value);
    cameraPitch  = parseFloat($('pitchInput').value);
    cameraDOF    = parseFloat($('dofInput').value);
    focus        = parseFloat($('focusInput').value);
    opacity      = parseFloat($('opacityInput').value);
    translucency = parseFloat($('glowInput').value);
    horizon      = parseFloat($('horizonInput').value);
    zres         = parseFloat($('zresInput').value);
    setCamera();
    parseFormula($('formulaInput').value);
}

function setConstants() {
    $('zoomInput').value      = zoom;
    $('xcenInput').value      = xcen;
    $('ycenInput').value      = ycen;
    $('perspInput').value     = cameraPersp;
    $('zposInput').value      = cameraZpos;
    $('yawInput').value       = cameraYaw;
    $('pitchInput').value     = cameraPitch;
    $('dofInput').value       = cameraDOF;
    $('focusInput').value     = focus;
    $('opacityInput').value   = opacity;
    $('glowInput').value      = translucency;
    $('horizonInput').value   = horizon;
    $('zresInput').value      = zres;
    $('formulaInput').value   = buildFormulaText();
}

// --------------------------------------------------------------------------
// Default formulas (from Flame.java loadDefault)
// --------------------------------------------------------------------------

var DEFAULT_JULIA = [
    "0.31 0.0 30.0 1.0 0.1 0.1 1.0 0.1 0.1",
    "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.2 0.0 0.19",
    "|",
    "5.1 0.4 0.0 0.452548 -0.452548 0.452548 0.452548 0.19 0.1",
    "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.1 0.0 0.0 0.0 0.0 0.0 0.9 0.0 0.0",
    "|",
    "5.8 0.1 255.0 1.0 0.0 0.0 1.0 0.0 0.0",
    "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.5 0.0 0.0",
    "|"
].join(" ");

var DEFAULT_SPIRAL = [
    "0.25 0.0 0.0 -0.0454 -0.606852 0.620299 -0.212502 0.784553 -0.258353",
    "0.0 0.5 0.0 0.0 0.01 0.0 0.0 0.0 0.0 0.0 0.0 0.49 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0",
    "|",
    "0.25 0.0 128.0 -0.238986 -0.147846 0.147846 -0.238986 -0.07213 -0.801111",
    "1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0",
    "|",
    "0.5 1.0 255.0 0.687624 -0.726067 0.726067 0.687624 0.038586 0.584278",
    "1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0",
    "|"
].join(" ");

function loadDefault() {
    if (Math.random() > 0.5)
        parseFormula(DEFAULT_JULIA);
    else
        parseFormula(DEFAULT_SPIRAL);
}

// --------------------------------------------------------------------------
// Render loop
// --------------------------------------------------------------------------

function render() {
    if (!running) return;

    if (reset === 1) {
        reset = 0;
        getConstants();
        clearScreenAndReset();
    } else if (reset === 2) {
        reset = 0;
        getConstants();
        updateHistogram();
    }

    if (!m_down) {
        nextpoints(150000);
        renderpass++;

        if (renderpass % 2 === 0) {
            updateHistogram();

            var now     = (new Date()).getTime();
            var elapsed = (now - t1) / 1000.0;
            $('renderTime').innerHTML  = elapsed.toFixed(1) + "s  pass: " + renderpass;
            $('renderSpeed').innerHTML = Math.floor(iterCount / elapsed) + " iter/s";
        }
    }

    setTimeout(render, 1);
}

// --------------------------------------------------------------------------
// Init
// --------------------------------------------------------------------------

function init() {
    canvas  = $('canvasFlame');
    ccanvas = $('canvasControls');
    canvas.width  = 1200;
    canvas.height = 1200;
    ccanvas.width  = 1200;
    ccanvas.height = 1200;
    ctx     = canvas.getContext('2d');
    ctx_img = ctx.createImageData(canvas.width, canvas.height);

    ximlen  = canvas.width;
    yimlen  = canvas.height;

    for (var i = 0; i < NVARS; i++) xforms[i] = createXForm();

    pallet = palettes["Geyser"];

    loadDefault();
    setCamera();
    clearScreenAndReset();
    setConstants();
}

// --------------------------------------------------------------------------
// Main
// --------------------------------------------------------------------------

function main() {
    init();

    $('viewPNG').onclick = function() {
        var link = document.createElement('a');
        link.download = 'flame.png';
        link.href = canvas.toDataURL('image/png');
        link.click();
    };

    $('randomBtn').onclick = function() {
        randomXform(-1);
        setConstants();
        reset = 1;
    };

    $('randomCBtn').onclick = function() {
        randomizeC();
        setConstants();
        reset = 1;
    };

    $('settingsForm').onsubmit = function(e) {
        e.preventDefault();
        reset = 1;
    };

    $('contrastSlider').oninput = function() {
        gradient = parseFloat((this.value / 100.0).toFixed(2));
        reset = 2;
    };

    $('brightnessSlider').oninput = function() {
        brightness = parseFloat((this.value / 100.0).toFixed(2));
        reset = 2;
    };

    $('colorPalette').onchange = function() {
        pallet = palettes[this.value];
        reset = 1;
    };

    $('zoomInput').onchange     = function() { zoom    = parseFloat(this.value); reset = 1; };
    $('xcenInput').onchange     = function() { xcen    = parseFloat(this.value); reset = 1; };
    $('ycenInput').onchange     = function() { ycen    = parseFloat(this.value); reset = 1; };
    $('perspInput').onchange    = function() { cameraPersp  = parseFloat(this.value); setCamera(); reset = 1; };
    $('zposInput').onchange     = function() { cameraZpos   = parseFloat(this.value); setCamera(); reset = 1; };
    $('yawInput').onchange      = function() { cameraYaw    = parseFloat(this.value); setCamera(); reset = 1; };
    $('pitchInput').onchange    = function() { cameraPitch  = parseFloat(this.value); setCamera(); reset = 1; };
    $('dofInput').onchange      = function() { cameraDOF    = parseFloat(this.value); reset = 1; };
    $('focusInput').onchange    = function() { focus        = parseFloat(this.value); reset = 1; };
    $('opacityInput').onchange  = function() { opacity      = parseFloat(this.value); reset = 1; };
    $('glowInput').onchange     = function() { translucency = parseFloat(this.value); reset = 1; };
    $('horizonInput').onchange  = function() { horizon      = parseFloat(this.value); reset = 1; };
    $('zresInput').onchange     = function() { zres         = parseFloat(this.value); reset = 1; };

    // Mouse: click+drag to zoom
    $('canvasControls').onmousedown = function(e) {
        m_down  = true;
        xanchor = e.clientX - this.offsetLeft;
        yanchor = e.clientY - this.offsetTop;
    };

    $('canvasControls').onmousemove = function(e) {
        if (!m_down) return;
        var c  = ccanvas.getContext('2d');
        c.clearRect(0, 0, ccanvas.width, ccanvas.height);
        xcurr  = e.clientX - this.offsetLeft;
        ycurr  = e.clientY - this.offsetTop;
        var dx = Math.abs(xcurr - xanchor);
        var dy = Math.abs(ycurr - yanchor);
        c.strokeStyle = '#FF3B03';
        c.lineWidth   = 1;
        c.strokeRect(xanchor, yanchor, dx, dy);
    };

    $('canvasControls').onmouseup = function(e) {
        var c = ccanvas.getContext('2d');
        c.clearRect(0, 0, ccanvas.width, ccanvas.height);
        m_down = false;
        xcurr  = e.clientX - this.offsetLeft;
        ycurr  = e.clientY - this.offsetTop;

        var dx = Math.abs(xcurr - xanchor);
        var dy = Math.abs(ycurr - yanchor);
        if (dy > dx) dx = dy;

        if (dx > 10) {
            var newxcen = ((xanchor + dx/2 - ximlen/2) / ximlen) * zoom;
            var newycen = ((yanchor + dx/2 - yimlen/2) / yimlen) * zoom;
            xcen = xcen - newxcen;
            ycen = ycen - newycen;
            zoom = (dx / ximlen) * zoom;
            $('zoomInput').value = zoom.toFixed(4);
            $('xcenInput').value = xcen.toFixed(4);
            $('ycenInput').value = ycen.toFixed(4);
            reset = 1;
        }
    };

    // Start rendering
    t1      = (new Date()).getTime();
    running = true;
    render();
}

main();


// portals
/*
0.31 0 30
1 0.1 0.1 1 0.1 0.1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0.2 0 0.19 
|
5.1 0.4 0
0.452548 -0.452548 0.452548 0.452548 0.19 0.1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
-0.1 0 0 0 0 0 0.9 0 0 
|
5.8 0.1 255
1 0 0 1 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0.5 0 0 
|
*/
