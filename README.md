## [Fractal Gallery and Explorers](http://rsweny.github.io/mandelbrot-js)

in HTML5 canvas, javascript and WebGPU.

[Mandelbrot](http://rsweny.github.io/mandelbrot-js/mandelbrot.html) | [Buddhabrot](http://rsweny.github.io/mandelbrot-js/buddhabrot.html) | [Mandelbulb](http://rsweny.github.io/mandelbrot-js/mandelbulb.html) | [Newton](http://rsweny.github.io/mandelbrot-js/newton.html) | [Newtonbrot](http://rsweny.github.io/mandelbrot-js/newtonbrot.html) | [Taylor](http://rsweny.github.io/mandelbrot-js/taylor.html)

### Theory

The famous [Mandelbrot set](http://en.wikipedia.org/wiki/Mandelbrot_set) is
a set of points in the complex plane.  In essence, what we want to find out
is if the iterative function C below will _converge_ to some constant or _diverge_
to infinity.

The function is

  `C_{n+1} = C_{n}^2 + C_{0}`

with the initial condition simply formed by taking the coordinates in the
complex plane,

  `C_{0} = x + iy`

Looking at the function, one can easily see that for big initial values, the
function should diverge.  But for values close to origo (i.e., for |x| and
|y| less than 1), we would expect the function to converge to zero, since
the product of two numbers less than one will always be less than either of
the factors (e.g., 0.5 x 0.4 = 0.2, which is less than both factors).

But if we actually plot it, what we get out isn't any nice plot.  Instead,
we get an amazingly complex and fractured plot.  This is the Mandelbrot set.

You can zoom forever into the plot, and it will present you with an unending
complex shape.  One can also calculate it's
so-called [Hausdorff dimension](http://en.wikipedia.org/wiki/Hausdorff_dimension),
which yields a noninteger number.  Thus, it's a fractal.

### Calculating the Mandelbrot Set

Calculating the Mandelbrot set is easy if you do it numerically.

Take any point `C_0 = (x, y)` and then calculate `C_1 = (x+iy)^2 + (x+iy)`
and continue doing this.  For practical purposes, let's predetermine a
_threshold_ value.  If the magnitude of `C` (defined for complex numbers
as being the distance to origo, or `sqrt(x^2+y^2)`) ever becomes larger than
this threshold value we will assume that it will diverge into infinity.  
If so, stop the calculation and plot a _black dot_ at the current location.

If `C` has not exceeded the threshold value after a predetermined number of
iterations, we will assume that the current parameters makes the function
converge.  In this case, plot a non-black dot at the current location.

### Colorizing the plot

I said above that if the function diverges, one should plot a non-black dot.
One could simply paint a white dot here.  But instead, maybe we want to get
an idea of _how fast_ the function is diverging to infinity at this point.

To do this, just take the current value of the number of steps performed
and _map_ that against a color spectrum, and paint that color.

So, functions diverging quickly will get about the same color.

### Smooth coloring

If you use the number of iterations to pick a color, you'll get ugly color
bands in the plot.  There is a really cool trick to get smooth, gradual
color changes.

So, you basically calculate `Z = Z^2` until it diverges and make a note of
the iteration count.  What we really want, though, is a _fractional_
iteration count, so we can multiply that with a color value to get smooth
colors.

The trick is to note that when you calculate `Z = Z^2` you'll get values `Z,
Z^2, Z^4, Z^8` and so on.  If you take the logarithm of this, you'll get the
values 1, 2, 4, 8 etc.  If you take the logarithm one more time, you'll get
1, 2, 3, 4, 5 etc.  So to get a fractional number of iterations, just do:

    log(log |Z|) / log 2

This is all explained over at http://linas.org/art-gallery/escape/smooth.html

In my code, I originally used the following smoothing equation:

    1 + n - Math.log(Math.log(Math.sqrt(Zr*Zr+Zi*Zi)))/Math.log(2.0);

With some elementary logarithm rules, we can simplify this to

    // Some constants
    var logBase = 1.0 / Math.log(2.0);
    var logHalfBase = Math.log(0.5)*logBase;
    // ...
    return 5 + n - logHalfBase - Math.log(Math.log(Tr+Ti))*logBase;

which is faster.  The constant `5` is another little trick, which should
be explained in the code itself.


### Older Versions

http://cabin.stasis.org/mandelbrot.html | http://cabin.stasis.org/buddhabrot.html | http://cabin.stasis.org/newton.html

