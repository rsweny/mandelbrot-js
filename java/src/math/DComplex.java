package math;

/**
*  This class represents double complex type. 
*  @author Herman Harjono.
*  @version Oct 5, 1998.
*/
public final class DComplex extends Number implements Cloneable//extends java.lang.Number
{
  /** 
  * Real part
  */
  private double re;

  /**
  *  Imaginary part
  */
  private double im;
  
  /**
  *  Returns the eucledian distance.
  */
  private static double hypot(double x, double y)
  {
    return Math.sqrt(x*x + y*y);
  }

  /**
  *  Rectangular to Polar coordinate conversion.
  */
  private DComplex rtop() 
  { return new DComplex(hypot(re, im), Math.atan2(im, re));}
  
  /**
  *  Polar to Rectangular coordinate conversion. 
  */
  private DComplex ptor() 
    {return new DComplex(re*Math.cos(im), re*Math.sin(im));}
  
  /**
  *   Create a DComplex from two doubles.
  */
  public DComplex(double r, double i)
  {re = r; im = i;}

  /**
  *   Create a DComplex (r,0)
  */  
  public DComplex(double r)
  {this(r,0);}
  
  /**
  *   Create a DComplex (0,0)
  */
  public DComplex()
  {this(0,0);}
  
  /**
  *  Create a DComplex whose members equal d
  */
  public DComplex(DComplex d)
  { re = d.re; im = d.im;}  

  /**
  *   Set this to (re, im).
  */
  public void set(double re, double im)
  {
    this.re = re; this.im = im;
  }

  /**
  *  Returns the real part.
  */
  public double real()
  {return re;}
  
  /**
  *  Returns the imaginary part.
  */
  public double aimag()
  {return im;}
 
  /**
  *  Returns the eucledian distance of the (real, imaginary).
  */
  public double abs()
  { return hypot(re, im);}
  
  /**
  *  Returns the atan2(imaginary, real)
  */
  public double arg()
  { return Math.atan2(im, re);}
  
  /**
  *   Returns the DComplex conjugate of this.
  */
  public DComplex conjg()
  {return new DComplex(re, -im);}
  
  /**
  *   Returns the exponent of this.
  */
  public DComplex exp( DComplex d)
  {
    double r = Math.exp(re);
    return new DComplex(r*Math.cos(im), r*Math.sin(im));
  }
  
  /**
  *  Returns the value of this raised to r.
  */
  public DComplex pow(double r)
  {
    DComplex polar = this.rtop();
    polar.re = Math.pow(polar.re, r);
    polar.im *= r;
    return polar.ptor();
  }

  /**
  *  Returns the value of this raised to r.
  */
  public DComplex power(int r)
  {
    switch(r)
    {
      
    case 0: return new DComplex(1);
    case 1: return new DComplex(this);
    case 2: return this.multiply(this);
    default: 
      DComplex polar = this.rtop();
      polar.re = Math.pow(polar.re, r);
      polar.im *= r;
      return polar.ptor();
    }
  }
  
  /**
  *  Returns sqrt(this)
  */
  public DComplex sqrt()
  {
    DComplex polar = this.rtop();
    polar.re = Math.sqrt(polar.re);
    polar.im = 0.5*polar.im;
    return polar.ptor();
  }
  
  
  /**
  *  Returns a DComplex whose value is this + d2
  */
  public DComplex add( DComplex d2)
  {
    return new DComplex(re+d2.re, im+d2.im);
  }
 
  /**
  *  Returns a DComplex whose value is (-real, -imaginary)
  */
  public DComplex negate()
  {
    return new DComplex(-re, -im);
  }
  
  /**
  *  Returns a DComplex whose value is this - d2.
  */
  public DComplex subtract(DComplex d2)
  {
    return new DComplex(re-d2.re, im-d2.im);
  }
  
  /**
  *  Returns a DComplex whose value is this * d2.
  */
  public DComplex multiply(DComplex d2)
  {
    return new DComplex(re*d2.re-im*d2.im, re*d2.im+im*d2.re);
  }
  
  /**
  *  Returns a DComplex whose value is this / d2.
  */
  public DComplex divide(DComplex d2)
  {
    double denom = d2.re*d2.re+d2.im*d2.im;
    return new DComplex((re*d2.re+im*d2.im)/denom,
      (im*d2.re-re*d2.im)/denom);
  }
 
  /**
  *  this += d.
  *  @return this
  */
  public DComplex addEqual(DComplex d)
  {
    re += d.re; im += d.im;
    return this;
  }

  /**
  *   this -= d.
  *  @return this
  */
  public DComplex subtractEqual( DComplex d)
  {
    re -= d.re; im -= d.im;
    return this;
  }
  
  /**
  *   this *= d
  *  @return this
  */
  public DComplex multiplyEqual( DComplex d)
  {
    double oldRe = re;
    double oldDRe = d.re; // important if this == &d
    re = re*d.re - im*d.im; 
    im = im*oldDRe + oldRe*d.im;
    return this;
   }

  /**
  *   this /= d
  *  @return this
  */
  public DComplex divideEqual( DComplex d)
  {
    double denom = d.re*d.re+d.im*d.im;
    double oldRe = re;
    double oldDRe = d.re; // important if this == &d
    re = (re*d.re + im*d.im)/denom;
    im = (im*oldDRe - oldRe*d.im)/denom;
    return this;
  }
  
  /**
  *  Compares the equality of this and d1.
  */
  public boolean equals(Object d1)
  {
    boolean result;
    try
    {
      result = (re==((DComplex)d1).re 
        && im==((DComplex)d1).im);
    }
    catch(ClassCastException e)
    {
      result = false;
    }
    return result;
  }
  
  /**
  *   Return the clone of this.
  */
  public Object clone() 
  {
    return new DComplex(re,im);
  }

  /**
  *   Return the string representation of this.
  */
  public String toString()
  {
    return  "(" + new Double(re) + "," + new Double(im) + ")";
  }

  public int hashCode()
  {
    return (new Double(re/2+im/2)).hashCode();
  }

  /**
  *  Returns the real as a byte. 
  *
  *  This may involves rounding or truncation.  
  */
  public byte byteValue() 
  {
    return (byte) re;
  }
  
  /**
  *  Returns the real as a double. 
  */    
  public double doubleValue() 
  {
    return re;
  }
  /**
  *  Returns the real as a float. 
  *
  *  This may involves rounding or truncation.  
  */
  public float floatValue() 
  {
    return (float) re;
  }

  /**
  *  Returns the real as an int. 
  *
  *  This may involves rounding or truncation.  
  */  
  public int intValue() 
  {
    return (int) re;
  }

 /**
  *  Returns the real as a long. 
  *
  *  This may involves rounding or truncation.  
  */
  public long longValue() 
  {
    return (long) re;
  }

 /**
  *  Returns the real as a short. 
  *
  *  This may involves rounding or truncation.  
  */
  public short shortValue() 
  {
    return (short) re;
  }
}

