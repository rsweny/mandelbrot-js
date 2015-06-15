package math;


/**  
 * This class represents multi-precision real and is intended to be directly
 * used by the user.
 *
 * While some operations (like comparisons) are inherited from MP,
 * this class implements other operations that are specific to MPReal.
 * In addition to the arithmetic operators, MPReal also provides a collection 
 * of transcedental functions, such as: cos, sin, etc.
 *
 * If a MPReal object needs a lower precision level than the current maximum 
 * precision level (mpipl), the user can pass a precision parameter to the 
 * constructor. Otherwise, the current maximum precision level is used.
 * In order to change the maximum precision level, the user should use the
 * accessor functions defined in MPGlobal
 *
 * All intermediate computations are done in maximum precision level. Only 
 * the final result is truncated (if necessary) to the precision level of the 
 * object.
 *
 * @author Herman Harjono
 * @version Oct 7, 1998.
 */ 
public final class MPReal extends MP
{
 /**
  * Log 2.
  */
  public static MPReal mpl02; 
  
  /**
  * Log 10.
  */  
  public static MPReal mpl10;
  
  /**
  * Pi.
  */
  public static MPReal mppic;
  
  /**
  * Current MP epsilon value. 10^MPIEP.
  */
  public static MPReal mpeps;

  static
  {
    // initializing pi
    mppic = new MPReal(mp21,false);
    MP.mppix(mppic, mpnw+1);

    // initialing l02
    mpl02 = new MPReal(mp21,false);
    MP t2 = new MP(6,false);
    MP.mpdmc(new MPDPE(2), t2);
    MP.mplogx(t2, mppic, mpl02, mpl02, mpnw+1);
  
    // initializing l10
    mpl10 = new MPReal(mp21,false);
    MP.mpdmc(new MPDPE(10), t2);
    MP.mplogx(t2, mppic, mpl02, mpl10, mpnw+1);
 
    // initializing eps
    mpeps = new MPReal(mp21,false);
    MP.mpdmc(new MPDPE(10), t2);
    MP.mpnpwx(t2, mpiep, mpeps, mpnw+1);
  
    // last words is not reliable
    mppic.nw--; mpl02.nw--;mpl10.nw--;mpeps.nw--;

  }

 /**
  *  Creates a MPReal (initialized to zero) with default precision level. 
  *  @see MPGlobal#mpipl
  */
  public MPReal()
  {super(true, mpipl);}
  
  /**
  *  Creates a MPReal (initialized to zero). 
  *  @param precision The precision level of the object, in digits. 
  */
  public MPReal(boolean b, int precision)
  {super(b, precision);}
  
  
  public MPReal(int size, boolean b)
  {super(size, b);}
  
  /**
  *  Copy constructor. Calls mp's copy constructor.
  */
  public MPReal(/*const*/ MPReal in) 
  {super((MP)in);}
  
  /**
  *  Creates a MPReal from double with default precision level, mpipl. 
  *  @see MPGlobal#mpipl
  */
  public MPReal(double d)
  { super(d, mpipl);}

 /**
  *  Creates a MPReal from double. 
  *  @param precision The precision level of the object, in digits. 
  *  @see MPGlobal#mpipl
  */
  public MPReal(double d, int precision)
  { super(d, precision);}
  
 /**
  *  Creates a MPReal from String with default precision level, mpipl.
  *  @see MPGlobal#mpipl
  */
  public MPReal(String str)
  {super(str, mpipl);}
  
 /**
  *  Creates a MPReal from String.
  *  @param precision The precision level of the object, in digits. 
  */
  public MPReal(String str, int precision)
  {
    super(str, precision);
  }
  
  /**
   *  Creates a MPReal from MPInt. Inherits the MPInt object's precision.
   *  Calls mp's copy constructor.
   */
  public MPReal(/*const*/ MPInt in)
  {super((MP)in);};
 
  /** 
  *   Non-public constructor
  */
  MPReal(int size)
  {super(size, false);}

  /**
  *  Creates a MPReal from MPComplex's real part with default precision, mpipl.
  *
  *  @see MPGlobal#mpipl
  */
  public MPReal(/*const*/MPComplex mpc)
  {
    this(mpc, mpipl);
  }
  
  /**
  *  Creates a MPReal from MPComplex's real part.
  *
  *  @param precision The precision level of the object, in digits. 
  */ 
  public MPReal(/*const*/MPComplex mpc, int precision)
  {
    super(true, precision);
    int lmpnw = Math.min(mpnw, maxnw - 1);
    mpeq(mpc.r, this, lmpnw);
  }  

  /**
  *  Copies the parameter object to this object. This function only copies
  *  up to the maximum precision level of this object. 
  *
  *  @param ja The object to be copied. 
  *  @return this object.
  */
  public MPReal assign(MP ja)
  {
    if(ja != this)
      MP.mpeq(ja, this, Math.min(mpnw, this.maxnw - 1));
  
    return this;
  }

 /**
  *  Returns a MPReal whose value is this + ja
  */
  public MPReal add(/*const*/ MPReal ja)
  {
    MPReal res = new MPReal();  
    mpadd(this, ja, res, mpnw);
    return res;
  }
  
  /**
  *  Returns a MPReal whose value is this - ja
  */
  public MPReal subtract(/*const*/ MPReal ja)
  {
    MPReal res = new MPReal();  
    mpsub(this, ja, res, mpnw);
    return res;
  }
  
  /**
  *  Returns a MPReal whose value is -this
  */
  public MPReal negate()
  {
    MPReal res = new MPReal();
    
    mpeq(this, res, mpnw);
    res.sign = !this.sign;
    return res;
  }
  
  /**
  *  Returns a MPReal whose value is this * ja
  */
  public MPReal multiply(/*const*/ MPReal ja)
  {
    MPReal res = new MPReal();  
    mpmulx(this, ja, res, mpnw);
    return res;
  }
  
  /**
  *  Returns a MPReal whose value is this / ja
  */
  public MPReal divide(/*const*/ MPReal ja)
  {
    MPReal res = new MPReal();  
    mpdivx(this, ja, res, mpnw);
    return res;
  }


  /**
  *  Returns a MPReal whose value is the absolute value of this number.
  */
  public MPReal abs()
  {
    MPReal res = new MPReal();  
    mpeq(this, res, mpnw);
    res.sign=true;
    return res;
  }

  /**
  *  Returns the MPReal whose value is the greater of this and val. 
  *  If the values are equal, either may be returned. 
  *  Note that no copy is made.
  */  
  public MPReal max(/*const*/ MPReal val)
  {
    int ic = mpcpr(this, val, mpnw);
    if (ic >= 0)
      return this;
    else
      return val;
  }
 
 /**
   *  Returns the MPReal whose value is the lesser of this and val. 
   *  If the values are equal, either may be returned.
   *  Note that no copy is made.
   */  
   public MPReal min(/*const*/ MPReal val)
   {
    int ic = mpcpr(this, val, mpnw);
    if (ic < 0)
      return this;
    else
      return val;
   }
  
  /**
   *  FORTRAN's sign operator.
   *  @return The absolute value of this with sign of val.
   */
   public MPReal sign(/*const*/ MP val)
   {
     MPReal res = new MPReal();
     mpeq(this, res, mpnw);
     res.sign=val.sign;
     return res;
   }
   
  /**
   *  Returns a MPReal whose value is (this ** exponent).
   *  Note that exponent can be MPInt or MPReal.
   */  
  public MPReal pow(/*const*/ MP exponent)
  {
    MPReal res = new MPReal();
    MP mpt1 = new MP(), mpt2 = new MP();
    mplogx(this, mppic, mpl02, mpt1, mpnw);
    mpmulx(mpt1, exponent, mpt2, mpnw);
    mpexpx(mpt2, mppic, mpl02, mpt1, mpnw);
    return res;    
  }
    
  /**
   *  Returns a MPReal whose value is (this ** exponent).
   */
  public MPReal pow(/*const*/ int exponent)
  {
    MPReal res = new MPReal();
    mpnpwx(this, exponent, res, mpnw);
    return res;
  }

 /**
  *  Returns a MPReal whose value is (this ** exponent).
  */
  public MPReal pow(double exponent)
  {
    MPReal res = new MPReal();
    MP mpt1 = new MP(), mpt2 = new MP();
    
    mplogx(this, mppic, mpl02, mpt1, mpnw);
    mpmuld(mpt1, new MPDPE(exponent), mpt2, mpnw);
    mpexpx(mpt2, mppic, mpl02, res, mpnw);
    return res;
  }

  /**
   *  Returns the arccosine of this.
   */
  public MPReal acos()
  {
    MPReal res  = new MPReal();
    MP mpt1 = new MP(), mpt2 = new MP(), mpt3 = new MP();
    
    mpdmc(new MPDPE(1), mpt1);
    mpmulx(this, this, mpt2, mpnw);
    mpsub(mpt1, mpt2, mpt3, mpnw);
    mpsqrx(mpt3, mpt1, mpnw);
    MPComplex.mpangx(this, mpt1, mppic, res, mpnw);
    return res;
  }

  
 /**
  *  Returns the truncation to whole number of this.
  */
  public MPReal aint()
  {
    MPReal res = new MPReal(); 
    MP mpt1 = new MP();
    mpinfr(this, res, mpt1, mpnw);
    return res;
  }
  
  /**
  *  Returns the nearest whole number of this.
  */
  public MPReal anint()
  {
    MPReal res = new MPReal();
    mpnint(this, res, mpnw);
    return res;
  }
  
  /**
  *  Returns the arcsine of this.
  */
  public MPReal asin()
  {
    MPReal res = new MPReal();
    MP mpt1 = new MP(), mpt2 = new MP(), mpt3 = new MP();
    
    mpdmc(new MPDPE(1), mpt1);
    mpmulx(this, this, mpt2, mpnw);
    mpsub(mpt1, mpt2, mpt3, mpnw);
    mpsqrx(mpt3, mpt1, mpnw);
    MPComplex.mpangx(mpt1, this, mppic, res, mpnw);
    return res;
  }

  /**
  *  Returns the arctangent of this.
  */
  public MPReal atan()
  {
    MPReal res = new MPReal();
    MP mpt1 = new MP(6,false);
    
    mpdmc(new MPDPE(1), mpt1);
    MPComplex.mpangx(this, mpt1, mppic, res, mpnw);
    return res;
  }

  /**
   *  Returns the arctangent of val/this.
   */
  public MPReal atan2(/*const*/ MPReal val)
  {
    MPReal res = new MPReal();
    
    MPComplex.mpangx(val, this, mppic, res, mpnw);
    return res;
  }

  /**
   *  Returns the cosine of this.
   */
  public MPReal cos()
  {
    MPReal res = new MPReal();
    MP mpt1 = new MP();
    
    MPComplex.mpcssx(this, mppic, res, mpt1, mpnw);
    return res;
  }
  
  /**
  *  Returns the hyperbolic cosine of this.
  */
  public MPReal cosh()
  {
    MPReal res = new MPReal();
    MP mpt1 = new MP();
    
    mpcshx(this, mppic, mpl02, res, mpt1, mpnw);
    return res;
  }

  /**
  *  Returns the exp(this).
  */
  public MPReal exp()
  {
    MPReal res = new MPReal();
    
    mpexpx(this, mppic, mpl02, res, mpnw);
    return res;
  }
  
  /**
  *  Retuns the natural logarithm of this.
  */
  public MPReal log()
  {
    MPReal res = new MPReal();
    
    mplogx(this, mppic, mpl02, res, mpnw);
    return res;
  }
 
  /**
  *  Return the common logarithm (base 10) of this.
  */
  public MPReal log10()
  {
    MPReal res = new MPReal();
    MP mpt1 = new MP();
    
    mplogx(this, mppic, mpl02, mpt1, mpnw);
    mpdivx(mpt1, mpl10, res, mpnw);
    return res;    
  }

  /**
  *  Hyperbolic cosine and sine of this.
  *  @param cosh out parameter for hyperbolic cosine.
  *  @param sinh out parameter for hyperbolic sine.
  */
  public void csshf(MPReal cosh, MPReal sinh)
  {
    mpcshx(this, mppic, mpl02, cosh, sinh, mpnw);
  }

  /**
  *  Cosine and sine of this.
  *  @param cosine out parameter for cosine.
  *  @param sine out parameter for sine.
  */
  public void cssnf(MPReal cosine, MPReal sine)
  {
    MPComplex.mpcssx(this, mppic, cosine, sine, mpnw);
   
  }
  
  /**
  *  Returns the N-th root of this.
  */
  public MPReal nrtf(int ib)
  {
    MPReal res = new MPReal();
    mpnrtx(this, ib, res, mpnw);
    return res;
  }
  
  /**
  *  Returns a pseudo-random MPReal number between 0 and 1.
  */
  public static MPReal rand()
  {
    MPReal res = new MPReal();
    
    mprand(res, mpnw);
    return res;
  }
  
  /**
  *  Returns the nearest integer (MPInt) of this.
  */
  public MPInt nint()
  {
    MPInt res  = new MPInt();
    
    mpnint(this, res, mpnw);
    return res;
  }
  
  /**
   *  Retuns the sine of this.
   */
  public MPReal sin()
  {
    MPReal res = new MPReal();
    MP mpt1 = new MP();
    
    MPComplex.mpcssx(this, mppic, mpt1, res, mpnw);
    return res;
  }
  
  /**
  *  Returns the hyperbolic sine of this.
  */  
  public MPReal sinh()
  {
    MPReal res = new MPReal();
    MP mpt1 = new MP();
    
    mpcshx(this, mppic, mpl02, mpt1, res, mpnw);
    return res;
  }

  /**
   *  Retuns the square root of this.
   */  
  public MPReal sqrt()
  {
    MPReal res = new MPReal();
    
    mpsqrx(this, res, mpnw);
    return res;
  }

  /**
  *  Returns the tangent of this.
  */
  public MPReal tan()
  {
    MPReal res = new MPReal();
    MP mpt1 = new MP(), mpt2 = new MP();
    
    MPComplex.mpcssx(this, mppic, mpt1, mpt2, mpnw);
    mpdivx(mpt1, mpt2, res, mpnw);
    return res;
  }
  
  /**
   *  Returns the hyperbolic tangent of this.
   */
  public MPReal tanh()
  {
    MPReal res = new MPReal();
    MP mpt1 = new MP(), mpt2 = new MP();
    
    mpcshx(this, mppic, mppic, mpt1, mpt2, mpnw);
    mpdivx(mpt1, mpt2, res, mpnw);
    return res;
  }
   
}

