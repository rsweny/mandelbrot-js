package math;


/**  
* This class represents multi-precision integer and is intended to be directly
* used by the user.
*
* While some operations (like comparisons) are inherited from MP,
* this class implements other operations that are specific to MPInt.
* The implementations of the operations are similar to those of MPReal except 
* that it truncates the fractional part before returning the result.
*
* If a MPInt object needs a lower precision level than the current maximum 
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
public final class MPInt extends MP
{
 /**
  *  Creates a MPInt (initialized to zero) with default precision level. 
  *  @see MPGlobal#mpipl
  */
  public MPInt()
  {super(true, mpipl);}
  
  /**
  *  Creates a MPInt (initialized to zero). 
  *  @param precision The precision level of the object, in digits. 
  */
  public MPInt(int precision)
  {super(true, precision);}
  
  /**
  *  Copy constructor. Calls mp's copy constructor.
  */
  public MPInt(/*const*/ MPInt in) 
  {super((MP)in);}
  
  /**
  *  Creates a MPInt from int with default precision level, mpipl.
  *  @param precision The precision level of the object, in digits. 
  */
  //public MPInt(int ia)
  //{ super(ia, new MPPrecision(mpipl));}
  
  /**
  *  Creates a MPInt from int. 
  *  @param precision The precision level of the object, in digits. 
  *  @see MPGlobal#mpipl
  */
  //public MPInt(int ia, /*const*/ MPPrecision precision)
  //{ super(ia, precision);}
  
  /**
  *  Creates a MPInt from double with default precision level, mpipl. 
  *  This may involves truncation or rounding.
  *  @see MPGlobal#mpipl
  */
  //public MPInt(double d)
  //{ this(d, new MPPrecision(mpipl));}
  
  /**
  *  Creates a MPInt from double. 
  *  This may involves truncation or rounding.
  *  @param precision The precision level of the object, in digits. 
  */
  public MPInt(double d, int precision)
  {
    super(true, precision);
    MP mpt1 = new MP(6,false), 
      mpt2 = new MP(8,false);
    mpdmc(new MPDPE(d), mpt1);
    mpinfr(mpt1, this, mpt2, Math.min(maxnw - 2, mpnw));
  }
  
  /**
  *  Creates a MPInt from String with default precision level, mpipl.
  *  This may involves truncation or rounding.
  *  @see MPGlobal#mpipl
  */
  public MPInt(String str)
  {this(str, mpipl);}
  
  /**
  *  Creates a MPInt from String.
  *  This may involves truncation or rounding.
  *  @param precision The precision level of the object, in digits. 
  */
  public MPInt(String str, int precision)
  {
    super(true, precision);
    int lmpnw = Math.min(mpnw, maxnw-2);
    int mpnw2 = lmpnw+2;
    MP mpt1 = new MP(mpnw2, false), mpt2 = new MP(mpnw2, false);
    mpdexc (str.toCharArray(), str.length(), mpt1, lmpnw);
    mpinfr (mpt1, this, mpt2, lmpnw);
  }
  
  /**
  *  Creates a MPInt from MPReal with default precision, mpipl.
  *
  *  This may involve truncation.
  *  @see MPGlobal#mpipl
  */
  public MPInt(/*const*/MPReal mpr)
  {
    this(mpr, mpipl);
  }
  
  /**
  *  Creates a MPInt from MPReal.
  *
  *  This may involve truncation.
  *  @param precision The precision level of the object, in digits. 
  */
  public MPInt(/*const*/MPReal mpr, int precision)
  {
    super(true, precision);
    int lmpnw = Math.min(mpnw, maxnw - 2);
    MP mpt1 = new MP(lmpnw+1,false), 
      mpt2 = new MP(lmpnw+2,false);
    mpeq(mpr, mpt1, lmpnw);
    mpinfr(mpt1, this, mpt2, lmpnw);
    
  }
  
  /**
  *  Creates a MPInt from MPComplex's real part with default precision, mpipl.
  *
  *  This may involve truncation.
  *  @see MPGlobal#mpipl
  */
  public MPInt(/*const*/MPComplex mpc)
  {
    this(mpc, mpipl);
  }
  
  /**
  *  Creates a MPInt from MPComplex's real part.
  *
  *  This may involve truncation.
  *  @param precision The precision level of the object, in digits. 
  */ 
  public MPInt(/*const*/MPComplex mpc, int precision)
  {
    super(true, precision);
    int lmpnw = Math.min(mpnw, maxnw - 2);
    MP mpt1 = new MP(lmpnw+1,false), 
      mpt2 = new MP(lmpnw+2,false);
    mpeq(mpc.r, mpt1, lmpnw);
    mpinfr(mpt1, this, mpt2, lmpnw);
  }
  
  /**
  *  Copies the parameter object to this object. This function only copies
  *  up to the maximum precision level of this object. 
  *
  *  @param ja The object to be copied. 
  *  @return this object.
  */
  public MPInt assign(MP ja)
  {
    if(ja != this)
    {
      if(ja.maxnw == this.maxnw && ja instanceof MPInt)
        MP.mpeq(ja, this, Math.min(mpnw, this.maxnw - 1));
      else
      {
        int mpnw1 = Math.min(mpnw, ja.maxnw-1);
        int mpnw2 = Math.min(mpnw, maxnw-2);
        MPReal mpt1 = new MPReal(),  
          mpt2 = new MPReal();
        MP.mpeq(ja, mpt1, mpnw1);
        MP.mpinfr(mpt1, this, mpt2, mpnw2);
      }
    }
    return this;
  }

  /**
  *  Returns a MPInt whose value is this + ja
  */
  public MPInt add(/*const*/ MPInt ja)
  {
    MPInt res = new MPInt(), 
      mpt1 = new MPInt(), mpt2 = new MPInt();
    
    mpadd(ja, this, mpt1, mpnw);
    mpinfr(mpt1, res, mpt2, mpnw);
    return res;
  }
  
  /**
  *  Returns a MPInt whose value is this - ja
  */
  public MPInt subtract(/*const*/ MPInt ja)
  {
    MPInt res = new MPInt(), 
      mpt1 = new MPInt(), mpt2 = new MPInt();
    
    mpsub(this, ja, mpt1, mpnw);
    mpinfr(mpt1, res, mpt2, mpnw);
    return res;
  }
  
  /**
  *  Returns a MPInt whose value is -this
  */
  public MPInt negate()
  {
    MPInt res = new MPInt();
    
    mpeq(this, res, mpnw);
    res.sign = ! this.sign;
    return res;
  }
  
  /**
  *  Returns a MPInt whose value is this * ja
  */
  public MPInt multiply(/*const*/ MPInt ja)
  {
    MPInt res = new MPInt(), 
      mpt1 = new MPInt(), mpt2 = new MPInt();
    mpmulx(ja, this, mpt1, mpnw);
    mpinfr(mpt1, res, mpt2, mpnw);
    return res;
  }
  
  /**
  *  Returns a MPInt whose value is this / ja
  */
  public MPInt divide(/*const*/ MPInt ja)
  {
    MPInt res = new MPInt(), 
      mpt1 = new MPInt(), mpt2 = new MPInt();
    mpdivx(this, ja, mpt1, mpnw);
    mpinfr(mpt1, res, mpt2, mpnw);
    return res;
  }
  
  /**
  *  Returns a MPInt whose value is this % ja
  */
  public MPInt mod(/*const*/ MPInt ja)
  {
    MPInt res = new MPInt(), 
      mpt1 = new MPInt(), mpt2 = new MPInt(), 
      mpt3 = new MPInt();
    mpdivx(this, ja, mpt1, mpnw);
    mpinfr(mpt1, mpt2, mpt3, mpnw);
    mpmulx(ja, mpt2, mpt1, mpnw);
    mpsub(this, mpt1, res, mpnw);
    return res;
  }
  

  /**
   *  Returns a MPInt whose value is the absolute value of this number.
   */
  public MPInt abs()
  {
    MPInt res = new MPInt();  
    mpeq(this, res, mpnw);
    res.sign=true;
    return res;
  }

  /**
  *  Returns the MPInt whose value is the greater of this and val. 
  *  If the values are equal, either may be returned. 
  *  Note that no copy is made.
  */  
  public MPInt max(/*const*/ MPInt val)
  {
    int ic = mpcpr(this, val, mpnw);
    if (ic >= 0)
      return this;
    else
      return val;
  }
 
 /**
   *  Returns the MPInt whose value is the lesser of this and val. 
   *  If the values are equal, either may be returned.
   *  Note that no copy is made.
   */  
   public MPInt min(/*const*/ MPInt val)
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
   public MPInt sign(/*const*/ MP val)
   {
     MPInt res = new MPInt();
     mpeq(this, res, mpnw);
     res.sign=val.sign;
     return res;
   }
  
  /**
   *  Returns a MPInt whose value is (this ** exponent).
   */
   public MPInt pow(/*const*/ MPInt exponent)
   {
     MPInt res = new MPInt(), 
       mpt1 = new MPInt(), mpt2 = new MPInt();
     mplogx(this, MPReal.mppic, MPReal.mpl02, mpt1, mpnw);
     mpmulx(mpt1, exponent, mpt2, mpnw);
     mpexpx(mpt2, MPReal.mppic, MPReal.mpl02, mpt1, mpnw);
     mpnint(mpt1, res, mpnw);
     return res;
   }
  
  /**
   *  Returns a MPInt whose value is (this ** exponent).
   */
   public MPInt pow(/*const*/ int exponent)
   {
     MPInt res= new MPInt(), mpt1= new MPInt();  
     mpnpwx(this, exponent, mpt1, mpnw);
     mpnint(mpt1, res, mpnw);
     return res;
   }
  

}

