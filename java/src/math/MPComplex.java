package math;

/**  
*  This class represents multi-precision complex.
*  A MPComplex object consists of two MP objects. One for the real part
*  and another for the imaginary part. 
*  If no parameter of type MPPrecision is passed, the default level of 
*  precision if used.
*
*  @author Herman Harjono
*  @version Oct 7, 1998.
*/ 
public final class MPComplex extends MPGlobal 
  //implements Cloneable//,Comparable 
{
  public void debug(String functionName)
  {
    r.debug(functionName);i.debug(functionName);
  }

  /** 
  *  Real part.
  */
  MP r;
  
  /**
  *  Imaginary part.
  */
  MP i;

  /**
  *  Create a MPComplex whose mantissa is of length size
  */
  private MPComplex(int size, boolean b)
  {
    r = new MP(size, b);
    i = new MP(size, b);
  }

  /**
  *  Returns a MPComplex whose mantissa is of length size
  */
  private static MPComplex createComplex(int size)
  {
    MPComplex res = new MPComplex(size,false);
    return res;
  }

  private int constHelper(int precision)
  {
    int maxnw = precisionToSize(precision);
    r = new MP(maxnw,false); 
    i = new MP(maxnw,false); 
    return maxnw;
  }

  /**
  *  Creates a MPComplex (initialized to zero) with default precision, mpipl. 
  *  @see mp#mpipl
  */
  public MPComplex()
  {
    constHelper(mpipl);
  }
 
  /**
  *  Creates a MPComplex (initialized to zero). 
  *  @param precision The precision level of the object, in digits. I
  */
  public MPComplex(boolean b, int precision)
  {
    constHelper(precision);
  }
  
  /**
  *  Creates a MPComplex from int with default precision, mpipl. 
  *  @see MPGlobal#mpipl
  */
  //public MPComplex(int real)
  //{
   // this(real, new MPPrecision(mpipl));
  //}

 /**
  *  Creates a MPComplex from int.
  *  @param precision The precision level of the object, in digits. 
  */
  //public MPComplex(int real, /*const*/ MPPrecision precision)
  //{
  //  this((double)real, precision);
  //}
  
  /**
  *  Creates a MPComplex from double with default precision, mpipl. 
  *  @see MPGlobal#mpipl
  */
  //public MPComplex(double real)
  //{ 
  //  this(real, new MPPrecision(mpipl));
  //}

 /**
  *  Creates a MPComplex from double. 
  *  @param precision The precision level of the object, in digits. 
  */
  public MPComplex(double real, int precision)
  {
    constHelper(precision);
    MP.mpdmc(new MPDPE(real), r);
  }

  /**
  *  Creates a MPComplex from DComplex with default precision, mpipl. 
  *  @see MPGlobal#mpipl
  */
  public MPComplex(/*const*/ DComplex dc)
  {
    this(dc, mpipl);
  }

   /**
  *  Creates a MPComplex from DComplex. 
  *  @param precision The precision level of the object, in digits. 
  */
  public MPComplex(/*const*/ DComplex dc, int precision)
  {
    constHelper(precision);
    MP.mpdmc(new MPDPE(dc.real()), r); // mpxzc
    MP.mpdmc(new MPDPE((double)dc.aimag()), i);
  }
  
  /**
  *  Creates a MPComplex from String with default precision level, mpipl. 
  *  @see MPGlobal#mpipl
  */
  public MPComplex(String real)
  { 
    this(real, mpipl);
  }
  
 /**
  *  Creates a MPComplex from String. 
  *  @param precision The precision level of the object, in digits. 
  */
  public MPComplex(String real, int precision)
  {
    int maxnw = constHelper(precision);
    MP.mpinpc(real.toCharArray(), real.length(), r, Math.min(mpnw, maxnw-2));
  }

  /**
  *  Creates a MPReal from MP. Inherits the MP object's precision.
  *  Note that MPInt and MPReal are concrete subclasses of MP.
  */
  public MPComplex(/*const*/ MP real)
  {
    int maxnw = real.maxnw;
    r = new MP(maxnw,false);
    i = new MP(maxnw,false);
    
    MP.mpeq(real, r, mpnw); // mpmzc
  }
  
  /**
  *  Copy constructor. 
  */
  public MPComplex(/*const*/ MPComplex in)
  {
    r = new MP(in.r.maxnw,false);
    i = new MP(in.i.maxnw,false);
    
    mpceq(in, this, mpnw);
  }

  /**
  *  Creates a MPComplex from int with default precision level, mpipl. 
  *  @see MPGlobal#mpipl
  */
  public MPComplex(int real, int imag)
  {
    this((double)real, (double)imag, mpipl);
  }
  
  /**
  *  Creates a MPComplex from int. 
  *  @param precision The precision level of the object, in digits. 
  */
  public MPComplex(int real, int imag, int precision)
  {
    this((double)real, (double)imag, precision);
  }
  
  /**
  *  Creates a MPComplex from double with default precision level, mpipl. 
  *  @see MPGlobal#mpipl
  */
  public MPComplex(double real, double imag)
  {
    this(real, imag, mpipl); 
  }
  
  /**
  *  Creates a MPComplex from double. 
  *  @param precision The precision level of the object, in digits. 
  */
  public MPComplex(double real, double imag, int precision)
  {
    constHelper(precision);
    
    MP.mpdmc(new MPDPE(real), r); 
    MP.mpdmc(new MPDPE(imag), i);
  }

  /**
  *  Creates a MPComplex from String with default precision level, mpipl.
  *  @see MPGlobal#mpipl
  */
  public MPComplex(String real, String imag)
  { 
    this (real, imag, mpipl);
  }
  
 /**
  *  Creates a MPComplex from String. 
  *  @param precision The precision level of the object, in digits. 
  */
  public MPComplex(String real, String imag, int precision)
  {
    int maxnw = constHelper(precision); 
    int lmpnw = Math.min(mpnw, maxnw -1);
    MP.mpinpc(real.toCharArray(), real.length(), r, lmpnw);
    MP.mpinpc(imag.toCharArray(), imag.length(), i, lmpnw);
  }

  /**
  *  Creates a MPComplex from mp. Inherits the mp objects precisions.
  */
  public MPComplex(/*const*/ MP real , /*const*/ MP imag)
  {
    r = new MP(real.maxnw,false);
    i = new MP(imag.maxnw,false);
    
    MP.mpeq(real, r, mpnw); 
    MP.mpeq(imag, i, mpnw);
  }

 /**
  *  Copies the parameter object to this object. This function only copies
  *  up to the maximum precision level of this object. 
  *
  *  @param zb The object to be copied. 
  *  @return this object.
  */
  public MPComplex assign(MPComplex zb)
  {
    mpceq(zb, this, Math.min(mpnw, r.maxnw));
    return this;
  }

 /**
  *  Returns a MPComplex whose value is this + ja
  */
  public MPComplex add(/*const*/ MPComplex ja)
  {
    MPComplex res = new MPComplex();  
    mpcadd(this, ja, res, mpnw);
    return res;
  }
  
  /**
  *  Returns a MPComplex whose value is this - ja
  */
  public MPComplex subtract(/*const*/ MPComplex ja)
  {
    MPComplex res = new MPComplex();  
    mpcsub(this, ja, res, mpnw);
    return res;
  }
  
  /**
  *  Returns a MPComplex whose value is -this
  */
  public MPComplex negate()
  {
    MPComplex res = new MPComplex();
    
    mpceq(this, res, mpnw);
    res.r.sign = !res.r.sign;
    res.i.sign = !res.i.sign;
    return res;
  }
  
  /**
  *  Returns a MPComplex whose value is this * ja
  */
  public MPComplex multiply(/*const*/ MPComplex ja)
  {
    MPComplex res = new MPComplex();  
    mpcmlx(this, ja, res, mpnw);
    return res;
  }
  
  /**
  *  Returns a MPComplex whose value is this / ja
  */
  public MPComplex divide(/*const*/ MPComplex ja)
  {
    MPComplex res = new MPComplex();  
    mpcdvx(this, ja, res, mpnw);
    return res;
  }

  /**
  *  Compares the equality of this Object and specified object.
  *
  *  Note that this method only compares up to max. precision level (mpipl).
  *  @param o the Object to be compared. 
  *  @return true if equal and false otherwise.
  */
  public boolean equals(Object o)
  {
    if(o==this) return true;
    try{
      return 0 == MP.mpcpr(r, ((MPComplex)o).r, mpnw) 
        && 0 == MP.mpcpr(i, ((MPComplex)o).i, mpnw);
    }
    catch(ClassCastException cce)
    { return false;}
  }

 /**
  *  Return the string representation (overwrite Object's)
  */
  public String toString()
  {
    return (r.toString() + i.toString());
  }

  /**
  *  Returns the absolute value of the real part.
  */
  public MPReal abs()
  {
    MPReal res = new MPReal();
    MP mpt1 = new MP(), mpt2 = new MP(), mpt3 = new MP();
    
    MP.mpmulx(r, r, mpt1, mpnw);
    MP.mpmulx(i, i, mpt2, mpnw);
    MP.mpadd(mpt1, mpt2, mpt3, mpnw);
    MP.mpsqrx(mpt3, res, mpnw);
    return res;
  }
  
  /**
  *  Returns (copy) the imaginary part of this.
  */
  public MPReal real()
  {
    MPReal res = new MPReal();
    
    MP.mpeq(r, res, mpnw);
    return res;
  }

  /**
  *  Returns (copy) the imaginary part of this.
  */
  public MPReal aimag()
  {
    MPReal res = new MPReal();
    
    MP.mpeq(i, res, mpnw);
    return res;
  }
  
  /**
  *  Returns the conjugate of this.
  */
  public MPComplex conjg()
  {
    MPComplex res = new MPComplex();
      
      mpceq(this, res, mpnw);
    res.i.sign = ! i.sign;
    return res;
  }

  /**
  *  Returns the value of the real part as double.
  *  This may involves truncation or rounding.
  */
  public double doubleValue()
  {
    return r.doubleValue();
  }
 
  
  /**
  *  Returns the value of the real part as int.
  *  This may involves truncation or rounding.
  */
  public int intValue()
  {
    return r.intValue();
  }
  
  
  /**
  *  Returns the value of the real part as float.
  *  This may involves truncation or rounding.
  */
  public float floatValue()
  {
    return r.floatValue();
  }
  
  /**
  *  Returns the value of the object as DComplex
  *  This may involves truncation or rounding.
  */
  public DComplex complexValue()
  {
    return new DComplex(r.toDPE().value(),
      i.toDPE().value());
  }
    
  
  /**
  *  Returns the square root of this.
  */
  public MPComplex sqrt()
  {
    MPComplex res = new MPComplex();
    
    mpcsqt(this, res, mpnw);
    return res;
  }
  
  /**
  *  Returns the value of (this ** exponent)
  */
  public MPComplex pow(/*const*/ int exponent)
  {
    MPComplex res = new MPComplex();  
    mpcpwx(this, exponent, res, mpnw);
    return res;
  }


 //************************* ARITHMETIC ROUTINES ****************************
  
  /**
  *  This computes the sum of the MPC numbers A and B and returns the MPC
  *  result in C.  
  *
  *  L is the offset between real and i parts in A, B, and c.
  *  Debug output starts with MPIDB = 9.
  *  <br>
  *  Dimension a(2*l), b(2*l), c(2*l).L must be at least MPNW + 2.
  */
  static void mpcadd (/*const*/ MPComplex a, /*const*/ MPComplex b, 
    MPComplex c, int lmpnw)
  {
    MP.mpadd (a.r, b.r, c.r, lmpnw);
    MP.mpadd (a.i, b.i, c.i, lmpnw);
  }

  
 /**
  *  This routine divides the MP complex numbers A and B to yield the MPC
  *  quotient c.  
  *
  *  For extra high levels of precision, use MPCDVX. 
  *  The last word of the result is not reliable.  
  *  Debug output starts with MPIDB = 7
  *  <br>
  *  Dimension a(2*l), b(2*l), c(2*l). L must be at least MPNW + 2. 
  *  <br>
  *  This routine employs the formula described in MPCMUL to save multiprecision
  *  multiplications.
  *  @see #mpcdvx
  *  @see #mpcmul
  */
  static void mpcdiv (/*const*/ MPComplex a, /*const*/ MPComplex b, 
    MPComplex c, int lmpnw)
  {
    int mpnw2 = lmpnw+2;
    MP f = new MP(6,false), sk0 = new MP(mpnw2,false), 
      sk1 = new MP(mpnw2,false), sk2 = new MP(mpnw2,false), 
      sk3 = new MP(mpnw2,false), sk4 = new MP(mpnw2,false); 
    MP ar=a.r, ai=a.i, br=b.r, bi=b.i;
    
    if (b.r.nw  == 0 && b.i.nw == 0) 
    {
      throw new ArithmeticException
        ("mpcdiv: Divisor is zero.");
    }
    
    f.nw=1; f.exponent=0; f.sign=true; f.mantissa[0]=1; f.mantissa[1]=0;
    
    MP.mpmul (ar, br, sk0, lmpnw);
    MP.mpmul (ai, bi, sk1, lmpnw);
    MP.mpadd (sk0, sk1, sk2, lmpnw);
    MP.mpsub (sk0, sk1, sk3, lmpnw);
    MP.mpadd (ar, ai, sk0, lmpnw);
    MP.mpsub (br, bi, sk1, lmpnw);
    MP.mpmul (sk0, sk1, sk4, lmpnw);
    MP.mpsub (sk4, sk3, sk1, lmpnw);
    MP.mpmul (br, br, sk0, lmpnw);
    MP.mpmul (bi, bi, sk3, lmpnw);
    MP.mpadd (sk0, sk3, sk4, lmpnw);
    MP.mpdiv (f, sk4, sk0, lmpnw);
    MP.mpmul (sk2, sk0, c.r, lmpnw);
    MP.mpmul (sk1, sk0, c.i, lmpnw);    
  }

  /**
  *  This sets the MPC number B equal to the MPC number A. 
  *
  *  Debug output starts with MPIDB = 10.
  *  <br>
  *  Dimension a(2*l), b(2*l).
  */
  static void mpceq (/*const*/ MPComplex a, MPComplex b, int lmpnw)
  {
    int i;
    
    int n1 = Math.min (a.r.nw, lmpnw);
    int n2 = Math.min (a.i.nw, lmpnw);
    b.r.sign = a.r.sign;
    b.r.nw = n1;
    b.r.exponent=a.r.exponent;
    b.i.sign=a.i.sign;
    b.i.nw=n2;
    b.i.exponent=a.i.exponent;
    
    for (i = 0; i<n1; i++)
      b.r.mantissa[i] = a.r.mantissa[i];
    for (i=0; i<n2; i++)
      b.i.mantissa[i] = a.i.mantissa[i];
  }
    
 /**
  *  This routine multiplies the MP complex numbers A and B to yield the MPC
  *  product c. 
  *
  *  For extra high levels of precision, use MPCMLX.  
  *  The last word of the result is not reliable.  
  *  Debug output starts with MPIDB = 7.
  *  <br>
  *  Dimension a(2*l), b(2*l), c(2*l). L must be at least MPNW + 2.
  *  <br>
  *  This routine employs the formula:
  *  <pre>
  *  (a_1 + a_2 i) (b_1 + b_2 i)  =  [a_1 b_1 - a_2 b_2]  +
  *   [(a_1 + a_2) (b_1 + b_2) - (a_1 b_1 + a_2 b_2)] i
  *  </pre>
  *  Note that this formula can be implemented with only three multiplications
  *  whereas the conventional formula requires four.
  *  @see #mpcmlx
  */
  static void mpcmul(/*const*/ MPComplex a, /*const*/ MPComplex b, 
    MPComplex c, int lmpnw)
  { 
    int mpnw2 = lmpnw+2;
    MP sk0 = new MP(mpnw2,false), sk1 = new MP(mpnw2,false), 
      sk2 = new MP(mpnw2,false), sk3 = new MP(mpnw2,false);
    
    MP.mpmul (a.r, b.r, sk0, lmpnw);
    MP.mpmul (a.i, b.i, sk1, lmpnw);
    MP.mpsub (sk0, sk1, c.r, lmpnw);
    MP.mpadd (sk0, sk1, sk2, lmpnw);
    MP.mpadd (a.r, a.i, sk0, lmpnw);
    MP.mpadd (b.r, b.i, sk1, lmpnw);
    MP.mpmul (sk0, sk1, sk3, lmpnw);
    MP.mpsub (sk3, sk2, c.i, lmpnw);    
  }
  
  /**
  *  This computes the N-th power of the MPC number A and returns the MPC
  *  result C in B.  
  *
  *  When N is zero, 1 is returned.  When N is negative, the
  *  reciprocal of A ^ |N| is returned. For extra high levels of precision, 
  *  use MPCPWX.  Debug output starts with MPIDB = 7.
  *  <br>
  *  Dimension a(2*l), b(2*l). L should be at least MPNW + 2. 
  *  <br>
  *  This routine employs the binary method for exponentiation.
  *  @see #mpcpwx
  */
  static void mpcpwr (/*const*/ MPComplex a, int n, 
    MPComplex b, int lmpnw)
  {
    int j;
    double t1;
    
    int na1 = Math.min (a.r.nw, lmpnw);
    int na2 = Math.min (a.i.nw, lmpnw);
    
    if (na1 == 0 && na2 == 0) 
    {
      if (n >= 0) 
      {
        zero(b);
        return;
      } 
      else 
      {
        throw new ArithmeticException
          ("mpcpwr: Argument is zero and N is negative or zero.");
      }
    }
    
    int nws = lmpnw;
    lmpnw++;
    int nn = Math.abs (n);
    
    if (nn == 0) 
    {
      b.r.nw=1; b.r.sign=true; b.r.exponent=0; (b.r).mantissa[0]=1; (b.r).mantissa[1]=0;
      MP.zero(b.i);
      return;
    } 
    
    MPComplex f = createComplex(6), sk0 = createComplex(lmpnw+3), 
      sk1 = createComplex(lmpnw+3), sk2 = createComplex(lmpnw+3);   
    mpceq (a, sk0, lmpnw); 
    f.r.nw=1; f.r.sign=true; f.r.exponent=0; (f.r).mantissa[0]=1; (f.r).mantissa[1]=0;
    boolean skip = false;
    
    if (nn == 1) 
    {
      mpceq (sk0, sk2, lmpnw); 
      skip = true;
    } 
    else if (nn == 2) 
    {
      mpcmul (sk0, sk0, sk2, lmpnw);
      skip = true;
    }
    
    //  Determine the least integer MN such that 2 ^ MN > NN.
    if (!skip)
    {
      int mn, kn, kk;
      t1 = nn;
      mn = (int)(CL2 * Math.log (t1) + 1.0 + MPRXX); // *cast*
      mpceq (f, sk2, lmpnw);
      kn = nn;
      
      //  Compute B ^ N using the binary rule for exponentiation.
      
      for (j = 1; j<=mn; j++)
      {
        kk = kn / 2;
        if (kn != 2 * kk) 
        {
          mpcmul (sk2, sk0, sk1, lmpnw);
          mpceq (sk1, sk2, lmpnw);
        }
        kn = kk;
        if (j < mn) 
        {
          mpcmul (sk0, sk0, sk1, lmpnw);
          mpceq (sk1, sk0, lmpnw);
        }
      }
    }
    //  Compute reciprocal if N is negative.
    
    if (n < 0) 
    {
      mpceq (f, sk1, lmpnw);
      mpcdiv (sk1, sk2, sk0, lmpnw);
      mpceq (sk0, sk2, lmpnw);
    }
    mpceq (sk2, b, lmpnw);
    
    mpcroun (b, nws);
  }
  
  /**
  * Rounding a complex number.
  */
  private static void mpcroun(MPComplex in, int lmpnw)
  { MP.mproun(in.r, lmpnw); MP.mproun(in.i, lmpnw);}

  /**
  *  This routine computes the complex square root of the MPC number A and 
  *  put the result in C. 
  *
  *  For extra high levels of precision, use MPCSQX.  
  *  The last word of the result is not reliable.  
  *  Debug output starts with MPIDB = 6.
  *  <br>
  *  Dimension a(2*l), b(2*l). L must be at least MPNW + 2.
  *  <br>
  *  This routine uses the following formula, where A1 and A2 are the real and
  *  imaginary parts of A, and where R = Sqrt [A1 ^ 2 + A2 ^2]:
  *  <pre>
  *  B = Sqrt [(R + A1) / 2] + I Sqrt [(R - A1) / 2]
  *  </pre>
  *  If the imaginary part of A is < 0, then the imaginary part of B is also
  *  set to be < 0.
  *  @see #mpcsqx
  */
  static void mpcsqt (/*const*/ MPComplex a, MPComplex b, int lmpnw)
  {
    int mpnw2 = lmpnw+2;
    MP sk0 = new MP(mpnw2,false), sk1 = new MP(mpnw2,false), 
      sk2 = new MP(mpnw2,false);
    
    if (a.r.nw == 0 && a.i.nw == 0) 
    {
      zero(b);
      return;
    }
    
    MP.mpmul (a.r, a.r, sk0, lmpnw);
    MP.mpmul (a.i, a.i, sk1, lmpnw);
    MP.mpadd (sk0, sk1, sk2, lmpnw);
    MP.mpsqrt (sk2, sk0, lmpnw);
    MP.mpeq (a.r, sk1, lmpnw);
    sk1.sign = true;
    MP.mpadd (sk0, sk1, sk2, lmpnw);
    MP.mpmuld (sk2, new MPDPE(0.50), sk1, lmpnw);
    MP.mpsqrt (sk1, sk0, lmpnw);
    MP.mpmuld (sk0, new MPDPE(2.0), sk1, lmpnw);
    if (a.r.nw >= 0) 
    {
      MP.mpeq (sk0, b.r, lmpnw);
      MP.mpdiv (a.i, sk1, b.i, lmpnw);
    } 
    else 
    {
      MP.mpdiv (a.i, sk1, b.r, lmpnw);
      b.r.sign = true;
      MP.mpeq (sk0, b.i, lmpnw);
      b.i.nw = a.i.nw;
    }  
  }

  /**
  *  This subracts the MPC numbers A and B and returns the MPC difference in
  *  C. 
  *
  *  Debug output starts with MPIDB = 9.
  *  <br>
  *  Dimension a(2*l), b(2*l), c(2*l). L must be at least MPNW + 2.
  */
  static void mpcsub (/*const*/ MPComplex a, /*const*/ MPComplex b, 
    MPComplex c, int lmpnw)
  {
    MP.mpsub (a.r, b.r, c.r, lmpnw);
    MP.mpsub (a.i, b.i, c.i, lmpnw);
  }
    
  /**
  * Initialize the complex number to zero.
  */
  private static void zero(MPComplex in)
  {MP.zero(in.r); MP.zero(in.i);}

  //*********************** ADVANCE ARITHMETIC ROUTINES ******************************
  
  /**
  *  This routine divides the MP complex numbers A and B to yield the MPC
  *  quotient c.  
  *
  *  Before calling MPCDVX, the arrays UU1 and UU2 must be initialized
  *  by calling MPINIX. For modest levels of precision, use MPCDIV.
  *  The last word of the result is not reliable.  
  *  Debug output starts with MPIDB = 7
  *  <br>
  *  Dimension a(2*l), b(2*l), c(2*l). L must be at least MPNW + 2. 
  *  <br>
  *  This routine employs the same scheme as MPCDIV.
  *  @see #mpcdiv
  */
  static void mpcdvx (/*const*/ MPComplex a, /*const*/ MPComplex b, 
    MPComplex c, int lmpnw)
  {   
    int mpnw2 = lmpnw+2;
    MP f = new MP(6,false), sk0 = new MP(mpnw2,false), 
      sk1 = new MP(mpnw2,false), sk2 = new MP(mpnw2,false), 
      sk3 = new MP(mpnw2,false), sk4 = new MP(mpnw2,false); 
    MP ar=a.r, ai=a.i, br=b.r, bi=b.i;
    if (b.r.nw  == 0 && b.i.nw == 0) 
    {
      throw new ArithmeticException
        ("mpcdvx: Divisor is zero.");
    }
    
    f.nw=1; f.exponent=0; f.sign=true; f.mantissa[0]=1; f.mantissa[1]=0;
    
    MP.mpmulx (ar, br, sk0, lmpnw);
    MP.mpmulx (ai, bi, sk1, lmpnw);
    MP.mpadd (sk0, sk1, sk2, lmpnw);
    MP.mpsub (sk0, sk1, sk3, lmpnw);
    MP.mpadd (ar, ai, sk0, lmpnw);
    MP.mpsub (br, bi, sk1, lmpnw);
    MP.mpmulx (sk0, sk1, sk4, lmpnw);
    MP.mpsub (sk4, sk3, sk1, lmpnw);
    MP.mpsqx (br,  sk0, lmpnw);
    MP.mpsqx (bi, sk3, lmpnw);
    MP.mpadd (sk0, sk3, sk4, lmpnw);
    MP.mpdivx (f, sk4, sk0, lmpnw);
    MP.mpmul (sk2, sk0, c.r, lmpnw);
    MP.mpmul (sk1, sk0, c.i, lmpnw);  
  }
  
  /**
  *  This routine multiplies the MP complex numbers A and B to yield the MPC
  *  product c. 
  *
  *  Before calling MPCMLX, the arrays UU1 and UU2 must be initialized by 
  *  calling MPINIX. For modest levels of precision, use MPCMUL. 
  *  The last word of the result is not reliable.  
  *  Debug output starts with MPIDB = 7.
  *  <br>
  *  Dimension a(2*l), b(2*l), c(2*l). L must be at least MPNW + 2.
  *  <br>
  *  This routine employs the same scheme as MPCMUL.
  *  @see #mpcmul
  */
  static void mpcmlx (/*const*/ MPComplex a, /*const*/ MPComplex b, 
    MPComplex c, int lmpnw)
  {
    int mpnw2 = lmpnw+2;
    MP sk0 = new MP(mpnw2,false), sk1 = new MP(mpnw2,false), 
      sk2 = new MP(mpnw2,false), sk3 = new MP(mpnw2,false);
    
    MP.mpmulx (a.r, b.r, sk0, lmpnw);
    MP.mpmulx (a.i, b.i, sk1, lmpnw);
    MP.mpsub (sk0, sk1, c.r, lmpnw);
    MP.mpadd (sk0, sk1, sk2, lmpnw);
    MP.mpadd (a.r, a.i, sk0, lmpnw);
    MP.mpadd (b.r, b.i, sk1, lmpnw);
    MP.mpmulx (sk0, sk1, sk3, lmpnw);
    MP.mpsub (sk3, sk2, c.i, lmpnw);
  }
  
  /**
  *  This computes the N-th power of the MPC number A and returns the MPC
  *  result c in B.  
  *
  *  When N is zero, 1 is returned.  When N is negative, the
  *  reciprocal of A ^ |N| is returned. Before calling MPCPWX, 
  *  the arrays UU1 and UU2 must be initialized by calling MPINIX.  
  *  For modest levels of precision, use MPCPWR.  The last word of the result 
  *  is not reliable.  Debug output starts with MPIDB = 6.
  *  <br>
  *  Dimension a(2*l), b(2*l). L should be at least MPNW + 2. 
  *  <br>
  *  This routine employs the binary method for exponentiation.
  *  @see #mpcpwr
  */
  static void mpcpwx (/*const*/ MPComplex a, int n, 
    MPComplex b, int lmpnw)
  {  
    int j;
    double t1;
    
    int na1 = Math.min (a.r.nw, lmpnw);
    int na2 = Math.min (a.i.nw, lmpnw);
    int ncr = (int)(Math.pow(2, mpmcr)); // *cast*
    
    //  Check if precision level of A is too low to justify advanced routine.
    if (na1 <= ncr && na2 <= ncr) 
    {
      mpcpwr (a, n, b, lmpnw);
      return;
    }
    
    if (na1 == 0 && na2 == 0) 
    {
      if (n >= 0) 
      {
        zero(b);
        return;
      } 
      else 
      {
        throw new ArithmeticException
          ("mpcpwx: Argument is zero and N is negative or zero.");
      }
    }
    int nn = Math.abs (n);
    
    if (nn == 0) 
    {
      b.r.nw=1; b.r.sign=true; b.r.exponent=0; (b.r).mantissa[0]=1; (b.r).mantissa[1]=0;
      MP.zero(b.i);
      return;
    } 
    
    MPComplex f = createComplex(6), sk0 = createComplex(lmpnw+3), 
      sk1 = createComplex(lmpnw+3), sk2 = createComplex(lmpnw+3);   
    mpceq (a, sk0, lmpnw); 
    f.r.nw=1; f.r.sign=true; f.r.exponent=0; (f.r).mantissa[0]=1; (f.r).mantissa[1]=0;
    boolean skip = false;
    
    if (nn == 1) 
    {
      mpceq (sk0, sk2, lmpnw); 
      skip = true;
    } 
    else if (nn == 2) 
    {
      mpcmul (sk0, sk0, sk2, lmpnw);
      skip = true;
    }
    
    //  Determine the least integer MN such that 2 ^ MN > NN.
    if (!skip)
    {
      int mn, kn, kk;
      t1 = nn;
      mn = (int)(CL2 * Math.log (t1) + 1.0 + MPRXX); // *cast*
      mpceq (f, sk2, lmpnw);
      kn = nn;
      
      //  Compute B ^ N using the binary rule for exponentiation.
      for (j = 1; j<=mn; j++)
      {
        kk = kn / 2;
        if (kn != 2 * kk) 
        {
          mpcmlx (sk2, sk0, sk1, lmpnw);
          mpceq (sk1, sk2, lmpnw);
        }
        kn = kk;
        if (j < mn) 
        {
          mpcmlx (sk0, sk0, sk1, lmpnw);
          mpceq (sk1, sk0, lmpnw);
        }
      }
    }
    //  Compute reciprocal if N is negative.
    if (n < 0) 
    {
      mpceq (f, sk1, lmpnw);
      mpcdvx (sk1, sk2, sk0, lmpnw);
      mpceq (sk0, sk2, lmpnw);
    }
    mpceq (sk2, b, lmpnw);  
  }
    
  /**
  *  This routine computes the complex square root of the MPC number A.  
  *
  *  For modest levels of precision, use MPCSQT.  The last
  *  word of the result is not reliable.  Debug output starts with MPIDB = 5.
  *  <br>
  *  Dimension a(2*l), b(2*l). L must be at least MPNW + 2.
  *  <br>
  *  This routine uses the same algorithm as MPCSQT.
  *  @see #mpcsqt
  */
  static void mpcsqx (/*const*/ MPComplex a, MPComplex b, int lmpnw)
  {
    int mpnw2 = lmpnw+2;
    MP sk0 = new MP(mpnw2,false), sk1 = new MP(mpnw2,false), sk2 = new MP(mpnw2,false); 
    
    if (a.r.nw == 0 && a.i.nw == 0) 
    {
      zero(b);
      return;
    }
    
    MP.mpsqx (a.r, sk0, lmpnw);
    MP.mpsqx (a.i, sk1, lmpnw);
    MP.mpadd (sk0, sk1, sk2, lmpnw);
    MP.mpsqrx (sk2, sk0, lmpnw);
    MP.mpeq (a.r, sk1, lmpnw);
    sk1.sign = true;
    MP.mpadd (sk0, sk1, sk2, lmpnw);
    MP.mpmuld (sk2, new MPDPE(0.50), sk1, lmpnw);
    MP.mpsqrx (sk1, sk0, lmpnw);
    MP.mpmuld (sk0, new MPDPE(2.0), sk1, lmpnw);
    
    if (a.r.sign) 
    {
      MP.mpeq (sk0, b.r, lmpnw);
      MP.mpdivx (a.i, sk1, b.i, lmpnw);
    } 
    else 
    {
      MP.mpdivx (a.i, sk1, b.r, lmpnw);
      b.r.sign = true;
      MP.mpeq (sk0, b.i, lmpnw);
      b.i.sign = a.i.sign;
    }
  }
  

  //**************** ADVANCE ALGEBRAIC AND TRANSCENDENTAL ROUTINES ******************
  
  /**
  *  This computes the MP angle A subtended by the MP pair (X, Y) considered as
  *  a point in the x-y plane.  
  *
  *  This is more useful than an arctan or arcsin
  *  routine, since it places the result correctly in the full circle, i.e. 
  *  -Pi < A <= Pi.  PI is the MP value of Pi computed by a previous to
  *  MPPI or MPPIX.  Before calling MPANGX, the arrays UU1 and UU2 must be
  *  initialized by calling MPINIX.  For modest levels of precision, use MPANG.
  *  The last word of the result is not reliable.  Debug output starts 
  *  with MPIDB = 6.
  *  <br>
  *  Dimension a(lmpnw+2), pi(lmpnw), x(lmpnw), y(lmpnw).
  *  <br>
  *  This routine employs a complex arithmetic version of the MPLOGX alogirthm.
  *  @see MP#mplogx
  *  @see MP#mppi
  *  @see MP#mppix
  *  @see MP#mpang
  */
  static void mpangx (/*const*/ MP x, /*const*/ MP y, /*const*/ MP pi, MP a, int lmpnw)
  {   
    int ix = x.sign ? 1 : -1;
    int nx = Math.min (x.nw, lmpnw);
    int iy = y.sign ? 1: -1;
    int ny = Math.min (y.nw, lmpnw);
    
    int ncr = (int)(Math.pow(2, mpmcr));
    
    //  Check if precision level is too low to justify the advanced routine.
    
    if (lmpnw <= ncr) 
    {
      MP.mpang (x, y, pi, a, lmpnw);
      return;
    }
    
    //  Check if both X and Y are zero.
    
    if (nx == 0 && ny == 0) 
    {
      throw new ArithmeticException
        ("mpangx: Both arguments are zero.");
    }
    
    //  Check if Pi has been precomputed.
    MPDPE t1 = new MPDPE();
    MP.mpmdc (pi, t1);
    if (t1.n != 0 || Math.abs (t1.a - CPI) > MPRX2) 
    {
      throw new ArithmeticException
        ("mpangx: PI must be precomputed.");
      
    }
    
    //  Check if one of X or Y is zero.
    if (nx == 0) 
    {
      if (iy > 0) 
        MP.mpmuld (pi, new MPDPE(0.5), a, lmpnw); 
      else 
        MP.mpmuld (pi, new MPDPE(-0.5), a, lmpnw); 
      
      return;
      
    } 
    else if (ny == 0) 
    {
      if (ix > 0) 
        MP.zero(a);
      else 
        MP.mpeq (pi, a, lmpnw);
      return;
    }
    
    //  Define scratch space.
    
    MPComplex sk0 = createComplex(lmpnw+2), 
      sk1 = createComplex(lmpnw+2), 
      sk2 = createComplex(lmpnw+2), 
      sk3 = createComplex(lmpnw+2);
    
    //  Multiply the input by a large power of two.
    MP.mpmdc (x, t1);
    int n2 = MPNBT * (lmpnw / 2 + 2) - t1.n;
    MPDPE dpe1 = new MPDPE(1.0, n2);
    MP.mpmuld (x, dpe1, sk0.r, lmpnw);
    MP.mpmuld (y, dpe1, sk0.i, lmpnw);
    
    //  Perform AGM iterations.
    sk1.r.nw=1; sk1.r.sign=true; sk1.r.exponent=0; (sk1.r).mantissa[0]=1; (sk1.r).mantissa[1]=0;
    MP.zero(sk1.i);
    sk3.r.nw=1; sk3.r.sign=true; sk3.r.exponent=0; (sk3.r).mantissa[0]=4; (sk3.r).mantissa[1]=0;
    MP.zero(sk3.i);
    
    mpcdvx (sk3, sk0, sk2, lmpnw);
    mpcagx (sk1, sk2, lmpnw);
    
    //  Compute A = Imag (Pi / (2 * Z)), where Z is the limit of the complex AGM.
    
    dpe1.a = 2.0; dpe1.n = 0;
    MP.mpmuld (sk1.r, dpe1, sk0.r, lmpnw);
    MP.mpmuld (sk1.i, dpe1, sk0.i, lmpnw);
    MP.mpeq(pi, sk2.r, lmpnw); MP.zero(sk2.i);
    mpcdvx (sk2, sk0, sk1, lmpnw);
    MP.mpeq (sk1.i, a, lmpnw);  
  }

  /**
  *  This performs the arithmetic-geometric mean (AGM) iterations.  
  *
  *  This routine is called by MPANGX.  It is not intended to be called directly by the user.
  *  <br>
  *  Dimension a(2*lmpnw+4), b(2*lmpnw+4).
  *  @see #mpangx
  */
  static void mpcagx (MPComplex a, MPComplex b, int lmpnw)
  {    
    int l1 = 0;
    MPComplex sk0 = createComplex(lmpnw+2), sk1 = createComplex(lmpnw+2);  
    int s1;
    //100  
    MPDPE dpe1 = new MPDPE(0.50,0);
    do
    { 
      l1++;
      if (l1 == 50) 
      {
        throw new ArithmeticException
          ("mpcagx: Iteration limit exceeded.");
      }
      
      s1 = sk0.r.exponent;
      mpcadd (a, b, sk0, lmpnw);
      MP.mpmuld (sk0.r, dpe1, sk1.r, lmpnw);
      MP.mpmuld (sk0.i, dpe1, sk1.i, lmpnw);
      mpcmlx (a, b, sk0, lmpnw);
      mpcsqx (sk0, b, lmpnw);
      mpceq (sk1, a, lmpnw);
      MP.mpsub (a.r, b.r, sk0.r, lmpnw);
      
    }
    //  Check for convergence.
    while (sk0.r.nw != 0 && (sk0.r.exponent < s1 || sk0.r.exponent >= -2)) ;
  }

  /**
  *  This computes the cosine and sine of the MP number A and returns the two MP
  *  results in X and Y, respectively.  
  *
  *  PI is the MP value of Pi computed by a previous  to MPPI or MPPIX. 
  *  Before calling MPCSSX, the arrays UU1 and UU2 must be initialized 
  *  by calling MPINIX.  For modest levels of precision, use MPCSSN.  
  *  The last word of the result is not reliable.
  *  Debug output starts with MPIDB = 5.
  *  <br>
  *  Dimension a(lmpnw), pi(lmpnw), x(lmpnw+2), y(lmpnw+2).
  *  <br>
  *  This routine employs a complex arithmetic version of the scheme used in
  *  MPEXPX.  See the comment about the parameter NIT in MPDIVX.
  *  @see MP#mppi
  *  @see MP#mppix
  *  @see MP#mpcssn
  *  @see MP#mpexpx
  *  @see MP#mpdivx
  */
  static void mpcssx (/*const*/ MP a, /*const*/ MP pi, MP x, MP y, int lmpnw)
  {
   
    int k;
    double t2;
    final int nit = 1;
    MP f1 = new MP(6,false); 
    
    int na = Math.min (a.nw, lmpnw);
    int ncr = (int)(Math.pow(2, mpmcr)); // *cast*
    
    //  Check if precision level is too low to justify advanced routine.
    
    if (lmpnw <= ncr) 
    {
      
      MP.mpcssn (a, pi, x, y, lmpnw);
      return;
    }
    
    //  Check if input is zero
    if (na == 0) 
    {
      x.sign=true; x.nw=1; x.exponent=0; x.mantissa[0]=1;
      MP.zero(y);
      return;
    }
    
    //  Check if Pi has been precomputed.
    
    MPDPE t1 = new MPDPE();
    MP.mpmdc (pi, t1);
    if (t1.n != 0 || Math.abs (t1.a - CPI) > MPRX2) 
    {
      throw new ArithmeticException
        ("mpcssx: PI must be precomputed.");    
    }
    f1.nw=1; f1.sign=true; f1.exponent=0;
    f1.mantissa[0]=1; f1.mantissa[1]=0;
    int nws = lmpnw;
    
    MPComplex sk0 = createComplex(lmpnw+2), 
      sk1 = createComplex(lmpnw+2), 
      sk2 = createComplex(lmpnw+2),
      sk3 = createComplex(lmpnw+2);
    
    //  Reduce argument to between - Pi and Pi.
    
    MP.mpmuld (pi, new MPDPE(2.0), sk0.r, lmpnw);
    MP.mpdivx (a, sk0.r, sk1.r, lmpnw);
    MP.mpnint (sk1.r, sk2.r, lmpnw);
    MP.mpmul (sk2.r, sk0.r, sk1.r, lmpnw);
    MP.mpsub (a, sk1.r, sk0.r,  lmpnw);
    
    //  Determine the least integer MQ such that 2 ^ MQ >= MPNW.
    
    t2 = nws;
    int mq = (int)(CL2 * Math.log (t2) + 1.0 - MPRXX); // *cast*
    
    MP.mpeq (f1, sk2.r, lmpnw);
    
    //  Compute initial approximation to [Cos (A), Sin (A)].
    
    lmpnw = ncr;
    MP.mpcssn (sk0.r, pi, sk3.r, sk3.i, lmpnw);
    int iq = 0;
    
    //  Perform the Newton-Raphson iteration with a dynamically changing precision
    //  level MPNW.
    for (k = mpmcr + 1; k<=mq; k++)
    {
      lmpnw = Math.min (2 * lmpnw, nws);
      boolean cont = true;
      while (cont)
      {
        mpangx (sk3.r, sk3.i, pi, sk1.r, lmpnw);
        MP.mpsub (sk0.r, sk1.r, sk2.i, lmpnw);
        mpcmlx (sk3, sk2, sk1, lmpnw);
        mpceq (sk1, sk3, lmpnw);
        if (k == mq - nit && iq == 0) 
          iq = 1;
        else
          cont = false;
      }
    }
    
    //  The final (cos, sin) result must be normalized to have magnitude 1.
    
    MP.mpsqx (sk3.r, sk0.r, lmpnw);
    MP.mpsqx (sk3.i, sk0.i, lmpnw);
    MP.mpadd (sk0.r, sk0.i, sk1.r, lmpnw);
    MP.mpsqrx (sk1.r, sk2.r, lmpnw);
    MP.mpdivx (sk3.r, sk2.r, sk0.r, lmpnw);
    MP.mpdivx (sk3.i, sk2.r, sk0.i, lmpnw);
    MP.mpeq(sk0.r, x, lmpnw);
    MP.mpeq(sk0.i, y, lmpnw);    
  }
  
}

