package math;

/**  
* This class represents MP objects.
* It also contains functions that  manipulate the mp objects.
* The core algorithms of the package is defined in this class.
* The format of an mp number:
*    sign (int): the sign. 
*    nw (int): the number of mantissa words 
*    exponent (int): the power of the radix. 
*                    The radix is 2^24 (1677216) for IEEE system. 
*    mantissa (float[]): floating point whole numbers between 0 and 
*                        one less than the radix.
*
* This class is not intended to be used directly by the user.
*
* The comments of the methods in this classes have some keywords:
*   Dimension: the bound of the mantissa array of the MP number.
*
* @author Herman Harjono
* @version Oct 5, 1998.
*/ 
public class MP extends MPGlobal implements Cloneable//Comparable 
{
  
  public void debug(String functionName)
  {
    for(int count=0;count< nw+2;count++)
      System.out.println(count + "  " + functionName + "  " 
        + (int)mantissa[count] + ".");
  }

  /********************** DATA MEMBERS ***********************/
  /**
  *  The maximum number of mantissa words. 
  *  
  *  This is the size of the mantissa array.
  *  Invariant: maxnw - 2 >= nw.
  */
  int maxnw;
  
  /**
  *  The sign. False if negative and true otherwise.
  */
  boolean sign;
  
  /**
  *  The number of significant words.
  *  
  *  Invariant: nw <= maxnw - 2.
  */
  int nw;
  
  /**
  *  Exponent of the radix
  */
  int exponent;
  
  /**
  *  Mantissa array. Size of the array is maxnw.
  *  @see #maxnw
  */
  public float mantissa[];   

  /**
  *  Variables used by mprand
  */
  static double t30, r30=0, sd = 314159265;
  
  /**
  *  Mutex variable for mprand
  */
  static Integer randMutex = new Integer(0);
  

  
  
 /*********************** CONSTRUCTORS **********************/
 /** 
  *  Constructor. Intiliaze number to 0 with default size.
  */
  MP()
  {this(mp2, false);}
 
 /**
  *  Copy constructor. Calls mpeq.
  */
  MP(/*const*/ MP in)
  {
    maxnw = in.maxnw;
    if(maxnw>0)
    {
      mantissa = new float[maxnw];
      mpeq(in, this, mpnw);
    }
    else 
    {
      exponent = in.exponent;
      sign = in.sign;
      nw = in.nw;
      mantissa = null;
    } 
  }

  /** 
  *  Constructor. Intiliaze number to 0.
  *  @param mpnw initial value of mpnw. If initialized to 
  *  a value < 1, mantissa array won't be allocated.
  *  @see #mantissa
  */
  MP(int maxNW, boolean b)
  {
    maxnw=maxNW;
    sign=true;
    nw=0;
    exponent=0;
    if(maxnw>=1)
    {
      mantissa = new float[maxnw];		
      mantissa[0]=0;
    }
    else
      mantissa = null;
  } 
  
  
 /**
  *  Creates a MP (initialized to zero). 
  *  @param precision The precision level of the object, in digits.
  */
  MP(boolean b, int precision)
  {
    maxnw = (precisionToSize(precision));
    mantissa = new float[maxnw];
    mantissa[0]=0;
    sign=true;
    nw=0;
    exponent=0;   
  }
  

  /**
  *  Creates a MP from float / double. 
  *  @param precision The precision level of the object, in digits. 
  */
  MP(double ia, int precision)
  {
    maxnw = precisionToSize(precision);
    mantissa = new float[maxnw];
    mpdmc(new MPDPE(ia), this);
  }
  
 /**
  *  Creates a MP from String. 
  *  @param precision The precision level of the object, in digits. If equals "", 
  *  it is assigned the value of mpipl.
  *  @see MP#mpipl
  */
  MP(String str)
  {
    this(str, mpipl);
  }
  
  /**
  *  Creates a MP from String. 
  *  @param precision The precision level of the object, in digits. If equals "", 
  *  it is assigned the value of mpipl.
  *  @see MP#mpipl
  */
  MP(String str, int precision)
  {
    char temp[] = str.toCharArray();
    maxnw = precisionToSize(precision);
    mantissa = new float[maxnw];
    mpdexc (temp, temp.length, this, Math.min(mpnw, maxnw-2));
  }

  /**
  *   Return the clone of this.
  */
  public Object clone() 
  {
    return new MP(this);
  }

 /**
  *  Return the string representation (overwrite Object's)
  *  10 ^        -4 x  3.14159265358979323846264338327950288419716939937510,
  */
  public String toString()
  {
    StringBuffer res = new StringBuffer();
    char az[] = new char[mpipl + 100];
    int nd = mpouts(this, mpoud, az, mpnw);
    int i;
    for(i=0; i<nd; i++)
    {
      res.append(az[i]);
   }
   
   String old = res.toString();
   int hat = old.indexOf("^");
   int xx = old.indexOf("x");
   String exponent = old.substring(hat+1, xx).trim();
   String number = old.substring(xx+1).trim();
   return  number +"E" + exponent;
  }
  
  public int getExponent()
  {
    StringBuffer res = new StringBuffer();
    char az[] = new char[mpipl + 100];
    int nd = mpouts(this, mpoud, az, mpnw);
    int i;
    for(i=0; i<nd; i++)
    {
      res.append(az[i]);
   }
   
   String old = res.toString();
   int hat = old.indexOf("^");
   int xx = old.indexOf("x");
   String exponent = old.substring(hat+1, xx).trim();
   return Integer.parseInt(exponent);
  }
  
 /**
  *  Compares the equality of this Object and specified object.
  *
  *  Note that this method only compares up to max. precision level (mpipl).
  *  @param o the Object to be compared. 
  *  @return true if equal and false otherwise.
  *  @see MPGlobal#mpipl
  */
  public boolean equals(Object o)
  {
    if(o==this) return true;
    try{
      return 0 == mpcpr(this, (MP)o, mpnw);
    }
    catch(ClassCastException cce)
    { return false;}
  }


  /**
  * Compares this Object with the specified Object for order.
  *
  * Note that this method only compares up to the maximum precision level (mpipl).
  * @param o the Object to be compared. 
  * @return Returns -1, 0, or 1 as this Object is less than, equal to, or 
  * greater than the given Object.
  * @exception  ClassCastException
  *  when the specified Object's type prevents it from being compared to this Object.
  * @see MPGlobal#mpipl
  */
  public int compareTo(Object o)
  {
    return mpcpr(this, (MP)o, mpnw);
  }
  
  public boolean isNegative()
  {
  	MP zero = new MP(true, 10);
  	return (compareTo(zero) < 0 );
  }

 /**
  *  Returns the value of this Object as a double. 
  *
  *  This may involves rounding or truncation.  
  */    
  public double doubleValue() 
  {
    return this.toDPE().value();
  }

  /**
  *  Returns the value of this Object as a float. 
  *
  *  This may involves rounding or truncation.  
  */    
  public float floatValue() 
  {
    return (float)this.toDPE().value();
  }

  /**
  *  Returns the value of this Object as an int. 
  *
  *  This may involves rounding or truncation.  
  */  
  public int intValue() 
  {
    return (int)this.toDPE().value();
  }

 /**
  *  Returns the value of this Object as a long. 
  *
  *  This may involves rounding or truncation.  
  */
  public long longValue() 
  {
    return (long)this.toDPE().value();
  }

 /**
  *  Returns the value of this Object as a short. 
  *
  *  This may involves rounding or truncation.  
  */
  public short shortValue() 
  {
    return (short)this.toDPE().value();
  }

  /**
  *  Return the value of mpipl
  */
  public static int getMpipl()
  {return mpipl;}  

  /***********************  MP ARITHMETIC FUNCTIONS ******************************/
  
  /** 
  *  This routine adds MP numbers a and b to yield the MP sum c.  
  *
  *  It attempts to include all significance of a and b in the result, up to the maximum
  *  mantissa length mpnw. This is a new simplified version. 
  *  <br>
  *  Dimension a(mpnw), b(mpnw), c(mpnw+2).
  */
  static void mpadd (/*const*/ MP a, /*const*/ MP b, 
    MP c, int lmpnw)
  {
    int i;
    
    int na = Math.min (a.nw, lmpnw);
    int nb = Math.min (b.nw, lmpnw);
    
    
    /* Sweny
    //  checks for zero inputs.
    if (na == 0) 
    {
      //  a is zero -- the result is b.
      c.sign=b.sign;
      c.nw=nb;
      c.exponent=b.exponent;
      
      for(i=0;i<nb;i++)
        c.mantissa[i] = b.mantissa[i];
      
      return;
    } 
    else if (nb == 0) 
    {
      
      //  b is zero -- the result is a.
      c.sign=a.sign;
      c.nw=na;
      c.exponent=a.exponent;
      
      for (i=0;i<na;i++)
        c.mantissa[i]=a.mantissa[i];
      
      return;
    }
    */
    
    double d[] = new double[lmpnw+5];  
    double db;
    if (a.sign == b.sign) 
      db = 1.0;  
    else
      db = -1.0;
    
    int ixa = a.exponent;  
    int ixb = b.exponent; 
    int ish = ixa - ixb;
    
    int nd,ixd;
    
    if (ish >= 0)
    {
      //  A has greater exponent than B, so B must be shifted to the right.
      int m1 = Math.min (na, ish);
      int m2 = Math.min (na, nb + ish);
      int m3 = na;
      int m4 = Math.min (Math.max (na, ish), lmpnw + 1); 
      int m5 = Math.min (Math.max (na, nb + ish), lmpnw + 1); 
      d[0]=0;
      d[1]=0;
      
      for (i = 0; i<m1; i++)
        d[i+2] = a.mantissa[i];
      
      for (i = m1; i<m2; i++)
        d[i+2] = a.mantissa[i] + db * b.mantissa[i-ish];
      
      for (i = m2; i<m3; i++)
        d[i+2] = a.mantissa[i];
      
      for (i = m3; i<m4; i++)
        d[i+2] = 0.0;
      
      for (i = m4; i <m5; i++)
        d[i+2] = db * b.mantissa[i-ish];
      
      nd = m5;
      ixd = ixa;
      d[nd+2] = 0.0;
      d[nd+3] = 0.0;
    } 
    else
    {
      //  B has greater exponent than A, so A must be shifted to the right.
      int nsh = - ish;
      int m1 = Math.min (nb, nsh);
      int m2 = Math.min (nb, na + nsh);
      int m3 = nb;
      int m4 = Math.min (Math.max (nb, nsh), lmpnw+1); 
      int m5 = Math.min (Math.max (nb, na + nsh), lmpnw+1); 
      d[0]=0;
      d[1]=0;
      
      for (i = 0; i<m1; i++)
        d[i+2] = db * b.mantissa[i]; 
      
      for (i = m1; i<m2; i++)
        d[i+2] = a.mantissa[i-nsh] + db * b.mantissa[i]; 
      
      for (i = m2; i<m3; i++)
        d[i+2] = db * b.mantissa[i]; 
      
      for (i = m3; i<m4; i++)
        d[i+2] = 0.0;
      
      for (i = m4; i< m5; i++)
        d[i+2] = a.mantissa[i-nsh];      
      
      nd = m5;
      ixd = ixb;
      d[nd+2] = 0.0;
      d[nd+3] = 0.0;
    }
    
    //  Call mpnorm to fix up result and store in c.
    d[0] = fSign(nd, a.sign ? 1 : -1);
    d[1]=ixd;
    
    mpnorm(d,c,lmpnw);
  }

  /**
  *  This computes the cube root of the MP number A and returns the MP result
  *  in B.  
  *
  *  For extra high levels of precision, use MPCBRX.  Debug output
  *  starts with MPIDB = 7.
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw+2).
  *  <br>
  *  This subroutine employs the following Newton-Raphson iteration, which
  *  converges to A ^ (-2/3):
  *  <pre>  
  *  X_{k+1} = X_k + (1 - X_k^3 * A^2) * X_k / 3
  *  </pre>
  *  where the muliplication () * X_k is performed with only half of the
  *  normal level of precision.  These iterations are performed with a
  *  maximum precision level MPNW that is dynamically changed, doubling with
  *  each iteration.  The final iteration is performed as follows (this is
  *  due to A. Karp):
  *  <pre>
  *  Cbrt(A) = (A * X_n) + [A - (A * X_n)^3] * X_n / 3 (approx.)
  *  </pre>
  *  where the multiplications A * X_n and [] * X_n are performed with only
  *  half of the final level of precision.  See the comment about the parameter
  *  NIT in MPDIVX.
  *  @see #mpcbrx
  *  @see #mpidvx
  */
  static void mpcbrt (/*const*/ MP a, MP b, int lmpnw)
  {
    int lmpnw3 = lmpnw+3;
    MP f =  new MP(6, false) ;
    MP sk0 = new MP(lmpnw3,false);
    MP sk1 = new MP(lmpnw3,false);  
    MP sk2 = new MP(lmpnw3,false);
    
    int na = Math.min (a.nw, lmpnw); 
    
    if (na == 0)
    {
      zero(b);
      return;
    }
    
    if (a.sign == false) 
    {
      throw new ArithmeticException(
        "mpcbrt: argument is negative --> " + 
        a);
    }
    
    // saving initial lmpnw;
    int nws=lmpnw;
    
    //  Determine the least integer MQ such that 2 ^ MQ >= MPNW.
    MPDPE t1 = new MPDPE();
    t1.a = lmpnw; 
    int mq = (int)(CL2 * Math.log (t1.a) + 1.0 - MPRXX); // *cast*
    
    //  Compute A^2 outside of the iteration loop.
    
    lmpnw++;
    mpmul (a, a, sk0, lmpnw);
    
    //  Compute the initial approximation of A ^ (-2/3).
    mpmdc (a, t1);
    MPDPE t2 = new MPDPE();
    t2.n = - 2 * t1.n / 3;
    t2.a = Math.pow((t1.a * Math.pow(2.0,(t1.n + 3.0 * t2.n / 2.0))), (-2.0 / 3.0));
    mpdmc (t2, b);
    f.sign=true;
    f.nw=1;
    f.exponent=0;
    f.mantissa[0]=1;
    f.mantissa[1]=0;
    
    lmpnw=3; 
    int iq = 0;
    
    //  Perform the Newton-Raphson iteration described above with a dynamically
    //  changing precision level MPNW (one greater than Math.powers of two).
    
    for (int k = 2; k<=mq - 1;k++)
    {
      int nw1 = lmpnw; 
      lmpnw=Math.min (2 * lmpnw - 2, nws) + 1;
      int nw2 = lmpnw; 
      
      boolean cont=true;
      while (cont)
      {
        mpmul (b, b, sk1, lmpnw);
        mpmul (b, sk1, sk2, lmpnw);
        mpmul (sk0, sk2, sk1, lmpnw);
        mpsub (f, sk1, sk2, lmpnw);
        
        lmpnw=nw1;
        
        mpmul (b, sk2, sk1, lmpnw);
        mpdivd (sk1, new MPDPE(3.0), sk2, lmpnw);
        
        lmpnw=nw2;
        
        mpadd (b, sk2, sk1, lmpnw);
        mpeq (sk1, b, lmpnw);
        
        if (k == mq - NIT && iq == 0) 
        {
          iq = 1;
        }
        else
          cont=false;
      }
    }
    
    //  Perform last iteration using Karp's trick.
    
    mpmul (a, b, sk0, lmpnw);
    int nw1 = lmpnw; 
    lmpnw = Math.min (2 * lmpnw - 2, nws) + 1;
    
    int nw2 = lmpnw; 
    
    mpmul (sk0, sk0, sk1, lmpnw);
    mpmul (sk0, sk1, sk2, lmpnw);
    mpsub (a, sk2, sk1, lmpnw);
    
    lmpnw=nw1;
    
    mpmul (sk1, b, sk2, lmpnw);
    mpdivd (sk2, new MPDPE(3.0), sk1, lmpnw);
    
    lmpnw=nw2;
    
    mpadd (sk0, sk1, sk2, lmpnw);
    mpeq (sk2, b, lmpnw);
    
    mproun (b, nws);    
  }

  /**
  *  This routine compares the MP numbers A and B and returns in IC the value
  *  -1, 0, or 1 depending on whether A < B, A = B, or A > B.  
  *  
  *  It is faster than merely subtracting A and B and looking at the fSign of the result.
  *  Debug output begins with MPIDB = 9.
  *  <br>
  *  Dimension a(lmpnw+2), b(lmpnw+2).
  */
  static int mpcpr (/*const*/ MP a, /*const*/ MP b, int lmpnw)
  {
      
    int i, ic=0;      
    int ia = (a.sign ? 1 : -1 );
    if (a.nw == 0.) ia = 0;
    int ib = (b.sign ? 1 : -1);
    if (b.nw == 0.) ib = 0;
    
    int na = Math.min (a.nw, lmpnw);
    int nb = Math.min (b.nw, lmpnw);
    
    //  Compare signs.
    if (ia != ib) 
      ic = (int)(fSign (1, ia - ib)); // *cast*
    
    //  The signs are the same.  Compare exponents.
    else if (a.exponent != b.exponent) 
      ic = (int)(ia * fSign (1, a.exponent - b.exponent)); // *cast*
    
    //  The signs are the exponents are the same. 
    //  Compare the actual (common) mantissas.
    else  
    { 
      boolean sameMantissas=true;
      for (i = 0; i<Math.min (na, nb);i++)
      {
        if (a.mantissa[i] != b.mantissa[i]) 
        {
          ic = (int)(ia * fSign (1., a.mantissa[i] - b.mantissa[i])); // *cast*
          sameMantissas=false;
          break;
        }
      }
      
      if(sameMantissas) 
        if (na != nb) // compares the length of mantissa array
        ic = (int)(ia * fSign (1, na - nb));// *cast*
      else
        //  The signs, exponents, mantissas and lengths are the same.  Thus A = B.
        ic = 0;
    }
    return ic;
  }

  /**
  *  This divides the MP number A by the MP number B to yield the MP quotient c.
  *
  *  For extra high levels of precision, use MPDIVX.  Debug output starts with
  *  MPIDB = 8.
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw), c(lmpnw+2).
  *  @see #mpdivx
  */
  static void mpdiv (/*const*/ MP a, /*const*/ MP b, MP c, int lmpnw)
  {   
    int i,j;
    double rb, ss, t0, t1, t2;
     
    int na = Math.min (a.nw, lmpnw);
    int nb = Math.min (b.nw, lmpnw);
    
    //  Check if dividend is zero.
    if (na == 0) 
    {
      zero(c);
      return;
    }
    
    if (nb == 1 && b.mantissa[0] == 1.) 
    {
      
      //  Divisor is 1 or -1 -- result is A or -A.
      c.nw=na;
      c.sign=!(a.sign^b.sign); // !XOR
      c.exponent=a.exponent-b.exponent;
      
      for (i = 0; i<na;i++)
        c.mantissa[i] = a.mantissa[i];
      
      return;
    }
    
    //  Check if divisor is zero.
    if (nb == 0) 
    {
      throw new ArithmeticException(
        "mpdiv: Divisor is zero.");
    }
    
    //  Initialize trial divisor and trial dividend.
    double d[] = new double[lmpnw+4];  
    t0 = MPBDX * b.mantissa[0];
    if (nb >= 2) t0 += b.mantissa[1];
    if (nb >= 3) t0 += MPRDX * b.mantissa[2];
    if (nb >= 4) t0 += MPRX2 * b.mantissa[3];
    rb = 1.0 / t0;
    int md = Math.min (na + nb, lmpnw);
    d[0]=0.0;
    
    for (i=1; i<=na; i++)
      d[i] = a.mantissa[i-1];
    
    for (i = na + 1; i<md + 4; i++)
      d[i] = 0.0;
    
    
    //  Perform ordinary long division algorithm.  First compute only the first
    //  NA words of the quotient.
    
    for(j = 2; j<=na + 1; j++)
    {
      t1 = MPBX2 * d[j-2] + MPBDX * d[j-1] + d[j] + MPRDX * d[j+1];
      t0 = (int) (rb * t1);
      int j3 = j - 3;
      int i2 = Math.min (nb, lmpnw + 2 - j3) + 2;
      int ij = i2 + j3;
      
      int i3;
      for (i = 2; i<i2; i++)
      {
        i3 = i + j3;
        d[i3] -= t0 * b.mantissa[i-2]; // *cast*
      }
      
      //  Release carries periodically to avoid overflowing the exact integer
      //  capacity of double floating point words in D.
      if (((j - 1)%MPNPR) == 0) 
      {
        //dir$ ivdep
        for (i = j; i<ij;i++)
        {
          t1 = d[i];
          t2 = (int) (MPRDX * t1);
          d[i] = t1 - MPBDX * t2; // *cast*
          d[i-1] += t2; // *cast*
        }
      }
      d[j-1] += MPBDX * d[j-2]; // *cast*
      d[j-2] = t0; // *cast* 
    }
    
    //  Compute additional words of the quotient, as long as the remainder
    //  is nonzero.
    
    
    boolean stopped=false;
    for (j = na + 2; j<=lmpnw + 3; j++)
    {
      t1 = MPBX2 * d[j-2] + MPBDX * d[j-1] + d[j];
      if (j <= lmpnw + 2) t1 = t1 + MPRDX * d[j+1];
      t0 = (int) (rb * t1);
      int j3 = j - 3;
      int i2 = Math.min (nb, lmpnw + 2 - j3) + 2;
      int ij = i2 + j3;
      ss = 0.0;
      
      int i3;
      for (i = 2; i<i2;i++)
      {
        i3 = i + j3;
        d[i3] -= t0 * b.mantissa[i-2]; // *cast*
        ss = ss + Math.abs (d[i3]);
      }
      
      if (((j - 1)%MPNPR) == 0) 
      {
        //dir$ ivdep
        for (i = j; i<ij; i++)
        {
          t1 = d[i];
          t2 = (int) (MPRDX * t1);
          d[i] = t1 - MPBDX * t2; // *cast*
          d[i-1] += t2; // *cast*
        }
      }
      d[j-1] += MPBDX * d[j-2]; // *cast*
      d[j-2] = t0; // *cast*
      if (ss == 0.0)
      {
        stopped=true;
        break;
      }
      if (ij <= lmpnw + 1) d[ij+2] = 0.0;
    }
    
    //  Set fSign and exponent, and fix up result.
    
    if(!stopped)
      j = lmpnw + 3;
    
    d[j-1] = 0.0;
    int is;
    if (d[0] == 0.0) 
      is=1;
    else 
      is=2;
    
    int nc = Math.min (j - 1, lmpnw);
    d[nc+2] = 0.0;
    d[nc+3] = 0.0;
    
    for(i=j; i>=2; i--)
      d[i] =  d[i-is];
    
    d[0] = fSign(nc, (!(a.sign^b.sign)) ? 1 : -1);
    d[1] = a.exponent - b.exponent + is - 2;  
    
    mpnorm (d, c, lmpnw);
  }


  /** 
  *  This routine divides the MP number A by the DPE number (B, N) to yield
  *  the MP quotient c.  
  *
  *  Debug output starts with MPIDB = 9.
  *  <br>
  *  Dimension a(lmpnw), c(lmpnw+2).
  */
  static void mpdivd(/*const*/ MP a,/*const*/ MPDPE b,  
    MP c, int lmpnw)
  {
    int j,k;
    double bb, br, dd, t1;
    MP f =  new MP(6, false) ;
      
    int na = Math.min (a.nw, lmpnw);
    int ib = (int)(fSign (1.0, b.a)); // *cast*
    
    //  Check if dividend is zero.
    
    if (na == 0) 
    {
      zero(c);
      return;
    }
    
    //  Check if divisor is zero.
    if (b.a == 0.0) 
    {
      throw new ArithmeticException(
        "mpdivd: Divisor is zero");
    }

    int n1 = b.n / MPNBT;
    int n2 = b.n - MPNBT * n1;
    bb = Math.abs (b.a) * Math.pow(2.0,n2);
    
    //  Reduce BB to within 1 and MPBDX.
    if (bb >= MPBDX) 
    {
      
      for (k = 1; k<=100;k++)
      {
        bb = MPRDX * bb;
        if (bb < MPBDX) 
        {
          n1 += k;
          break;
        }
      }
    }
    else if (bb < 1.0) 
    {
      for (k = 1; k<=100; k++)
      {
        bb = MPBDX * bb;
        if (bb >= 1.0)
        {
          n1 -= k;
          break;
        }
      }
    }
    
    //  If B cannot be represented exactly in a single mantissa word, use MPDIV.
    
    if (bb != (int)(bb)) 
    {
      MPDPE dpb = new MPDPE(fSign(bb, b.a), n1 * MPNBT);
      mpdmc (dpb , f);
      mpdiv (a, f, c, lmpnw);
      return;
    }
    
    double d[] = new double[lmpnw+4];  
    
    br = 1.0 / bb;
    dd = a.mantissa[0];
    
    //  Perform short division (not vectorizable at present).  Continue as long as
    //  the remainder remains nonzero.
    
    boolean skipJ=false;
    
    for(j = 2; j<= lmpnw + 3;j++)
    {
      t1 = (int)(br * dd);
      d[j] = t1; 
      dd = MPBDX * (dd - t1 * bb);
      if (j <= na) 
        dd = dd + a.mantissa[j-1];
      
      else 
      if (dd == 0.0)
      {
        skipJ=true;
        break;
      }
    }
    //  Set fSign and exponent of result.
    
    if(!skipJ)
      j = lmpnw + 3;
    
    int nc = Math.min (j - 1, lmpnw);
    d[0] = fSign(nc, (a.sign ? 1 : -1)*ib);
    d[1] = a.exponent - n1;
    
    if (j <= lmpnw + 2) d[j+1] = 0.0;
    if (j <= lmpnw + 1) d[j+2] = 0.0;
    
    mpnorm (d, c, lmpnw);
  }

      
  /**
  *  This routine converts the DPE number (A, N) to MP form in B.  
  *  
  *  All bits of A are recovered in B.  However, note for example that if A = 0.1 and N
  *  is 0, then B will NOT be the multiprecision equivalent of 1/10.  Debug
  *  output starts with MPIDB = 9.
  *  <br>
  *  Dimension b(6).
  */
  static void mpdmc(/*const*/ MPDPE a, MP b)
  {
    int i,k;
    double aa;
      
    //  Check for zero.
    if (a.a == 0.0) 
    {
      zero(b);
      return;
    }
    
    int n1 = a.n / MPNBT;
    int n2 = a.n - MPNBT * n1;
    aa = Math.abs (a.a) * Math.pow(2.0, n2);
    
    //  Reduce AA to within 1 and MPBDX.
    if (aa >= MPBDX) 
    {
      for(k = 1; k<=100; k++)
      {
        aa = MPRDX * aa;
        if (aa < MPBDX) 
        {
          n1 = n1 + k;
          break;
        }
      }
    }
    else if (aa < 1.0) 
    {
      for (k = 1; k<=100; k++)
      {
        aa = MPBDX * aa;
        if (aa >= 1.0) 
        {
          n1 = n1 - k;
          break;
        }
      } 
    }
    
    //  Store successive sections of AA into B.
    b.exponent=n1;
    b.mantissa[0] = (float)((int)(aa));  // *cast* float
    aa = MPBDX * (aa - b.mantissa[0]);
    b.mantissa[1] = (float)((int)(aa)); // *cast* float
    aa = MPBDX * (aa - b.mantissa[1]);
    b.mantissa[2] = (float)((int)(aa)); // *cast* float
    aa = MPBDX * (aa - b.mantissa[2]);
    b.mantissa[3] = (float)((int)(aa)); // *cast* float
    b.mantissa[4] = 0;
    b.mantissa[5] = 0;
    
    for (i = 3; i>=0; i--)
      if (b.mantissa[i] != 0.) break;
    
    aa = i + 1;
    b.sign=(a.a>=0); 
    b.nw=(int)(aa); // *cast*
    
  }

  /**
  *  This routine sets the MP number B equal to the MP number A. 
  *
  *  Debug output starts with MPIDB = 10.
  *  <br>
  *  Dimension a(lmpnw + 1), b(lmpnw + 1).
  *  <br>
  *  The fact that only MPNW + 1 cells, and not MPNW + 2 cells, are copied is
  *  important in some routines that increase the precision level by one.
  */
  static void mpeq(/*const*/ MP a, MP b, int lmpnw)
  {
    int i;
    
    int na = Math.min (a.nw, lmpnw);
    if (na == 0) 
    {
      zero(b);
      return;
    }
    b.nw=na;
    b.sign=a.sign;
    b.exponent=a.exponent;
    for (i = 0; i<=na;i++)
      b.mantissa[i] = a.mantissa[i];
  }

 /**
  *  Sets B to the integer part of the MP number A and sets c equal to the
  *  fractional part of A.  
  *
  *  Note that if A = -3.3, { B = -3 and c = -0.3.
  *  Debug output starts with MPIDB = 9.
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw+2), c(lmpnw+2).
  */
  static void mpinfr (/*const*/ MP a, MP b, MP c, int lmpnw)
  {
    int i;
    
    //  Check if  A  is zero.
    int na = Math.min (a.nw, lmpnw);
    int ma = a.exponent;
    if (na == 0)  
    {
      zero(b); zero(c);
    }
    
    if (ma >= lmpnw - 1) 
    {
      throw new ArithmeticException(
        "mpinfr: Argument is too large -->" + a);
    }
    
    //  Place integer part in  B.
    int nb = Math.min (Math.max (ma + 1, 0), na);
    if (nb == 0) 
      zero(b); 
    else 
    {
      b.nw=nb; 
      b.sign=a.sign; 
      b.exponent=ma;
      b.mantissa[nb] = 0;
      b.mantissa[nb+1] = 0;
      
      for(i = 0; i<nb; i++)
        b.mantissa[i] = a.mantissa[i];
    }
    
    //  Place fractional part in c.
    int nc = na - nb;
    if (nc <= 0) 
      zero(c);
    else 
    {
      c.nw=nc;
      c.sign=a.sign;
      c.exponent=ma-nb;
      c.mantissa[nc] = 0;
      c.mantissa[nc+1] = 0;
      
      for(i = 0; i<nc; i++)
        c.mantissa[i] = a.mantissa[i+nb];   
    }
    
    //  Fix up results.  B may have trailing zeros and c may have leading zeros.
    mproun (b, lmpnw);
    mproun (c, lmpnw);  
  }

  /**
  *  This converts the MP number A to the DPE number b, accurate to between
  *  14 and 17 digits, depending on system. 
  *
  *  B will be between 1 and MPBDX.
  *  Debug output starts with MPIDB = 9.
  *  dimension a(lmpnw)
  */
  static void mpmdc(/*const*/ MP a, MPDPE b)
  {
    double aa;     
    boolean isAZero = false;
    if (a.nw==0)  
    {
      b.a = 0.0;
      b.n = 0;
      isAZero=true;
    }
    
    if(!isAZero)
    {
      int na = a.nw;
      aa = a.mantissa[0];
      if (na >= 2) aa += MPRDX * a.mantissa[1];
      if (na >= 3) aa += MPRX2 * a.mantissa[2];
      if (na >= 4) aa += MPRDX * MPRX2 * a.mantissa[3];
      
      b.n = MPNBT * a.exponent;
      b.a = fSign (aa, a.sign ? 1.0 : -1.0);
    }
  }
  
 /**
  *  This routine multiplies MP numbers A and B to yield the MP product c.
  *
  *  When one of the arguments has a much higher level of precision than the
  *  other, this routine is slightly more efficient if A has the lower level of
  *  precision.  For extra high levels of precision, use MPMULX.  Debug output
  *  starts with MPIDB = 8. 
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw), c(lmpnw+2).
  *  <br>
  *  This routine returns up to MPNW mantissa words of the product.  If the
  *  complete double-long product of A and B is desired (for example in large
  *  integer applications), then MPNW must be at least as large as the sum of
  *  the mantissa lengths of A and B.  In other words, if the precision levels
  *  of A and B are both 64 words, then MPNW must be at least 128 words to
  *  obtain the complete double-long product in c.
  *  @see #mpmulx
  */  
  static  void mpmul(/*const*/ MP a, /*const*/ MP b, MP c, int lmpnw)
  {
    int i, j;
    double t1, t2;
      
    int na = Math.min (a.nw, lmpnw);
    int nb = Math.min (b.nw, lmpnw);
    
    if (na == 0 || nb == 0) 
    {
      //  One of the inputs is zero -- result is zero.
      zero(c);
      return;
    }

    if (na == 1 && a.mantissa[0] == 1)
    {
      //  A is 1 or -1 -- result is B or -B.
      c.sign=!(a.sign^b.sign);
      c.nw=nb;
      c.exponent=a.exponent+b.exponent;
      
      for(i = 0; i<nb;i++)
        c.mantissa[i] = b.mantissa[i];
      return;
    }
    else if (nb == 1 && b.mantissa[0] == 1.) 
    {     
      //  B is 1 or -1 -- result is A or -A.
      c.sign = !(a.sign^b.sign);
      c.nw = na;
      c.exponent = a.exponent + b.exponent;
      
      for (i = 0; i<na; i++)
        c.mantissa[i] = a.mantissa[i];
      return;
    }
    
    double d[] = new double[lmpnw+4];
    
    int nc = Math.min (na + nb, lmpnw);
    double d2 = a.exponent + b.exponent; // *cast*
    
    for(i = 0; i<nc + 4; i++)
      d[i] = 0.0;
    
    //  Perform ordinary long multiplication algorithm.  Accumulate at most MPNW+4
    //  mantissa words of the product.
    
    for (j = 3; j<= na + 2; j++)
    {
      t1 = a.mantissa[j-3];
      int j3 = j - 3;
      int n2 = Math.min (nb + 2, lmpnw + 4 - j3);
      
      for(i = 2; i<n2; i++) 
        d[i+j3] += t1 * b.mantissa[i-2];
      
      //  Release carries periodically to avoid overflowing the exact integer
      //  capacity of double floating point words in D.
      
      if (((j - 2)% MPNPR) == 0) 
      {
        int i1 = Math.max (3, j - MPNPR);
        int i2 = n2 + j3;
        
        //dir$ ivdep
        for(i = i1-1; i< i2; i++)
        {
          t1 = d[i];
          t2 = (int)(MPRDX * t1);
          d[i] = t1 - MPBDX * t2;
          d[i-1] +=  t2;
        }
      }
    }
    
    //  If D(1) is nonzero, shift the result one cell right.
    
    if (d[1] != 0.0) 
    {
      d2 += 1.0;
      
      //dir$ ivdep
      for (i = nc + 3; i>=2; i--)
        d[i] = d[i-1];
      
    }
    d[0] = fSign (nc, (!(a.sign^b.sign)) ? 1.0 : -1.0);
    d[1] = d2;
    
    //  Fix up result, since some words may be negative or exceed MPBDX.
    
    mpnorm (d, c, lmpnw);
  }

 /**
  *  This routine multiplies the MP number A by the DPE number (B, N) to yield
  *  the MP product c. 
  *
  *  Debug output starts with MPIDB = 9.  
  *  <br>
  *  Dimension a(lmpnw), c(lmpnw+2).
  */
  static void mpmuld (/*const*/ MP a, /*const*/ MPDPE b, MP c, int lmpnw)
  {
    int i,k;
    double bb;
    MP f =  new MP(6, false) ;
    
    //  Check for zero inputs.
    
    int na = Math.min (a.nw, lmpnw);
    int ib = (int)(fSign(1.0, b.a)); // *cast*
    if (na == 0 || b.a == 0.0) 
    {
      zero(c);
      return;
    }
    int n1 = b.n / MPNBT;
    int n2 = b.n - MPNBT * n1;
    bb = Math.abs (b.a) * Math.pow(2.0, n2);
    
    //  Reduce BB to within 1 and MPBDX.
    
    if (bb >= MPBDX) 
    {
      
      for (k = 1; k<=100; k++)
      {
        bb = MPRDX * bb;
        if (bb < MPBDX) 
        {
          n1 += k;
          break;
        }
      }
    }
    else if (bb < 1.0) 
    {
      for (k = 1; k<= 100; k++)
      {
        bb = MPBDX * bb;
        if (bb >= 1.0) 
        {
          n1 -= k;
          break;
        }
      }
    }
    
    //  If B cannot be represented exactly in a single mantissa word, use MPMUL.
    
    if (bb != (int)(bb)) 
    {
      MPDPE dpb = new MPDPE(fSign (bb, b.a), n1 * MPNBT); 
      mpdmc (dpb, f);
      mpmul (f, a, c, lmpnw);
      return; 
    }
    double d[] = new double[lmpnw+4];

    //  Perform short multiply operation.
    
    //dir$ ivdep
    for(i = 2; i<na + 2; i++)
      d[i] = bb * a.mantissa[i-2];
    
    //  Set the exponent and fix up the result.
    
    d[0] = fSign (na, (a.sign ? 1 : -1) * ib);
    d[1] = a.exponent + n1;
    d[na+2] = 0.0;
    d[na+3] = 0.0;
    
    mpnorm (d, c, lmpnw);
  }
 
 /**
  *  This sets B equal to the integer nearest to the MP number A.  
  *
  *  Debug output starts with MPIDB = 8.
  *  <br>
  *  Dimension a(lmpnw+2), b(lmpnw+2).
  */
  static void mpnint (/*const*/ MP a, MP b, int lmpnw)
  {
    boolean ini = a.mantissa[0] == 524288;
    /*    if(ini)
       {
	  a.debug("a");
	  b.debug("b");
       }*/

    int i;
    MP f =  new MP(6, false) ;
    MP s = new MP(lmpnw+2,false);
    
    int na = Math.min (a.nw, lmpnw);
    if (na == 0)  
    {
      //  A is zero -- result is zero.
      zero(b);
      return;
    }
    if (a.exponent >= lmpnw) 
    {
      //  A cannot be represented exactly as an integer.
      throw new ArithmeticException(
        "mpnint: Argument is too large --> " + a);
   }
    
    f.nw=1;
    f.sign=true;
    f.exponent=-1;
    f.mantissa[0]=(float)(0.5 * MPBDX); // *cast*
    f.mantissa[1]=0;
    
    //  Add or subtract 1/2 from the input, depending on its fSign.
    if (a.sign) 
      mpadd (a, f, s, lmpnw);
    else 
      mpsub (a, f, s, lmpnw);
    

    int ic = s.sign ? 1 : -1; 
    int nc = s.nw;
    int mc = s.exponent;

    //  Place integer part of S in B.
    
    int nb = Math.min (Math.max (mc + 1, 0), nc);
    if (nb == 0) 
      zero(b);
    else 
    {
      b.nw=nb;
      b.sign=(ic>=0);
      b.exponent=mc;
      b.mantissa[nb] = 0;
      b.mantissa[nb+1] = 0;
      
      for(i = 0; i<nb; i++)
        b.mantissa[i] = s.mantissa[i];
    }
  }  

  /**
  *  This converts the MP number in array D of MPCOM4 to the standard
  *  normalized form in A.  
  *
  *  The MP routines often leave negative numbers or values exceeding the 
  *  radix MPBDX in result arrays, and this fixes them.
  *  MPNORM assumes that two extra mantissa words are input at the end of D.
  *  This reduces precision loss when it is necessary to shift the result to
  *  the left.  This routine is not intended to be called directly by the user.
  *  Debug output starts with MPIDB = 10. 
  *  <br>
  *  Dimension d(lmpnw+4), a(lmpnw+2).
  *  d is a double array with David's MP number format. 
  */
  private static void mpnorm(double d[], MP a, int lmpnw)
  {
    final boolean risc = true; 
    double t1, t2, t3;
    int i;       
    int ia = (int)(fSign (1.0, d[0])); // sign
    int na = Math.min ((int)(Math.abs (d[0])), lmpnw);
    if (na == 0)  
    {
      zero(a);
      return;
    }
    int n4 = na + 4;
    double a2 = d[1]; 
    d[1] = 0.0;
    
    boolean needToNormalize=true;
    while(needToNormalize)
    {
      /****  110 ****/
      // !>
      //  Try a vectorized fixup loop three times, unless A is very short.  This
      //  should handle 99% of the inputs.  On RISC computers, it is more
      //  efficient to completely bypass this loop, by defining RISC = 1 
      
      boolean breakLoop = false;
      if(!risc)
      {
        double s1;
        int k;
        if (na > 8) 
        { 
          for(k = 1; k<=3; k++)
          {
            s1 = 0.0;
            
            //dir$ ivdep
            for (i = 2; i<n4; i++)
            {
              t2 = MPRDX * d[i];
              t1 = (int)(t2);
              if (t2 < 0.0 && t1 != t2) t1--;
              d[i] -= t1 * MPBDX;
              d[i-1] += t1;
              s1 += Math.abs (t1);
            }
            
            if (s1 == 0.0) 
            {
              breakLoop=true;
              break;
            }
          }
        }
      }
      
      
      //  Still not fixed - use recursive loop.  This loop is not vectorizable,
      //  but it is guaranteed to complete the job in one pass.
      if(!breakLoop) //here that changes
      { 
        t1=0;
        for(i = n4-1; i>=2; i--)
        {
          t3 = t1 + d[i];
          t2 = MPRDX * (t3);
          t1 = (int)(t2);
          if (t2 < 0.0 && t1 != t2) t1--;
          d[i] = t3 - t1 * MPBDX;
        } 
        d[1] += t1;
      }
      
      if (d[1] < 0.0) 
      { 
        //  D(1) is negative -- negate all words and re-normalize.
        ia = - ia;
        d[2] += MPBDX * d[1];
        d[1] = 0.0;
        
        for(i = 1; i< n4; i++)
          d[i] = - d[i];
      }
      else if (d[1] > 0.0) 
      {
        
        //  The fixup loops above "spilled" a nonzero number into D(1).  Shift the
        //  entire number right one cell.  The exponent and length of the result
        //  are increased by one.
        
        for (i = n4-2; i>=1; i--)
          a.mantissa[i-1] = (float)d[i]; // *cast*
        
        na = Math.min (na + 1, lmpnw);
        
        a2++;
        needToNormalize=false;
      } 
      else  //here the last one
      {
        for (i = 2; i< n4; i++)
          a.mantissa[i-2] = (float)d[i]; // *cast*
        needToNormalize=false;
      }
    }
    
    //  Perform rounding and truncation.
    a.nw=na;
    a.sign=(ia>=0);
    a.exponent = (int)(a2);
    
    mproun (a, lmpnw);
  }

  /**
  *  This computes the N-th power of the MP number A and returns the MP result
  *  in B. 
  *
  *  When N is zero, 1 is returned.  When N is negative, the reciprocal
  *  of A ^ |N| is returned.  For extra high levels of precision, use MPNPWX.
  *  Debug output starts with MPIDB = 7.  
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw+2).
  *  <br>
  *  This routine employs the binary method for exponentiation.
  *  @see #mpnpwx
  */
  static void mpnpwr (/*const*/ MP a, int n, MP b, int lmpnw)
  {
    int j;
    double t1;
    int lmpnw3 = lmpnw+3;
    MP f1 =  new MP(6, false) ;
    MP sk0 = new MP(lmpnw3,false);
    MP sk1 = new MP(lmpnw3,false);
    
    int na = Math.min (a.nw, lmpnw);
    if (na == 0) 
    {
      if (n >= 0) 
      {
        zero(b);
        return;
      } 
      else 
      {
        throw new ArithmeticException(
          "mpnpwr: Argument is zero and n is negative or zero --> " +
          a + "\n" + n);
      }
    }
    
    int nws = lmpnw;
    lmpnw++;
    int nn = Math.abs (n);
    f1.sign=true;
    f1.nw=1;
    f1.exponent=0;
    f1.mantissa[0]=1;
    f1.mantissa[1]=0;
    
    boolean skip=false;
    
    switch(nn)
    {
    case 0:
      mpeq (f1, b, lmpnw);
      return;
      
    case 1:
      mpeq (a, b, lmpnw);
      skip=true;
      break;
      
    case 2:
      mpmul (a, a, sk0, lmpnw);
      mpeq (sk0, b, lmpnw);
      skip=true;
      break;
    }
    
    if(!skip)
    {
      //  Determine the least integer MN such that 2 ^ MN > NN.
      t1 = nn;
      int mn = (int)(CL2 * Math.log (t1) + 1.0 + MPRXX); // cast
      mpeq (f1, b, lmpnw);
      mpeq (a, sk0, lmpnw);
      int kn = nn;
      
      //  Compute B ^ N using the binary rule for exponentiation.
      
      for (j = 1; j<= mn; j++)
      {
        int kk = kn / 2;
        if (kn != 2 * kk) 
        {
          mpmul (b, sk0, sk1, lmpnw);
          mpeq (sk1, b, lmpnw);
        }
        kn = kk;
        if (j < mn) 
        {
          mpmul (sk0, sk0, sk1, lmpnw);
          mpeq (sk1, sk0, lmpnw);
        }
      }
    }
    
    //  Compute reciprocal if N is negative.
    if (n < 0) 
    {
      mpdiv (f1, b, sk0, lmpnw);
      mpeq (sk0, b, lmpnw);
    }
    mproun (b, nws);    
  }

  /**
  *  This computes the N-th root of the MP number A and returns the MP result
  *  in B.  
  *  
  *  N must be at least one and must not exceed 2 ^ 30.  For extra high
  *  levels of precision, use MPNRTX.  Debug output starts with MPIDB = 7.  
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw+2).
  *  <br>
  *  This subroutine employs the following Newton-Raphson iteration, which
  *  converges to A ^ (-1/N): 
  *  <pre>
  *   X_{k+1} = X_k + (X_k / N) * (1 - A * X_k^N)
  *  </pre>
  *  The reciprocal of the final approximation to A ^ (-1/N) is the N-th root.
  *  These iterations are performed with a maximum precision level MPNW that
  *  is dynamically changed, approximately doubling with each iteration.
  *  See the comment about the parameter NIT in MPDIVX. 
  *  <br>
  *  When N is large and A is very near one, the following binomial series is
  *  employed instead of the Newton scheme: 
  *  <pre>
  *  (1 + x)^(1/N)  =  1  +  x / N  +  x^2 * (1 - N) / (2! N^2)  +  ...
  *  </pre>
  *  See the comment about the parameter NIT in MPDIVX.
  *  @see #mpnrtx
  *  @see #mpdivx
  */  
  static void mpnrt (/*const*/ MP a, int n, MP b, int lmpnw)
  {
    int k;
    double t2, tn;
    int lmpnw3 = lmpnw+3;
    MP f1 =  new MP(6, false) ;
    MP f2 =  new MP(6, false) ;
    MP sk0 = new MP(lmpnw3,false), sk1 = new MP(lmpnw3,false), 
      sk2 = new MP(lmpnw3,false), sk3 = new MP(lmpnw3,false);

    int na = Math.min (a.nw, lmpnw);
    
    if (na == 0) 
    {
      zero(b);
      return;
    }
    
    if (!a.sign) 
    {
      throw new ArithmeticException(
        "mpnrt: Argument is negative -->" + a);
    }
    
    if (n <= 0 || n > N30) 
    {
      throw new ArithmeticException(
        "mpnrt: Improper value of n -->" + n);
    }
    
    //  If N = 1, 2 or 3,  MPEQ, MPSQRT or MPCBRT.  These are faster.
    switch (n)
    {
    case 1:
      mpeq (a, b, lmpnw);
      return;
    case 2:
      mpsqrt (a, b, lmpnw);
      return;
    case 3:
      mpcbrt (a, b, lmpnw);
      return;
    }
    
    int nws = lmpnw;
    f1.nw=1;
    f1.sign=true;
    f1.exponent=0;
    f1.mantissa[0]=1;
    f1.mantissa[1]=0;
    
    //  Determine the least integer MQ such that 2 ^ MQ >= MPNW.
    MPDPE t1 = new MPDPE();
    t1.a = lmpnw;
    int mq = (int)(CL2 * Math.log (t1.a) + 1.0 - MPRXX); // *cast*
    
    //  Check how close A is to 1.
    mpsub (a, f1, sk0, lmpnw);
    if (sk0.nw == 0) 
    {
      mpeq (f1, b, lmpnw);
      return;
    }
    
    mpmdc (sk0, t1);
    int n2 = (int)(CL2 * Math.log (Math.abs (t1.a))); // *cast*
    t1.a *= Math.pow(0.5, n2);
    t1.n += n2;
    if (t1.n <= -30) 
    {
      t2 = n;
      n2 = (int)(CL2 * Math.log (t2) + 1.0 + MPRXX); // *cast*
      int n3 = - MPNBT * lmpnw / t1.n;
      if (n3 < 1.25 * n2) 
      {
        //  A is so close to 1 that it is cheaper to use the binomial series.
        lmpnw++;
        mpdivd (sk0, new MPDPE(t2), sk1, lmpnw);
        mpadd (f1, sk1, sk2, lmpnw);
        k = 0;
        
        int temp = t1.n;
        t1.n=0;
        do 
        {
          k++;
          t1.a = 1 - k * n;
          t2 = (k + 1) * n;
          mpmuld (sk1, t1, sk3, lmpnw);
          mpdivd (sk3, new MPDPE(t2), sk1, lmpnw);
          mpmul (sk0, sk1, sk3, lmpnw);
          mpeq (sk3, sk1, lmpnw);
          mpadd (sk1, sk2, sk3, lmpnw);
          mpeq (sk3, sk2, lmpnw);
        }
        while (sk1.nw != 0 && sk1.exponent >= - lmpnw); 
        
        t1.n = temp;
        mpeq (sk2, b, lmpnw);
        mpdiv (f1, sk2, sk0, lmpnw);
        mproun (b, nws);
        return;
      }
    }
    
    //  Compute the initial approximation of A ^ (-1/N).
    tn = n;
    MPDPE dpn = new MPDPE(n, 0);
    mpmdc (a, t1);
    MPDPE dp1 = new MPDPE();
    dp1.n = (int)(- t1.n / tn); // *cast*
    dp1.a = Math.exp (-1.0 / tn * (Math.log (t1.a) + (t1.n + tn * dp1.n) * ALT));
    mpdmc (dp1, b);
    mpdmc (dpn, f2);
    lmpnw = 3;
    int iq = 0;
    
    //  Perform the Newton-Raphson iteration described above with a dynamically
    //  changing precision level MPNW (one greater than Math.powers of two).
    for (k = 2; k<= mq;k++)
    {
      lmpnw = Math.min (2 * lmpnw - 2, nws) + 1;
      boolean loop=true;
      while (loop)
      {
        mpnpwr (b, n, sk0, lmpnw);
        mpmul (a, sk0, sk1, lmpnw);
        mpsub (f1, sk1, sk0, lmpnw);
        mpmul (b, sk0, sk1, lmpnw);
        mpdivd (sk1, new MPDPE(tn), sk0, lmpnw);
        mpadd (b, sk0, sk1, lmpnw);
        mpeq (sk1, b, lmpnw);
        if (k == mq - NIT && iq == 0) 
          iq = 1;
        else loop=false;
      }
    }
    //  Take the reciprocal to give final result.
    mpdiv (f1, b, sk1, lmpnw);
    mpeq (sk1, b, lmpnw);
    mproun (b, nws);
  }

  /**
  *  This returns a pseudo-random MP number A between 0 and 1. 
  *
  *  Debug output starts with MPIDB = 9.
  *  <br>
  *  Dimension a(lmpnw+2).
  */
  static void mprand (MP a, int lmpnw)
  {
    int i;
    final double f7 = 78125;
    double t1, t2;
    
    synchronized(randMutex)
    {
      if(r30 == 0) 
      {
        r30 = 1.0;
        t30 = 1.0;
        
        for(i = 1; i<=30; i++)
        {
          r30 = 0.50 * r30;
          t30 = 2.0 * t30;
        }
      }
      
      a.nw=lmpnw;
      a.sign=true;
      a.exponent=-1;
      
      for(i = 0; i<=lmpnw + 1; i++)
      {
        t1 = f7 * sd;
        t2 = (int)(r30 * t1);
        sd = t1 - t30 * t2;
        a.mantissa[i] = (float)((int)(MPBDX * r30 * sd)); // *cast*
      }
    }
    mproun (a, lmpnw);  
  }

  
  /**
  *  This performs rounding and truncation of the MP number A.  
  *  
  *  It is called by MPNORM, and also by other subroutines when the precision level is
  *  reduced by one.  It is not intended to be directly called by the user.
  *  <br>
  *  Dimension a(lmpnw+2).
  *  <br>
  *  The parameter AMX is the absolute value of the largest exponent word
  *  allowed for MP numbers.
  *  @see #mpnorm
  */
  static void mproun(MP a, int lmpnw)
  {
    final double amx = 2.e6;
    int i;
    
    //  Check for initial zeroes.  
    int a2 = a.exponent;
    a.exponent = 0;
    int na = Math.min (a.nw, lmpnw);
    int n1 = na + 1;
    if (a.mantissa[0] == 0) 
    {
      //  Find the first nonzero word and shift the entire number left.  The length
      //  of the result is reduced by the length of the shift. 
      boolean allZero=true;
      for (i = 1; i<=n1; i++)
      {
        if (a.mantissa[i] != 0) 
        {
          allZero=false;
          break;          // bug fixed here
        }
      }
      if(allZero)
      {
        zero(a);
        return;
      }
      int k = i; // bug fixed here
      
      //dir$ ivdep
      for (i = 0; i<= n1 - k; i++)
        a.mantissa[i] = a.mantissa[i+k];
      
      
      a2 -= k;
      na -= Math.max (k - 2, 0);
    }
    
    //  Perform rounding depending on mpird.
    
    if (na == lmpnw && mpird >= 1) 
    {
       if ((mpird == 1) && (a.mantissa[na] >= 0.5 * MPBDX) || 
        (mpird == 2) && (a.mantissa[na] >= 1))
        a.mantissa[na-1]++;
      
      
      //  Release carries as far as necessary due to rounding.
      boolean loopBreak = false;
      for(i = na - 1; i>=0; i--) // bug fixed here
      {
        if (a.mantissa[i] < MPBDX) 
        {
          loopBreak=true;
          break;
        }
        a.mantissa[i] -= (float)MPBDX;  // *cast*
        if(i != 0)
          a.mantissa[i-1]++;
        else
          a.exponent++;
      }
      
      //  Release of carries due to rounding continued all the way to the start --
      //  i.e. number was entirely 9's.
      if(!loopBreak)
      {
        a.mantissa[0] = (float)a.exponent;//a[0] = float(a.exponent + 1); // *cast*
        na = 1;
        a2++;
      }
    }
    
    try 
    {
       if (a.mantissa[na-1] == 0) 
       {
	  
	  //  At least the last mantissa word is zero.  Find the last nonzero word
	  //  and adjust the length of the result accordingly.
	  
	  boolean allZero=true;
	  for (i = na - 1; i>=0; i--) // bug fixed here
	  {
	     if (a.mantissa[i] != 0)
	     {
		allZero=false;
		break;
	     }
	  } 
	  if(allZero)
	  {
	     zero(a);
	     return;
	  }
	  na = i + 1; // bug fixed here
       }
    }
    catch(ArrayIndexOutOfBoundsException e) {}

    //  Check for overflow and underflow.
    
    if (a2 < - amx) 
    {
      throw new ArithmeticException(
        "mproun: Exponent underflow.");
    } 
    else if (a2 > amx) 
    {
      throw new ArithmeticException(
        "mproun: Exponent overflow.");
    }
    
    //  Check for zero.
    if (a.mantissa[0] == 0) 
      zero(a);
    else 
    {
      a.nw=na;
      a.exponent=a2;
      a.mantissa[na] = 0;
      a.mantissa[na+1] = 0;
    }
  }
  
  /**
  *  This computes the square root of the MP number A and returns the MP result
  *  in B.
  *  
  *  For extra high levels of precision, use MPSQRX.  Debug output
  *  starts with MPIDB = 7. 
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw+2).
  *  <br>
  *  This subroutine employs the following Newton-Raphson iteration, which
  *  converges to 1 / Sqrt(A):  
  *  <pre>
  *   X_{k+1} = X_k + 0.5 * (1 - X_k^2 * A) * X_k
  *  </pre>   
  *  where the muliplication () * X_k is performed with only half of the
  *  normal level of precision.  These iterations are performed with a
  *  maximum precision level MPNW that is dynamically changed, doubling with
  *  each iteration.  The final iteration is performed as follows (this is
  *  due to A. Karp):
  *  <pre>
  *   Sqrt(A) = (A * X_n) + 0.5 * [A - (A * X_n)^2] * X_n  (approx.)
  *  </pre>
  *  where the multiplications A * X_n and [] * X_n are performed with only
  *  half of the final level of precision.  See the comment about the parameter
  *  NIT is MPDIVX.
  *  @see #mpsqrx
  *  @see #mpdivx
  */
  static void mpsqrt(/*const*/ MP a, MP b, int lmpnw)
  {
    int k;
    double t2;
    int lmpnw3 = lmpnw+3;
    MP f =  new MP(6, false) , 
      sk0 = new MP(lmpnw3,false), sk1 = new MP(lmpnw3,false), 
      sk2 = new MP(lmpnw3,false);
    
    int na = Math.min (a.nw, lmpnw);
    
    if (na == 0) 
    {
      zero(b);
      return;
    }
    
    if (!a.sign) 
    {
      throw new ArithmeticException(
        "mpsqrt: Argument is negative --> " + a);
    }
    int nws = lmpnw;
    
    //  Determine the least integer MQ such that 2 ^ MQ >= MPNW.
    MPDPE t1 = new MPDPE(lmpnw,0);
    int mq = (int)(CL2 * Math.log (t1.a) + 1.0 - MPRXX); // *cast* 
    
    //  Compute the initial approximation of 1 / Sqrt(A).
    int iq;
    
    mpmdc (a, t1);
    MPDPE dp1 = new MPDPE();
    
    dp1.n = - t1.n / 2;
    t2 = Math.sqrt (t1.a * Math.pow(2.0, (t1.n + 2 * dp1.n)));
    dp1.a = 1.0/t2; 
    mpdmc (dp1, b);
    
    f.nw=1;
    f.sign=true;
    f.exponent=0;
    f.mantissa[0]=1;
    f.mantissa[1]=0;
    lmpnw = 3;
    iq = 0;
    
    //  Perform the Newton-Raphson iteration described above with a dynamically
    //  changing precision level MPNW (one greater than Math.powers of two). 
    int nw1, nw2;
    for(k = 2; k<= mq - 1; k++)
    {
      nw1 = lmpnw;
      lmpnw = Math.min (2 * lmpnw - 2, nws) + 1;
      nw2 = lmpnw;
      boolean stop=false;
      while(!stop)
      {
        mpmul (b, b, sk0, lmpnw);
        mpmul (a, sk0, sk1, lmpnw);
        mpsub (f, sk1, sk0, lmpnw);
        lmpnw = nw1;
        mpmul (b, sk0, sk1, lmpnw);
        mpmuld (sk1, new MPDPE(0.50), sk0, lmpnw);
        lmpnw = nw2;
        mpadd (b, sk0, sk1, lmpnw);
        mpeq (sk1, b, lmpnw);
        if (k == mq - NIT && iq == 0) 
          iq = 1;
        else
          stop=true;
      }
    }
    
    //  Perform last iteration using Karp's trick.  
    mpmul (a, b, sk0, lmpnw);
    nw1 = lmpnw;
    lmpnw = Math.min (2 * lmpnw - 2, nws) + 1;
    nw2 = lmpnw;
    mpmul (sk0, sk0, sk1, lmpnw);
    mpsub (a, sk1, sk2, lmpnw);
    lmpnw = nw1;
    mpmul (sk2, b, sk1, lmpnw);
    mpmuld (sk1, new MPDPE(0.50), sk2, lmpnw);
    lmpnw = nw2;
    mpadd (sk0, sk2, sk1, lmpnw);
    mpeq (sk1, b, lmpnw);
    
    //  Restore original precision level.
    mproun (b, nws);  
  }

  /**
  *  This routine subtracts MP numbers A and B to yield the MP difference c,
  *  by negating B and adding.  
  *
  *  Debug output starts with MPIDB = 9.
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw), c(lmpnw+2).
  */
  static void mpsub (/*const*/ MP a, /*const*/ MP b, MP c, int lmpnw)
  {
    int i;
    
    //  Check if A = B.  This is necessary because A and B might be same array,
    //  in which case negating B below won't work.
    
    
    
    
    /* Sweny
    // check if A == B points to the same object 
    if(a == b) 
    {
      zero(c);
      return;
    }
    
    // check if their mantissas are the same
    if(a.sign == b.sign && a.nw == b.nw && a.exponent == b.exponent)
    {
      for(i=0; i<a.nw;i++)
      {
        if(a.mantissa[i] != b.mantissa[i]) break;
      }
      if (i==a.nw)
      {
        zero(c);
        return;
      }
    }
    */
    
    //  Make bb shares b's mantissa and negate bb's sign.
    MP bb = new MP(0,false);
    bb.sign=!b.sign;
    bb.nw=b.nw;
    bb.exponent=b.exponent;
    bb.mantissa=b.mantissa;
    
    mpadd (a, bb, c, lmpnw);
    bb.mantissa = null;
  }

 /**
  *  Assign 0 to MP number in
  */
 static void zero(MP in)
 {in.nw=0; in.sign=true; in.exponent=0;}
  
 /**
  *   Return the DPE value of this
  */
  MPDPE toDPE()
  { MPDPE res = new MPDPE(); mpmdc(this, res); return res;}


 /******************************  IO FUNCTIONS ******************************/

 /**
  *  Converts the char array A of length N into the MP number B.  
  *
  *  The string A must be in the format '10^s a x tb.c' where a, b and c are digit
  *  strings; s and t are '-', '+' or blank; x is either 'x' or '*'.  Blanks may
  *  be embedded anywhere.  The digit string a is limited to nine digits and
  *  80 total characters, including blanks.  The exponent portion (i.e. the
  *  portion up to and including x) and the period may optionally be omitted.
  *  Debug output starts with MPIDB = 7.
  *  <br>
  *  Dimension a(n), b(lmpnw+2).
  */
  static void mpinpc (char a[], int n, MP b, int lmpnw)
  { 
    int i, j, k, is, id; 
    double bi;
    char ai;
    
    char ca[] = new char[81];
    
    int mpnw3 = lmpnw+3;
    MP f =  new MP(6, false) , 
      sk0 = new MP(mpnw3,false), 
      sk1 = new MP(mpnw3,false), 
      sk2 = new MP(mpnw3,false);
    int nws = lmpnw++;
    int i1 = 0;
    int nn = 0;
    
    //  Find the carat, period, plus or minus sign, whichever comes first.
    
    boolean caretFound=false;
    
    for (i = 0; i<n; i++)
    {
      ai = a[i];
      if (ai == '^')
      {
        caretFound=true;
        break;
      }
      if (ai == '.' || ai == '+' || ai == '-') 
        break;
    }
    
    for(j=0;j<81;j++)
      ca[j]='\0';
    
    if (caretFound)
    {
      //  Make sure number preceding the caret is 10.
      int i2 = i-1;
      if (i2 > 79)
      {
        throw new NumberFormatException(
          "mpinpc: Syntax error in literal string.");
      }
      
      j=0;
      for(i = 0; i<=i2; i++)
      {
        ai = a[i];
        if (ai == ' ')
          continue;
        else if (!Character.isDigit(ai))
        {
          throw new NumberFormatException(
           "mpinpc: Syntax error in literal string.");
        }
        ca[j++] = ai;
      }
      
      if (ca[0]!='1' || ca[1]!='0')
      {
        throw new NumberFormatException(
           "mpinpc: Syntax error in literal string.");
      }
      i1 = i2 + 2;
      
      
      boolean exit=true;
      //  Find the x or *.
      for(i = i1; i< n; i++)
      {
        ai = a[i];
        if (ai == 'x' || ai == '*')
        {
          exit=false;
          break;
        }
      }
      if(exit)
      {
         throw new NumberFormatException(
           "mpinpc: Syntax error in literal string.");
      }
      
      
      //  Convert the exponent.
      
      i2 = i - 1;
      int l1 = i2 - i1;
      if (l1 > 79) 
      {
        throw new NumberFormatException(
           "mpinpc: Syntax error in literal string.");
      }
      
      id = 0;
      is = 1;
      
      j=0;
      for (i = 0; i<=l1; i++)
      {
        ai = a[i+i1];
        if (ai == ' ' || ai == '+') 
          continue;
        else if (ai == '-' && id == 0) 
        {
          id = 1;
          is = -1;
        } 
        else 
        {
          if (!Character.isDigit(ai))
          {
            throw new NumberFormatException(
              "mpinpc: Syntax error in literal string.");
          }
          id = 1; 
          ca[j++] = ai;
        }
      }
      
      ca[j]='\0';
      nn=Integer.parseInt(new String(ca,0,j));
      
      nn = is * nn;
      i1 = i2 + 2;
    }
    //  Find the next nonblank character.
    boolean exit=true;
    
    for(i = i1; i< n; i++)
    {
      if (a[i] != ' ')
      {
        exit=false;
        break;
      }
    }
    
    if(exit)
    {
      throw new NumberFormatException(
        "mpinpc: Syntax error in literal string.");
    }

    //  Check if the nonblank character is a plus or minus SIGN.
    
    i1 = i;
    if (a[i1] == '+') 
    {
      i1 = i1 + 1;
      is = 1;
    } 
    else if (a[i1] == '-') 
    {
      i1 = i1 + 1;
      is = -1;
    } 
    else 
      is = 1;
    
    int nb = 0;
    int ib = 0;
    id = 0;
    int ip = 0;
    zero(sk2);
    f.nw=1;
    f.sign=true;
    f.exponent=0;
    int it = 0;
    
    int mm;
    boolean cont=true;
    while(cont)
    {
      ip = 0;
      
      for(mm=0; mm<6; mm++)
        ca[mm]='0';
      
      //  Scan for digits, looking for the period also.  On the first pass we just
      //  count, so that on the second pass it will come out right.
      
      for(i = i1; i< n; i++)
      {
        ai = a[i];
        if (ai == ' ') 
        {
        } 
        else if (ai == '.') 
        {
          if (ip != 0)
          {
            throw new NumberFormatException(
              "mpinpc: Syntax error in literal string.");
          }
          ip = id;
        } 
        else if (!Character.isDigit(ai))
        {
          throw new NumberFormatException(
              "mpinpc: Syntax error in literal string.");
        } 
        else 
        {
          id++;
          ca[ib++] = ai;
        }
        if (ib == 6 || i == (n-1) && ib != 0) 
        {
          if (it != 0) 
          {
            nb++;
            ca[ib]='\0';
            
            bi=Integer.parseInt(new String(ca,0,ib));
            
            mpmuld (sk2, new MPDPE(1e6), sk0, lmpnw);
            
            if (bi != 0) 
            {
              f.nw=1; f.sign=true;
              f.mantissa[0]=(float)bi; // *cast*
            } 
            else 
            {
              f.nw=0; f.sign=true;
            }
            mpadd (sk0, f, sk2, lmpnw);
            for(mm=0; mm<6; mm++)
              ca[mm]='0';
          }
          if ((i+1) != n) ib = 0;
        }
      }
      
      if (it == 0) 
      {
        ib = 6 - ib;
        if (ib == 6) 
          ib = 0;
        it = 1;
      }
      else
        cont=false;
    }
    
    
    if (is == -1) 
      sk2.sign = ! sk2.sign;
    
    if (ip == 0) 
      ip = id;
    
    nn += ip - id;
    f.nw=1; f.sign=true;
    f.mantissa[0]=10;
    mpnpwr (f, nn, sk0, lmpnw);
    mpmul (sk2, sk0, sk1, lmpnw);
    mpeq (sk1, b, lmpnw);
    mproun (b, nws);    
  }
  
  /**
  *  Converts the MP number A into character form in the char array B.
  *
  *  In other words, B is contained in B[0], ..., B[N-1].  The format is 
  *  analogous to the Fortran exponential format (E format), except that 
  *  the exponent is placed first.
  *  <br>
  *  Dimension a(lmpnw), b(7.225 * MPNW + 30).
  *  <br>
  *  This routine is called by MPOUT, but it may be directly called by the user
  *  if desired for custom output. 
  *
  *  @return  the length of the output.  
  *  @see #mpout
  *  @see #mpoutx
  */
  static int mpoutc(/*const*/ MP a, char b[], int lmpnw)
  {
    int i, j,k, nn, n; 
    double aa, t1;
    
    int mpnw3 = lmpnw+3;
    MP f =  new MP(6, false) , 
      sk0 = new MP(mpnw3,false), sk1 = new MP(mpnw3,false);
    
    int na = Math.min (a.nw, lmpnw);
    lmpnw++;
    f.sign=true; f.nw=1; f.exponent=0;
    f.mantissa[0]=10;
    
    //  Determine exact power of ten for exponent.
    int nx;
    if (na != 0) 
    {
      aa = a.mantissa[0];
      if (na >= 2) aa += MPRDX * a.mantissa[1];
      if (na >= 3) aa += MPRX2 * a.mantissa[2];
      if (na >= 4) aa += MPRDX * MPRX2 * a.mantissa[3];
      t1 = AL2 * MPNBT * a.exponent + log10 (aa);
      if (t1 >= 0.0) 
        nx = (int)t1; // *cast*
      else 
        nx = (int)(t1 - 1.0); // *cast*
      
      mpnpwr (f, nx, sk0, lmpnw);
      mpdiv (a, sk0, sk1, lmpnw);
      
      //  If we didn't quite get it exactly right, multiply or divide by 10 to fix.
      boolean cont = true;
      while (cont) 
      {
        if (sk1.exponent < 0)
        {
          nx = nx - 1;
          mpmuld (sk1, new MPDPE(10.0), sk0, lmpnw);
          mpeq (sk0, sk1, lmpnw);
        }
        else if(sk1.mantissa[0] >= 10) 
        {
          nx++;
          mpdivd (sk1, new MPDPE(10.0), sk0, lmpnw);
          mpeq (sk0, sk1, lmpnw);
        }
        else
          cont = false;
      }
      
      sk1.sign = true; 
    }
    else 
      nx = 0;
    
    
    //  Place exponent first instead of at the very end as in Fortran.
    
    b[0] = '1';
    b[1] = '0';
    b[2] = ' ';
    b[3] = '^';
    
    char ca[]=String.valueOf(nx).toCharArray();
 
    int len = ca.length;
    int blank = 14-len;
    for(i=4;i<blank;i++)
      b[i]=' ';
    
    for(i = 0;i< len; i++)
      b[blank+i] = ca[i]; 
    
    b[14] = ' ';
    b[15] = 'x';
    b[16] = ' ';
    
    
    //  Insert sign and first digit.
    
    if (a.sign == false) 
      b[17] = '-';
    else  
      b[17] = ' ';
    
    if (na != 0) 
      nn = (int)sk1.mantissa[0]; // *cast*
    else 
      nn = 0;
    
   ca = String.valueOf(nn).toCharArray();
    
    b[18] = ca[0];
    b[19] = '.';
    
    int ix = 20;
    if (na == 0) 
    {     
      b[ix]='\0';
      return ix;
    }
    
    f.mantissa[0] = (float)nn; // *cast*
    mpsub (sk1, f, sk0, lmpnw);
    
    if (sk0.nw == 0) 
    {     
      b[ix]='\0';
      return ix;
    }
    
    mpmuld (sk0, new MPDPE(1e6), sk1, lmpnw);
    
    int nl = (int)(Math.max (lmpnw * log10 (MPBDX) / 6.0 - 1.0, 1.0)); // *cast*
    
    //  Insert the digits of the remaining words.
    boolean skip = false;
    for (j = 1; j<=nl; j++)
    {
      if (sk1.exponent == 0.) 
      {
        nn = (int)(sk1.mantissa[0]); // *cast*
        f.nw=1; f.sign=true; 
        f.mantissa[0] = (float) nn; // *cast*
      } 
      else 
      {
        f.nw=0; f.sign=true;
        nn = 0;
      }
      
      //sprintf(ca, "%06d", nn);
      ca = String.valueOf(nn).toCharArray();
      for(i = 0; i<6-ca.length;i++)
        b[i+ix] = '0';
      k=0;
      for(;i<6;i++)
        b[i+ix] = ca[k++];
      
      ix  += 6;
      mpsub (sk1, f, sk0, lmpnw);
      mpmuld (sk0, new MPDPE(1e6), sk1, lmpnw);
      if (sk1.nw == 0) 
      {
        skip = true;
        break;
      }
    }
    //  Check if trailing zeroes should be trimmed.
    
    if(!skip)
      j = nl + 1;
    
    int l = --ix;
    if (b[l] == '0' || (j > nl && b[l-1] == '0' && 
      b[l-2] == '0' && b[l-3] == '0')) 
    {
      b[l] = '\0';
      boolean loopbreak = false;
      for (i = l - 1; i>=20; i--)
      {
        if (b[i] != '0') 
        {
          ix = i;
          loopbreak = true;
          break;
        }
        b[i] = '\0';
      }
      if(!loopbreak)
        ix = 20;
      
      //  Check if trailing nines should be rounded up.
      
    } 
    else if (j > nl && b[l-1] == '9' && b[l-2] == '9' 
      && b[l-3] == '9') 
    {
      b[l] = '\0';
      
      skip = false;
      for (i = l - 1; i>=20; i--)
      {
        if (b[i] != '9') 
        {
          skip = true;
          break;
        }
        b[i] = '\0';
      }
      
      //  We have rounded away all digits to the right of the decimal point, and the
      //  digit to the left of the digit is a 9.  Set the digit to 1 and increase
      //  the exponent by one.
      
      if(!skip)
      {
        ix = 20;
        if (b[18] == '9') 
        {
          b[18] = '1';
          ca = String.valueOf(nx+1).toCharArray();
          k=0;
          for (i = 0; i<10-ca.length; i++)
            b[i+4] = ' ';
          
          for(;i<10;i++)
            b[i+4] = ca[k++];   
        } 
        else 
        {
          ca[0] = b[18]; ca[1]='\0';
          nn = Integer.parseInt(new String(ca,0,1));   
          ca = String.valueOf(nn+1).toCharArray();      
          b[18] = ca[0];
        }
      }
      
      else
      {
        ca[0] = b[i]; ca[1]='\0';
        nn = Integer.parseInt(new String(ca,0,1));
        ca = String.valueOf(nn+1).toCharArray();      
        
        b[i] = ca[0];
        ix = i;
      }
    }
    
    n = ix;
    b[++n]='\0';
    return n;
  }
  
  
 /**
  *  This routine converts the char array  A, which
  *  represents a multiprecision number in C++ style, i.e.
  *  "1234567890" or "1.23456789e-21", into standard MP binary format.
  *
  *  This routine is not intended to be called directly by the user.
  *  <br>
  *  Dimension b(lmpnw+2)
  */
  static void mpdexc (/*const*/ char a[], int l, MP b, int lmpnw)
  {
    int i;
    boolean foundExponent = false;
    for( i = 0; i< l; i++)
    {
      if (a[i] == 'D' || a[i] == 'E' || a[i] == 'd' ||
        a[i] == 'e') 
      {
        foundExponent = true;
        break;
      }
    }
    if(!foundExponent)
    {
      mpinpc (a, l, b, lmpnw);
      return;
    }
    
    char c[] = new char[mpipl+101];
    int i1 = i + 1; 
    int l1 = i; 
    int l2 = l - i1; 
    c[0] = '1';
    c[1] = '0';
    c[2] = '^';
    
    for(i = 0; i<l2; i++)
      c[i+3] = a[i+i1];
    
    c[l2+3] = 'x';
    
    for(i = 0; i<l1; i++)
      c[i+l2+4] = a[i];
    c[i+l2+4]='\0';

    mpinpc (c, l1 + l2 + 4, b, lmpnw);
  }

 /**
  *  This function is a wrapper for mpoutc. 
  *  
  *  cs is the output char array.  cs must be dimensioned at least la + 25. 
  *  A comma is placed at the end of the last line to denote the end of
  *  the mp number. Here is an example of the output:
  *  <pre>
  *   10 ^        -4 x  3.14159265358979323846264338327950288419716939937510,
  *  </pre>
  *  @see #mpoutc
  */
  static int mpouts (/*const*/ MP a, int la, char[] cs, int lmpnw)
  {
    int l;
    int ll = (int)(la / log10 (MPBDX) + 2.0);
    lmpnw = Math.min (lmpnw, ll);
    if(mpipl<10000)
      l = mpoutc (a, cs, lmpnw);
    else
      l = mpoutx(a,cs, lmpnw);
    l = Math.min (l, la + 20) + 1;
    //cs[l-1] = ',';
    //cs[l]='\0';
    return l;
  }

 /***************** ALGEBRAIC & TRANSCEDENTAL ROUTINES ***********************/
  
 /**
  *  This computes the MP angle A subtended by the MP pair (X, Y) considered as
  *  a point in the x-y plane.  
  *
  *  This is more useful than an arctan or arcsin
  *  routine, since it places the result correctly in the full circle, i.e.
  *  -Pi < A <= Pi.  PI is the MP value of Pi computed by a previous call to
  *  MPPI.  For extra high levels of precision, use MPANGX.  The last word of
  *  the result is not reliable.  Debug output starts with MPIDB = 5.
  *  <br>
  *  Dimension a(lmpnw+2), pi(lmpnw), x(lmpnw), y(lmpnw).
  *  <br>
  *  The Taylor series for Sin converges much more slowly than that of Arcsin.
  *  Thus this routine does not employ Taylor series, but instead computes
  *  Arccos or Arcsin by solving Cos (a) = x or Sin (a) = y using one of the
  *  following Newton iterations, both of which converge to a:
  *  <pre>
  *    z_{k+1} = z_k - [x - Cos (z_k)] / Sin (z_k)
  *    z_{k+1} = z_k + [y - Sin (z_k)] / Cos (z_k)
  *  </pre>
  *  The first is selected if Abs (x) <= Abs (y); otherwise the second is used.
  *  These iterations are performed with a maximum precision level MPNW that
  *  is dynamically changed, approximately doubling with each iteration.
  *  See the comment about the parameter NIT in MPDIVX.
  *  @see #mppi
  *  @see mp_complex#mpangx
  *  @see #mpidvx
  */
  static void mpang (/*const*/ MP x, /*const*/ MP y, /*const*/ MP pi, MP a, int lmpnw)
  {
    int k;
    int mpnw3 = lmpnw+3;
    MP sk0 = new MP(mpnw3,false), sk1 = new MP(mpnw3,false),
      sk2 = new MP(mpnw3,false), 
      sk3 = new MP(mpnw3,false), sk4 = new MP(mpnw3,false);
    
    int ix = x.sign ? 1 : -1;
    int nx = Math.min (x.nw, lmpnw);
    int iy = y.sign ? 1 : -1;
    int ny = Math.min (y.nw, lmpnw);
    
    //  Check if both X and Y are zero.
    
    if (nx == 0 && ny == 0) 
    {
      throw new ArithmeticException(
        "mpang: Both arguments are zero.");
    }
    
    //  Check if Pi has been precomputed.
    MPDPE t1 = new MPDPE();
    mpmdc (pi, t1);
    if (t1.n != 0 || Math.abs (t1.a - CPI) > MPRX2) 
    {
      throw new ArithmeticException(
        "mpang: PI must be precomputed");
    }
    
    //  Check if one of X or Y is zero.
    if (nx == 0) 
    {
      if (iy > 0) 
        mpmuld (pi, new MPDPE(0.5), a, lmpnw); 
      else 
        mpmuld (pi, new MPDPE(-0.5), a, lmpnw); 
      
      return;
      
    } 
    else if (ny == 0) 
    {
      if (ix > 0) 
        zero(a);
      else 
        mpeq (pi, a, lmpnw);
      
      return;
    }
    
    int nws = lmpnw;
    lmpnw++;
    
    //  Determine the least integer MQ such that 2 ^ MQ >= MPNW.
    
    t1.a = nws;
    int mq = (int)(CL2 * Math.log (t1.a) + 1.0 - MPRXX); // *cast*
    
    //  Normalize x and y so that x^2 + y^2 = 1.
    
    mpmul (x, x, sk0, lmpnw);
    mpmul (y, y, sk1, lmpnw);
    mpadd (sk0, sk1, sk2, lmpnw);
    mpsqrt (sk2, sk3, lmpnw);
    mpdiv (x, sk3, sk1, lmpnw);
    mpdiv (y, sk3, sk2, lmpnw);
    
    //  Compute initial approximation of the angle.
    
    mpmdc (sk1, t1);
    MPDPE t2 = new MPDPE();
    mpmdc (sk2, t2);
    t1.n = Math.max (t1.n, -66);
    t2.n = Math.max (t2.n, -66);
    t1.a = t1.a * Math.pow(2.0, t1.n);
    t2.a = t2.a * Math.pow(2.0, t2.n);
    MPDPE t3 = new MPDPE();
    t3.a = Math.atan2 (t2.a, t1.a);
    mpdmc (t3, a);
    
    //  The smaller of x or y will be used from now on to measure convergence.
    //  This selects the Newton iteration (of the two listed above) that has the
    //  largest denominator.
    
    int kk;
    if (Math.abs (t1.a) <= Math.abs (t2.a)) 
    {
      kk = 1;
      mpeq (sk1, sk0, lmpnw);
    } 
    else {
      kk = 2;
      mpeq (sk2, sk0, lmpnw);
    }
    
    lmpnw = 3;
    int iq = 0;
    
    //  Perform the Newton-Raphson iteration described above with a dynamically
    //  changing precision level MPNW (one greater than powers of two).
    
    for (k = 2;  k<=mq; k++)
    {
      lmpnw = Math.min (2 * lmpnw - 2, nws) + 1;
      boolean cont = true;
      while(cont)
      {
        mpcssn (a, pi, sk1, sk2, lmpnw);
        if (kk == 1)
        {
          mpsub (sk0, sk1, sk3, lmpnw);
          mpdiv (sk3, sk2, sk4, lmpnw);
          mpsub (a, sk4, sk1, lmpnw);
        } 
        else 
        {
          mpsub (sk0, sk2, sk3, lmpnw);
          mpdiv (sk3, sk1, sk4, lmpnw);
          mpadd (a, sk4, sk1, lmpnw);
        }
        mpeq (sk1, a, lmpnw);
        if (k == mq - NIT && iq == 0)
          iq = 1;
        else
          cont = false;
      }
    }
    mproun (a, nws);	
  }
  
  /**
  *  This computes the hyperbolic cosine and sine of the MP number A and
  *  returns the two MP results in X and Y, respectively.  
  *  
  *  AL2 is the MP value of Log (10) computed by a previous call to MPLOG.  
  *  For extra high levels of precision, use MPCSHX.  The last word of the result 
  *  is not reliable. Debug output starts with MPIDB = 5.
  *  <br>
  *  Dimension a(lmpnw), al2(lmpnw), x(lmpnw+2), y(lmpnw+2).
  *  <br>
  *  @see #mplog
  *  @see #mpcshx
  */
  static void mpcssh (/*const*/ MP a, /*const*/ MP al2, MP x, MP y, int lmpnw)
  {
    int mpnw3 = lmpnw+3;
    MP f =  new MP(6, false) ,
      sk0 = new MP(mpnw3,false), sk1 = new MP(mpnw3,false), 
      sk2 = new MP(mpnw3,false), sk3 = new MP(mpnw3,false); 
    
    int nws = lmpnw;
    lmpnw++;
    f.sign=true; f.nw=1; f.exponent=0;
    f.mantissa[0]=1; f.mantissa[1]=0;
    
    mpexp (a, al2, sk0, lmpnw);
    mpdiv (f, sk0, sk1, lmpnw);
    mpadd (sk0, sk1, sk2, lmpnw);
    mpmuld (sk2, new MPDPE(0.5), sk3, lmpnw);
    mpeq (sk3, x, lmpnw);
    mpsub (sk0, sk1, sk2, lmpnw);
    mpmuld (sk2, new MPDPE(0.5), sk3, lmpnw);
    mpeq (sk3, y, lmpnw);  
    
    mproun (x, nws);
    mproun (y, nws);
  }

/**
  *  This computes the cosine and sine of the MP number A and returns the two MP
  *  results in X and Y, respectively. 
  *
  *  PI is the MP value of Pi computed by a previous  to MPPI.  
  *  For extra high levels of precision, use MPCSSX.
  *  The last word of the result is not reliable.  Debug output starts with
  *  MPIDB = 6.
  *  <br>
  *  Dimension a(lmpnw), pi(lmpnw), x(lmpnw+2), y(lmpnw+2).
  *  <br>
  *  This routine uses the conventional Taylor's series for Sin (s):
  *  <pre>
  *  Sin (s) =  s - s^3 / 3! + s^5 / 5! - s^7 / 7! ...
  *  </pre>
  *  where s = t - a * pi / 2 - b * pi / 16 and the integers a and b are chosen
  *  to minimize the absolute value of s.  We can then compute
  *  <pre>
  *  Sin (t) = Sin (s + a * pi / 2 + b * pi / 16)
  *  Cos (t) = Cos (s + a * pi / 2 + b * pi / 16)
  *  </pre>
  *  by applying elementary trig identities for sums.  The sine and cosine of
  *  b * pi / 16 are of the form 1/2 * Sqrt {2 +- Sqrt [2 +- Sqrt(2)]}.
  *  Reducing t in this manner insures that -Pi / 32 < s <= Pi / 32, which
  *  accelerates convergence in the above series.
  *  @see #mppi
  *  @see mp_complex#mpcssx
  */
  static void mpcssn (/*const*/ MP a, /*const*/ MP pi, MP x, MP y, int lmpnw)
  {
   boolean ini = a.nw == 14 && a.mantissa[13]==3145728;

   /* if(ini)
      {
	 a.debug("a");
	 pi.debug("pi");
	 x.debug("x");
	 y.debug("y");

      }*/

    double t2;
    int mpnw3 = lmpnw+3;
    MP f =  new MP(6, false) , 
      sk0 = new MP(mpnw3,false), sk1 = new MP(mpnw3,false), 
      sk2 = new MP(mpnw3,false), sk3 = new MP(mpnw3,false), sk4 = new MP(mpnw3,false),
      sk5 = new MP(mpnw3,false), sk6 = new MP(mpnw3,false);
    
    int na = Math.min (a.nw, lmpnw);
    int l1; 
    if (na == 0) 
    {
      x.sign=true; x.nw=1; x.exponent=0; x.mantissa[0]=1;
      zero(y);
      l1 = 0;
      return;   
    }
    
    //  Check if Pi has been precomputed.
    MPDPE t1 = new MPDPE();
    mpmdc (pi, t1);
    if (t1.n != 0 || Math.abs (t1.a - CPI) > MPRX2) 
    {
      throw new ArithmeticException(
        "mpccsn: pi must be precomputed.");
    }
    
    int nws = lmpnw;
    lmpnw++;
    f.nw=1; f.sign=true; f.exponent=0;
    f.mantissa[0]=1; f.mantissa[1]=0;
    
    //  Reduce to between - Pi and Pi.
    
    mpmuld (pi, new MPDPE(2.0), sk0, lmpnw);
    mpdiv (a, sk0, sk1, lmpnw);
    mpnint (sk1, sk2, lmpnw);
    mpsub (sk1, sk2, sk3, lmpnw);
    /*    if(ini)
       {
	  sk0.debug("sk0");
	  sk1.debug("sk1");
	  sk2.debug("sk2"); 
	  sk3.debug("sk3");
       }*/
    
    //  Determine nearest multiple of Pi / 2, and within a quadrant, the nearest
    //  multiple of Pi / 16.  Through most of the rest of this subroutine, KA and
    //  KB are the integers a and b of the algorithm above.
    
    mpmdc (sk3, t1);
    int ka, kb;
    if (t1.n >= - MPNBT) 
    {
      t1.a *= Math.pow(2.0, t1.n);
      t2 = 4.0 * t1.a;  
      
      ka = (int)(nint (t2)); // *cast* 
      kb = (int)(nint (8.0 * (t2 - ka))); // *cast*
    } 
    else 
    {
      ka = 0;
      kb = 0;
    }
    
    t1.a = (8 * ka + kb) / 32.0;
    t1.n = 0;
    mpdmc (t1, sk1);
    mpsub (sk3, sk1, sk2, lmpnw);
    mpmul (sk0, sk2, sk1, lmpnw);
    
    //  Compute cosine and sine of the reduced argument s using Taylor's series.
    
    if (sk1.nw == 0) 
    {
      zero(sk0);
      l1 = 0;
    }
    else
    {
      mpeq (sk1, sk0, lmpnw);
      mpmul (sk0, sk0, sk2, lmpnw);
      l1 = 0;
      
      // 100
      do 
      {
        l1 = l1 + 1;
        if (l1 == 10000) 
        {
          throw new ArithmeticException(
            "mpcssn: Iteration limit exceeded.");
        }   
        t2 = - (2.0 * l1) * (2.0 * l1 + 1.0);
        mpmul (sk2, sk1, sk3, lmpnw);
        mpdivd (sk3, new MPDPE(t2), sk1, lmpnw);
        mpadd (sk1, sk0, sk3, lmpnw);
        mpeq (sk3, sk0, lmpnw);
      }
      //  Check for convergence of the series.
      while (sk1.nw != 0 && sk1.exponent >= sk0.exponent - lmpnw); 
    }
       
    //  Compute Cos (s) = Sqrt [1 - Sin^2 (s)].
    
    mpeq (sk0, sk1, lmpnw);
    mpmul (sk0, sk0, sk2, lmpnw);
    mpsub (f, sk2, sk3, lmpnw);
    mpsqrt (sk3, sk0, lmpnw);
   
    //  Compute cosine and sine of b * Pi / 16.
    int kc;
    kc = Math.abs (kb);
    f.mantissa[0] = 2;
    if (kc == 0) 
    {
      sk2.sign=true; sk2.nw=1; sk2.exponent=0;
      sk2.mantissa[0]=1;
      zero(sk3);
    } 
    else 
    {
      switch (kc)
      {
      case 1: 
        mpsqrt (f, sk4, lmpnw);
        mpadd (f, sk4, sk5, lmpnw);
        mpsqrt (sk5, sk4, lmpnw);
        break;
      case 2: 
        mpsqrt (f, sk4, lmpnw);
        break;
      case 3: 
        mpsqrt (f, sk4, lmpnw);
        mpsub (f, sk4, sk5, lmpnw);
        mpsqrt (sk5, sk4, lmpnw);  
        break;
      case 4: 
        zero(sk4);
        break;
      }
      mpadd (f, sk4, sk5, lmpnw);
      mpsqrt (sk5, sk3, lmpnw);
      mpmuld (sk3, new MPDPE(0.5), sk2, lmpnw);
      mpsub (f, sk4, sk5, lmpnw);
      mpsqrt (sk5, sk4, lmpnw);
      mpmuld (sk4, new MPDPE(0.5), sk3, lmpnw);
    }
    if (kb < 0) sk3.sign = ! sk3.sign;
        
    //  Apply the trigonometric summation identities to compute cosine and sine of
    //  s + b * Pi / 16.
    
    mpmul (sk0, sk2, sk4, lmpnw);
    mpmul (sk1, sk3, sk5, lmpnw);
    mpsub (sk4, sk5, sk6, lmpnw);
    mpmul (sk1, sk2, sk4, lmpnw);
    mpmul (sk0, sk3, sk5, lmpnw);
    mpadd (sk4, sk5, sk1, lmpnw);
    mpeq (sk6, sk0, lmpnw);
       
       //  This code in effect applies the trigonometric summation identities for
    //  (s + b * Pi / 16) + a * Pi / 2.
    
    switch(ka)
    {
    case 0: 
      mpeq (sk0, x, lmpnw);
      mpeq (sk1, y, lmpnw);
      break;
    case 1: 
      mpeq (sk1, x, lmpnw);
      x.sign = ! x.sign;
      mpeq (sk0, y, lmpnw);
      break;
    case -1: 
      mpeq (sk1, x, lmpnw);
      mpeq (sk0, y, lmpnw);
      y.sign = ! y.sign;
      break; 
    case 2: case -2: 
      mpeq (sk0, x, lmpnw);
      x.sign = ! x.sign;
      mpeq (sk1, y, lmpnw);
      y.sign = ! y.sign;
      break;
    }
    

    mproun (x, nws);
    mproun (y, nws);
    /* if(ini)
	{
	 x.debug("x");
	 y.debug("y");
	}*/

  }
/**
  *  This computes the exponential function of the MP number A and returns the
  *  MP result in B.  
  *
  *  AL2 is the MP value of Log(2) produced by a prior to MPLOG.
  *  For extra high levels of precision, use MPEXPX.  The last
  *  word of the result is not reliable.  Debug output starts with MPIDB = 7.
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw+2), al2(lmpnw).
  *  <br>
  *  This routine uses a modification of the Taylor's series for Exp (t):
  *  <pre>
  *  Exp (t) =  (1 + r + r^2 / 2! + r^3 / 3! + r^4 / 4! ...) ^ q * 2 ^ n
  *  </pre>
  *  where q = 256, r = t' / q, t' = t - n Log(2) and where n is chosen so
  *  that -0.5 Log(2) < t' <= 0.5 Log(2).  Reducing t mod Log(2) and
  *  dividing by 256 insures that -0.001 < r <= 0.001, which accelerates
  *  convergence in the above series.
  *  @see #mpexpx
  *  @see #mplog
  */
  static void mpexp (/*const*/ MP a, /*const*/ MP al2, MP b, int lmpnw)
  {
    int i, l1;
    final int nq = 8;
    int mpnw3 = lmpnw+3;
    MP f =  new MP(6, false) , 
      sk0 = new MP(mpnw3,false), sk1 = new MP(mpnw3,false), 
      sk2 = new MP(mpnw3,false), sk3 = new MP(mpnw3,false);
    
    MPDPE t1 = new MPDPE();
    mpmdc (a, t1);
    t1.a = t1.value();
    
    //  Unless the argument is near Log (2), Log(2) must be precomputed.  This
    //  exception is necessary because MPLOG calls MPEXP to initialize Log (2).
    
    MPDPE t2 = new MPDPE();
    if (Math.abs (t1.a - ALT) > MPRDX) 
    {
      mpmdc (al2, t2);
      if (t2.n != - MPNBT || Math.abs (t2.a * Math.pow(0.50, MPNBT) - ALT) > MPRX2) 
      {
        throw new ArithmeticException(
          "mpexp: LOG (2) must be precomputed.");
      }
    }  

    //  Check for overflows and underflows.
    if (t1.a >= 1e9) 
    {
      if (t1.a > 0.0) 
      {
        throw new ArithmeticException(
          "MPEXP: Argument is too large --> "
          + t1.a + " x 10 ^" + t1.n); 
      } 
      else 
      {
        zero(b);
        l1 = 0;
        return;
      }
    }
    
    int nws = lmpnw;
    lmpnw++;
    f.sign=true; f.nw=1; f.exponent=0;
    f.mantissa[0]=1; f.mantissa[1]=0;
    
    //  Compute the reduced argument A' = A - Log(2) * Nint [A / Log(2)].  Save
    //  NZ = Nint [A / Log(2)] for correcting the exponent of the final result.
    int nz;
    if (Math.abs (t1.a - ALT) > MPRDX) 
    {
      mpdiv (a, al2, sk0, lmpnw);
      mpnint (sk0, sk1, lmpnw);
      mpmdc (sk1, t1);
      nz = (int)(t1.value() + fSign (MPRXX, t1.a)); // *cast*
      mpmul (al2, sk1, sk2, lmpnw);
      mpsub (a, sk2, sk0, lmpnw);
    } 
    else 
    {
      mpeq (a, sk0, lmpnw);
      nz = 0;
    }
    
    double tl = sk0.exponent - lmpnw;
    
    //  Check if the reduced argument is zero.
    
    boolean skip = false;
    if (sk0.nw == 0) 
    {
      sk0.nw=1; sk0.sign=true;
      sk0.exponent=0;
      l1 = 0;
      skip = true;
    }
    
    if (!skip)
    {
      //  Divide the reduced argument by 2 ^ NQ.  
      mpdivd (sk0, new MPDPE(1.0, nq), sk1, lmpnw);
      
      //  Compute Exp using the usual Taylor series.    
      mpeq (f, sk2, lmpnw);
      mpeq (f, sk3, lmpnw);
      l1 = 0;
      t2.n = 0; 
      do 
      {
        l1 = l1 + 1;
        if (l1 == 10000) 
        {
          throw new ArithmeticException(
            "mpexp: Iteration limit exceeded.");
        }
        
        t2.a = l1;
        mpmul (sk2, sk1, sk0, lmpnw);
        mpdivd (sk0, t2, sk2, lmpnw);
        mpadd (sk3, sk2, sk0, lmpnw);
        mpeq (sk0, sk3, lmpnw);
      }
      //  Check for convergence of the series.
      while (sk2.nw != 0. && sk2.exponent >= tl); 
      
      //  Raise to the (2 ^ NQ)-th power.   
      for (i = 1; i<=nq; i++)
      {
        mpmul (sk0, sk0, sk1, lmpnw);
        mpeq (sk1, sk0, lmpnw);
      }   
    }
    //  Multiply by 2 ^ NZ.
    mpmuld (sk0, new MPDPE(1.0, nz), sk1, lmpnw);
    mpeq (sk1, b, lmpnw);  
    mproun (b, nws);
  }

/**
  *  This computes the natural logarithm of the MP number A and returns the MP
  *  result in B.  
  *
  *  AL2 is the MP value of Log(2) produced by a prior  to MPLOG. 
  *  For extra high levels of precision, use MPLOGX.  The last word of
  *  the result is not reliable.  Debug output starts with MPIDB = 6.
  *  <br>
  *  Dimension a(lmpnw), al2(lmpnw), b(lmpnw+2).
  *  <br>
  *  The Taylor series for Log converges much more slowly than that of Exp.
  *  Thus this routine does not employ Taylor series, but instead computes
  *  logarithms by solving Exp (b) = a using the following Newton iteration,
  *  which converges to b:
  *  <pre>
  *    x_{k+1} = x_k + [a - Exp (x_k)] / Exp (x_k)
  *  </pre>
  *  These iterations are performed with a maximum precision level MPNW that
  *  is dynamically changed, approximately doubling with each iteration.
  *  See the comment about the parameter NIT in MPDIVX.
  *  @see #mpdivx
  *  @see #mplogx
  */
  static void mplog (/*const*/ MP a, /*const*/ MP al2, MP b, int lmpnw)
  {
    int k;
    int mpnw3 = lmpnw+3;
    MP sk0 = new MP(mpnw3,false), sk1 = new MP(mpnw3,false), 
      sk2 = new MP(mpnw3,false);
    
    int na = Math.min (a.nw, lmpnw);
    
    if (a.sign == false || na == 0) 
    {
      throw new ArithmeticException(
        "mplog: Argument is less than or equal to zero -->"
        + a);
    }
    
    //  Unless the input is close to 2, Log (2) must have been precomputed.
    MPDPE t1 = new MPDPE(), t2 = new MPDPE();
    mpmdc (a, t1);
    if (Math.abs (t1.a - 2.0) > 1e-3 || t1.n != 0) 
    {
      mpmdc (al2, t2);
      if (t2.n != - MPNBT || Math.abs (t2.a * Math.pow(0.50, MPNBT) - ALT) > MPRX2) 
      {
        throw new ArithmeticException(
          "mplog: LOG (2) must be precomputed.");
      }
    }
    
    //  Check if input is exactly one.
    if (a.nw == 1 && a.sign && a.exponent == 0 && a.mantissa[0] == 1.) 
    {
      b.nw = 0; b.exponent=0; b.sign=true;
      return;
    }
    int nws = lmpnw;
    
    //  Determine the least integer MQ such that 2 ^ MQ >= MPNW.
    t2.a = nws;
    int mq = (int)(CL2 * Math.log (t2.a) + 1.0 - MPRXX); // *cast*
    
    
    //  Compute initial approximation of Log (A).
    t1.a = Math.log (t1.a) + t1.n * ALT;
    t1.n = 0;
    mpdmc (t1, b);
    lmpnw = 3;
    int iq = 0;
    
    //  Perform the Newton-Raphson iteration described above with a dynamically
    //  changing precision level MPNW (one greater than powers of two).
    
    for(k = 2; k<=mq; k++)
    {
      lmpnw = Math.min (2 * lmpnw - 2, nws) + 1 ;
      boolean cont = true;
      while(cont)
      { 
        mpexp (b, al2, sk0, lmpnw);
        mpsub (a, sk0, sk1, lmpnw);
        mpdiv (sk1, sk0, sk2, lmpnw);
        
        mpadd (b, sk2, sk1, lmpnw);
        
        mpeq (sk1, b, lmpnw);
        if (k == mq - NIT && iq == 0) 
          iq = 1;
        else
          cont = false;
      }
    }
    mproun (b, nws);
  }

 /**
  *  This computes Pi to available precision (MPNW mantissa words).  
  *
  *  For extra high levels of precision, use MPPIX.  The last word of 
  *  the result is not reliable.  Debug output starts with MPIDB = 7.
  *  <br>
  *  Dimension pi(lmpnw+2).
  *  <br>
  *  The algorithm that is used for computing Pi, which is due to Salamin
  *  and Brent, is as follows:
  *  <pre>
  *  Set  A_0 = 1,  B_0 = 1/Sqrt(2)  and  D_0 = Sqrt(2) - 1/2.
  *  </pre>
  *  then from k = 1 iterate the following operations:
  *  <pre>
  *  A_k = 0.5 * (A_{k-1} + B_{k-1})
  *  B_k = Sqrt (A_{k-1} * B_{k-1})
  *  D_k = D_{k-1} - 2^k * (A_k - B_k) ^ 2
  *  </pre>
  *  Then  P_k = (A_k + B_k) ^ 2 / D_k  converges quadratically to Pi.
  *  In other words, each iteration approximately doubles the number of correct
  *  digits, providing all iterations are done with the maximum precision.
  *  @see #mppix
  */
  static void mppi(MP pi, int lmpnw)
  {
    int k;
    double t1;
    int mpnw3 = lmpnw+3;
    MP f =  new MP(6, false) , 
      sk0 = new MP(mpnw3,false), sk1 = new MP(mpnw3,false), 
      sk2 = new MP(mpnw3,false), sk3 = new MP(mpnw3,false), sk4 = new MP(mpnw3,false);
    
    //  Perform calculations to one extra word accuracy.
    int nws = lmpnw;
    lmpnw++;
    
    //  Determine the number of iterations required for the given precision level.
    //  This formula is good only for this Pi algorithm.
    
    t1 = nws * log10 (MPBDX);
    int mq = (int)(CL2 * (Math.log (t1) - 1.0) + 1.0); // *cast*
    
    //  Initialize as above.
    
    sk0.nw=1; sk0.sign=true; sk0.exponent=0; sk0.mantissa[0]=1;
    f.nw=1; f.sign=true; f.exponent=0; f.mantissa[0]=2; f.mantissa[1]=0;
    
    mpsqrt (f, sk2, lmpnw);
    mpmuld (sk2, new MPDPE(0.50), sk1, lmpnw);
    f.exponent = -1;
    f.mantissa[0] = (float)(0.50 * MPBDX); // *cast*
    mpsub (sk2, f, sk4, lmpnw);
    
    //  Perform iterations as described above.
    
    for (k = 1; k<= mq; k++)
    {
      mpadd (sk0, sk1, sk2, lmpnw);
      mpmul (sk0, sk1, sk3, lmpnw);
      mpsqrt (sk3, sk1, lmpnw);
      mpmuld (sk2, new MPDPE(0.50), sk0, lmpnw);
      mpsub (sk0, sk1, sk2, lmpnw);
      mpmul (sk2, sk2, sk3, lmpnw);
      t1 = Math.pow(2.0, k);
      mpmuld (sk3, new MPDPE(t1), sk2, lmpnw);
      mpsub (sk4, sk2, sk3, lmpnw);
      mpeq (sk3, sk4, lmpnw);
    }
    
    //  Complete computation.
    mpadd (sk0, sk1, sk2, lmpnw);
    mpmul (sk2, sk2, sk3, lmpnw);
    mpdiv (sk3, sk4, sk2, lmpnw);
    mpeq (sk2, pi, lmpnw);
    
    mproun (pi, nws);
  }


//************************ FFT CORE ROUTINES  *******************************
  
  /**
  *  This performs an N-point complex-to-real FFT, where N = 2^M.  X is the
  *  double complex input array, and Y is the double output array.
  * 
  *  The array X is used as a scratch array in MPFFT1, and so is overwritten.
  *  X must be dimensioned with N/2+N1*NSP1+1 DC cells, and Y with N DP cells, 
  *  where N = 2^M and N1 = 2^int(M/2).   This Dimension requirement for X is 
  *  somewhat greater than shown in the Dimension statement below, because 
  *  MPFFT1, which is called by this routine, requires more.  IS is the SIGN of
  *  the transform.  Before calling MPFFTCR, the UU1 and UU2 arrays must be 
  *  initialized by calling MPINIX.  This routine is not intended to be called 
  *  directly by the user.
  *  @see #init_mpuu1
  *  @see #mpfft1
  */
  static void mpfftcr (int is, int m, int n, DComplex x[], double y[])
  {
    int k;
    final DComplex pointFive = new DComplex(0.5);
    final DComplex zeroOne = new DComplex(0.0, 1.0);
  
    int mx = (int)mpuu1[0].real(); // *cast*
    
    //  Check if input parameters are invalid.
    
    if ((is != 1 && is != -1) || m < 3 || m > mx) 
    {
      throw new ArithmeticException(
        "mpfftcr: Either the UU arrays have not been initialized " +
        "or one of the input parameters is invalid: " +
        is + "\t" + m + "\t" + mx);
    }
    
    DComplex dc1[] = new DComplex[n/2], 
      a1, a2, x1, x2;
    
    int n1 = (int)(Math.pow(2, (m / 2))); // *cast*
    int n2 = n / 2;
    int n4 = n / 4;
    
    //  Construct the input to MPFFT1.
    dc1[0] = pointFive.multiply(
      new DComplex ((x[0].add(x[n2])).real(),  
                     (x[0].subtract(x[n2])).real()
                    ));
    
    if (is == 1) 
      dc1[n4] = x[n4].conjg();
    else
   //   try
   //     {
          dc1[n4] = (DComplex)x[n4].clone();
    //    }
    //  catch(CloneNotSupportedException e){}

    int ku = n2;
    
    if (is == 1) 
    {
      //dir$ ivdep
      for (k = 1; k< n4; k++)
      {
        x1 = x[k];
        x2 = x[n2-k].conjg();
        a1 = x1.add(x2);
        a2 = zeroOne.multiply(
          mpuu1[k+ku]).multiply(x1.subtract(x2));
        dc1[k] = pointFive.multiply(a1.add(a2));
        dc1[n2-k] = pointFive.multiply((a1.subtract(a2)).conjg());
      }
    } 
    else 
    {
      //dir$ ivdep
      for (k = 1; k< n4; k++)
      {
        x1 = x[k];
        x2 = (x[n2-k]).conjg();
        a1 = x1.add(x2);
        a2 = zeroOne.multiply(
          mpuu1[k+ku].conjg()).multiply(x1.subtract(x2));

        dc1[k] = pointFive.multiply(a1.add(a2));
        dc1[n2-k] = pointFive.multiply((a1.subtract(a2)).conjg());
      }
    }
    
    //  Perform a normal N/2-point FFT on DC1.
    
    mpfft1 (is, m - 1, n1, n2 / n1, dc1, x);
    
    //  Copy DC1 to Y such that DC1(k) = Y(2k-1) + i Y(2k).
    
    for(k = 0; k< n / 2; k++)
    {
      y[2*k] = dc1[k].real();
      y[2*k+1] = dc1[k].aimag();
    }
  }

  /**
  *  This performs an N-point real-to-complex FFT, where N = 2^M.  X is the
  *  double input array, and Y is the double complex output array.
  *
  *  X must be dimensioned with N DP cells, and Y with N/2+N1*NSP1+1 DC cells,
  *  where N = 2^M and N1 = 2^int(M/2).  This Dimension requirement for Y is 
  *  somewhat greater than that shown in the Dimension statement below, because
  *  MPFFT1, which is called by this routine, requires more.  IS is the SIGN of
  *  the transform.  Before calling MPFFTRC, the UU1 and UU2 arrays must be 
  *  initialized by calling MPINIX.  This routine is not intended to be called 
  *  directly by the user.
  *  @see #init_mpuu1
  */
  static void mpfftrc (int is, int m, int n, /*const*/ double x[], DComplex y[])
  {
    int k;
    int mx = (int)mpuu1[0].real(); // *cast*
    
    //  Check if input parameters are invalid.
    
    if ((is != 1 && is != -1) || m < 3 || m > mx) 
    {
      throw new ArithmeticException(
        "mpfftrc: Either the UU arrays have not been initialized " +
        "or one of the input parameters is invalid: " +
        is + "\t" + m + "\t" + mx);
    }
    
    DComplex dc1[] = new DComplex[n/2], a1, a2, z1, z2;
    
    int n1 = (int) Math.pow(2,(m / 2)); // *cast*
    int n2 = n / 2;
    int n4 = n / 4;
    
    //  Copy X to DC1 such that DC1(k) = X(2k-1) + i X(2k).
    
    for(k = 0; k< n2; k++)
      dc1[k] = new DComplex (x[2*k], x[2*k+1]);
    
    //  Perform a normal N/2-point FFT on DC1.
    
    mpfft1 (is, m - 1, n1, n2 / n1, dc1, y);
    
    //  Reconstruct the FFT of X.
    
    y[0] = new DComplex (2.0 * (dc1[0].real() + dc1[0].aimag()), 0.0);
    
    final DComplex two = new DComplex(2.0);
    if (is == 1) 
      y[n4] = dc1[n4].multiply(two);
    else 
      y[n4] = dc1[n4].conjg().multiply(two);
    
    y[n2] = new DComplex (2.0 * (dc1[0].real() - dc1[0].aimag()), 0.0);
    int ku = n2;

    final DComplex zeroMinOne = new DComplex(0.0, -1.0);
    if (is == 1) 
    {
      //!dir$ ivdep
      for(k = 1; k< n4; k++)
      {
        z1 = dc1[k];
        z2 = dc1[n2-k].conjg();
        a1 = z1.add(z2);
        a2 = zeroMinOne.multiply(mpuu1[k+ku]).multiply(z1.subtract(z2));
        y[k] = a1.add(a2);
        y[n2-k] = a1.subtract(a2).conjg();
      }
    } 
    else 
    {
      //!dir$ ivdep
      for (k = 1; k< n4; k++)
      {
        z1 = dc1[k];
        z2 = dc1[n2-k].conjg();
        a1 = z1.add(z2);
        a2 = zeroMinOne.multiply(mpuu1[k+ku].conjg()).multiply(z1.subtract(z2));
        y[k] = a1.add(a2);
        y[n2-k] = a1.subtract(a2).conjg();
      }
    }
  }


/**
  *  This routine performs a complex-to-complex FFT.  
  *  
  *  IS is the SIGN of the transform, 
  *  N = 2^M is the size of the transform.  N1 = 2^M1 and N2 = 2^M2,
  *  where M1 and M2 are defined as below.  X is the input and output array,
  *  and Y is a scratch array.  X must have at N, and Y at least N + N1*MPNSP1,
  *  double complex cells.  The arrays MPUU1 and MPUU2 must have been 
  *  initialized by calling MPINIX.  This routine is not intended to be called 
  *  directly by the user.
  *  <br>
  *  This employs the two-pass variant of the "four-step" FFT.  See the
  *  article by David H. Bailey in J. of Supercomputing, March 1990, p. 23-35.
  *  @see #init_mpuu1
  */
  static void mpfft1 (int is, int m, int n1, int n2, DComplex x[], DComplex y[])
  {
    int i, j, k;
    DComplex q1[] = new DComplex[MPNSP2], q2[] = new DComplex[MPNSP2], 
      z1[][] = new DComplex[MPNROW+MPNSP1][n1], 
      z2[][] = new DComplex[MPNROW+MPNSP1][n1];
    
    int yrow = n2 + MPNSP1;
    
    int m1 = (m + 1) / 2;
    int m2 = m - m1;
    int nr1 = Math.min (n1, MPNROW);
    int nr2 = Math.min (n2, MPNROW);
    int ku = (int)mpuu2[m-1].real(); // *cast*
    
    for (i = 0; i<n1; i += nr1)
    {
      //  Copy NR1 rows of X (treated as a N1 x N2 complex array) into Z1.
      for(k = 0; k<nr1; k++)
      {
        for(j = 0; j<n2; j++)
          z1[k][j] = x[j*n1 + i+k];
                 
      }
      
      //  Perform NR1 FFTs, each of length N2.
      
      mpfft2 (is, nr1, m2, n2, z1, z2);
      
      //  Multiply the resulting NR1 x N2 complex block by roots of unity and
      //  store transposed into the appropriate section of Y.
      
      int iu = i + ku - n1 - 1;
      
      if (is == 1) 
      {
        for(k = 0; k<nr1; k++)
        {
          for (j = 0; j< n2; j++)
            y[(i+k)*yrow + j] = mpuu2[iu+k+ (j+1)*n1].multiply(z1[k][j]);
        } 
      } 
      else 
      {
        for(k = 0; k<nr1; k++)
        {
          for (j = 0; j<n2; j++)
            y[(i+k)*yrow + j] = mpuu2[iu+k+ (j+1)*n1].conjg().multiply(z1[k][j]);
        }
      }
    }
    
    for (i = 0; i< n2; i += nr2)
    {
      
      //  Copy NR2 rows of the Y array into Z2.
      
      for(k = 0; k< nr2; k++)
      {
        for(j = 0; j< n1; j++)
          z2[k][j] = y[j*yrow + i+k];
      }
      
      //  Perform NR2 FFTs, each of length N1.
      
      mpfft2 (is, nr2, m1, n1, z2, z1);
      
      //  Copy NR2 x N1 complex block back into X array.  It's a little more
      //  complicated if M is odd.
      
      if ((m%2) == 0) 
      {
        for(k = 0; k< nr2; k++)
        {
          // don't need to clone since z2 is temporary
          for (j = 0; j< n1; j++)
            x[i+k + j*n1] = z2[k][j];
        }
      } 
      else 
      {
        int j2;
        for (j = 0; j<n1 / 2; j++)
        {
          j2 = 2 * j;
          //dir$ ivdep
          for (k = 0; k< nr2; k++)
          {
            // don't need to clone sin z2 is temporary
            x[i+k + j*n1] = z2[k][j2]; 
            x[i+k+n2 + n1*j] = z2[k][j2+1]; 
          }   
        }
      }
    }    
  }

  /**
  *  This performs NS simultaneous N-point complex-to-complex FFTs, where
  *  N = 2^M.  
  *
  *  X is the input and output array, UU1 is the root of unity array,
  *  and Y is a scratch array.  X, Y and UU1 are double complex.   This routine
  *  is not intended to be called directly by the user.
  */
  static void mpfft2 (int is, int ns, int m, int n, DComplex x[][], DComplex y[][])
  {
    int l, j, i;
    
    //  Perform the second variant of the Stockham FFT.
    
    for(l = 1; l<=m; l+=2)
    {
      mpfft3 (is, l, ns, m, n, x, y);
      if (l == m)
      {
        for(i = 0; i<ns; i++) 
        {
          for (j = 0; j<n; j++)
            x[i][j] = y[i][j];
        }
        return;
      }
      mpfft3 (is, l + 1, ns, m, n, y, x);
    }  
  }

  /**
  *  This performs the L-th iteration of the second variant of the Stockham FFT
  *  on the NS vectors in X.  
  *
  *  Y is a scratch array, and UU1 is the root of unity array.  
  *  X, Y and UU1 are double complex.  This routine is not
  *  intended to be called directly by the user.
  */
  static void mpfft3 (int is, int l, int ns, int m, int n, DComplex x[][], DComplex y[][])
  {  
    DComplex u1, x1, x2;
    int i,j, k;
    
    //  Set initial parameters.
    int n1 = n / 2;
    int lk = (int)Math.pow(2,(l - 1)); // *cast*
    int li = (int)Math.pow(2, (m - l)); // *cast*
    int lj = 2 * lk;
    int ku = li;
    
    int i11, i12, i21, i22;
    for(i = 0; i<= li - 1; i++)
    {
      i11 = i * lk + 1;
      i12 = i11 + n1;
      i21 = i * lj + 1;
      i22 = i21 + lk;
      if (is == 1) 
        u1 = mpuu1[i+ku];
      else 
        u1 = mpuu1[i+ku].conjg();
      
      for (k = -1; k< lk - 1; k++)
      {
        //!dir$ ivdep
        for (j = 0; j< ns; j++)
        {
          x1 = x[j][i11+k];
          x2 = x[j][i12+k];
          y[j][i21+k] = x1.add(x2);
          y[j][i22+k] = u1.multiply(x1.subtract(x2));
        }
      }
    }
  }
  
  /**
  * This computes the linear convolution of the N-long inputs A and B.  
  *
  * |IQ| is the number of arguments (i.e., if IQ = 1, then B is ignored).  
  * If IQ is negative (and N < 64) then only the second half of the result vector is
  * required (i.e. this is a  by itself -- see below). NSQ = int(sqrt(3*N))
  * is an input required for the Dimension of DC1 and DC2 (see below). 
  * <br>
  * This routine employs an advanced FFT-based scheme, except for small n.
  * This routine is not intended to be called directly by the user.
  * <br> 
  *  Two machine-dependent parameters are set in this routine:
  *  <br>ERM = Maximum tolerated FFT roundoff error.  On IEEE systems ERM =
  *    0.438D0.  It is not necessary to specify ERM for modest levels of
  *    precision -- see comments below.
  *  <br>MBT = Number of mantissa bits in double data.  MBT = 53 on
  *    IEEE systems, and MBT = 48 (i.e. single precision) on Crays.
  *    It is not necessary to specify MBT for modest levels of precision.
  */
  static void mplconv (int iq, int n, int nsq, double a[], double b[], double c[])
  {
    int i,j,k;
    final double erm = 0.438;
    final int mbt = 53;
    double an, t1, t2;
    
    //  Handle the case where N is less than NCR1 = 2 ** (mpmcr-1).  If IQ < 0, 
    //  only the second half of the result vector is returned, since the first 
    //  half won't be used.
    
    int ncr1 = (int)Math.pow(2, (mpmcr - 1)); // *cast*
    int n1, n2;
    if (n < ncr1) 
    {
      switch(iq)
      {
      case 1:
        for(k = 0; k<2 * n; k++)
        {
          t1 = 0.0;
          n1 = Math.max (k - n + 2, 1);
          n2 = Math.min (k+1, n);
          
          for(j = n1-1; j<n2; j++)
            t1 += a[j] * a[k-j];
          
          c[k] = t1;
        }
        break;
        
      case 2: 
        for(k = 0; k< 2 * n; k++)
        {
          t1 = 0.0;
          n1 = Math.max (k - n + 2, 1);
          n2 = Math.min (k+1, n);
          
          for(j = n1-1; j< n2; j++)
            t1 += a[j] * b[k-j];
          
          c[k] = t1;
        }
        break;
        
      case -1:
        for (k = 0; k< n - 1; k++)
          c[k] = 0.0;
        
        for(k = n-1; k< 2 * n; k++)
        {
          t1 = 0.0;
          n1 = k - n + 2;
          n2 = n;
          
          for(j = n1-1; j< n2; j++)
            t1 += a[j] * a[k-j];
          
          c[k] = t1;
        }
        break;
        
      case -2:
        for (k = 0; k< n - 1; k++)
          c[k] = 0.0;
        
        for(k = n-1; k<2 * n; k++)
        {
          t1 = 0.0;
          n1 = k - n + 2;
          n2 = n;
          
          for (j = n1 - 1; j< n2; j++)
            t1 += a[j] * b[k-j];
          
          c[k] = t1;
        }
        break;
      }
      return;
    }
    
    double d1[] = new double[3*n+2], d2[] = new double[3*n+2], 
      d3[] = new double[3*n+2];
    DComplex dc1[] = new DComplex[3*n/2+nsq*MPNSP1+3], 
      dc2[] = new DComplex[3*n/2+nsq*MPNSP1+3];
    
    //  Determine M1 and N1.  Note that by this reckoning, N1 <= 1.5 N.  This is
    //  the reason for the 3*n/2 dimensions above.
    
    t1 = 0.750 * n;
    int m1 = (int)(CL2 * Math.log (t1) + 1.0 - MPRXX); // *cast*
    n1 = (int) (Math.pow(2, m1)); // *cast*
    int m2 = m1 + 1;
    n2 = 2 * n1;
    int n4 = 2 * n2;
    int nm = Math.min (2 * n, n2);
    
    if (Math.abs (iq) == 1) 
    {
      for(i = 0; i<n; i++)
        d1[i] = a[i];
      
      for(i = n ; i< n2; i++)
        d1[i] = 0.0;
      
      //  Perform a forward real-to-complex FFT on the vector in a.
      
      mpfftrc (1, m2, n2, d1, dc1);
      
      //  Square the resulting complex vector.
      
      for(i = 0; i< n1 + 1; i++)
        dc1[i]  = dc1[i].multiply(dc1[i]);
    } 
    else 
    {
      for(i = 0; i<n; i++)
      {
        d1[i] = a[i];
        d2[i] = b[i];
      }
      
      for(i = n; i< n2; i++)
      {
        d1[i] = 0.0;
        d2[i] = 0.0;
      }
      
      //  Perform forward real-to-complex FFTs on the vectors in a and b.
      
      mpfftrc (1, m2, n2, d1, dc1);
      mpfftrc (1, m2, n2, d2, dc2);
      
      //  Multiply the resulting complex vectors.
      
      for(i = 0; i<= n1; i++)
        dc1[i] = dc1[i].multiply(dc2[i]);
    }
    
    //  Perform an inverse complex-to-real FFT on the resulting data.
    
    mpfftcr (-1, m2, n2, dc1, d3);
    
    //  Divide by N4.
    // here!!!!!!!!!
    an = 1.0 / n4;
    
    for(i = 0; i< nm; i++)
    {
      t1 = an * d3[i];
      t2 = nint(t1);
      // D1(I) = ABS (T2 - T1)
      c[i] = t2;
    }
    
    //  Find the largest FFT roundoff error.  Roundoff error is minimal unless
    //  exceedingly high precision (i.e. over one million digits) is used.  Thus
    //  this test may be disabled in normal use.  To disable this test, set the
    //  skip flag to false and comment out the line of the previous loop
    //  that begins D1(I) =.  To enable, for  the reverse.
    
    //  This code can be used as a rigorous system integrity test.  First set
    //  MBT according to the system being used, and { set ERM to be fairly
    //  small, say 0.001 or whatever is somewhat larger than the largest FFT
    //  roundoff error typically encountered for a given precision level on the
    //  computer being used.  Enable this test as explained in the previous
    //  paragraph.  Then if an anomalously large roundoff error is detected, a
    //  hardware or compiler error has likely occurred.
    
    boolean skip = true;
    if(!skip)
    {
      int i1=0;
      t1 = 0.0;
      
      for(i = 0; i< nm; i++)
      {
        if (d1[i] > t1) 
        {
          i1 = i;
          t1 = d1[i];
        }
      }
      //  Check if maximum roundoff error exceeds the limit ERM, which is set above.
      //  Also determine the number of fractional bits and how large the error is in
      //  terms of units in the last place (ulp).
      
      if (t1 > erm)  
      {
        int i2, i3, i4, i5;
        t2 = an * d1[i1];
        i2 = (int)(CL2 * Math.log (t1) + 1.0 + MPRXX); // *cast*
        i3 = (int)(CL2 * Math.log (t2) + 1.0 + MPRXX); // *cast*
        i4 = mbt + i2 - i3;
        i5 = (int)(t1 * Math.pow(2, i4) + MPRXX); // *cast
        throw new ArithmeticException(
          "mplconv: Excessive FFT roundoff error --> " 
          + "\t" + i1 + "\t" + t1 + "\t" + i4 + "\t" + i5);        
      }
      
    }
    //  Handle case where n > n1.
    int m, m21, ms;
    if (n > n1) 
    {
      m = n - n1;
      m2 = 2 * m;
      m21 = 2 * m - 1;
      ms = (int)(Math.sqrt (3.0 * m21) + MPRXX); // *cast*
      k = n1 - m + 1;
      
      if (Math.abs (iq) == 1) 
      {
        for(i = 0; i< m21; i++)
          d1[i] = a[k+i];
        
        mplconv (-1, m21, ms, d1, d2, d3);
      } 
      else 
      {
        for (i = 0; i<m21; i++)
        {
          d1[i] = a[k+i];
          d2[i] = b[k+i];
        }
        
        mplconv (-2, m21, ms, d1, d2, d3);
      }
      
      int ii;
      for(i = 0; i< m2; i++)
      {
        ii = i + m2 - 2;
        c[i] -= d3[ii];
        c[i+n2] = d3[ii];
      }
    }
  }
  

  //************************ADVANCE ARITHMETIC ROUTINES ************************

  /**
  *  This computes the cube root of the MP number A and returns the MP result
  *  in B.  
  *  
  *  Before calling MPCBRX, the arrays UU1 and UU2 must be initialized by
  *  calling MPINIX.  For modest levels of precision, use MPCBRT.  Debug output
  *  starts with MPIDB = 6.
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw+2).
  *  <br>
  *  This routine uses basically the same Newton iteration algorithm as MPCBRT.
  *  In fact, this routine calls MPCBRT to obtain an initial approximation.
  *  See the comment about the parameter NIT in MPDIVX.
  *  @see #mpcbrt
  *  @see #mpdivx
  *  @see #mpinix
  */
  static void mpcbrx (/*const*/ MP a, MP b, int lmpnw)
  {
    int k;
    double t1;
    int mpnw3 = lmpnw+3;
    MP f =  new MP(6, false) , sk0 = new MP(mpnw3,false), 
      sk1 = new MP(mpnw3,false), sk2 = new MP(mpnw3,false);
    
    int na = Math.min (a.nw, lmpnw);
    int ncr = (int)(Math.pow(2, mpmcr)); // *cast*
    
    if (na == 0) 
    {
      zero(b);
      return;
    }
    if (a.sign == false) 
    {
      throw new ArithmeticException(
        "mpcbrx: Argument is negative --> " + a);
    }
    
    //  Check if precision level is too low to justify the advanced routine.
    
    if (lmpnw <= ncr) 
    {
      mpcbrt (a, b, lmpnw);
      return;
    }
    int nws = lmpnw;
    
    //  Determine the least integer MQ such that 2 ^ MQ >= MPNW.
    
    t1 = lmpnw;
    int mq = (int)(CL2 * Math.log (t1) + 1.0 - MPRXX); // *Cast*
    
    //  Compute A^2 outside of the iteration loop.
    
    mpsqx (a, sk0, lmpnw);
    
    //  Compute the initial approximation of A ^ (-2/3).
    
    lmpnw = ncr + 1;
    mpcbrt (a, sk1, lmpnw);
    mpdiv (sk1, a, b, lmpnw);
    f.nw=1; f.sign=true; f.exponent=0; f.mantissa[0]=1; f.mantissa[1]=0;
    int iq = 0, nw1, nw2;
    
    //  Perform the Newton-Raphson iteration described above with a dynamically
    //  changing precision level MPNW (powers of two).
    
    for(k = mpmcr + 1; k<=mq - 1; k++)
    {
      nw1 = lmpnw; 
      lmpnw=Math.min (2 * lmpnw - 2, nws) + 1;
      nw2 = lmpnw; 
      
      boolean cont=true;
      while (cont)
      {
        mpsqx (b, sk1, lmpnw);
        mpmulx (b, sk1, sk2, lmpnw);
        mpmulx (sk0, sk2, sk1, lmpnw);
        mpsub (f, sk1, sk2, lmpnw);
        
        lmpnw = nw1;
        
        mpmulx (b, sk2, sk1, lmpnw);
        mpdivd (sk1, new MPDPE(3.0), sk2, lmpnw);
        
        lmpnw = nw2;
        
        mpadd (b, sk2, sk1, lmpnw);
        mpeq (sk1, b, lmpnw);
        
        if (k == mq - NIT && iq == 0) 
          iq = 1;
        
        else
          cont=false;
      }
    }
    
    //  Perform last iteration using Karp's trick.
    
    mpmulx(a, b, sk0, lmpnw);
    nw1 = lmpnw; 
    lmpnw = Math.min (2 * lmpnw - 2, nws) + 1;
    
    nw2 = lmpnw;
    
    mpsqx (sk0, sk1, lmpnw);
    mpmulx (sk0, sk1, sk2, lmpnw);
    mpsub (a, sk2, sk1, lmpnw);
    
    lmpnw = nw1;
    
    mpmulx (sk1, b, sk2, lmpnw);
    mpdivd (sk2, new MPDPE(3.0), sk1, lmpnw);
    
    lmpnw = nw2;
    
    mpadd (sk0, sk1, sk2, lmpnw);
    mpeq (sk2, b, lmpnw);  
    
    mproun (b, nws);
  }

  /**
  *  This divides the MP number A by the MP number B and returns the MP result
  *  in c.  
  *
  *  Before calling MPDIVX, the arrays UU1 and UU2 must be initialized by
  *  calling MPINIX.  For modest levels of precision, use MPDIV.  Debug output 
  *  starts with MPIDB = 7.
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw), c(lmpnw+2). 
  *  <br>
  *  This subroutine employs the following Newton-Raphson iteration, which
  *  converges to 1 / B:
  *  <pre>
  *   X_{k+1} = X_k + (1 - X_k * B) * X_k
  *  </pre>
  *  where the muliplication () * X_k is performed with only half of the
  *  normal level of precision.  These iterations are performed with a
  *  maximum precision level MPNW that is dynamically changed, doubling with
  *  each iteration.  The final iteration is performed as follows (this is
  *  due to A. Karp):
  *  <pre>
  *   A / B = (A * X_n) + [A - (A * X_n) * B] * X_n  (approx.)
  *  </pre>
  *  where the multiplications A * X_n and [] * X_n are performed with only
  *  half of the final level of precision.
  *  <br>
  *  One difficulty with this procedure is that errors often accumulate in the
  *  trailing mantissa words.  This error can be controlled by repeating one of
  *  the iterations.  The iteration that is repeated is controlled by setting
  *  the parameter NIT below:  If NIT = 0, the last iteration is repeated (this
  *  is most effective but most expensive).  If NIT = 1, { the next-to-last
  *  iteration is repeated, etc.
  *  @see #mpdiv
  *  @see #mpinix
  */
  static void mpdivx (/*const*/ MP a, /*const*/ MP b, MP c, int lmpnw)
  {
    int k;
    double t1;
    
    int mpnw3 = lmpnw+3;
    MP f =  new MP(6, false) , sk0 = new MP(mpnw3,false), 
      sk1 = new MP(mpnw3,false), sk2 = new MP(mpnw3,false);
    
    int na = Math.min (a.nw, lmpnw);
    int nb = Math.min (b.nw, lmpnw);
    int ncr = (int)(Math.pow(2, mpmcr)); // *cast*
    
    //  Check if dividend is zero.
    /* Sweny
    if (na == 0) 
    {
      zero(c);
      return;
    }
    */
    
    //  Check if divisor is zero. 
    if (nb == 0)  
    {
      throw new ArithmeticException(
        "mpdivx: Divisor is zero");
    }
    
    //  Check if precision level of divisor is too low to justify the advanced
    //  routine.
    if (nb <= ncr) 
    {
      mpdiv (a, b, c, lmpnw);
      return;
    }
    int nws = lmpnw;
    
    //  Determine the least integer MQ such that 2 ^ MQ >= MPNW.
    
    t1 = lmpnw;
    int mq = (int)(CL2 * Math.log (t1) + 1.0 - MPRXX); // *cast*
    
    //  Compute the initial approximation of 1 / B to a precision of NCR words.
    lmpnw = ncr + 1;
    
    f.nw=1; f.sign=true; f.exponent=0; f.mantissa[0]=1; f.mantissa[1]=0;
    mpdiv(f,b,c, lmpnw);
    int iq = 0;
    int nw1, nw2;
    
    //  Perform the Newton-Raphson iterations described above.
    for(k = mpmcr + 1; k<= mq - 1; k++)
    {
      nw1 = lmpnw;
      lmpnw = Math.min (2 * lmpnw - 2, nws) + 1;
      nw2 = lmpnw;
      boolean cont = true;
      while (cont)
      {
        mpmulx (b, c, sk0, lmpnw);
        mpsub (f, sk0, sk1, lmpnw);
        lmpnw = nw1;
        mpmulx (c, sk1, sk0, lmpnw);
        lmpnw = nw2;
        mpadd (c, sk0, sk1, lmpnw);
        mpeq (sk1, c, lmpnw);
        if (k == mq - NIT && iq == 0) 
          iq = 1;
        else
          cont=false;
      }
    }
    
    //  Perform last iteration using Karp's trick.
    mpmulx (a, c, sk0, lmpnw);
    nw1 = lmpnw;
    lmpnw = Math.min (2 * lmpnw - 2, nws) + 1;
    nw2 = lmpnw;
    mpmulx (sk0, b, sk1, lmpnw);
    mpsub (a, sk1, sk2, lmpnw);
    lmpnw = nw1;
    mpmulx (sk2, c, sk1, lmpnw);
    lmpnw = nw2;
    mpadd (sk0, sk1, sk2, lmpnw);
    mpeq (sk2, c, lmpnw);
    
    //  Restore original precision level.
    
    mproun (c, nws);  
  }

 /**
  *  This routine multiplies MP numbers A and B to yield the MP product c.
  *
  *  Before calling MPMULX, the arrays UU1 and UU2 must be initialized by
  *  calling MPINIX.  For modest levels of precision, use MPMUL.  Debug output 
  *  starts with MPIDB = 8.
  *  <br>
  *  Dimension a(lmpnw+), b(lmpnw), c(lmpnw+2).
  *  <br>
  *  This routine returns up to MPNW mantissa words of the product.  If the
  *  complete double-long product of A and B is desired (for example in large
  *  integer applications), then MPNW must be at least as large as the sum of
  *  the mantissa lengths of A and B.  In other words, if the precision levels
  *  of A and B are both 256 words, { MPNW must be at least 512 words to
  *  obtain the complete double-long product in c.
  *  <br>
  *  This subroutine uses an advanced technique involving the fast Fourier
  *  transform (FFT).  For high precision it is significantly faster than the
  *  conventional scheme used in MPMUL.
  *  @see #mpmul
  *  @see #mpinix
  */
  static void mpmulx(/*const*/ MP a, /*const*/ MP b, MP c, int lmpnw)
  {
    double t1, t2, t3, t4;
    int i;
    
    int na = Math.min (a.nw, lmpnw);
    int nb = Math.min (b.nw, lmpnw);
    int ncr =(int)( Math.pow(2, mpmcr)); // *cast*
    
    /* Sweny
    if (na == 0 || nb == 0) 
    {
      //  One of the inputs is zero -- result is zero.
      zero(c);
      return;
    }
    */
    
    //  Check if precision level of one of the arguments is too low to justify the
    //  advanced routine.
    
    if (na <= ncr || nb <= ncr) 
    {
      mpmul (a, b, c, lmpnw);
      return;
    }
    double d1[] = new double[2*lmpnw+4], d2[] = new double[2*lmpnw+4], 
      d3[] = new double[4*lmpnw+8];
    
    //  Place the input data in A and B into the scratch arrays DD1 and DD2.
    //  This code also splits the input data into half-sized words.
    
    //!dir$ ivdep
    int i2;
    for (i = 0; i< na; i++)
    {
      i2 = 2 * i;
      t1 = a.mantissa[i];
      t2 = (int) (MPRBX * t1);
      d1[i2] = t2;
      d1[i2+1] = t1 - MPBBX * t2;
    }
    
    for(i = 2 * na; i<2 * nb; i++)
      d1[i] = 0.0;
    
    // !dir$ ivdep
    for(i = 0; i< nb; i++)
    {
      i2 = 2 * i;
      t1 = b.mantissa[i];
      t2 = (int) (MPRBX * t1);
      d2[i2] = t2;
      d2[i2+1] = t1 - MPBBX * t2;
    }
    
    for (i = 2 * nb; i< 2 * na; i++)
      d2[i] = 0.0;
    
    int nn = 2 * Math.max (na, nb);
    int nx = (int)(Math.sqrt (3.0 * nn) + MPRXX); // *cast*
    
    mplconv (2, nn, nx, d1, d2, d3);
    
    //  Recombine words and release carries.
    
    int nc = Math.min (na + nb, lmpnw);
    int nc1 = Math.min (lmpnw + 1, na + nb - 1);
    d1[0] = fSign (nc, (!(a.sign^b.sign)) ? 1.0 : -1.0);
    d1[1] = a.exponent + b.exponent + 1;
    d1[2] = d3[0];
    d1[nc+2] = 0.0;
    d1[nc+3] = 0.0;
    
    //!dir$ ivdep
    for(i = 1; i<= nc1; i++)
    {
      i2 = 2 * i;
      t3 = d3[i2-1];
      t4 = d3[i2];
      t1 = (int) (MPRDX * t3);
      t2 = t3 - MPBDX * t1;
      t3 = (int) (MPRDX * t4);
      t4 -= MPBDX * t3;
      d1[i+2] = MPBBX * t2 + t4;
      d1[i+1] += MPBBX * t1 + t3;
    }
    //  Fix up the result.
    mpnorm (d1, c, lmpnw);  
  }
  
  /**
  *  This computes the N-th power of the MP number A and returns the MP result
  *  in B. 
  *
  *  When N is zero, 1 is returned.  When N is negative, the reciprocal
  *  of A ^ |N| is returned.  Before calling MPNPWX, the arrays UU1 and UU2
  *  must be initialized by calling MPINIX.  For modest levels of precision, use
  *  MPNPWR.  Debug output starts with MPIDB = 6.
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw+2).
  *  <br>
  *  This routine employs the binary method for exponentiation.
  *  @see #mpnpwr
  *  @see #mpinix
  */
  static void mpnpwx (/*const*/ MP a, int n, MP b, int lmpnw)
  {
    int j;
    double t1;
    int mpnw2 =lmpnw+2;
    MP f1 =  new MP(6, false) , 
      sk0 = new MP(mpnw2,false), sk1 = new MP(mpnw2,false);
    
    
    int ncr = (int)(Math.pow(2, mpmcr));
    int na = Math.min (a.nw, lmpnw);
    
    //  Check if precision level of A is too low to justify the advanced routine.
    if (na <= ncr && n >= 0 && n <= 4) 
    {
      mpnpwr (a, n, b, lmpnw);
      return;
    }
    
    if (na == 0) 
    {
      if (n >= 0) 
      {
        zero(b);
        return;
      } 
      else 
      {
        throw new ArithmeticException(
          "mpnpwx: Argument is zero and n is negative or zero");
      }
    }
    
    int nn = Math.abs(n);
    f1.sign=true;
    f1.nw=1;
    f1.exponent=0;
    f1.mantissa[0]=1;
    f1.mantissa[1]=0;
    
    boolean skip=false;
    switch(nn)
    {
    case 0:
      mpeq (f1, b, lmpnw);
      return;
      
    case 1:
      mpeq (a, b, lmpnw);
      skip=true;
      break;
      
    case 2:
      mpsqx (a, b, lmpnw);
      skip=true;
      break;
    }
    
    if(!skip)
    {
      //  Determine the least integer MN such that 2 ^ MN > NN.
      t1 = nn;
      int mn = (int)(CL2 * Math.log (t1) + 1.0 + MPRXX); // cast
      mpeq (f1, b, lmpnw);
      mpeq (a, sk0, lmpnw);
      int kn = nn;    
      
      //  Compute B ^ N using the binary rule for exponentiation.
      for (j = 1; j<= mn; j++)
      {
        int kk = kn / 2;
        if (kn != 2 * kk) 
        {
          mpmulx (b, sk0, sk1, lmpnw);
          mpeq (sk1, b, lmpnw);
        }
        
        kn = kk;
        if (j < mn) 
        {
          mpsqx (sk0, sk1, lmpnw);
          mpeq (sk1, sk0, lmpnw);
        }
      }
    }
    
    //  Compute reciprocal if N is negative.
    
    if (n < 0) 
    {
      mpdivx (f1, b, sk0, lmpnw);
      mpeq (sk0, b, lmpnw);
    }  
  }

 
  /**
  *  This computes the N-th root of the MP number A and returns the MP result
  *  in B.
  *
  *  N must be at least one and must not exceed 2 ^ 30.  Before calling
  *  MPNRTX, the arrays UU1 and UU2 must be initialized by calling MPINIX.  For 
  *  modest levels of precision, use MPNRT.  Debug output starts with MPIDB = 6.
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw+2).
  *  <br>
  *  This routine uses basically the same Newton iteration algorithm as MPNRT.
  *  In fact, this routine calls MPNRT to obtain an initial approximation.
  *  See the comment about the parameter NIT in MPDIVX.
  *  @see #mpdivx.
  *  @see #mpnrt
  *  @see #mpinix
  */
  static void mpnrtx (/*const*/ MP a, int n, MP b, int lmpnw)
  {
    int k;
    double t2, tn;
    
    int mpnw3 = lmpnw+3;
    MP f1 =  new MP(6, false) , f2 =  new MP(6, false) , 
      sk0 = new MP(mpnw3,false), sk1 = new MP(mpnw3,false), 
      sk2 = new MP(mpnw3,false), sk3 = new MP(mpnw3,false);
    
    int ncr = (int)(Math.pow(2, mpmcr)); // *cast*
    
    int na = Math.min (a.nw, lmpnw);
    
    if (na == 0) 
    {
      zero(b);
      return;
    }
    if (!a.sign) 
    {
      throw new ArithmeticException(
        "mpnrtx: Argument is negative --> " + a);
    }
    
    if (n <= 0 || n > N30) 
    {
      throw new ArithmeticException(
        "mpnrtx: Improper value of N --> " + n);
    }
    
    //  Check if precision level is too low to justify the advanced routine.
    if (lmpnw <= ncr) 
    {
      mpnrt (a, n, b, lmpnw);
      return; 
    }
    
    //  If N = 1, 2 or 3,  MPEQ, MPSQRX or MPCBRX.  These are faster.
    switch (n)
    {
    case 1:
      mpeq (a, b, lmpnw);
      return;
    case 2:
      mpsqrx (a, b, lmpnw);
      return;
    case 3:
      mpcbrx (a, b, lmpnw);
      return;
    }
    
    int nws = lmpnw;
    f1.nw=1;
    f1.sign=true;
    f1.exponent=0;
    f1.mantissa[0]=1;
    f1.mantissa[1]=0;
    
    //  Determine the least integer MQ such that 2 ^ MQ >= MPNW.
    MPDPE t1 = new MPDPE();
    t1.a = lmpnw;
    int mq = (int)(CL2 * Math.log (t1.a) + 1.0 - MPRXX); // *cast*
    
    //  Check how close A is to 1.
    mpsub (a, f1, sk0, lmpnw);
    if (sk0.nw == 0) 
    {
      mpeq (f1, b, lmpnw);
      return;
    }
    
    mpmdc (sk0, t1);
    int n2 = (int)(CL2 * Math.log (Math.abs (t1.a))); // *cast*
    t1.a *=  Math.pow(0.5, n2);
    t1.n += n2;
    if (t1.n <= -30) 
    {
      t2 = n;
      n2 = (int)(CL2 * Math.log (t2) + 1.0 + MPRXX); // *cast*
      int n3 = - MPNBT * lmpnw / t1.n;
      if (n3 < 1.25 * n2) 
      {
        //  A is so close to 1 that it is cheaper to use the binomial series.
        mpdivd (sk0, new MPDPE(t2), sk1, lmpnw);
        mpadd (f1, sk1, sk2, lmpnw);
        k = 0;
        do 
        {
          k++;
          t1.a = 1 - k * n;
          t2 = (k + 1) * n;
          mpmuld (sk1, new MPDPE(t1.a), sk3, lmpnw);
          mpdivd (sk3, new MPDPE(t2), sk1, lmpnw);
          mpmulx (sk0, sk1, sk3, lmpnw);
          mpeq (sk3, sk1, lmpnw);
          mpadd (sk1, sk2, sk3, lmpnw);
          mpeq (sk3, sk2, lmpnw);
        }
        while (sk1.nw != 0 && sk1.exponent >= - lmpnw);        
        
        mpeq (sk2, b, lmpnw);
        mproun (b, nws);
        
        return;
      }
    }
    
    //  Compute the initial approximation of A ^ (-1/N).
    lmpnw = ncr + 1;
    mpnrt (a, n, sk0, lmpnw);
    mpdiv (f1, sk0, b, lmpnw);
    
    tn = n;
    MPDPE dpn = new MPDPE(n, 0);
    mpdmc (dpn, f2);
    int iq = 0;
    
    //  Perform the Newton-Raphson iteration described above with a dynamically
    //  changing precision level MPNW (powers of two).  
    int nw1, nw2;
    for (k = mpmcr+1; k<= mq;k++)
    {
      nw1 = lmpnw;
      
      lmpnw = Math.min (2 * lmpnw - 2, nws) + 1;
      nw2 = lmpnw;
      boolean loop=true;
      while (loop)
      {
        
        mpnpwx (b, n, sk0, lmpnw);
        mpmulx (a, sk0, sk1, lmpnw);
        mpsub (f1, sk1, sk0, lmpnw);
        lmpnw = nw1;
        mpmulx (b, sk0, sk1, lmpnw);
        mpdivd (sk1, new MPDPE(tn), sk0, lmpnw);
        lmpnw = nw2;
        mpadd (b, sk0, sk1, lmpnw);
        mpeq (sk1, b, lmpnw);
        if (k == mq - NIT && iq == 0) 
          iq = 1;
        else loop=false;
      }
    }
    
    //  Take the reciprocal to give final result.
    mpdivx (f1, b, sk0, lmpnw);
    mpeq (sk0, b, lmpnw);
    
    mproun (b, nws);
  }

  /**
  *  This computes the square root of the MP number A and returns the MP result
  *  in B.
  *
  *  Before calling MPSQRX, the arrays UU1 and UU2 must be initialized by
  *  calling MPINIX.  For modest levels of precision, use MPSQRT.  Debug output
  *  starts with MPIDB = 6.
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw+2).
  *  <br>
  *  This routine uses basically the same Newton iteration algorithm as MPSQRT.
  *  In fact, this routine calls MPSQRT to obtain an initial approximation.
  *  See the comment about the parameter NIT in MPDIVX.
  *  @see #mpsqrt
  *  @see #mpdivx
  *  @see #mpinix
  */
  static void mpsqrx (/*const*/ MP a, MP b, int lmpnw)
  {
    int k;
    
    int mpnw3 = lmpnw+3;
    MP f =  new MP(6, false) , sk0 = new MP(mpnw3,false), 
      sk1 = new MP(mpnw3,false), sk2 = new MP(mpnw3,false);
    
    int na = Math.min (a.nw, lmpnw);
    int ncr = (int)(Math.pow(2, mpmcr)); // *cast*
    
    if (na == 0) 
    {
      zero(b);
      return;
    }
    
    if (!a.sign) 
    {
      throw new ArithmeticException(
        "mpsqrx: Argument is negative --> " + a);
    }
    
    //  Check if precision level is too low to justify the advanced routine.  
    if (lmpnw <= ncr) 
    {
      mpsqrt (a, b, lmpnw);
      return;
    }
    int nws = lmpnw;
    
    //  Determine the least integer MQ such that 2 ^ MQ >= MPNW.
    double t1=lmpnw;
    int mq = (int)(CL2 * Math.log (t1) + 1.0 - MPRXX); // *cast* 
    
    //  Compute the initial approximation of 1 / Sqrt(A).
    
    lmpnw = ncr + 1;
    mpsqrt (a, sk0, lmpnw);
    mpdiv (sk0, a, b, lmpnw);
    
    f.nw=1;
    f.sign=true;
    f.exponent=0;
    f.mantissa[0]=1;
    f.mantissa[1]=0;
    int iq = 0;
    
    //  Perform the Newton-Raphson iteration described above with a dynamically
    //  changing precision level MPNW (one greater than powers of two).
    
    int nw1, nw2;
    for(k = mpmcr + 1; k<= mq - 1; k++)
    {
      nw1 = lmpnw;
      lmpnw = Math.min (2 * lmpnw - 2, nws) + 1;
      nw2 = lmpnw;
      boolean stop=false;
      while(!stop)
      {   
        mpsqx (b, sk0, lmpnw);
        mpmulx (a, sk0, sk1, lmpnw);
        mpsub (f, sk1, sk0, lmpnw);
        lmpnw = nw1;
        mpmulx (b, sk0, sk1, lmpnw);
        mpmuld (sk1, new MPDPE(0.50), sk0, lmpnw);
        lmpnw = nw2;
        mpadd (b, sk0, sk1, lmpnw);
        mpeq (sk1, b, lmpnw);
        if (k == mq - NIT && iq == 0) 
          iq = 1;
        else
          stop=true;
      }
    }
    //  Perform last iteration using Karp's trick.
    
    mpmulx (a, b, sk0, lmpnw);
    nw1 = lmpnw;
    lmpnw = Math.min (2 * lmpnw - 2, nws) + 1;
    
    nw2 = lmpnw;
    
    mpsqx (sk0, sk1, lmpnw);
    mpsub (a, sk1, sk2, lmpnw);
    lmpnw = nw1;
    mpmulx (sk2, b, sk1, lmpnw);
    mpmuld (sk1, new MPDPE(0.50), sk2, lmpnw);
    lmpnw = nw2;
    mpadd (sk0, sk2, sk1, lmpnw);
    mpeq (sk1, b, lmpnw);
    
    mproun (b, nws);
  }
  
  /**
  *  This routine squares the MP number A to yield the MP product B.
  *
  *  Before calling MPSQX, the arrays UU1 and UU2 must be initialized by calling
  *  MPINIX.  For modest levels of precision, use MPMUL.  MPNW should be a power
  *  of two.  Debug output starts with MPIDB = 8.
  *  <br>
  *  Dimension a(lmpnw), b(lmpnw+2).
  *  <br>
  *  This subroutine uses the same FFT technique as MPMULX.  It is faster
  *  because only one forward FFT has to be computed.  See the comments in
  *  MPMULX about obtaining the complete double-long result.
  *  @see #mpinix
  *  @see #mpmul
  */
  static void mpsqx(/*const*/ MP a, MP b, int lmpnw)
  {
    double t1, t2, t3, t4;
    int i;
    
    int na = Math.min (a.nw, lmpnw);
    int ncr = (int)(Math.pow(2, mpmcr));
    
    if (na == 0) 
    {
      zero(b);
      return;
    }
    
    //  Check if precision level of the argument is too low to justify the
    //  advanced routine.
    
    if (na <= ncr) 
    {
      mpmul (a, a, b, lmpnw);
      return;
    }
    
    double d1[] = new double[2*lmpnw+4], d2[] = new double[4*lmpnw+8];
    
    //  Place the input data in A into the scratch array DD1.
    //  This code also splits the input data into half-sized words.
    
    //!dir$ ivdep
    int i2;
    for(i = 0; i<na; i++)
    {
      i2 = 2 * i ;
      t1 = a.mantissa[i];
      t2 = (int) (MPRBX * t1);
      d1[i2] = t2;
      d1[i2+1] = t1 - MPBBX * t2;
    }
    
    int nn = 2 * na;
    int nx = (int)(Math.sqrt (3.0 * nn) + MPRXX); // *cast*
    mplconv (1, nn, nx, d1, d1, d2);
    
    //  Recombine words and release carries.
    
    int nc = Math.min (2 * na, lmpnw);
    int nc1 = Math.min (lmpnw + 1, 2 * na - 1);
    d1[0] = nc;
    d1[1] = 2 * a.exponent + 1;
    d1[2] = d2[0];
    d1[nc+2] = 0.0;
    d1[nc+3] = 0.0;
    
    //!dir$ ivdep
    for(i = 1; i<= nc1; i++)
    {
      i2 = 2 * i;
      t3 = d2[i2-1];
      t4 = d2[i2];
      t1 = (int) (MPRDX * t3);
      t2 = t3 - MPBDX * t1;
      t3 = (int) (MPRDX * t4);
      t4 -= MPBDX * t3;
      d1[i+2] = MPBBX * t2 + t4;
      d1[i+1] += MPBBX * t1 + t3;
    }
    
    //  Fix up the result.
    mpnorm (d1, b, lmpnw);
  }
 
  //************* ADVANCE ALGEBRAIC & TRANSCEDENTAL ROUTINES ******************
  
  
  /**
  *  This performs the arithmetic-geometric mean (AGM) iterations. 
  *
  *  This routine is called by MPLOGX.  It is not intended to be called directly by the user.
  *  <br>
  *  Dimension a(lmpnw+2), b(lmpnw+2).
  *  @see #mplogx
  */
  private static void mpagmx (MP a, MP b, int lmpnw)
  {
    int mpnw2 = lmpnw+2;
    MP sk0 = new MP(mpnw2,false), sk1 = new MP(mpnw2,false);  
    
    // sk0 is initialized to zero by constructor
    int l1 = 0;
    int s1;
    // 100 
    MPDPE dpe1 = new MPDPE(0.50, 0);
    do
    { 
      l1++;
      if (l1 == 50) 
      {
        throw new ArithmeticException(
          "mpagmx: Iteration limit exceeded.");
      }
      s1 = sk0.exponent;
      mpadd (a, b, sk0, lmpnw);
      mpmuld (sk0, dpe1, sk1, lmpnw);
      mpmulx (a, b, sk0, lmpnw);
      mpsqrx (sk0, b, lmpnw);
      mpeq (sk1, a, lmpnw);
      mpsub (a, b, sk0, lmpnw);
    }
    //  Check for convergence.
    while (sk0.nw != 0. && (sk0.exponent < s1 || sk0.exponent >= -2));
  }
  
  /**
  *  This computes the hyperbolic cosine and sine of the MP number A and
  *  returns the two MP results in X and Y, respectively. 
  *  
  *  PI is the MP value of Pi computed by a previous  to MPPI or MPPIX.  
  *  AL2 is the MP value of Log (10) computed by a previous  to MPLOG or MPLOGX.  
  *  Before calling MPCSHX, the arrays UU1 and Uu2 must be initialized by calling
  *  MPINIX.  For modest levels of precision, use MPCSSH.  The last word 
  *  of the result is not reliable.  Debug output starts with MPIDB = 5.
  *  <br>
  *  Dimension a(lmpnw), al2(lmpnw), x(lmpnw+2), y(lmpnw+2), pi(lmpnw).
  *  @see #mppi
  *  @see #mppix
  *  @see #mplog
  *  @see #mplogx
  *  @see #mpcssh
  *  @see #mpinix
  */
  static void mpcshx (/*const*/ MP a, /*const*/ MP pi, /*const*/ MP al2, 
    MP x, MP y, int lmpnw)
  {
    int mpnw2 = lmpnw+2;
    MP f =  new MP(6, false) , sk0 = new MP(mpnw2,false), 
      sk1 = new MP(mpnw2,false), sk2 = new MP(mpnw2,false);
    
    
    f.sign=true; f.nw=1; f.exponent=0;
    f.mantissa[0]=1; f.mantissa[1]=0;
    
    MPDPE dpe1 = new MPDPE(0.5);
    mpexpx (a, pi, al2, sk0, lmpnw);
    mpdivx (f, sk0, sk1, lmpnw);
    mpadd (sk0, sk1, sk2, lmpnw);
    mpmuld (sk2, dpe1, x, lmpnw);
    mpsub (sk0, sk1, sk2, lmpnw);
    mpmuld (sk2, dpe1, y, lmpnw);
  }
 
  /**
  *  This computes the exponential function of the MP number A and returns the
  *  MP result in B.  
  *
  *  PI is the MP value of Pi produced by a prior  to MPPI or MPPIX.  
  *  AL2 is the MP value of Log(2) produced by a prior  to
  *  MPLOG  or MPLOGX.  Before calling MPEXPX, the arrays UU1 and UU2 must be
  *  initialized by calling MPINIX.  For modest levels of precision, use MPEXP.
  *  The last word of the result is not reliable.  Debug output starts 
  *  with MPIDB = 5.  
  *  <br>
  *  Dimension a(lmpnw), al2(lmpnw), b(lmpnw+2), pi(lmpnw).
  *  <br>
  *  This routine uses the Newton iteration
  *  <pre>
  *    b_{k+1} = b_k [a + 1 - log b_k]
  *  </pre>
  *  with a dynamically changing level of precision.  Logs are performed using
  *  MPLOGX.  See the comment about the parameter NIT in MPDIVX.
  *  @see #mppi
  *  @see #mppix
  *  @see #mplog
  *  @see #mplogx
  *  @see #mpinix
  *  @see #mpepx
  */
  static void mpexpx (/*const*/ MP a, /*const*/ MP pi, /*const*/ MP al2, MP b, int lmpnw)
  {
    int k;
    final int nit = 1;
    int mpnw2 =lmpnw+2;
    MP f1 =  new MP(6, false) , sk0 = new MP(mpnw2,false), 
      sk1 = new MP(mpnw2,false), sk2 = new MP(mpnw2,false);
    
    int ncr = (int)(Math.pow(2, mpmcr));
    
    MPDPE t1 = new MPDPE();
    mpmdc (a, t1);
    t1.a *=  Math.pow(2.0, t1.n);
    
    //  Check if precision level is too low to justify the advanced routine.
    
    if (lmpnw <= ncr) 
    {
      mpexp (a, al2, b, lmpnw);
      return;
    }
    
    //  Check if Log(2) has been precomputed.
    
    MPDPE t2 = new MPDPE();
    mpmdc (al2, t2);
    if (t2.n != - MPNBT || Math.abs (t2.a * Math.pow(0.50, MPNBT) - ALT) > MPRX2)
    {
      throw new ArithmeticException(
        "mpexpx: LOG(2) must be precomputed.");
    }
    
    
    //    Check if Pi has been precomputed.
    
    mpmdc (pi, t2);
    if (t2.n != 0 || Math.abs (t2.a - PI) > MPRX2) 
    {
      throw new ArithmeticException(
        "mpexpx: PI must be precomputed.");
    }
    
    //  Check for overflows and underflows.
    
    if (t1.a >= 1e9) 
    {
      if (t1.a > 0.0) 
      {
        throw new ArithmeticException
          ("mpexpx: Argument is too large --> " + 
          t1.a + " x 10 ^" + t1.n);
      } 
      else 
      {
        zero(b);
        return;
      }
    }
    
    int nws = lmpnw;
    f1.sign=true; f1.nw=1; f1.exponent=0;
    f1.mantissa[0]=1; f1.mantissa[1]=0;
    
    //  Determine the least integer MQ such that 2 ^ MQ >= MPNW.
    t2.a = nws;
    int mq = (int)(CL2 * Math.log (t2.a) + 1.0 - MPRXX);
    
    mpadd (a, f1, sk0, lmpnw);
    
    //  Compute initial approximation to Exp (A).
    
    lmpnw = ncr;
    mpexp (a, al2, b, lmpnw);
    int iq = 0;
    
    //  Perform the Newton-Raphson iteration described above with a dynamically
    //  changing precision level MPNW.
    for(k = mpmcr + 1; k<=mq; k++)
    {
      lmpnw = Math.min (2 * lmpnw, nws);
      boolean cont = true;
      while(cont)
      {
        mplogx (b, pi, al2, sk1, lmpnw);
        mpsub (sk0, sk1, sk2, lmpnw);
        mpmulx (b, sk2, sk1, lmpnw);
        mpeq (sk1, b, lmpnw);
        if (k == mq - nit && iq == 0) 
          iq = 1;
        else
          cont=false;
      }
    }
  }
  
  /**
  *  This computes the natural logarithm of the MP number A and returns the MP
  *  result in B.  
  *
  *  PI is the MP value of Pi produced by a prior  to MPPI or MPPIX.  
  *  AL2 is the MP value of Log(2) produced by a prior  to MPLOG
  *  or MPLOGX.  Before calling MPLOGX, the arrays UU1 and UU2 must be
  *  initialized by calling MPINIX.  For modest levels of precision, use MPLOG.
  *  The last word of the result is not reliable.  Debug output starts 
  *  with MPIDB = 6.
  *  <br>
  *  Dimension al2(lmpnw), pi(lmpnw), a(lmpnw), b(lmpnw+2).
  *  <br>
  *  This uses the following algorithm, which is due to Salamin.  If a is
  *  extremely close to 1, use a Taylor series.  Otherwise select n such that
  *  z = x 2^n is at least 2^m, where m is the number of bits of desired
  *  precision in the result.  Then
  *  <pre>
  *  Log(x) = Pi / [2 AGM (1, 4/x)]
  *  </pre>
  *  @see #mppi
  *  @see #mppix
  *  @see #mplog
  *  @see #mplogx
  *  @see #mpinix
  */
  static void mplogx (/*const*/ MP a, /*const*/ MP pi, /*const*/ MP al2, MP b, int lmpnw)
  {
    final int mzl = -5;
    
    double st, tn;
    int mpnw2 = lmpnw+2;
    MP f1 =  new MP(6, false) , f4 =  new MP(6, false) , 
      sk0 = new MP(mpnw2,false), sk1 = new MP(mpnw2,false), 
      sk2 = new MP(mpnw2,false), sk3 = new MP(mpnw2,false);
    
    int na = Math.min (a.nw, lmpnw);
    int ncr = (int)(Math.pow(2, mpmcr)); // *cast*
    
    //  Check if precision level is too low to justify the advanced routine.
    
    if (lmpnw <= ncr) 
    {
      mplog (a, al2, b, lmpnw);
      return;    
    }
    
    if (!a.sign || na == 0) 
    {
      throw new ArithmeticException
        ("mplogx: Argument is less than or equal to zero -->" 
        + a);
    }
    
    //  Check if Pi has been precomputed.
    MPDPE t1 = new MPDPE();
    mpmdc (pi, t1);
    if (t1.n != 0 || Math.abs (t1.a - CPI) > MPRX2) 
    {
      throw new ArithmeticException
        ("mplogx: PI must be precomputed.");
    }
    
    //  Unless the input is 2, Log (2) must have been precomputed.
    MPDPE t2 = new MPDPE();
    int it2;
    if (a.nw != 1 || a.exponent != 0 || a.mantissa[0] != 2 || !a.sign) 
    {
      it2 = 0;
      mpmdc (al2, t2);
      if (t2.n != - MPNBT || Math.abs (t2.a * Math.pow(0.50, MPNBT) - ALT) > MPRX2) 
      {
        throw new ArithmeticException
          ("mplogx: LOG(2) must be precomputed.");
      }
    }
    else  
      it2 = 1;
    
    
    f1.nw=1; f1.sign=true; f1.exponent=0; f1.mantissa[0]=1; f1.mantissa[1]=0;
    f4.nw=1; f4.sign=true; f4.exponent=0; f4.mantissa[0]=4; f4.mantissa[1]=0;
    
    
    //  If argument is 1, the result is zero.  If the argument is extremely close
    //  to 1.  If so, employ a Taylor's series instead.
    
    mpsub (a, f1, sk0, lmpnw);
    
    if (sk0.nw == 0) 
    {
      zero(b);
      return;
    } 
    else if (sk0.exponent <= mzl) 
    {
      mpeq (sk0, sk1, lmpnw);
      mpeq (sk1, sk2, lmpnw);
      int i1 = 1;
      int is = 1;
      int tl = sk0.exponent - lmpnw - 1;
      
      do
      {
        i1++;
        is = - is;
        st = is * i1;
        mpmulx (sk1, sk2, sk3, lmpnw);
        mpdivd (sk3, new MPDPE(st), sk2, lmpnw);
        mpadd (sk0, sk2, sk3, lmpnw);
        mpeq (sk3, sk0, lmpnw);
      }
      while (sk2.exponent >= tl); 
      
      mpeq (sk0, b, lmpnw);
      return;
    }
    
    //  If input is exactly 2, set the exponent to a large value.  Otherwise
    //  multiply the input by a large power of two.
    
    mpmdc (a, t1);
    t2.n = MPNBT * (lmpnw / 2 + 2) - t1.n;
    tn = t2.n;
    t2.a = 1.0;
    if (it2 == 1) 
      mpdmc (t2, sk0);
    else 
      mpmuld (a, t2, sk0, lmpnw);
    
    
    //  Perform AGM iterations. 
    mpeq (f1, sk1, lmpnw);
    mpdivx (f4, sk0, sk2, lmpnw);
    mpagmx (sk1, sk2, lmpnw);
    
    //  Compute B = Pi / (2 * A), where A is the limit of the AGM iterations.
    
    mpmuld (sk1, new MPDPE(2.0), sk0, lmpnw);
    mpdivx (pi, sk0, sk1, lmpnw);
    
    //  If the input was exactly 2, divide by TN.  Otherwise subtract TN * Log(2).
    
    if (it2 == 1) 
      mpdivd (sk1, new MPDPE(tn), sk0, lmpnw);
    else 
    {
      mpmuld (al2, new MPDPE(tn), sk2, lmpnw);
      mpsub (sk1, sk2, sk0, lmpnw);
    }
    mpeq (sk0, b, lmpnw);
  }

  /**
  *  This computes Pi to available precision (MPNW mantissa words). 
  *
  *  Before calling MPPIX, the arrays UU1 and UU2 must be initialized by calling
  *  MPINIX.  For modest levels of precision, use MPPI.  MPNW should be a power
  *  of two.  The last word of the result is not reliable.  Debug
  *  output starts with MPIDB = 7.
  *  <br>
  *  Dimension pi(lmpnw+2).
  *  <br>
  *  This routine uses basically the same algorithm as MPPI.
  *  @see #mppi
  *  @see #mpinix
  */
  static void mppix (MP pi, int lmpnw)
  {
    int k;
    double t1;
    int mpnw2 = lmpnw+2;
    MP f =  new MP(6, false) , sk0 = new MP(mpnw2,false), 
      sk1 = new MP(mpnw2,false), sk2 = new MP(mpnw2,false), 
      sk3 = new MP(mpnw2,false), sk4 = new MP(mpnw2,false);
    
    int ncr = (int)(Math.pow(2, mpmcr)); // *cast*
    
    //  Check if precision level is too low to justify the advanced routine.
    
    if (lmpnw <= ncr) 
    {
      mppi (pi, lmpnw);
      return;
    }
    
    //  Determine the number of iterations required for the given precision level.
    //  This formula is good only for this Pi algorithm.
    
    t1 = lmpnw * log10 (MPBDX);
    int mq = (int)(CL2 * (Math.log (t1) - 1.0) + 1.0); // *cast*
    
    //  Initialize as above.
    sk0.nw=1; sk0.sign=true; sk0.exponent=0; sk0.mantissa[0]=1;
    f.nw=1; f.sign=true; f.exponent=0; f.mantissa[0]=2; f.mantissa[1]=0;
    
    mpsqrx (f, sk2, lmpnw);
    mpmuld (sk2, new MPDPE(0.50), sk1, lmpnw);
    f.exponent = -1;
    f.mantissa[0] = (float)(0.50 * MPBDX); // *cast*
    mpsub (sk2, f, sk4, lmpnw);
    
    //  Perform iterations as described above.
    
    MPDPE dpe1 = new MPDPE(0.50);
    for (k = 1; k<= mq; k++)
    {
      mpadd (sk0, sk1, sk2, lmpnw);
      mpmulx (sk0, sk1, sk3, lmpnw);
      mpsqrx (sk3, sk1, lmpnw);
      mpmuld (sk2, dpe1, sk0, lmpnw);
      mpsub (sk0, sk1, sk2, lmpnw);
      
      mpsqx (sk2, sk3, lmpnw);
      t1 = Math.pow(2.0, k);
      
      mpmuld (sk3, new MPDPE(t1), sk2, lmpnw);
      
      mpsub (sk4, sk2, sk3, lmpnw);
      mpeq (sk3, sk4, lmpnw);
    }
    
    //  Complete computation.
    mpadd (sk0, sk1, sk2, lmpnw);
    mpsqx (sk2, sk3, lmpnw);
    mpdivx (sk3, sk4, sk2, lmpnw);
    mpeq (sk2, pi, lmpnw);
  }

  //******************** ADVANCE IO / DEBUG ROUTINES ***************************
  
  /**
  *  Converts the MP number A into character form in the char array B.
  *
  *  N (an output parameter) is the length of the output.  In other words, B is 
  *  contained in B[0], ..., B[N-1].  The format is analogous to the Fortran
  *  exponential format (E format), except that the exponent is placed first.
  *  Before calling MPOUTX, the arrays UU1 and UU2 must be initialized by
  *  calling MPINIX.  For modest levels of precision, use MPOUTC.  Debug output
  *  starts with MPIDB = 7.
  *  <br>
  *  Dimension a(lmpnw), b(7.225 * MPNW + 30 + 1).
  *  @see #mpinix
  *  @see #mpoutc
  */
  static int mpoutx (/*const*/ MP a, char b[], int lmpnw)
  {
    int i, n;  
    int na = Math.min (a.nw, lmpnw);
    int ncr = (int)(Math.pow(2, mpmcr));
    
    //  Check if actual precision level of argument is too low to justify the
    //  advanced routine.
    if (na <= ncr) 
      return mpoutc (a, b, lmpnw);
    
    double t1, t2;
    char c1[] = new char[17], c2[] = new char[17];
    char b1[] = new char[8*lmpnw+31], b2[] = new char[8*lmpnw+31];
    
    int mpnw2 = lmpnw+2;
    MP sk0 = new MP(mpnw2,false), sk1 = new MP(mpnw2,false),
      sk2 = new MP(mpnw2,false), sk3 = new MP(mpnw2,false), 
      sk4 = new MP(mpnw2,false);
    
    //  Normalize input to an integer by multiplying by a suitable power of 10.
    
    t1 = a.mantissa[0] + MPRDX * a.mantissa[1] + MPRX2 * a.mantissa[2];
    t2 = log10 (t1);
    int m1 = (int) (Math.max (AL2 * MPNBT * (a.nw - a.exponent) - t2, 0.0)); // *cast*
    MPDPE dpe1 = new MPDPE(10,0); 
    mpdmc (dpe1, sk0);
    mpnpwx (sk0, m1, sk2, lmpnw);
    mpmulx (a, sk2, sk1, lmpnw);
    sk1.sign = true;
    
    //  Split large integer into two approximately equal decimal sections.
    MPDPE dpe2 = new MPDPE(), dpe3 = new MPDPE();
    mpmdc(sk1, dpe2);
    
    MPDPE.dpdec (dpe2, dpe3);
    
    int m2 = dpe3.n / 2;
    mpnpwx (sk0, m2, sk3, lmpnw);
    mpdivx (sk1, sk3, sk0, lmpnw);
    mpinfr (sk0, sk2, sk4, lmpnw);
    mpmulx (sk2, sk3, sk0, lmpnw);
    mpsub (sk1, sk0, sk3, lmpnw);
    
    //  Recursively convert each section.
    
    int mpnws = lmpnw;
    lmpnw = sk2.nw + 1;
    int nb1, nb2;
    nb1 = mpoutx (sk2, b1, lmpnw);
    lmpnw = sk3.nw + 1;
    nb2 = mpoutx (sk3, b2, lmpnw);
    lmpnw = mpnws;
    
    //  Obtain decimal exponents from each section.
    for(i = 0; i<17; i++)
    { 
      c1[i] = '\0';
      c2[i] = '\0';
    }
    
    for(i = 0; i<10; i++)
    {
      c1[i] = b1[i+4];
      c2[i] = b2[i+4];
    }
    
    String tempStr = new String(c1);
    tempStr = tempStr.substring(0, tempStr.indexOf('\0'));
    int ie1 = Integer.parseInt(tempStr);
    tempStr = new String(c2);
    tempStr = tempStr.substring(0, tempStr.indexOf('\0'));
    int ie2 = Integer.parseInt(tempStr);
    
    //  Set exponent of result.
    
    int ie = ie1 + m2 - m1;
    
    c1 = (new Integer(ie)).toString().toCharArray();
    
    for (i = 0; i< 4; i++)
      b[i] = b1[i];
    
    for(;i < 14 - c1.length; i++)
      b[i] = ' ';
    
    int ii=0;
    for(;i<14;i++)
      b[i] = c1[ii++];
    
    //  Copy mantissa of first section.
    
    for(i = 14; i<nb1; i++)
      b[i] = b1[i];
    
    int i2, i1 = ie1 + m2 - ie2 + 19;
    
    //  If first section is too long, then round trailing digits (probably 9s).
    if (nb1 > i1) 
    {
      i2 = Integer.parseInt(new String(b,i1,1));
      if (i2 >= 5) 
      {
        boolean skip = false;
        for(i = i1-1; i>=20; i--)
        {
          if (b[i] != '9')
          {
            skip = true;
            break;
          }
          b[i] = '0';
        }
        
        if(!skip)
        {
          throw new ArithmeticException
            ("mpoutx: Exceptional case -- contact author.");
        }
        
        i2 = Integer.parseInt(new String(b,i,1));
        c1 = (new Integer(i2+1)).toString().toCharArray();
        b[i] = c1[0];
      }
    } 
    else if (nb1 < i1) 
    {
      
      //  If first section is too short, then insert zeroes in gap.
      
      for(i = nb1; i< i1; i++)
        b[i] = '0';
    }
    
    //  Copy mantissa of second section.
    
    b[i1] = b2[18];
    n = Math.min (i1 + nb2 - 19, (int) (7.225 * lmpnw + 30));
    
    for(i = i1 + 1; i< n; i++)
      b[i] = b2[i-i1+19];
    
    //  Fix fSign.
    
    if (!a.sign) b[17] = '-';
    
    // Put null terminating character;
    b[n] = '\0';
    return n;
  }
}

