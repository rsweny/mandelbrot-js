package math;

/**  
* This class is the ancestor class for all other classes in the package.
* It contains constants, parameters, short mathematical methods, and
* initializer methods.
*
* @author Herman Harjono
* @version Oct 5, 1998.
*/ 

public abstract class MPGlobal
{ 
  static final double AL2 = 0.301029995663981195;
  static final double CL2 = 1.4426950408889633;
  static final double PI = 3.141592653589793238;
  static final double CPI = 3.141592653589793;
  static final double ALT = 0.693147180559945309;
  static final double LOGE10 = Math.log(10);
  static final int NIT = 3;
  static final int N30 = (int)(Math.pow(2,30));

  /**  
  *  Depends on the system word size.
  * 
  *  On IEEE systems and most other 32 bit systems, set to 4096.D0
  *  On Cray system, set to 2048.D0.
  */
  static final double MPBBX = 4096.0;
  static final double MPBDX=MPBBX*MPBBX;  
  static final double MPBX2=MPBDX*MPBDX; 
  static final double MPRBX=1.0/MPBBX;
  static final double MPRDX=MPRBX*MPRBX; 
  static final double MPRX2=MPRDX*MPRDX; 
  static final double MPRXX=16.0*MPRX2;
  
  /**  
  *  Depends on the system word size.
  * 
  *  On IEEE systems and most other 32 bit systems, set to 24.
  *  On Cray system, set to 22.
  */
  static final int MPNBT=24;
  
  /**  
  *  Depends on the system word size.
  *
  *  On IEEE systems and most other 32 bit systems, set to 32.
  *  On Cray system, set to 16.
  */
  static final int MPNPR=32;
  
  /**  
  *  Depends on the system word size.
  *
  *  On IEEE systems and most other 32 bit systems, set to 7.
  *  On Cray system, set to 8.
  */
  static final int MPMCRX=7;
  
  /**
  *  This is spacing parameters to avoid bank and cache conflicts in the FFT routines.
  *
  *  Value of 32 appears to work well on most systems.
  *  Set MPNROW = 64 on Crays.
  */
  static final int MPNROW=16;
  
  /**
  *  This is spacing parameters to avoid bank and cache conflicts in the FFT routines.
  *
  *  Value of 3 appears to work well on most systems. 
  */
  static final int MPNSP1=2; 
  
  /**
  *  This is spacing parameters to avoid bank and cache conflicts in the FFT routines.
  *
  *  Value of 17 appears to work well on most systems.
  */
  static final int MPNSP2=9;
  
  /**
  *  Cross-over point for advanced routines. 
  *  Initial value: MPMCRX.
  *  @see #MPMCRX
  */
  static int mpmcr;
  
  /**
  *  MPFUN rounding option. Initial value: 1.
  */
  static int mpird;
  
  /**
  *  Precision level, in digits
  */
  static int mpipl;
  
  /**
  *  Output precision, in digits.
  */
  static int mpiou; 
  
  /**
  *  Log_10 of MP epsilon level
  */
  static int mpiep; 
  
  /**
  *  Precision level, in words.
  */
  static int mpwds; 
  
  /**
  * Current output precision (in digits).
  */
  static int mpoud;
  
  static int mp2; 
  
  static int mp21;
  
  static int mpnw; 
  
  //static MPSize DEFAULTSIZE;
  
  
  static DComplex mpuu1[]; 
  static DComplex mpuu2[];
  
  /** 
  * Static Initializer
  */
  static
  {
     setMaximumPrecision(28);
  }

  private static int currMpipl;
  /**
   *  Sets the maximum precision level for all threads (in decimal digits).
   */
  public static final void setMaximumPrecision(int p)
  {
    if(p == currMpipl)
       return;

    currMpipl = p;
    mpmcr=MPMCRX;
    mpird=1;
    
    mpipl = currMpipl; 
    mpiou = 56; 
    mpiep = 10 - mpipl; 
    mpwds = (int)((mpipl /  7.224719896) + 1);
    mp2 = mpwds + 2;
    mp21 = mp2 + 1;
    mpoud = mpiou; 
    mpnw = mpwds;
    //const MPSize DEFAULTSIZE = MPSize(int((mpipl / 7.224719896) + 4));
    
    // initializing mpuu1
    int n = mpnw+1;
    double t1, ti;
    int j, i;
    
    t1 = 0.75 * n;
    int m = (int)(CL2 * Math.log (t1) + 1.0 - MPRXX); // *cast*
    int mq = m + 2;
    int nq = (int)(Math.pow(2, mq)); // *cast*
    mpuu1 = new DComplex[nq];
    mpuu1[0] = new DComplex(mq,0);
    
    int ku = 1;
    int ln = 1;
    
    for(j = 1; j<= mq; j++)
    {
      t1 = PI / ln;      
      for (i = 0; i<= ln - 1; i++)
      {
        ti = i * t1;
        mpuu1[i+ku] = new DComplex(Math.cos (ti), Math.sin (ti));
      }
      
      ku +=  ln;
      ln *= 2;
    }
    
    // initializing mpuu2
    double tpn;
    int k;
    
    t1 = 0.75 * n;
    mpuu2 = new DComplex [mq+nq];
    
    ku = mq + 1;
    mpuu2[0] = new DComplex(mq,0);
    
    int mm, nn, mm1, mm2, nn1, nn2, iu;
    for (k = 2; k<= mq - 1; k++)
    {
      mpuu2[k-1] = new DComplex(ku,0);
      mm = k;
      nn = (int)(Math.pow(2, mm)); // *cast*
      mm1 = (mm + 1) / 2;
      mm2 = mm - mm1;
      nn1 = (int)(Math.pow(2, mm1)); // *cast*
      nn2 = (int)(Math.pow(2, mm2)); // *cast*
      tpn = 2.0 * PI / nn;
      
      for(j = 0; j<nn2; j++)
      {
        for(i = 0; i<nn1; i++)
        {
          iu = ku + i + j * nn1;
          t1 = tpn * i * j;
          mpuu2[iu-1] = new DComplex(Math.cos (t1), Math.sin (t1));
        }
      }
      ku += nn;
    }
    
  }
   


  public static final double nint(double x)
  {
     if(x < 0)
	return Math.ceil(x - .5);
     else
	return Math.floor(x + .5);
  }
  /** 
  *  Fortran's sign function.
  *  @return the absolute value of a with sign of b.
  */
  public static final double fSign(double a, double b)
  { return (b>=0 ? Math.abs(a) : -Math.abs(a));}

  /** 
  *   Log base 10.
  */
  public static final double log10(double val)
  {
    return (Math.log(val)/LOGE10);
  }

  static final int precisionToSize(int precision)
  {
    int maxnw=(int)((precision / 7.224719896) + 4);
    if(maxnw<8) maxnw=8;  // since at least we need 6 (such as in mpdmc), therefore 
	  // it's safe to allocate at least 8.
    return maxnw;
  }
  
public static int[] HSBtoRGB(float hue, float saturation, float brightness)
{
	if (saturation == 0)
	return convert(brightness, brightness, brightness, 0);
	
	hue = hue - (float) Math.floor(hue);
	int i = (int) (6 * hue);
	float f = 6 * hue - i;
	float p = brightness * (1 - saturation);
	float q = brightness * (1 - saturation * f);
	float t = brightness * (1 - saturation * (1 - f));
	switch (i)
	{
		case 0:
			return convert(brightness, t, p, 0);
		case 1:
			return convert(q, brightness, p, 0);
		case 2:
			return convert(p, brightness, t, 0);
		case 3:
			return convert(p, q, brightness, 0);
		case 4:
			return convert(t, p, brightness, 0);
		default:
			return convert(brightness, p, q, 0);
		
	}
}

private static int[] convert(float red, float green, float blue, float alpha)
{
	int redval = Math.round(255 * red);
	int greenval = Math.round(255 * green);
	int blueval = Math.round(255 * blue);
	return new int[]{redval,greenval,blueval};
}

}

