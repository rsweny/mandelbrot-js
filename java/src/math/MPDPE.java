package math;

/**
*  This class represents DPE numbers.
*  A DPE number has two components, a and n.
*  The double-value of a DPE number x is (x.a * 2 ^ x.n).
*  @author Herman Harjono.
*  @version Oct 5, 1998.
*/
final class MPDPE extends MPGlobal
{
  double a;
  int n;
  
  MPDPE() {}
  MPDPE(double A, int N) {a = A; n = N;}
  MPDPE(double A) {this(A,0);}
  
  /** 
  *  Returns the value of a * 2^n. 
  */
  double value() 
  {return a * Math.pow(2,n);}
  
  /** 
  * This converts the dpe number a to decimal form
  * i.e. b * 10^nb, where |b| is between 1 and 10.
  */
  static void dpdec(/*const*/ MPDPE a, MPDPE b)
  {
    final double xlt = 0.3010299956639812;
    
    if (a.a != 0.0) 
    {
      double t1 = xlt * a.n + log10(Math.abs(a.a));
      b.n = (int)t1; 
      if (t1 < 0.0) b.n -= 1;
      b.a = fSign (Math.pow(10.0, (t1 - b.n)), a.a);
    }
    else
    {
      b.a = 0.0;
      b.n = 0;
    }   
  }
}

