package fractal;

//**************************************************************************
// Newton Fractor Generator
// Ryan Sweny
// rsweny@alumni.uwaterloo.ca
// 2004
//**************************************************************************

import java.awt.*;
import java.util.*;
import math.Complex;
import math.MPGlobal;
import java.awt.event.*;
import java.awt.image.MemoryImageSource;


public class Newton extends java.applet.Applet implements Runnable, KeyListener, MouseListener, MouseMotionListener
{
	MemoryImageSource screenMem;
    Graphics offscreenGraphics;
    Graphics g_g;
    Image offscreenImage;
   
	Thread[] runner;
	boolean running = false;
    
    TextField cError;
    TextField c;
	TextField e;
	TextField tz;
	Panel p = new Panel(new GridLayout(9,2));
	TextField tX;
	TextField tY;
	TextField tContrast;
	TextField tBright;
	TextField tDepth;
  
    double zoom = 3;
    double xcen = 0.0;
    double ycen = 0.0;
    
    float brightness = 6.0f;
    float intensity = 0.35f;
    int alg = 1;
    float rootBoundry = 0.0000000000001f;
    float initColor = 0.6f;
    
    long runtotal = 0;
    float quad = 25;
    
 
    boolean mDown = false; 
    int depth;
    int points;
    
    float order = 3;
    float img_order = 0;
    
    int yimlen;
    int ximlen;
    
    int xcurr,ycurr,xanchor,yanchor;

    double p1 = Math.random();
    double p2 = Math.random();
    
    boolean mandelbrotAddition = false;
	long lastPassTime;
	long lastRunTime = 0;
	int reset = -1;
	boolean fixHoles = false;
	
	//constants
	Complex two = new Complex(2, 0);
	Complex one = new Complex(1, 0);
	Complex pointone = new Complex(0.1, 0);
	
	double complex_error = 0;
	Complex oneError = new Complex(1,complex_error);
	
	int[] pixels;
   	int[][] img_red = null;
	int[][] img_green = null;
	int[][] img_blue = null;
	int[][] img_alpha = null;
	LinkedList roots = new LinkedList();

	public void start()
	{
		running = true;
		for (int i = 0; i < runner.length; i++)
		{
			if (runner[i] == null);
			{
				runner[i] = new Thread(this);
				runner[i].start();
			}
		}
	}


	public void stop()
	{
		running = false;
		
		for (int i = 0; i < runner.length; i++)
			runner[i] = null;
			
		System.out.println("stop()");
	}
	
	public void destroy()
	{
		super.destroy();
		stop();
	}


    public void init()
    {
		Runtime runtime = Runtime.getRuntime();
        int nrOfProcessors = runtime.availableProcessors();
		runner = new Thread[nrOfProcessors];
		System.out.println(nrOfProcessors + " processors detected.");
		
	   	String szheight = getParameter("height");
	    if (szheight == null) 
	    {
	      szheight ="500";
	    }
	    yimlen = Integer.parseInt(szheight);
	    
	    String szwidth = getParameter("width");
	    if (szwidth == null) 
	    {
	      szwidth ="500";
	    }
    	ximlen = Integer.parseInt(szwidth);
    	
    	String szdepth = getParameter("depth");
	    if (szdepth == null) 
	    {
	      szdepth ="1000";
	    }
    	depth = Integer.parseInt(szdepth);
    	
    	points = ximlen;
    	
		setBackground(Color.black);
		setLayout(new FlowLayout());
		
		cError = new TextField(5);
		c = new TextField(5);
		e = new TextField(5);
		tz = new TextField(3);
		tContrast = new TextField(4);
		tBright = new TextField(4);
		tDepth = new TextField(4);
		tX = new TextField(8);
		tY = new TextField(8);
		
		tX.addKeyListener(this);
		tY.addKeyListener(this);
		
		cError.addKeyListener(this);
		c.addKeyListener(this);
		e.addKeyListener(this);
		tz.addKeyListener(this);
		tContrast.addKeyListener(this);
		tBright.addKeyListener(this);
		tDepth.addKeyListener(this);
		
		Label complexError = new Label("Complex Error:");
		complexError.setForeground(Color.WHITE);
		
		Label complex = new Label("Complex Exponent:");
		complex.setForeground(Color.WHITE);
		Label real = new Label("Real Exponent:");
		real.setForeground(Color.WHITE);
		Label lzoom = new Label("Zoom:");
		lzoom.setForeground(Color.WHITE);
		Label lblX = new Label("X:");
		lblX.setForeground(Color.WHITE);
		Label lblY = new Label("Y:");
		lblY.setForeground(Color.WHITE);
		Label lblContrast = new Label("Contrast:");
		lblContrast.setForeground(Color.WHITE);
		Label lblBright = new Label("Brightness:");
		lblBright.setForeground(Color.WHITE);
		Label lblDepth = new Label("Depth:");
		lblDepth.setForeground(Color.WHITE);

		cError.setText("0");
		c.setText("0");
		e.setText("3");
		tz.setText("5");

		p.add(lzoom);
		p.add(tz);
		p.add(real);
		p.add(e);
		p.add(complex);
		p.add(c);
		p.add(complexError);
		p.add(cError);
		
		p.add(lblX);
		p.add(tX);
		p.add(lblY);
		p.add(tY);
		
		p.add(lblContrast);
		p.add(tContrast);
		
		p.add(lblBright);
		p.add(tBright);
		
		p.add(lblDepth);
		p.add(tDepth);
		
		add(p);
		addMouseListener(this);
		addMouseMotionListener(this);
		
		setConstants();
		getConstants();
		clearScreenAndReset(1);
    }
    
    public void newFractal()
    {
    	roots = new LinkedList();
    	p1 = Math.random();
		runtotal = 0;
		quad = 5;
    }
    

	public void run()
	{
		Thread thisThread = Thread.currentThread();
		while ( running )
		{
			try 
			{
				if (!mDown) 
				{
					if (reset > -1)
					{ 
						clearScreenAndReset(reset);
					}
					offscreenGraphics.setPaintMode();
					nextpoints(alg);
		        }
		        Thread.sleep(100);
		        
				long diff = (System.currentTimeMillis() - lastRunTime);
				if (diff > 5000) 
				{
					lastRunTime = System.currentTimeMillis();
					updateHistogram();
					repaint();
				}
			} 
			catch (Exception e) 
			{ 
				e.printStackTrace();
			}
		}
		System.out.println("run() exited: " + thisThread.toString());
	}
	
	public synchronized void updateHistogram()
	{
		//System.out.println("quad is " + quad);
		if (fixHoles)
		{
			fixHoles = false;
			int cnt = 0;
			System.out.println("Hitting remaining pixels...");
			for (int i=0; i<ximlen; i++) 
	        {
	            for (int j=0; j<yimlen; j++) 
	            {
	            	if (img_alpha[i][j] == 0)
	            	{
			            int xi = (i%ximlen);
						double xf = (double)xi / (double)ximlen;
						double yf = (double)j / (double)yimlen;
						doCalculation(xf,yf,xi,j);
						cnt++;
					}
	            }
	        }
	        
	        System.out.println("Hitting remaining done: " + cnt);
		}
		else if (runtotal > pixels.length*2) 
		{
    		offscreenImage = createImage(screenMem);
    	} 
    		
	}


    public void update(Graphics g) 
    {
        paint(g);		
    }
   
    
	public void paint(Graphics g)
	{
		g.drawImage(offscreenImage, 0, 0, this);
    	if (mDown) 
        {
			g.setColor(Color.white);
			g.drawRect(xanchor,yanchor,xcurr-xanchor,ycurr-yanchor);
        } 
	}
   
    public void clearScreenAndReset(int r) 
    {
    	offscreenImage = createImage(ximlen, yimlen);
		offscreenGraphics = offscreenImage.getGraphics();
		oneError = new Complex(1,complex_error);
		
    	if (r == 1) newFractal();
    	
    	runtotal = 0;
        quad = 10;
        reset = -1;
        
    	img_red = new int[ximlen][yimlen];
		img_green = new int[ximlen][yimlen];
		img_blue = new int[ximlen][yimlen];
		img_alpha = new int[ximlen][yimlen];
		pixels = new int[ximlen*yimlen];
		screenMem = new MemoryImageSource(ximlen, yimlen, pixels, 0, ximlen);
		g_g = getGraphics();
		
		lastPassTime = System.currentTimeMillis();
		lastRunTime = 0;
    }
    
    private void doCalculation(double xf, double yf, int xi, int yi)
    {
    	double a,b,aa,bb,aanew,bbnew,mag2;
		float hue, huesat, new_gamma;
		double[] iter;
		
		a = 2*(xf - 0.5)*zoom + xcen;
		b = 2*(yf - 0.5)*zoom + ycen;

		aa = a;
		bb = b;
		iter = NewtonIterate(aa, bb, alg);

		if (mandelbrotAddition)
		{
			//get Color based on angle of z	
						
			//double angle = (float)Math.atan(bb/aa);
			double angle = (float)Math.atan(iter[2]/iter[3]);
			angle = angle + (Math.PI / 2);
			angle = angle / Math.PI;
			hue = (float)angle;
		}
		else
		{
			//limit the colors to 30% of the spectrum, unless a complex root
			float limitfactor = 0.3f;
			if (img_order != 0) limitfactor = 1.0f;	
			hue = initColor + (float)(iter[1]*limitfactor)/roots.size(); 
		}

		//normal shading
		huesat = (float)iter[0]/(depth/brightness);
		if (huesat >= 1.0f) huesat = 0.9999f;

		huesat = (iter[0] == depth) ? 1 : huesat * 16000000;
		if (huesat > 16000000) { huesat = Math.abs(16000000 + (16000000-huesat)); }
		huesat = ((int)huesat % 16000000)/16000000.0f;

		new_gamma = (float) Math.pow(huesat, intensity);

		int[] rgb = MPGlobal.HSBtoRGB(hue, 1.0f - huesat, new_gamma);
		addPixel(xi, yi, rgb);
		
        //rough painting
		if (runtotal < pixels.length*2)
		{
			offscreenGraphics.setColor(new Color( pixels[yi*ximlen + xi] ));
		    offscreenGraphics.fillRect(xi,yi,(int)quad,(int)quad);
		    g_g.setColor(new Color( pixels[yi*ximlen + xi] ));
		    g_g.fillRect(xi,yi,(int)quad,(int)quad);
		}
	}
	
	private synchronized void addPixel(int xi, int yi, int[] rgb)
	{
		int red, green, blue;
		if ( quad < 5 )
		{
			//once roots are calculated, start recording pixel data
			img_red[xi][yi] += rgb[0];
			img_green[xi][yi] += rgb[1];
			img_blue[xi][yi] += rgb[2];
			img_alpha[xi][yi]++;
			
			red = img_red[xi][yi]/img_alpha[xi][yi];
        	green = img_green[xi][yi]/img_alpha[xi][yi];
        	blue = img_blue[xi][yi]/img_alpha[xi][yi];
		}
		else
		{
			//rough colors
			red = rgb[0];
        	green = rgb[1];
        	blue = rgb[2];
		}
		
        red = red << 16;
		green = green << 8;
        pixels[yi*ximlen + xi] = 0xff000000 | (red & 0xffff0000) | (green & 0x0000ff00) | (blue & 0xff);
	}


	private void nextpoints(int alg)
	{
		for (int i=0; i<points; i++)
		{
			double xf = Math.random();
			double yf = Math.random();
			int xi = (int)Math.floor(ximlen*xf);
			int yi = (int)Math.floor(yimlen*yf);
			doCalculation(xf,yf,xi,yi);
		}
	
		runtotal += points;
		if (quad < 1.001)
		{
			quad = 1;
		}
		else
		{
			quad = (float)(50*Math.sqrt((float)points/runtotal));	
		}
		
		if (runtotal%pixels.length == 0)
		{
			if (runtotal == pixels.length*2) fixHoles = true;
			
			long diff = (System.currentTimeMillis() - lastPassTime) / 1000;
			
			String strStatus = "Pass: " + (runtotal / pixels.length) + " - " + diff + "s";
			System.out.println(strStatus);
			showStatus(strStatus);
			
			lastPassTime = System.currentTimeMillis();
		}
		
	}
	
	public void mouseClicked(MouseEvent e) {
	}
	
	public void mouseEntered ( MouseEvent e ) {
	}
	
	public void mouseExited ( MouseEvent e ) {
	}
	
	public void mouseMoved(MouseEvent e) {
	}

	public void mousePressed( MouseEvent e )
	{
	 	xanchor = e.getX();
		yanchor = e.getY();
	
		if (e.getButton() == MouseEvent.BUTTON1 && e.isControlDown())
		{
			alg = (alg+1)%7;
			System.out.println("alg is " + alg);
			reset = 1;
		}
		else if (e.getButton() == MouseEvent.BUTTON1 && e.isShiftDown())
		{
			p2 = Math.random();
			if (img_order == 0)
			{
				img_order = ((float)p2 * order);
			}
			else
			{
				img_order = 0;
			}
			
			reset = 1;
		}
		else if (e.getButton() == MouseEvent.BUTTON3 && e.isShiftDown())
		{
			mandelbrotAddition = !mandelbrotAddition;
			reset = 1;
		}
		else if (e.getButton() == MouseEvent.BUTTON3)
		{
			p.setVisible(!p.isVisible());
		}
		else if (e.getButton() == MouseEvent.BUTTON2)
		{
			initColor = (float)Math.random();
		    p2 = Math.random();
		    	
	   		order = 2 + (float)Math.random()*11;
	   		if (p2 > 0.5) order = (float)Math.floor(order);
	   		
	   		xcen = 0.0;
	   		ycen = 0.0;
	   		zoom = 3;
	   		
	   		reset = 1;
		}
	}
	
	
	public void mouseDragged(MouseEvent e)
	{
		mDown = true;
		xcurr = e.getX();
		ycurr = e.getY();
		repaint();
	}
	
	public void mouseReleased ( MouseEvent e ) 
	{
		xcurr = e.getX();
		ycurr = e.getY();
	    mDown = false;
	    repaint();
	    
	   	int dx = Math.abs(xcurr - xanchor);
	    int dy = Math.abs(ycurr - yanchor);
	    if (dy > dx)  dx = dy;
	   	if (dx > 10) 
	    { 
	       	double newxcen = (xanchor + xcurr)/2.0;
	        newxcen = (newxcen/ximlen-0.5)*2*zoom + xcen;
	        xcen = newxcen;
	        double newycen = (yanchor + ycurr)/2.0;
	        newycen = (newycen/yimlen-0.5)*2*zoom + ycen;
	        ycen = newycen;
	        zoom = (double)dx/ximlen*zoom;
	        
	        reset = 0;
	    }
	    
	    setConstants();
		System.gc();
	}
    
    public double addRoot(Complex root)
    {
    	ListIterator i = (ListIterator)roots.iterator();
    	double colorFactor = 10/rootBoundry;
    	while (i.hasNext())
    	{
    		int index = i.nextIndex();
    		Complex c = (Complex)i.next();
    		if ( distance(root, c) < rootBoundry )
    		{
    			double retVal = index + colorFactor*distance(root, c);
				return Math.min(roots.size() - 1, retVal);
    		}
    	}
 
		if (roots.size() < 20 && quad > 4)
		{
			roots.add(root);
			System.out.println("add root:" + roots.size());
		}
		return roots.size() - 1;
    }
    
    public double distance(Complex a, Complex b) 
    {
		return Math.sqrt((a.re-b.re)*(a.re-b.re) + (a.im-b.im)*(a.im-b.im));
	}
	
	//z^n - 1 / z
	public Complex DivZFunction(Complex z, double i, double j)
	{
		Complex exponent = new Complex(order, img_order);
		return (z.pow(exponent)).sub( one.div(z) );
	}

	//n*z^(n-1) + z^(-2)
	public Complex DivZDerivative(Complex z, double i, double j)
	{
		Complex minusTwo = new Complex(-2,0);
		Complex exponent = new Complex(order, img_order);
		Complex exponentLessOne = exponent.sub(oneError);
		return (exponent.mul( z.pow(exponentLessOne) )).add ( z.pow(minusTwo) );
	}
	
	//z^10 + 0.2 i * z^5 - 1.
	public Complex Poly2Function(Complex z, double i, double j) 
	{
		Complex exponent = new Complex(order, img_order);
		Complex point2i = new Complex(0, 0.2);
		Complex ten = new Complex(10, 0);
		return (z.pow(ten)).add (point2i.mul(z.pow(exponent))).sub ( one );
	}
	
	//10z^9 + 0.2i*5*z^4
	public Complex Poly2Derivative(Complex z, double i, double j) 
	{
		Complex exponent = new Complex(order, img_order);
		Complex exponentLessOne = exponent.sub(oneError);
		Complex ten = new Complex(10, 0);
		Complex nine = new Complex(9, 0);
		Complex point2i = new Complex(0, 0.2);

		return (ten.mul(z.pow(nine))).add(  point2i.mul(exponent).mul(z.pow(exponentLessOne)) );
	}
	

	//2z^3 - c + 1
	public Complex PolyMFunction(Complex z, double i, double j) 
	{
		Complex c = new Complex(i, j);
		Complex exponent = new Complex(order, img_order);
		return two.mul( z.pow(exponent) ).sub(c).add(one);
	}
	
	//6z^2 - 1
	public Complex PolyMDerivative(Complex z, double i, double j) 
	{
		Complex exponent = new Complex(order, img_order);
		Complex exponentLessOne = exponent.sub(oneError);
		return two.mul(exponent).mul( z.pow(exponentLessOne) ).sub(one);
	}
	
	
	
	//z^c - z + 0.1
	public Complex Poly3Function(Complex z, double i, double j) 
	{
		Complex c = new Complex(1+order, img_order);
		return z.pow(c).sub(z).add(pointone);
	}
	
	//c*z^(c-1) - 1
	public Complex Poly3Derivative(Complex z, double i, double j) 
	{
		Complex cminus1 = new Complex(1+order-1, (img_order == 0) ? 0 : img_order-1);
		Complex c = new Complex(1+order, img_order);
		return c.mul(z.pow(cminus1)).sub(one);
	}
	
	
	//z^n - 3z^5 + 6z^3 - 3z + 3
	public Complex PolyFunction(Complex z, double i, double j) 
	{
		Complex exponent = new Complex(order, img_order);
		Complex three = new Complex(3, 0);
		Complex six = new Complex(6, 0);
		Complex five = new Complex(5, 0);
		
		return (z.pow(exponent)).sub (three.mul(z.pow(five))).add ( six.mul(z.pow(three)) ).sub (three.mul(z)).add(three);
	}
	
	public Complex PolyDerivative(Complex z, double i, double j) 
	{
		Complex exponent = new Complex(order, img_order);
		Complex exponentLessOne = exponent.sub(oneError);
		Complex three = new Complex(3, 0);
		Complex eighteen = new Complex(18, 0);
		Complex fifteen = new Complex(15, 0);
		Complex four = new Complex(4, 0);
		return (exponent.mul(z.pow(exponentLessOne))).sub (fifteen.mul(z.pow(four))).add ( eighteen.mul(z.pow(two)) ).sub(three);
	}


	// z^n - 1 = 0
	public Complex UnityFunction(Complex z, double i, double j) 
	{
		Complex exponent = new Complex(order, img_order);
		return z.pow(exponent).sub( one );
	}
	
	public Complex UnityDerivative(Complex z, double i, double j) 
	{
		Complex exponent = new Complex(order, img_order);
		Complex exponentLessOne = exponent.sub(oneError);
		return exponent.mul( z.pow(exponentLessOne) );
	}
	
	
	//z^z - cz
	public Complex ZZFunction(Complex z, double i, double j) 
	{
		Complex con = new Complex(order, img_order);
		return (z.pow(z)).sub( con.mul(z) );
	}
	
	//z^z * (1 + lnz) - c
	public Complex ZZDerivative(Complex z, double i, double j) 
	{
		Complex con = new Complex(order, img_order);
		Complex lnz = z.log().add(one);
		return ((z.pow(z)).mul( lnz )).sub(con);
	}
	
	public Complex Newt(Complex z, double i, double j, int mode) 
	{
		if (mode == 1)
		{
			return z.sub( UnityFunction(z,i,j).div(UnityDerivative(z,i,j)) );
		}
		else if (mode == 2)
		{
			return z.sub( Poly2Function(z,i,j).div(Poly2Derivative(z,i,j)) );
		}
		else if (mode == 3)
		{
			return z.sub( DivZFunction(z,i,j).div(DivZDerivative(z,i,j)) );
		}
		else if (mode == 4)
		{
			return z.sub( Poly3Function(z,i,j).div(Poly3Derivative(z,i,j)) );
		}
		else if (mode == 5)
		{
			return z.sub( PolyMFunction(z,i,j).div(PolyMDerivative(z,i,j)) );
		}
		else if (mode == 6)
		{
			return z.sub( ZZFunction(z,i,j).div(ZZDerivative(z,i,j)) );
		}
		else
		{
			return z.sub( PolyFunction(z,i,j).div(PolyDerivative(z,i,j)) );
		}
	}
	
	public double[] NewtonIterate(double i, double j, int mode) 
	{
		int n = 0;
		Complex z = new Complex(i,j);
		Complex old = new Complex(i, j);
		Complex older = new Complex(i, j);
		
		z = Newt(z, i, j, mode);
		
		double hue = 0;
		double w = 0;
		while(n < depth && distance(old,z) > rootBoundry) 
		{
			older = old;
			old = z;
			
			z = Newt(z, i, j, mode);
			
			if (mandelbrotAddition)
				z = z.add(new Complex(i/2,j/2));
				
			n++;
			
			//normal smoothing
			w = 1.0 / z.sub(old).abs();
			hue += Math.pow(1.05, -w);
		}
		
		double[] vals = new double[4];
		if (n != depth) 
		{
			vals[0] = hue;
			vals[1] = addRoot(z);
		}
		else
		{
			vals[0] = n;
			vals[1] = 0;
		}
		
		vals[2] = z.re;
		vals[3] = z.im;
		
		return vals;
	}
	
	
	private void getConstants()
	{
		try
		{
			xcen = Double.parseDouble(tX.getText());
			ycen = Double.parseDouble(tY.getText());
			intensity = Float.parseFloat(tContrast.getText());
			brightness = Float.parseFloat(tBright.getText());
			zoom = Double.parseDouble(tz.getText());
			depth = Integer.parseInt(tDepth.getText());
			
			img_order = Float.parseFloat(c.getText());
			order = Float.parseFloat(e.getText());
			complex_error = Double.parseDouble(cError.getText());
		}
		catch(Exception e)
		{
		}
	}
	
	private void setConstants()
	{
		try
		{
			tContrast.setText("" + intensity);
			tBright.setText("" + brightness);
			tDepth.setText("" + depth);
			tX.setText("" + xcen);
			tY.setText("" + ycen);
			tz.setText("" + zoom);
			c.setText(img_order + ""); 
			e.setText(order + ""); 
			cError.setText(complex_error+"");
		}
		catch(Exception e)
		{
		}
	}
	
	public void keyPressed(KeyEvent e)
	{
	}
	
	public void keyReleased(KeyEvent ev)
	{
		getConstants();
		
		if (ev.getSource() == c || ev.getSource() == e || ev.getSource() == cError)
			reset = 1;
		else
			reset = 0;
	}
	
	public void keyTyped(KeyEvent e)
	{
	}
}







