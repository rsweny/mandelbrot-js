package fractal;

//**************************************************************************
// Mandelbrot Fractor Generator
// Ryan Sweny
// rsweny@alumni.uwaterloo.ca
// 2002
// v2.2
//**************************************************************************

import java.awt.*;
import java.awt.event.*;
import math.MPGlobal;
import math.MPReal;
import math.Complex;
import java.awt.image.MemoryImageSource;

public class Mandelbrot extends java.applet.Applet implements Runnable, MouseListener, MouseMotionListener, KeyListener, ItemListener
{
	Label lblColor = new Label("Color Bands:");
	Label lblSat = new Label("Saturation Bands:");
	Label lblDepth = new Label("Depth:");
	Label lblZoom = new Label("Zoom:");
	Label lblX = new Label("X:");
	Label lblY = new Label("Y:");
	Label lblDiameter = new Label("D:");
	Label lblStatus = new Label("");
	Label lblInitColor = new Label("Initial Color:");
	
	Checkbox chkOrbitTrap = new Checkbox("Orbit Trap", null, false);
	Label lblOrbitType = new Label("Type:");
	Checkbox chkOrbitAverage = new Checkbox("Average Dist", null, false);
	
	TextField tOrbitType;
	TextField tOrbitX;
	TextField tOrbitY;
	TextField tOrbitDiameter;
	
	TextField tcolor;
	TextField tsat;
	
	TextField tdepth;
	
	TextField tzoom;
	TextField tX;
	TextField tY;
	
	TextField tinitColor;
	
	Panel p = new Panel(new BorderLayout());
	Panel pan1 = new Panel(new GridLayout(7,2));
	Panel pan2 = new Panel(new GridLayout(1,1));
	Panel pan3 = new Panel(new FlowLayout(FlowLayout.LEFT));
	
	MemoryImageSource screenMem;
	Graphics g_g;
	Graphics offscreenGraphics;
	Image offscreenImage;
	
	Thread[] runner;
	boolean running = false;

	double zoom_limit = 9.0E-14;
	
	double zoom = 1.7;
	double xcen = -0.5;
	double ycen = 0.01;
	
	boolean arbitraryPrecision = false;
	MPReal big_xcen = new MPReal("-0.5");
	MPReal big_ycen = new MPReal("0.01");
	MPReal big_zoom = new MPReal("1.7");

	MPReal two = new MPReal(2.0d);
	MPReal point5 = new MPReal(0.5d);

	float initialColor = -0.057f;
	float gradient = 0.45f;
	float glow = 1.7f;
	float band = 1;
	float spectrum = 0.3f;
	boolean colorBands = false;
	boolean smoothColor = false;

	int quad = 1;
	boolean mDown = false;
	long depth;
	int points = 0;

	double cx = 0;
	double cy = 0;
	float order = 3;
	float img_order = 3;
	float tol = 0.05f;

	int yimlen = 500;
	int ximlen = 500;

	int xcurr,ycurr,xanchor,yanchor;
	int alg = 9;
	
	double orbit_trap = 0.5;
	double aapoint = 0.5;
	double bbpoint = -0.25;
	
	boolean doOrbit = false;
	int orbitType = 0;
	boolean doOrbitAverage = false;
	
	long lastPassTime;
	long lastRunTime;
	double reset = -10;
	
	int[] pixels;
	int[][] img_red = null;
	int[][] img_green = null;
	int[][] img_blue = null;
	int[][] img_alpha = null;
	
	int lasti;
	
	public void start()
	{
		running = true;
		if (runner[0] == null);
		{
			runner[0] = new Thread(this);
			runner[0].start();
		}
	}


	public void stop()
	{
		running = false;
		lasti = 0;
		
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

		String szx = getParameter("x");
		if (szx != null)
		{
			xcen = Double.parseDouble(szx);
		}


		String szy = getParameter("y");
		if (szy != null)
		{
			ycen = Double.parseDouble(szy);
		}

		String szzoom = getParameter("zoom");
		if (szzoom != null)
		{
			zoom = Double.parseDouble(szzoom);
		}

		points = ximlen;
		setBackground(Color.black);
		setLayout(new FlowLayout(FlowLayout.LEFT));
		
		tcolor = new TextField(6);
		tsat = new TextField(6);
		tdepth = new TextField(8);
		tzoom = new TextField(53);
		tX = new TextField(53);
		tY = new TextField(53);
		tinitColor = new TextField(6);
		tOrbitType = new TextField(1);
		tOrbitType.setText("0");
		tOrbitX = new TextField(2);
		tOrbitY = new TextField(2);
		tOrbitDiameter = new TextField(2);

		tcolor.addKeyListener(this);
		tsat.addKeyListener(this);
		tdepth.addKeyListener(this);
		tzoom.addKeyListener(this);
		tX.addKeyListener(this);
		tY.addKeyListener(this);
		tinitColor.addKeyListener(this);
		tOrbitType.addKeyListener(this);
		tOrbitX.addKeyListener(this);
		tOrbitY.addKeyListener(this);
		tOrbitDiameter.addKeyListener(this);
		
		chkOrbitTrap.addItemListener(this);
		chkOrbitAverage.addItemListener(this);
		
		setConstants();

		p.setBackground(Color.WHITE);
		pan1.add(lblColor);
		pan1.add (tcolor);
		
		pan1.add(lblSat);
		pan1.add (tsat);
		
		pan1.add(lblDepth);
		pan1.add (tdepth);
		
		pan1.add(lblZoom);
		pan1.add (tzoom);
		
		pan1.add(lblX);
		pan1.add (tX);
		
		pan1.add(lblY);
		pan1.add (tY);
		
		pan1.add(lblInitColor);
		pan1.add(tinitColor);
		
		pan2.add(lblStatus);
		
		pan3.add(chkOrbitTrap);
		
		pan3.add(lblOrbitType);
		pan3.add(tOrbitType);
		
		pan3.add(lblX);
		pan3.add(tOrbitX);
		pan3.add(lblY);
		pan3.add(tOrbitY);
		
		pan3.add(lblDiameter);
		pan3.add(tOrbitDiameter);
		
		pan3.add(chkOrbitAverage);

		p.add(pan1, BorderLayout.NORTH);
		p.add(pan2, BorderLayout.SOUTH);
		
		pan3.setBackground(new Color(128,128,128));
		p.add(pan3, BorderLayout.CENTER);
		
		add(p);
		addMouseListener(this);
		addMouseMotionListener(this);
		
		clearScreenAndReset(0);
	}

	public void newFractal()
	{
		initialColor += 0.043;
		spectrum -= (0.03 + spectrum/10);
		if (spectrum < 0.1) {spectrum += 3.0f;}
		alg = (alg+1)%10;
		System.out.println("Algorithm: " + alg);
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
					if (reset > -10)
					{ 
						clearScreenAndReset(reset);
					}
					nextpoints(alg);
				}
				Thread.sleep(10);

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

	private void clearScreenAndReset(double newcx)
	{	
		lasti = 0;
		
		try { Thread.sleep(500); } catch(Exception e) {}
		offscreenImage = createImage(ximlen, yimlen);
		offscreenGraphics = offscreenImage.getGraphics();
		
		quad = 10;
		lasti = 0;
		cx = newcx;
		reset = -10;

		pixels = new int[ximlen * yimlen];
		img_red = new int[ximlen][yimlen];
		img_green = new int[ximlen][yimlen];
		img_blue = new int[ximlen][yimlen];
		img_alpha = new int[ximlen][yimlen];
		
		screenMem = new MemoryImageSource(ximlen, yimlen, pixels, 0, ximlen);
		g_g = getGraphics();
		
		lastPassTime = System.currentTimeMillis();
		lastRunTime = 0;
	}

	private double[] doCalcBig(MPReal xf, MPReal yf)
	{
		long iter = 0;
		MPReal aaaa, bbbb, aabb;
		
		MPReal atemp = two.multiply( xf.subtract(point5) );
		MPReal btemp = two.multiply( yf.subtract(point5) );

		MPReal at = atemp.multiply(big_zoom);
		MPReal bt = btemp.multiply(big_zoom);

		MPReal a = at.add(big_xcen);
		MPReal b = bt.add(big_ycen);

		MPReal large_aa = a;
		MPReal large_bb = b;

		iter = 0;
		if (cx != 0)
		{
			a = new MPReal(cx);
			b = new MPReal(cy);
		}

		do
		{
			aaaa = large_aa.multiply(large_aa);
			bbbb = large_bb.multiply(large_bb);
			aabb = large_aa.multiply(large_bb);
			large_aa = ((aaaa).subtract( bbbb )).add(a);
			large_bb = (aabb.add(aabb)).add(b);
			iter++;
		}
		while ( iter < depth && aaaa.add(bbbb).doubleValue() < 64 );

		
		
		double[] vals = new double[4];
		vals[0] = iter;
		vals[1] = large_aa.doubleValue();
		vals[2] = large_bb.doubleValue();
		vals[3] = 0;
		return vals;
	}

	private double[] doCalc(double xf,double yf)
	{
		long iter = 0;
		double aaold, bbold, a, b, aa, bb,aanew, bbnew;
		aaold = bbold = 0;

		a = (2*(xf - 0.5))*zoom + xcen;
		b = (2*(yf - 0.5))*zoom + ycen;

		aa = a;
		bb = b;

		iter = 0;
		if (cx != 0)
		{
			a = cx;
			b = cy;
		}

		double distance = 0;
		
		do
		{
			if (alg == 0)
			{
				aanew = (aa*aa*aa - 3*aa*bb*bb) + a;
				bbnew = (3*aa*aa*bb - bb*bb*bb) + b;
			}
			else if (alg == 1)
			{
				aanew = Math.tan(aa*aa - bb*bb) + a;
				bbnew = 2*aa*bb + b;
			}
			else if (alg == 2)
			{
				aanew = Math.cos(aa*aa - bb*bb) + a;
				bbnew = 2*aa*bb + b;
			}
			else if (alg == 3)
			{
				aanew = Math.sin(aa*aa - bb*bb + a);
				bbnew = 2*aa*bb + b;
			}
			else if (alg == 4)
			{
				aanew = Math.atan(aa*aa - bb*bb + a);
				bbnew = 2*aa*bb + b;
			}
			else if (alg == 5)
			{
				aanew = aa*aa - bb*bb + a + aaold;
				bbnew = 2*aa*bb + b + bbold;
				aaold = aa;
				bbold = bb;
			}
			else if (alg == 6)
			{
				aanew = aa*aa - bb*bb + a;
				bbnew = Math.sin(2*aa*bb) + b;
			}
			else if (alg == 7)
			{
				//z=2zc*cos(pi/z) 
				Complex Z = new Complex(aa,bb);
				Complex C = new Complex(a,b);
				Complex cTwo = new Complex(2,0);
				Complex cPi = new Complex(Math.PI,0);
				
				Complex t1 = cTwo.mul(Z).mul(C);
				Complex t2 = cPi.div(Z);
				
				Complex result = t1.mul( t2.cos() );
				aanew = result.re;
				bbnew = result.im;
			}
			else if (alg == 8)
			{
				aanew = aa*aa - bb*bb + a + aa;
				bbnew = 2*aa*bb + b + bb;
			}
			else
			{
				aanew = aa*aa - bb*bb + a;
				bbnew = 2*aa*bb + b;
				
				/*Complex Z = new Complex(aa,bb);
				Complex C = new Complex(a,b);
				Complex I = new Complex(iter,0);
				Complex cOne = new Complex(1,0);
				Complex result = Z.add( cOne.div(Z.pow(I)) );
				result = (result).mul(C);
				aanew = result.re;
				bbnew = result.im;*/
			}
			
			aa = aanew;
			bb = bbnew;
			iter++;
			
			
			if (doOrbit)
			{
				double dist;
				if (orbitType == 0)
				{
					dist = Math.abs( (bb + bbpoint) * (aa + aapoint) );
			  	}
				else if (orbitType == 1)
				{
					//pointy rect
					dist = Math.abs( (bb - bbpoint) )	;	
					double dist2 = Math.abs((aa - aapoint)*(bb - bbpoint));	
					if (dist2 > dist)			
			  			dist = dist2;
			  	}
			  	else if (orbitType == 2)
			  	{
					//hyperbola
					dist = Math.abs( (bb - bbpoint) * (aa - aapoint) );
				}
				else if (orbitType == 3)
			  	{
					//circle
					dist = distance(aapoint, bbpoint, aa, bb);
				}
				else
				{
					//spiral
					dist = Math.atan( Math.abs( (bb - bbpoint) / (aa - aapoint) ) );
				}
				
				if (dist < orbit_trap) 
				{
					if (doOrbitAverage)
					{
						if (distance == 0) distance = dist;
						else distance = dist*0.25+distance*0.75;
					}
					else
					{		
						if (distance == 0 || dist < distance) 
							distance = dist;
					}
				}
			}
		}
		while ( iter < depth && (aa*aa + bb*bb) < 64 );
		
		
		double[] vals = new double[4];
		vals[0] = iter;
		vals[1] = aa;
		vals[2] = bb;
		vals[3] = distance;
		return vals;
	}
	
	private void getConstants()
	{
		try
		{
			if (!arbitraryPrecision)
			{
				zoom = Double.parseDouble(tzoom.getText());
				xcen = Double.parseDouble(tX.getText());
				ycen = Double.parseDouble(tY.getText());
			}
			else
			{
				big_zoom = new MPReal(tzoom.getText());
				big_xcen = new MPReal(tX.getText());
				big_ycen = new MPReal(tY.getText());
			}
			depth = Long.parseLong(tdepth.getText());
			spectrum = Float.parseFloat(tcolor.getText());
			band = Float.parseFloat(tsat.getText());
			initialColor = Float.parseFloat(tinitColor.getText());
			aapoint = Double.parseDouble(tOrbitX.getText());
			bbpoint = Double.parseDouble(tOrbitY.getText());
			orbit_trap = Double.parseDouble(tOrbitDiameter.getText());
			
			orbitType = Integer.parseInt(tOrbitType.getText());
			doOrbitAverage = chkOrbitAverage.getState();
			doOrbit = chkOrbitTrap.getState();
		}
		catch(Exception e)
		{
		}
		checkPrecision();
	}
	
	private void setConstants()
	{
		try
		{
			tcolor.setText("" + spectrum);
			tsat.setText("" + band);
			tdepth.setText("" + depth);
			tinitColor.setText("" + initialColor);
			
			tOrbitX.setText("" + aapoint);
			tOrbitY.setText("" + bbpoint);
			tOrbitDiameter.setText("" + orbit_trap);
			
			if (!arbitraryPrecision)
			{
				tzoom.setText("" + zoom);
				tX.setText("" + xcen);
				tY.setText("" + ycen);
			}
			else
			{
				tzoom.setText(big_zoom.toString());
				tX.setText(big_xcen.toString());
				tY.setText(big_ycen.toString());
			}
				
		}
		catch(Exception e)
		{
		}
	}
	
	
	private double distance(double x1, double y1, double x2, double y2)
	{
		return Math.sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
	}


	private void nextpoints(int alg)
	{
		Thread thisThread = Thread.currentThread();
		double xf,yf;
		MPReal big_xf, big_yf;
		xf = yf = 0;

		float hue, hue2, numColors, huesat, new_gamma, new_sat;
		long iter;
		int xi,yi;
		
		//pass in the power of the function, for mandelbrot its 2 
		double il = 1.0/Math.log(alg == 0 ? 3 : 2); 

		int i;
		double thread_rnd = Math.random()*0.00001;
		for (i=lasti; i<lasti+points-quad; i+=quad)
		{
			double[] vals;
			if (i >= pixels.length) //if we're done the first pass, do random points to anti-alias
			{
				xf = Math.random();
				yf = Math.random();
				xi = (int)(ximlen*xf);
				yi = (int)(yimlen*yf);
			}
			else //do every pixel in order
			{
				xi = (i%ximlen);
				yi = (i / ximlen);
				xf = (double)xi / (double)ximlen + thread_rnd;
				yf = (double)yi / (double)yimlen + thread_rnd;
			}

			if (arbitraryPrecision)
			{
				big_xf = new MPReal(xf);
				big_yf = new MPReal(yf);
				vals = doCalcBig(big_xf, big_yf);
				iter = (long)vals[0];
				
				//exponential smoothing
				double cabsZ = Math.sqrt(vals[1]*vals[1] + vals[2]*vals[2]);
				hue = (float)(iter - il*Math.log( Math.log(cabsZ) )) / depth;
			}
			else
			{
				vals = doCalc(xf,yf);
				iter = (long)vals[0];
				
				//exponential smoothing
				double cabsZ = Math.sqrt(vals[1]*vals[1] + vals[2]*vals[2]);
				hue = (float)(iter - il*Math.log( Math.log(cabsZ) )) / depth;
				
				//binary decomposition
				//if (bb > 0) hue = hue/2;
			}
			
			//orbit trap
			if (doOrbit)
			{
				float dx = (float)(orbit_trap - vals[3] + 0.03);
				hue = hue / dx;
			}
			
			//angle (metallic effect)
			//double d = Math.atan2(aa, bb); //get angle of z
			//if (d < 0) d = -d ;
			//hue = (float)( hue*d / (Math.PI*2));
			
			
			//scale the color gradients logarithmically for deep zooms
			if (depth > 10000)
			{
				hue = (hue + 1)*10;
				hue = (float)Math.log(1.0 + Math.log(hue));
			}

			huesat = hue;
			hue = hue * 1600000 * spectrum;
			hue = ((int)hue % 1600000)/1600000.0f;
			hue2 = hue;

			if (spectrum < 0.7) {hue2 = hue * (1 / spectrum);}

			if (colorBands)
			{
				huesat = huesat * 1600000 * band;
				huesat = ((int)huesat % 1600000)/1600000.0f;
				new_sat = 1.0f - huesat;

				if (smoothColor)
				{			
					if (huesat > 0.5) huesat = 1 - huesat;
					huesat *= 2;
				}

				if (iter == depth)
					new_gamma = 0;
				else
					new_gamma = (float) Math.min(Math.pow(huesat, gradient), 1);

				if (smoothColor)
				{					
					if (new_gamma > 0.5) new_gamma = 1 - new_gamma;
					new_gamma *= 2;
				}
			}
			else
			{
				if (iter == depth)
					new_gamma = 0;
				else
					new_gamma = (float) Math.min(glow*Math.pow(hue2, gradient), 1);

				if (smoothColor)
				{
					if (new_gamma > 0.5) new_gamma = 1 - new_gamma;
					new_gamma *= 2;
				}
				
				
				new_sat = (spectrum < 1.3 && ((int)(spectrum*10)) % 2 == 0 ) ? new_gamma : 1-new_gamma;
				new_sat = (float)Math.pow(new_sat/band, 0.75);
				if (new_sat > 1.0f) new_sat = 1.0f;
				
				
				if (smoothColor)
				{
					if (new_sat > 0.5) new_sat = 1 - new_sat;
					new_sat *= 2;
				}
			}
			

			int[] rgb = MPGlobal.HSBtoRGB(hue + initialColor, new_sat, new_gamma);
			addPixel(xi, yi, rgb);
			
			//rough painting
			if (i < pixels.length)
			{
				if (thisThread == runner[0])
				{
					offscreenGraphics.setColor( new Color(pixels[yi*ximlen + xi]) );
				    offscreenGraphics.fillRect(xi,yi,quad,quad);
				    g_g.setColor( new Color(pixels[yi*ximlen + xi]) );
				    g_g.fillRect(xi,yi,quad,quad);
				}
			}
		}
		
		if (thisThread == runner[0])
		{
			if (i <= pixels.length - points*quad)
			{
				lasti += points*quad;
			}
			else if (quad > 1)
			{
				lasti = 0;
				quad = 1;
				System.out.println("starting remaining threads");
				for (int k = 1; k < runner.length; k++)
				{
					if (runner[k] == null)
					{
						runner[k] = new Thread(this);
						runner[k].start();
					}
				}
			}
			else //start of random points
			{
				lasti += points;
	
				if (lasti%pixels.length == 0)
				{
					long diff = (System.currentTimeMillis() - lastPassTime)/1000;
					
					String strStatus = "Pass: " + (lasti/pixels.length)*runner.length + " - " + diff + "s";
					System.out.println(strStatus);
					showStatus(strStatus);
					
					lastPassTime = System.currentTimeMillis();  	
				}
			}
		}
		
	}
	
	private synchronized void addPixel(int xi, int yi, int[] rgb)
	{
		img_red[xi][yi] += rgb[0];
		img_green[xi][yi] += rgb[1];
		img_blue[xi][yi] += rgb[2];
		img_alpha[xi][yi]++;
		
        int red = img_red[xi][yi]/img_alpha[xi][yi];
        int green = img_green[xi][yi]/img_alpha[xi][yi];
        int blue = img_blue[xi][yi]/img_alpha[xi][yi];

        red = red << 16;
		green = green << 8;
        pixels[yi*ximlen + xi] = 0xff000000 | (red & 0xffff0000) | (green & 0x0000ff00) | (blue & 0xff);
	}
	
	public void updateHistogram()
	{
		if (lasti > pixels.length) 
    		offscreenImage = createImage(screenMem);   	
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
	
		if (e.getButton() == MouseEvent.BUTTON3 && e.isControlDown())
		{
			smoothColor = !smoothColor;
			glow=1;
			reset = cx;
		}
		else if ( e.getClickCount() == 2 )
		{
			if (!arbitraryPrecision)
			{
				findCenter();
			}
			else
			{
				findBigCenter();
			}
			setConstants();
		}	
		else if (e.getButton() == MouseEvent.BUTTON3)
		{
			p.setVisible(!p.isVisible());
		}
		else if (e.getButton() == MouseEvent.BUTTON1 && e.isControlDown())
		{
			colorBands = !colorBands;
			glow=1;
			reset = cx;
		}
		else if (e.getButton() == MouseEvent.BUTTON2)
		{
			xcen = 0;
			ycen = 0.01;
			zoom = 3;
			reset = 0;
			newFractal();
		}
		else if (e.getButton() == MouseEvent.BUTTON1 && e.isShiftDown())
		{
			cx = (float)  ( ((float)e.getX()/(float)ximlen-0.5f)*2*zoom + (float)xcen  );
			cy = (float) -( ((float)e.getY()/(float)yimlen-0.5f)*2*zoom + (float)ycen  );
			reset = cx;
			
			xcen = 0;
			ycen = 0;
			zoom = 2;
		}
		
		setConstants();
		e.consume(); 	
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
		if (dx < 10)
		{
			return;
		}
	
		checkPrecision();
	
		if (!arbitraryPrecision)
		{
			double newxcen = (xanchor + xcurr)/2.0;
			newxcen = ((newxcen/ximlen)-0.5)*2*zoom + xcen;
			xcen = newxcen;
	
			double newycen = (yanchor + ycurr)/2.0;
			newycen = ((newycen/yimlen)-0.5)*2*zoom + ycen;
			ycen = newycen;
	
			zoom = (dx*zoom)/ximlen;
			System.out.println("Zoom is " + zoom);
			System.out.println("Xcen " + xcen);
			System.out.println("Ycen " + ycen);
		}
		else
		{
			double newxcen = (xanchor + xcurr)/2.0;
			double newycen = (yanchor + ycurr)/2.0;
	
			MPReal big_ximlen = new MPReal((double)ximlen);
			MPReal big_yimlen = new MPReal((double)yimlen);
			MPReal big_dx = new MPReal((double)dx);
	
			MPReal tmp_xanchor = new MPReal( newxcen );
			MPReal tmp_yanchor = new MPReal( newycen );
	
			MPReal xtemp = tmp_xanchor.divide(big_ximlen).subtract(point5);
			MPReal ytemp = tmp_yanchor.divide(big_yimlen).subtract(point5);
	
			big_xcen = (xtemp.multiply(two).multiply(big_zoom)).add(big_xcen);
			big_ycen = (ytemp.multiply(two).multiply(big_zoom)).add(big_ycen);
	
			big_zoom = (big_dx.divide( big_ximlen )).multiply(big_zoom);
	
			System.out.println("BigZoom is " + big_zoom.toString());
			System.out.println("Xcen " + big_xcen.toString());
			System.out.println("Ycen " + big_ycen.toString());
		}
	
		reset = cx;
		setConstants();
		System.gc();
		e.consume();
	} 	
	
	public void mouseDragged(MouseEvent e)
	{
		mDown = true;
		xcurr = e.getX();
		ycurr = e.getY();
		repaint();
		e.consume();
	}

	private void checkPrecision()
	{
		int exp = - big_zoom.getExponent();
		
		
		if (zoom < zoom_limit && !arbitraryPrecision)
		{
			arbitraryPrecision = true;
			big_zoom = new MPReal(zoom);
			big_xcen = new MPReal(xcen);
			big_ycen = new MPReal(ycen);
			System.out.println("mantissa is " + big_zoom.mantissa.length);
		}
		else if (exp >= 20)
		{
			//-38, 76:14
			arbitraryPrecision = true;
			int newprecision = exp+10;
			math.MPGlobal.setMaximumPrecision(newprecision);
			System.out.println("Precision now " + newprecision);
			System.out.println("mantissa is " + big_zoom.mantissa.length);
		}
		else if (exp <= 13)
		{
			arbitraryPrecision = false;
			zoom = Double.parseDouble(tzoom.getText());
			xcen = Double.parseDouble(tX.getText());
			ycen = Double.parseDouble(tY.getText());
		}
	}
	
	
   public void findBigCenter()
   {
   		MPReal searchDistance = new MPReal(0.01d);
   		MPReal xdistance = new MPReal(0.01d);
   		MPReal ydistance = new MPReal(0.01d);

		long search = 0;
		long old_iter = 0;
		
		while (search < 1000)
		{
			search++;
			MPReal newx = big_xcen.add( (xdistance.multiply(big_zoom)).multiply(two) );
			MPReal newy = big_ycen.add( (ydistance.multiply(big_zoom)).multiply(two) );
			
			if (search %100 == 1) 
			{
				lblStatus.setText(search + "th search: " + newx.toString() + "," + newy.toString() );
				big_zoom = big_zoom.multiply(two);
			}
			searchDistance = searchDistance.add(new MPReal(0.001d));
			
			
			if (newx.abs().compareTo(two) > 0 || newy.abs().compareTo(two) > 0)
	   		{
	   			lblStatus.setText("Out of bounds, stopping search.");
	   			return;
	   		}
			

			MPReal xf = point5.add(xdistance);
			MPReal yf = point5.add(ydistance);
			
			double[] vals = doCalcBig(xf,yf);
			long iter = (long)vals[0];

			if ( iter == depth )
			{
				big_xcen = newx;
				big_ycen = newy;
				
				reset = 0;

				lblStatus.setText("Finished: " + big_xcen.toString() + "," + big_ycen.toString());
				return;
			}
			else if ( iter > old_iter )
			{
				big_xcen = newx;
				big_ycen = newy;
				
				if (old_iter > 0) 
					lblStatus.setText("Improved: " + big_xcen.toString() + "," + big_ycen.toString());
   				
   				searchDistance = new MPReal(0.01d);
   				xdistance = new MPReal(0.01d);
   				ydistance = new MPReal(0.01d);
   				old_iter = iter;
			}
			
			if (!xdistance.isNegative() && !ydistance.isNegative())
			{
				xdistance = xdistance.add(searchDistance);
				ydistance = ydistance.add(searchDistance).negate();
			}
			else if (!xdistance.isNegative() && ydistance.isNegative())
			{
				xdistance =  xdistance.add(searchDistance).negate() ;
				ydistance = ydistance.subtract(searchDistance);
			}
			else if (xdistance.isNegative() && ydistance.isNegative())
			{
				xdistance = xdistance.subtract(searchDistance);
				ydistance =  ydistance.subtract(searchDistance).negate() ;
			}
			else //if (xdistance.isNegative() && !ydistance.isNegative())
			{
				xdistance =  xdistance.subtract(searchDistance).negate() ;
				ydistance = ydistance.add(searchDistance);
			}
		}
	}
   
   
   public void findCenter()
   {
   		double searchDistance = 0.01;
   		double xdistance = searchDistance;
   		double ydistance = searchDistance;

		long search = 0;
		long old_iter = 0;
		
		while (search < 10000)
		{
			search++;
			if (search %1003 == 1) lblStatus.setText(search + "th search, now at " + (xcen+ (xdistance*zoom*2)) + "," + (ycen+ (ydistance*zoom*2) ) );
			searchDistance += 0.001;
			
			if (Math.abs(xcen+ (xdistance*zoom*2)) > 2 || Math.abs(ycen+ (ydistance*zoom*2)) > 2)
	   		{
	   			lblStatus.setText("Out of bounds, stopping search.");
	   			return;
	   		}
   		
   			double[] vals = doCalc(0.5+xdistance,0.5+ydistance);
			long iter = (long)vals[0];

			if ( iter == depth )
			{
				xcen += (xdistance*zoom*2);
				ycen += (ydistance*zoom*2);
				
				reset = 0;
				
				lblStatus.setText("Done: " + xcen + "," + ycen);
				return;
			}
			else if ( iter > old_iter )
			{
				xcen += (xdistance*zoom*2);
				ycen += (ydistance*zoom*2);
				lblStatus.setText("Improved: " + xcen + "," + ycen);
				searchDistance = 0.01;
   				xdistance = searchDistance;
   				ydistance = searchDistance;
   				old_iter = iter;
			}
			
			if (xdistance >= 0 && ydistance >= 0)
			{
				xdistance = xdistance+searchDistance;
				ydistance = -(ydistance+searchDistance);
			}
			else if (xdistance >= 0 && ydistance < 0)
			{
				xdistance = -(xdistance+searchDistance);
				ydistance = ydistance-searchDistance;
			}
			else if (xdistance < 0 && ydistance < 0)
			{
				xdistance = xdistance-searchDistance;
				ydistance = -(ydistance-searchDistance);
			}
			else //if (xdistance < 0 && ydistance >= 0)
			{
				xdistance = -(xdistance-searchDistance);
				ydistance = ydistance+searchDistance;
			}
		}
	}
	
	public void itemStateChanged(ItemEvent e) 
	{
		reset = cx;
		getConstants();
	}
	
	public void keyPressed(KeyEvent e)
	{
	}
	
	public void keyReleased(KeyEvent e)
	{
		getConstants();
		checkPrecision();
		reset = cx;
	}
	
	public void keyTyped(KeyEvent e)
	{
	}
}