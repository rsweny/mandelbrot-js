package fractal;

//**************************************************************************
// Mandelbrot Fractor Generator
// Ryan Sweny
// rsweny@alumni.uwaterloo.ca
// 2002
// v2.1
//**************************************************************************

import java.awt.*;
import java.awt.event.*;
import java.awt.image.MemoryImageSource;
import math.Complex;

public class Taylor extends java.applet.Applet implements Runnable, MouseListener, MouseMotionListener, KeyListener, ItemListener
{
	Label lblDepth = new Label("Depth:");
	Label lblZoom = new Label("Zoom:");
	Label lblX = new Label("X:");
	Label lblY = new Label("Y:");
	Label lblGlow = new Label("Glow:");
	Label lblBrightness = new Label("Brightness:");
	Label lblComplex = new Label("Complex Scale:");

	TextField tdepth;
	
	TextField tzoom;
	TextField tX;
	TextField tY;
	TextField tGlow;
	TextField tBrightness;
	TextField tComplex;
	
	Panel p = new Panel(new BorderLayout());
	Panel pan1 = new Panel(new GridLayout(7,2));

	
	int lasti;
	MemoryImageSource screenMem;
	Graphics g_g;
	Graphics offscreenGraphics;
	Image offscreenImage;
	Thread runner;
	
	double zoom = 9;
	double xcen = 0.0;
	double ycen = 0.01;

	int quad = 1;
	boolean mDown = false;
	long depth;
	int points = 0;


	boolean scalePoint;
	boolean compoundPoint;
	boolean recurse;
	boolean mandelbrot;
	double complexScale = 0;
	float cx = 0.0f;
	float cy = 0.0f;
	double glow = 0.01;
	
	double brightness = 1.0;
	double brightness_factor = 1.0;

	int yimlen = 500;
	int ximlen = 500;

	int[][] img_red = null;
	int[][] img_green = null;
	int[][] img_blue = null;
	int[][] img_alpha = null;
	int[] pixels;

	int xcurr,ycurr,xanchor,yanchor;
	int alg = 0;
	
	long lastPassTime;
	long lastRunTime;
	double reset = -1;
	
	public void start()
	{
		if (runner == null);
		{
			runner = new Thread(this);
			runner.start();
		}
	}

	public void stop()
	{
		runner = null;
		System.out.println("stop()");
	}
	
	public void destroy()
	{
		super.destroy();
		stop();
	}


	public void init()
	{
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
			szdepth ="20";
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
		
		//xarr = new int[points];
		//yarr = new int[points];

		offscreenImage = createImage(ximlen, yimlen);
		offscreenGraphics = offscreenImage.getGraphics();

		setBackground(Color.black);
		setLayout(new FlowLayout(FlowLayout.LEFT));
		

		tdepth = new TextField(8);
		
		tzoom = new TextField(20);
		tX = new TextField(9);
		tY = new TextField(9);
		tGlow = new TextField(9);
		tBrightness = new TextField(9);
		tComplex = new TextField(9);
		
	
		tdepth.addKeyListener(this);
		tzoom.addKeyListener(this);
		tX.addKeyListener(this);
		tY.addKeyListener(this);
		tGlow.addKeyListener(this);
		tBrightness.addKeyListener(this);
		tComplex.addKeyListener(this);
		

	
		setConstants();

		p.setBackground(Color.WHITE);

		
		pan1.add(lblDepth);
		pan1.add (tdepth);
		
		pan1.add(lblZoom);
		pan1.add (tzoom);
		
		pan1.add(lblX);
		pan1.add (tX);
		
		pan1.add(lblY);
		pan1.add (tY);
		
		pan1.add(lblGlow);
		pan1.add (tGlow);
		
		pan1.add(lblBrightness);
		pan1.add (tBrightness);
		
		pan1.add(lblComplex);
		pan1.add (tComplex);
	
	
		p.add(pan1, BorderLayout.NORTH);

		
		add(p);
		
		addMouseListener(this);
		addMouseMotionListener(this);
		
		clearScreenAndReset();
	}

	public void newFractal()
	{
		alg = (alg +1)%3;	
	}


	public void run()
	{
		Thread thisThread = Thread.currentThread();
		while ( runner == thisThread )
		{
			try
			{
				if (!mDown) 
				{
					if (reset > -1)
					{ 
						clearScreenAndReset();
					}
					nextpoints(alg);
				}
				
				long diff = (System.currentTimeMillis() - lastRunTime);
				if (diff > 5000) 
				{
					lastRunTime = System.currentTimeMillis();
					updateHistogram();
					repaint();
					Thread.sleep(100);
				}
				
			}
			catch (Exception e)
			{
				e.printStackTrace();
			}
		}
		System.out.println("run() exited");
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

	private void clearScreenAndReset()
	{
		offscreenImage = createImage(ximlen, yimlen);
		offscreenGraphics = offscreenImage.getGraphics();
		brightness_factor = Math.pow(depth,brightness);

		quad = 10;	
		lasti = 0;
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



	
	private void getConstants()
	{
		try
		{

			zoom = Double.parseDouble(tzoom.getText());
			xcen = Double.parseDouble(tX.getText());
			ycen = Double.parseDouble(tY.getText());
			glow = Double.parseDouble(tGlow.getText());
			brightness = Double.parseDouble(tBrightness.getText());
			
			brightness_factor = Math.pow(depth,brightness);
			
			complexScale = Double.parseDouble(tComplex.getText());

			depth = Long.parseLong(tdepth.getText());
		}
		catch(Exception e)
		{
		}
		
	}
	
	private void setConstants()
	{
		try
		{
			tdepth.setText("" + depth);
		
			tzoom.setText("" + zoom);
			tX.setText("" + xcen);
			tY.setText("" + ycen);
			tGlow.setText("" + glow);
			tBrightness.setText("" + brightness);
			tComplex.setText("" + complexScale);
		}
		catch(Exception e)
		{
		}
	}
	
	public synchronized void updateHistogram()
	{
		if (lasti > pixels.length) 
    		offscreenImage = createImage(screenMem);   	
	}
	
	private synchronized void nextpoints(int alg)
	{
		double xf,yf;
		int xi,yi,i;
		
		Complex one = new Complex(1,0);
		Complex two = new Complex(2,0);
	
		for (i=lasti; i<lasti+points-quad; i+=quad)
		{
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
				xf = (double)xi / (double)ximlen;
				yf = (double)yi / (double)yimlen;
			}
	
			//based on the algorithm published here: http://fractorama.com/doc/taylor.html
			//////////////////////////////////////
			double a = 2*(xf - 0.5)*zoom + xcen;
			double b = 2*(yf - 0.5)*zoom + ycen;
	
			//the current pixel we are calculating
			Complex current = new Complex(a,b);
			
			
			//the real value for this pixel
			Complex zvalue;
			if (alg == 0)
			{
				zvalue = Complex.pow(current,Math.E);
			}
			else if (alg == 1)
			{
				zvalue = ( one.add(current).div(one.sub(current)) ).log();
			}
			else
			{
				zvalue = current.sinh();
			}
			
			double xavg = 0;
			double yavg = 0;
			double mavg = 0;
	
			double xmax = 0;
			double ymax = 0;
			double mmax = 0;
	
			double xtot = 0;
			double ytot = 0;
			double mtot = 0;
	
			
			double den = 1;
			double den2 = 1;
			Complex denominator = new Complex(den);
			Complex denominator2 = new Complex(den);
			
			//expand at f(x - point)
			Complex point = Complex.pow(current,Math.E).add( new Complex(cx,cy) );
			
			Complex value = new Complex(0,0);
			Complex zSum = new Complex(0,0);
			
			// Loop until we've been through the loop 'depth' times
			long count = 0;
			while(count < depth)
			{
				Complex scale = new Complex(1,complexScale);
				Complex oldcurrent = current;
				if (scalePoint) //multiply by f(a)
				{
					if (alg == 0)
					{
						scale = Complex.pow(point,Math.E);	
					}
					else if (alg == 1)
					{
						scale = ( one.add(point).div(one.sub(point)) ).log();
					}
					else
					{
						scale = point.sinh();
					}
				}
				
				//f(x-a): we will approximate the series about a
				if (compoundPoint)
				{
					current = current.sub(point);
				}
				
				
				//factorial denominators
				if(count > 0) 
				{ 
					den *= count; 
					den2*= (count*2 + 1);
					denominator = new Complex(den);
					denominator2 = new Complex(den2);
				}
						

				if (alg == 0)
				{	// f(x) = x^n / n!
					Complex numerator = scale.mul( Complex.pow(current,count) );
					zSum = zSum.add( numerator.div(denominator) );
				}
				else if (alg == 1)
				{	// f(x) = 2x^(2n-1) / (2n-1)
					Complex twoN1 = new Complex( 2*count-1 );
					Complex numerator = scale.mul( two.mul(Complex.pow(current,twoN1)) );
					zSum = zSum.add( numerator.div(twoN1) );
				}
				else
				{	// f(x) = x^(2n+1) / (2n+1)!
					Complex numerator = scale.mul( Complex.pow(current,den2) );
					zSum = zSum.add( numerator.div(denominator2) );
				}
				
				if (mandelbrot)
					zSum = zSum.mul( zvalue.sub(zSum) );
				
				if (!recurse)
					current = oldcurrent;
				

				value = zSum.cot();
				//Complex cotz = zvalue.cot();
			
				//System.out.println("val re " + value.re + " val im " + value.im);
				if (Double.isNaN(value.re) || Double.isNaN(value.im))
				{
					//count = depth;
				}
				else
				{
					double x = Complex.pow(value, 0.1).abs();
					double y = Complex.pow(value, 0.5).abs();
					double m = Complex.pow(value, 0.9).abs();
					
					//double rx = Complex.pow(cotz, 0.1).abs();
					//double ry = Complex.pow(cotz, 0.5).abs();
					//double rm = Complex.pow(cotz, 0.9).abs();
		
					if (xavg == 0)
					{
						xavg = x;
						yavg = y;
						mavg = m;
					}

					double xdev = Math.abs(x - xavg);
					double ydev = Math.abs(y - yavg);
					double mdev = Math.abs(m - mavg);
					
					//double xdev = Math.abs(x - rx);
					//double ydev = Math.abs(y - ry);
					//double mdev = Math.abs(m - rm);
					
					xavg = (x*0.5 + xavg*(0.5 + 0.01*glow));
					yavg = (y*0.5 + yavg*(0.5 + 0.1*glow));
					mavg = (m*0.5 + mavg*(0.5 + 1*glow));
		
					xmax = Math.max(xmax, xdev);
					ymax = Math.max(ymax, ydev);
					mmax = Math.max(mmax, mdev);
		
					xtot += xdev;
					ytot += ydev;
					mtot += mdev;
				}
				count++;
			}
	
	
			int red = 0;
			int green = 0;
			int blue = 0;
	
			if(count > 0)
			{
				red = (int)(((xtot*255)/brightness_factor) / xmax);
				green = (int)(((ytot*255)/brightness_factor) / ymax);
				blue = (int)(((mtot*255)/brightness_factor) / mmax);

				//if (red > 255) red = 255;
				//if (green > 255) green = 255;
				//if (blue > 255) blue = 255;
			}
			//////////////////////////////////////
	
			img_red[xi][yi] += red;
			img_green[xi][yi] += green;
			img_blue[xi][yi] += blue;
			img_alpha[xi][yi]++;
			

	        red = img_red[xi][yi]/img_alpha[xi][yi];
	        green = img_green[xi][yi]/img_alpha[xi][yi];
	        blue = img_blue[xi][yi]/img_alpha[xi][yi];
	        
	        if (red > 255) red = 255;
			if (green > 255) green = 255;
			if (blue > 255) blue = 255;

            red = red << 16;
			green = green << 8;
            pixels[yi*ximlen + xi] = 0xff000000 | (red & 0xffff0000) | (green & 0x0000ff00) | (blue & 0xff);
			
			//rough painting
			if (i < pixels.length)
			{
				offscreenGraphics.setColor( new Color(pixels[yi*ximlen + xi]) );
			    offscreenGraphics.fillRect(xi,yi,quad,quad);
			    
			    g_g.setColor( new Color(pixels[yi*ximlen + xi]) );
			    g_g.fillRect(xi,yi,quad,quad);
			}
		}
	
		if (i <= pixels.length - points*quad)
		{
			lasti += points*quad;
		}
		else if (quad > 1)
		{
			lasti = 0;
			quad = 1;
		}
		else //start of random points
		{
			lasti += points;
			
			if (lasti%pixels.length == 0)
			{
				long diff = (System.currentTimeMillis() - lastPassTime) / 1000;
				
				String strStatus = "Pass: " + (lasti/pixels.length) + " - " + diff + "s";
				System.out.println(strStatus);
				showStatus(strStatus);

				lastPassTime = System.currentTimeMillis();  	
			}
		}
		
	}

	
   public void mouseClicked(MouseEvent e) {
   }

   public void mouseEntered ( MouseEvent e ) {
   }

   public void mouseExited ( MouseEvent e ) {
   }
   
	
   	

   public void mousePressed ( MouseEvent e )
   {
   	 	xanchor = e.getX();
		yanchor = e.getY();

		if (e.getButton() == MouseEvent.BUTTON3 && e.isControlDown())
		{
			recurse = !recurse;
			reset = 1;
		}
		else if (e.getButton() == MouseEvent.BUTTON3 && e.isShiftDown())
		{
			mandelbrot = !mandelbrot;
			reset = 1;
		}
		else if (e.getButton() == MouseEvent.BUTTON3)
		{
			p.setVisible(!p.isVisible());
		}
		else if (e.getButton() == MouseEvent.BUTTON2)
		{
			xcen = 0;
			ycen = 0;
			zoom = 4;
			compoundPoint = false;
			reset = 1;
			newFractal();
		}
		else if (e.getButton() == MouseEvent.BUTTON1 && e.isShiftDown())
		{
			compoundPoint = true;
			cx = (float)  ( ((float)e.getX()/(float)ximlen-0.5f)*2*zoom + (float)xcen  );
			cy = (float) -( ((float)e.getY()/(float)yimlen-0.5f)*2*zoom + (float)ycen  );
			System.out.println("Expand point: " + cx + " " + cy);
			reset = 1;
		}
		else if (e.getButton() == MouseEvent.BUTTON1 && e.isControlDown())
		{
			scalePoint = !scalePoint;
			reset = 1;
		}

		setConstants();
		e.consume(); 	
   }

   public void mouseReleased ( MouseEvent e ) 
   {
   		xcurr = e.getX();
		ycurr = e.getY();
		mDown = false;
		
		int dx = Math.abs(xcurr - xanchor);
		int dy = Math.abs(ycurr - yanchor);
		if (dy > dx)  dx = dy;
		if (dx < 10)
		{
			return;
		}

		double newxcen = (xanchor + xcurr)/2.0;
		newxcen = ((newxcen/ximlen)-0.5)*2*zoom + xcen;
		xcen = newxcen;

		double newycen = (yanchor + ycurr)/2.0;
		newycen = ((newycen/yimlen)-0.5)*2*zoom + ycen;
		ycen = newycen;

		zoom = ((double)dx/ximlen)*zoom;
		System.out.println("Zoom is " + zoom);
		System.out.println("Xcen " + xcen);
		System.out.println("Ycen " + ycen);

		setConstants();
		System.gc();
		e.consume();
		reset = 1;
   } 	
	
	public void mouseDragged(MouseEvent e)
	{
		mDown = true;
		xcurr = e.getX();
		ycurr = e.getY();
		repaint();
		e.consume();
	}
	
	public void mouseMoved(MouseEvent e) 
	{
	}

	public void itemStateChanged(ItemEvent e) 
	{
		reset = 1;
	}
	
	public void keyPressed(KeyEvent e)
	{
	}
	
	public void keyReleased(KeyEvent e)
	{
		getConstants();
		reset = 1;
	}
	
	public void keyTyped(KeyEvent e)
	{
	}
	
}