package fractal;

//**************************************************************************
// Newtonbrot Fractor Generator
// Ryan Sweny
// rsweny@alumni.uwaterloo.ca
// 2005
//**************************************************************************

import java.awt.*;
import java.awt.image.PixelGrabber;
import java.awt.image.MemoryImageSource;
import java.util.*;
import math.Complex;
import java.awt.event.*;


public class Newtonbrot extends java.applet.Applet implements Runnable, KeyListener
{
	LinkedList roots = new LinkedList();

	Image img;
	MemoryImageSource screenMem;
	int[] pixels;

	Thread runner;

	double zoom = 6;
	double xcen = 0.1;
	double ycen = 0;

	int runtotal = 1;
	boolean mDown = false;
	int depth_blue = 50;
	int depth_green = 300;
	int depth_red = 1000;
	int points = 1000;

	float rootBoundry = 0.0000000000001f;
	float order = 4.7242613f;
	float img_order = 0.47967252f;
	float img_error = 1.0f;

	TextField cError;
	TextField c;
	TextField e;
	TextField tz;
	TextField tr;
	TextField tg;
	TextField tb;
	Panel p = new Panel(new GridLayout(9,2));
	TextField tX;
	TextField tY;
	TextField tContrast;
	TextField tBright;

	double gradient = 0.11;
	double brightness = 2.0;

	int yimlen = 500;
	int ximlen = 500;

	long[][] img_alpha = null;
	long[][] img_red = null;
	long[][] img_green = null;
	long[][] img_blue = null;
	double[][] img_cache = null;
	int replacePoint[][] = new int[3][10];
	
	long accepted = 0;
	long rejected = 0;
	int percent_accepted = 0;
	double[][][] smart_random_points;
	double[][] smart_random_points_score;
	double[] avg_score;
	int border_buffer = 40;

	int xcurr,ycurr,xanchor,yanchor;
	int alg = 0;

	boolean doInverse = false;
	boolean byStructure = true;
	boolean mandelbrotAddition = false;
	
	int pal = 0;
	long oldtime = System.currentTimeMillis();

    private int LOW_RES_PASSES = 540;
    int quad = 36;
	int limit = 2;
	
	//constants
	Complex two = new Complex(2, 0);
	Complex one = new Complex(1, 0);
	Complex pointone = new Complex(0.1, 0);
	
	float complex_error = 1;
	Complex oneError = new Complex(1,complex_error);
   
	public void start()
	{
		if (runner == null);
		{
			runner = new Thread(this);
			runner.setPriority(Thread.MIN_PRIORITY);
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

		String dred = getParameter("depth_red");
		if (dred != null)
		{
			depth_red = Integer.parseInt(dred);
		}

		String dgreen = getParameter("depth_green");
		if (dgreen != null)
		{
			depth_green = Integer.parseInt(dgreen);
		}

		String dblue = getParameter("depth_blue");
		if (dblue != null)
		{
			depth_blue = Integer.parseInt(dblue);
		}

		pixels = new int[ ximlen * yimlen ];
		screenMem = new MemoryImageSource(ximlen, yimlen, pixels, 0, ximlen);
		setBackground(Color.black);
		setLayout(new FlowLayout(FlowLayout.RIGHT));

		cError = new TextField(5);
		c = new TextField(5);
		e = new TextField(5);
		tz = new TextField(3);
		tr = new TextField(4);
		tg = new TextField(4);
		tb = new TextField(4);
		tContrast = new TextField(4);
		tBright = new TextField(4);
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
		tr.addKeyListener(this);
		tg.addKeyListener(this);
		tb.addKeyListener(this);

		
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

		cError.setText("1");
		c.setText("" + 0.4841751596941611);
		e.setText("" + 4.9242613);
		tz.setText("" + 6);

		tb.setText("40");
		tg.setText("140");
		tr.setText("200");
		tr.setForeground(Color.RED);
		tg.setForeground(Color.GREEN);
		tb.setForeground(Color.BLUE);


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
		
		Panel depthPanel = new Panel();
		depthPanel.add (tr);
		depthPanel.add (tg);
		depthPanel.add (tb);
		
		p.add(lblDepth);
		p.add(depthPanel);

		add(p);
		
		setConstants();
		getConstants();


		clearScreenAndReset(true);
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
					nextpoints();
					runtotal++;
					
					//do more calculations per paint later on
					if (runtotal > 400)
					{
						nextpoints();
						nextpoints();
					}
				}
				
				long denom = rejected+accepted;
				if (denom > 0)
    				percent_accepted = (int)(accepted*100/denom);
				else
					percent_accepted = 0;
	    				
	    		accepted = rejected = 0;
				
				if (runtotal%5 == 0) 
				{
					long curtime = System.currentTimeMillis();
	    			
					String strStatus = "Pass: " + runtotal + " - " + ((curtime-oldtime)/1000) + "s Accept: " + percent_accepted + " - " + ((int)(avg_score[0]*100)) + " " + ((int)(avg_score[1]*100)) + " " + ((int)(avg_score[2]*100));
					System.out.println(strStatus);
					showStatus(strStatus);
	    			
	    			oldtime = curtime;

					updateHistogram();
					repaint();	
					calcAvgScore();
					Thread.sleep(100);
				}
				Thread.sleep(10);
				
				//fading for LOW_RES_PASSES
				if (runtotal > 100 && runtotal <= 700 && runtotal%50==40)
				{
			        ensureBlack();
			        System.out.println("dark pass, limit is: " + limit + " quad: " + (quad-14));
			    }
			}
			catch (Exception e)
			{
				e.printStackTrace();
				int max_depth = Math.max(depth_green, Math.max(depth_blue, depth_red));
				img_cache = new double[max_depth][2];
			}
		}
		System.out.println("run() exited");
	}

	public void update(Graphics g)
	{
		paint(g);
	}
	
	private void updateHistogram()
	{
		int red,green,blue;
		double max_red = FindPeak(img_red);
		double max_green = FindPeak(img_green);
		double max_blue = FindPeak(img_blue);
		double max_alpha = FindPeak(img_alpha);

		max_red = Math.pow(max_red, gradient);
		max_green = Math.pow(max_green, gradient);
		max_blue = Math.pow(max_blue, gradient);
		max_alpha = Math.pow(max_alpha, gradient);

		for (int i=0; i<ximlen; i++)
		{
			for (int j=0; j<yimlen; j++)
			{
				
				if (byStructure)
				{
					if (img_alpha[i][j] > 0) 
					{
						double z = Math.pow(img_alpha[i][j], gradient)/max_alpha;
						red = (int)(img_red[i][j]*z*brightness/img_alpha[i][j] );
						green = (int)(img_green[i][j]*z*brightness/img_alpha[i][j] );
						blue = (int)(img_blue[i][j]*z*brightness/img_alpha[i][j]);
					}
					else
					{
						red = 0;
						green = 0;
						blue = 0;
					}
				}
				else
				{
					red = (int)(( brightness*Math.pow(img_red[i][j], gradient)/max_red)*255);
					green = (int)(( brightness*Math.pow(img_green[i][j], gradient)/max_green)*255);
					blue = (int)(( brightness*Math.pow(img_blue[i][j], gradient)/max_blue)*255);
				}

				if (red > 255) red = 255;
				if (green > 255) green = 255;
				if (blue > 255) blue = 255;

				red = red << 16;
				green = green << 8;
				
				int newcolor = 0xff000000 |(red & 0x00ff0000) | (green & 0x0000ff00) | (blue & 0xff);;
				pixels[j*ximlen + i] = newcolor;
			}
		}

		img = createImage(screenMem);
	}

	public void paint(Graphics g)
	{
		g.drawImage(img, 0, 0, this);
    	if (mDown) 
        {
			g.setColor(Color.white);
			g.drawRect(xanchor,yanchor,xcurr-xanchor,ycurr-yanchor);
        } 
	}
	
	public void ensureBlack()
	{
		for (int i = 0; i < ximlen; i++)
		{
			for (int j = 0; j < yimlen; j++)
			{
				if (img_alpha[i][j] > 0)
				{
					img_red[i][j] -= (long)((double)img_red[i][j]/img_alpha[i][j]);
					img_green[i][j] -= (long)((double)img_green[i][j]/img_alpha[i][j]);
					img_blue[i][j] -= (long)((double)img_blue[i][j]/img_alpha[i][j]);
				}
				
				if (runtotal % 100 == 40)
				{
					img_red[i][j] /= 2;
					img_green[i][j] /= 2;
					img_blue[i][j] /= 2;
					img_alpha[i][j] /= 2;
				}
			}
		}
		limit += 2;
		quad -= 2;
	}

	public double FindPeak(long[][] arr)
	{
		long max = 0;
		for (int i = 0; i < ximlen; i++)
		{
			for (int j = 0; j < yimlen; j++)
			{
				if (arr[i][j] > max)
				{
					max = arr[i][j];
				}
			}
		}
		
		//this fractal has extreme values at certain points, cut down the max value
		if (gradient <= 0.3)
			return Math.pow(max,0.9);
		else if (gradient <= 0.2)
			return Math.pow(max,0.7);
		else if (gradient <= 0.1)
			return Math.pow(max,0.4);
		else if (gradient <= 0.05)
			return Math.pow(max,0.3);
		else
			return max;
	}
	
	private synchronized void resetScore()
	{
		System.out.println("reset score");
		for (int j = 0; j < 3; j++)
    	{
    		avg_score[j] = 0;
	    	for (int i = 0; i < points; i++)
	    	{
	    		smart_random_points[j][i][0] = Math.random();
	    		smart_random_points[j][i][1] = Math.random();
	    		smart_random_points_score[j][i] = 0;
	    	}
    	}
    }
    
	public void clearScreenAndReset(boolean newRoots)
	{
		if (newRoots) roots = new LinkedList();
		oneError = new Complex(1,complex_error);
		
		img_red = new long[ximlen][yimlen+40];
		img_green = new long[ximlen][yimlen+40];
		img_blue = new long[ximlen][yimlen+40];
		img_alpha = new long[ximlen][yimlen+40];

		int max_depth = Math.max(depth_green, Math.max(depth_blue, depth_red));
		img_cache = new double[max_depth][2];
		
		smart_random_points = new double[3][points][2];
    	smart_random_points_score = new double[3][points];
    	avg_score = new double[3];

		resetScore();

		runtotal = 1;
		quad = 36;
		limit = 2;
	}

	private void nextpoints()
	{
		int s2 = quad/2;
		
		int iter;
		int pointx;
		int pointy;
		double close = 0.001;

		for (int i=0; i<points; i++)
		{
			int counter = 0;
			Complex z = new Complex(0,0);

			if (byStructure)
			{
				//// color by palette  /////////////////////////////////////////////////////////////////////
				float[] data = calcOrbit(i, depth_red, 0);
				iter = (int)data[0];
				int basecolor = (int)Math.floor(data[1]);
				float fractioncolor = data[1] - basecolor;

				int amtred = (int)((1.0f-fractioncolor)*(float)Pallet.fpalette[pal][basecolor][0] + fractioncolor*(float)Pallet.fpalette[pal][basecolor+1][0]);
				int amtgreen = (int)((1.0f-fractioncolor)*(float)Pallet.fpalette[pal][basecolor][1] + fractioncolor*(float)Pallet.fpalette[pal][basecolor+1][1]);
				int amtblue = (int)((1.0f-fractioncolor)*(float)Pallet.fpalette[pal][basecolor][2] + fractioncolor*(float)Pallet.fpalette[pal][basecolor+1][2]);

				if ( (iter == depth_red && doInverse) || (iter < depth_red && !doInverse) )
				{
					for (int j = 0; j < iter; j++)
					{
						pointx = (int)Math.floor ( ((img_cache[j][0]+xcen)*ximlen)/zoom + ximlen/2 );
						pointy = (int)Math.floor ( ((img_cache[j][1]+ycen)*yimlen)/zoom + yimlen/2 );
						
						if (pointx >= -border_buffer && pointy >= -border_buffer && pointx < ximlen + border_buffer && pointy < yimlen + border_buffer)
						{
							counter++;
							if (pointx >= 0 && pointy >= 0 && pointx < ximlen && pointy < yimlen)
							{
								//tag this counter since it was actually drawn to the screen
								counter += 0.001;
							
								if ( j > 5 && doInverse &&
									(  (Math.abs(img_cache[j][0] - img_cache[j-1][0]) < close && Math.abs(img_cache[j][1] - img_cache[j-1][1]) < close ) 
									|| (Math.abs(img_cache[j][0] - img_cache[j-2][0]) < close && Math.abs(img_cache[j][1] - img_cache[j-2][1]) < close ) 
									|| (Math.abs(img_cache[j][0] - img_cache[j-3][0]) < close && Math.abs(img_cache[j][1] - img_cache[j-3][1]) < close )
									)) 
									j = iter;
						
								img_alpha[pointx][pointy]+=10;
								img_red[pointx][pointy]+=amtred*10;
								img_green[pointx][pointy]+=amtgreen*10;
								img_blue[pointx][pointy]+=amtblue*10;
								
								//fill in rough blocks after roots are calculated and until LOW_RES_PASSES
								if (runtotal > 50 && runtotal < LOW_RES_PASSES)
								{
									try
									{
										for (int k = 7; k < quad - 7; k++)
										{
											for (int m = 7; m < quad - 7; m++)
											{
												long diff = limit - img_alpha[pointx+k-s2][pointy+m-s2];
												if (diff > 0)
												{
													double center_bright = 1.0d - (Math.max(Math.abs(k-s2),Math.abs(m-s2))*2) / (double)quad;
													center_bright *= center_bright;
													
													//center_bright = Math.min(diff, center_bright*9)+1;
													center_bright *= 10;
													
													img_alpha[pointx+k-s2][pointy+m-s2] += center_bright;
													img_red[pointx+k-s2][pointy+m-s2] += amtred*center_bright;
													img_green[pointx+k-s2][pointy+m-s2] += amtgreen*center_bright;
													img_blue[pointx+k-s2][pointy+m-s2] += amtblue*center_bright;
												}
											}
										}
									}
									catch(Exception e)
									{
									}
								}
							}
						}
					}
				}
				else
				{
					counter = 0;
				}
				updateScore(counter, 0, i);
			} //end if byStructure
			else
			{
				//// red  /////////////////////////////////////////////////////////////////////
				float[] data = calcOrbit(i, depth_red, 0);
				iter = (int)data[0];
				if ( (iter == depth_red && doInverse) || (iter < depth_red && !doInverse) )
					counter = drawOrbit(i, iter, Color.red, img_red);
				else 
					counter = 0;
					
				updateScore(counter, 0, i);

				//// green  /////////////////////////////////////////////////////////////////////
				data = calcOrbit(i, depth_green, 1);
				iter = (int)data[0];
				if ( (iter == depth_green && doInverse) || (iter < depth_green && !doInverse) )
					counter = drawOrbit(i, iter, Color.green, img_green);
				else 
					counter = 0;

				updateScore(counter, 1, i);

				//// blue  /////////////////////////////////////////////////////////////////////
				data = calcOrbit(i, depth_blue, 2);
				iter = (int)data[0];
				if ( (iter == depth_blue && doInverse) || (iter < depth_blue && !doInverse) )
					counter = drawOrbit(i, iter, Color.blue, img_blue);
				else 
					counter = 0;

				updateScore(counter, 2, i);
			} //end byDepth
		}
		
		//get all new random points after a while
		if ( runtotal%5000 == 0 )
			resetScore();
	}
	
	private void calcAvgScore()
	{
		long total_c = 0;
		for (int i = 0; i < 3; i++)
		{
			double c = 0;
			for (int j = 0; j < points; j++)
			{
				c += smart_random_points_score[i][j];
			}
			avg_score[i] = c / points;
			total_c += (long)c;
			
			if (avg_score[i] > 1 && border_buffer == 40)
			{
				//propogate these good points to the other 2 color channels
				int nextColor = (i+1)%3;
				int lastColor = (i+2)%3;
				int start = (int)(Math.random() * (points/2));
				
				if (avg_score[nextColor] < 0.4)
				{
					System.out.println("propogate " + i + " to " + nextColor);
					for (int j = start; j < start + points/4; j++)
					{
						smart_random_points[nextColor][j][0] = smart_random_points[i][j][0];
						smart_random_points[nextColor][j][1] = smart_random_points[i][j][1];
						smart_random_points_score[nextColor][j] = smart_random_points_score[i][j];
					}
				}
				
				if (avg_score[lastColor] < 0.4)
				{
					System.out.println("propogate " + i + " to " + lastColor);
					for (int j = start; j < start + points/4; j++)
					{
						smart_random_points[lastColor][j][0] = smart_random_points[i][j][0];
						smart_random_points[lastColor][j][1] = smart_random_points[i][j][1];
						smart_random_points_score[lastColor][j] = smart_random_points_score[i][j];
					}
				}
			}
		}
		
		if (total_c < 2) 
		{
			//we have no points to draw, probably because this is a deep zoom. Increase our capture size
			border_buffer = border_buffer*2;
			border_buffer = Math.min(80000, border_buffer);
		}
		else
		{
			//we have some good points, gradually cut the border back to the normal 40
			if (border_buffer == 40 || (percent_accepted > 15 && zoom > 0.01)) 
			{
				border_buffer = 40;
			}
			else
			{
				border_buffer = (int)(border_buffer / 1.2);
				
				if (border_buffer < 40)
				{
					System.out.println("Border fixed, clearing scores.");
					for (int j = 0; j < 3; j++)
			    	{
			    		avg_score[j] = 0;
			    	}
				
					border_buffer = 40;
				}
			}	
		}
		if (border_buffer != 40) System.out.println("border_buffer " + border_buffer + " - " + total_c);
    }
	
	private int drawOrbit(int i, int iter, Color colorObj, long[][] img_color)
	{
		int pointx,pointy;
		int counter = 0;
		double close = 0.001;
		
		for (int j = 0; j < iter; j++)
		{
			pointx = (int)Math.floor ( ((img_cache[j][0]+xcen)*ximlen)/zoom + ximlen/2 );
			pointy = (int)Math.floor ( ((img_cache[j][1]+ycen)*yimlen)/zoom + yimlen/2 );
				
			if (pointx >= -border_buffer && pointy >= -border_buffer && pointx < ximlen + border_buffer && pointy < yimlen + border_buffer)
			{
				counter++;
				if (pointx >= 0 && pointy >= 0 && pointx < ximlen && pointy < yimlen)
				{
					//tag this counter since it was actually drawn to the screen
					counter += 0.001;
					
					img_color[pointx][pointy]++;
					if ( j > 5 && doInverse &&
						(  (Math.abs(img_cache[j][0] - img_cache[j-1][0]) < close && Math.abs(img_cache[j][1] - img_cache[j-1][1]) < close ) 
						|| (Math.abs(img_cache[j][0] - img_cache[j-2][0]) < close && Math.abs(img_cache[j][1] - img_cache[j-2][1]) < close ) 
						|| (Math.abs(img_cache[j][0] - img_cache[j-3][0]) < close && Math.abs(img_cache[j][1] - img_cache[j-3][1]) < close )
						)) 
						j = iter;
				}
			}
		}
		return counter;
	}
	
	
	private float[] calcOrbit(int i, int depth, int color)
	{
		double xf,yf,a,b,aa,bb,aanew,bbnew,mag2;
		int iter = 0;

		xf = smart_random_points[color][i][0];
		yf = smart_random_points[color][i][1];
		
		a = 2*(xf - 0.5f);
		b = 2*(yf - 0.5f);
		aa = a;
		bb = b;
		
		Complex z = new Complex(aa,bb);
		Complex old = new Complex(z.re, z.im);

		z = Newt(z, aa, bb, alg);

		while(iter < depth && distance(old,z) > rootBoundry)
		{
			old = new Complex(z.re,z.im);

			z = Newt(z, aa, bb, alg);
			if (mandelbrotAddition) z = z.add(new Complex(aa/2,bb/2));

			img_cache[iter][0] = z.re;
			img_cache[iter][1] = z.im;
			iter++;
		}
		
		float solutionColor;
		if (iter != depth_red || (doInverse && runtotal > 25) )
		{
			if (mandelbrotAddition) 
			{
				//get Color based on angle of z
				z = Newt(z, aa, bb, alg);
				double angle = Math.atan(z.im/z.re);
				angle = angle + (Math.PI / 2);
				angle = angle / Math.PI;
				solutionColor = (float)(angle*254.0);
			}
			else
			{
				float curroot = (float)addRoot(z);
				solutionColor = ( (curroot*254) / roots.size() );
			}
		}
		else 
		{
			solutionColor = 0;
		}
		
		float[] data = new float[2];
		data[0] = iter;
		data[1] = solutionColor;
		
		return data;
	}
	
	private synchronized void updateScore(double counter, int color, int i)
	{
		//score how interesting this point is
		if (smart_random_points_score[color][i] > 0)
			smart_random_points_score[color][i] = 0.5*counter + 0.5*smart_random_points_score[color][i];
		else
			smart_random_points_score[color][i] = counter;
			
		//if we are struggling and expanded the search zone, give points that actually display a +1 bonus
		if (border_buffer > 40)
		{	
			int roundCounter = (int)counter;
			double countDiff = ((double)roundCounter) - counter;
			if (countDiff != 0) smart_random_points_score[color][i] += 1;
		}
			
		//if the point had more than one hit, keep it
		if (smart_random_points_score[color][i] > 1.1)
		{
			accepted++;
			double factor = 0.01;
			
			if (smart_random_points_score[color][i] > avg_score[color]*(1.2+zoom))
			{
				//procreate this extra good point
				smart_random_points[color][replacePoint[color][0]][0] = smart_random_points[color][i][0] + (Math.random() - 0.5)*factor*zoom;
				smart_random_points[color][replacePoint[color][0]][1] = smart_random_points[color][i][1] + (Math.random() - 0.5)*factor*zoom;
				
				//overwrite a previously identified weak point
				int k = 0;
				for (k = 0; k < replacePoint[color].length - 1; k++)
				{
					replacePoint[color][k] = replacePoint[color][k+1];
				}
				replacePoint[color][k-1] = i;
			}
			
			//cull some points
			if (avg_score[color] > 3 && i%3 == 0)
				factor = factor*5;
		
			//this is an interesting point, try another one near it next time
			smart_random_points[color][i][0] += (Math.random() - 0.5)*factor*zoom;
			smart_random_points[color][i][1] += (Math.random() - 0.5)*factor*zoom;
		}
		else
		{
			//try another random point
			rejected++;
			smart_random_points_score[color][i] = 0;
			smart_random_points[color][i][0] = Math.random();
	    	smart_random_points[color][i][1] = Math.random();
			
			//add this weak point to the front of the list
			int k = 1;
			for (k = 1; k < replacePoint[color].length; k++)
			{
				replacePoint[color][k] = replacePoint[color][k-1];
			}
			replacePoint[color][0] = i;
		}
	}

	public boolean mouseDown(Event evt, int x, int y)
	{
		xanchor = x;
		yanchor = y;
		
		if (evt.modifiers == Event.META_MASK)
		{
			p.setVisible(!p.isVisible());
		}
		else if (evt.modifiers == Event.ALT_MASK)
		{
			alg++;
			if (alg > 5) {alg = 0;}
			System.out.println("Alg is " + alg);
			clearScreenAndReset(true);
		}
		else if (evt.modifiers == (Event.META_MASK + Event.CTRL_MASK))
		{
			doInverse = !doInverse;
			runtotal = 1;
			clearScreenAndReset(false);
		}
		else if (evt.modifiers == (Event.META_MASK + Event.SHIFT_MASK))
		{
			mandelbrotAddition = !mandelbrotAddition;
			clearScreenAndReset(true);
		}
		else if (evt.modifiers == Event.CTRL_MASK)
		{
			byStructure = !byStructure;
			clearScreenAndReset(false);
		}
		else if (evt.modifiers == Event.SHIFT_MASK)
		{
			pal = (pal +1)%12;
			runtotal = 1;
			clearScreenAndReset(false);
		}

		return true;
	}


	public boolean mouseDrag(Event evt, int x, int y)
	{
		mDown = true;
		xcurr = x;
		ycurr = y;
		repaint();
		return true;
	}

	public boolean mouseUp(Event evt, int x, int y)
	{
		xcurr = x;
		ycurr = y;
		mDown = false;
		
		int dx = Math.abs(xcurr - xanchor);
		int dy = Math.abs(ycurr - yanchor);
		if (dy > dx)  dx = dy;
		
        //zoom is too small
       	if (dx > 10)
       	{
			double newxcen = (xanchor + xcurr)/2.0f;
			newxcen = xcen - ((newxcen - ximlen/2)/(ximlen))*zoom;
			xcen = newxcen;
	
			double newycen = (yanchor + ycurr)/2.0f;
			newycen = ycen - ((newycen - yimlen/2)/(yimlen))*zoom;
			ycen = newycen;
	
			zoom = (double)dx/ximlen*zoom;
			clearScreenAndReset(false);
		}

		setConstants();
		return true;
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
			return z.sub( DivZ2Function(z,i,j).div(DivZ2Derivative(z,i,j)) );
		}
		else if (mode == 4)
		{
			return z.sub( Poly3Function(z,i,j).div(Poly3Derivative(z,i,j)) );
		}
		else if (mode == 5)
		{
			return z.sub( DivZFunction(z,i,j).div(DivZDerivative(z,i,j)) );
		}
		else
		{
			return z.sub( PolyFunction(z,i,j).div(PolyDerivative(z,i,j)) );
		}
	}

	public double distance(Complex a, Complex b)
	{
		return Math.sqrt((a.re-b.re)*(a.re-b.re) + (a.im-b.im)*(a.im-b.im));
	}
	


	// (z^n - 1) / z
	public Complex DivZ2Function(Complex z, double i, double j)
	{
		Complex exponent = new Complex(order, img_order);
		return (z.pow(exponent).sub( one )).div (z);
	}
 
 	// (n-1)z^(n-2) + z^(-2)
	public Complex DivZ2Derivative(Complex z, double i, double j)
	{
		Complex minusTwo = new Complex(-2,0);
		Complex exponent = new Complex(order, img_order);
		Complex exponentLessOne = exponent.sub(oneError);
		Complex exponentLessTwo = exponent.sub(two);
		return (exponentLessOne.mul( z.pow(exponentLessTwo) )).add ( z.pow(minusTwo) );
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
	
	//z^c - z + 1
	public Complex Poly3Function(Complex z, double i, double j) 
	{
		Complex c = new Complex(10+order, img_order);
		return z.pow(c).sub(z).add(pointone);
	}
	
	//c*z^(c-1) - 1
	public Complex Poly3Derivative(Complex z, double i, double j) 
	{
		Complex cminus1 = new Complex(10+order-1, (img_order == 0) ? 0 : img_order-1);
		Complex c = new Complex(10+order, img_order);
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
	

	public double addRoot(Complex root)
	{
		ListIterator i = (ListIterator)roots.iterator();
		double colorFactor = 2/rootBoundry;
		
		
		int closest = 0;
		double minDistance = 100;
		
		while (i.hasNext())
		{
			int index = i.nextIndex();
			Complex c = (Complex)i.next();
			
			double dist = distance(root, c);
			if (dist < minDistance) 
			{	
				minDistance = dist;
				closest = index;
			}
			
			if ( dist <= rootBoundry*1.1 )
			{
    			double retVal = index + colorFactor*distance(root, c);
				return Math.min(roots.size() - 1, retVal);
			}
		}

		if (roots.size() < 200 && runtotal < 25)
		{
			roots.add(root);
			System.out.println("("+runtotal+") add root:" + roots.size());
		}

		if (doInverse)
		{
			return closest;
		}
		else
		{
			return roots.size() - 1;
		}
	}
	
	private void getConstants()
	{
		try
		{
			xcen = Double.parseDouble(tX.getText());
			ycen = Double.parseDouble(tY.getText());
			gradient = Double.parseDouble(tContrast.getText());
			brightness = Double.parseDouble(tBright.getText());
			zoom = Double.parseDouble(tz.getText());
			
			img_order = Float.parseFloat(c.getText());
			order = Float.parseFloat(e.getText());
			complex_error = Float.parseFloat(cError.getText());
			
			depth_red = Integer.parseInt(tr.getText());
			depth_green = Integer.parseInt(tg.getText());
			depth_blue = Integer.parseInt(tb.getText());
		}
		catch(Exception e)
		{
		}
	}
	
	private void setConstants()
	{
		try
		{
			tContrast.setText("" + gradient);
			tBright.setText("" + brightness);
			tX.setText("" + xcen);
			tY.setText("" + ycen);
			tz.setText("" + zoom);
			cError.setText(complex_error + "");
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
		if (ev.getSource() != tContrast && ev.getSource() != tBright)
			clearScreenAndReset(ev.getSource() == c || ev.getSource() == e);

	}
	
	public void keyTyped(KeyEvent e)
	{
	}



}







