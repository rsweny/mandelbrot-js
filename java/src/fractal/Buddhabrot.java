package fractal;

//**************************************************************************
// Buddhabrot Fractor Generator
// Ryan Sweny
// rsweny@alumni.uwaterloo.ca
// 2002
//**************************************************************************

import java.awt.*;
import java.awt.image.PixelGrabber;
import java.awt.image.MemoryImageSource;
import java.awt.event.*;


public class Buddhabrot extends java.applet.Applet implements Runnable, MouseListener, MouseMotionListener, KeyListener
{
	Label lblZoom = new Label("Zoom:");
	Label lblX = new Label("X:");
	Label lblY = new Label("Y:");
	Label lblContrast = new Label("Contrast:");
	Label lblBrightness = new Label("Brightness:");
	Label lblJuliaX = new Label("Julia X:");
	Label lblJuliaY = new Label("Julia Y:");
	
	TextField tzoom;
	TextField tX;
	TextField tY;
	TextField tJuliaX;
	TextField tJuliaY;
	TextField tContrast;
	TextField tBrightness;
	TextField tr;
	TextField tg;
	TextField tb;
	Checkbox zoomType = new Checkbox("Magnify Zoom", null, true);
	Panel p = new Panel(new GridLayout(9,2));
	
    Graphics g_g;
    Image img;
    MemoryImageSource screenMem;
    int[] pixels;
    
	Thread[] runner;
	boolean running = false;

    double zoom = 5;
    double xcen = 0.0;
    double ycen = 0.0;

    long runtotal = 0;
    boolean mDown = false; 
    int depth_blue;
    int depth_green;
    int depth_red;
    int points = 2000;
    
    double cx = 0;
    double cy = 0;
    double gradient = 0.8;
    double brightness = 1.1;
    
    int yimlen;
    int ximlen;
    int numPixels;
    
    int xcurr,ycurr,xanchor,yanchor;
    int alg = 0;
    boolean doHolo = false;
    boolean doInverse = false;
    
    long oldtime = System.currentTimeMillis();
    
    long accepted = 0;
	long rejected = 0;
	int percent_accepted = 0;
	double[][][] smart_random_points;
	double[][] smart_random_points_score;
	double[] avg_score;
	int border_buffer = 40;
    
    long[][] img_red = null;
	long[][] img_green = null;
	long[][] img_blue = null;

	int replacePoint[][] = new int[3][10];
	int lasti;

	public void start()
	{
		running = true;
		if (runner[0] == null);
		{
			runner[0] = new Thread(this);
			runner[0].setPriority(Thread.MIN_PRIORITY);
			runner[0].start();
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
    	
    	numPixels = ximlen * (4*yimlen/5);
    	pixels = new int[ ximlen * yimlen ];
    	
    	smart_random_points = new double[3][points][2];
    	smart_random_points_score = new double[3][points];
    	avg_score = new double[3];
 
 		resetScore();
    	
		screenMem = new MemoryImageSource(ximlen, yimlen, pixels, 0, ximlen);
		g_g = getGraphics();
		g_g.setPaintMode();
		
		setBackground(Color.black);
		setLayout(new FlowLayout(FlowLayout.LEFT));
		
		tzoom = new TextField(20);
		tX = new TextField(20);
		tY = new TextField(20);
		tJuliaX = new TextField(20);
		tJuliaY = new TextField(20);
		tContrast = new TextField(20);
		tBrightness = new TextField(20);
		tr = new TextField(4);
		tg = new TextField(4);
		tb = new TextField(4);
		
		tzoom.addKeyListener(this);
		tContrast.addKeyListener(this);
		tBrightness.addKeyListener(this);
		tX.addKeyListener(this);
		tY.addKeyListener(this);
		tJuliaX.addKeyListener(this);
		tJuliaY.addKeyListener(this);
		tr.addKeyListener(this);
		tg.addKeyListener(this);
		tb.addKeyListener(this);
		
		setConstants();

		p.add(lblZoom);
		p.add (tzoom);
		
		p.add(lblX);
		p.add (tX);
		
		p.add(lblY);
		p.add (tY);
		
		p.add(lblJuliaX);
		p.add (tJuliaX);
		
		p.add(lblJuliaY);
		p.add (tJuliaY);
		
		p.add(lblContrast);
		p.add(tContrast);
		
		p.add(lblBrightness);
		p.add(tBrightness);
		
		Label lblDepth = new Label("Depth:");
		tb.setText("" + depth_blue);
		tg.setText("" + depth_green);
		tr.setText("" + depth_red);
		tr.setForeground(Color.RED);
		tg.setForeground(Color.GREEN);
		tb.setForeground(Color.BLUE);
		
		Panel depthPanel = new Panel();
		depthPanel.add (tr);
		depthPanel.add (tg);
		depthPanel.add (tb);
		
		p.add(lblDepth);
		p.add(depthPanel);
		
		zoomType.setState(true);
		p.add(zoomType);
		
		p.setBackground(Color.WHITE);
		add(p);
		
		addMouseListener(this);
		addMouseMotionListener(this);
		enableEvents(AWTEvent.MOUSE_EVENT_MASK | AWTEvent.MOUSE_MOTION_EVENT_MASK ); 
		
		clearScreenAndReset(0);
    }
   
	public void run()
	{
		Thread thisThread = Thread.currentThread();
		while ( running )
		{
			try 
			{
				if (!mDown) nextpoints();
				
				long denom = rejected+accepted;
				if (denom > 0)
    				percent_accepted = (int)(accepted*100/denom);
				else
					percent_accepted = 0;
				
    			accepted = rejected = 0;
				
				//only paint every 20th pass
				if (runtotal%20 == 0 || mDown) 
				{
					long curtime = System.currentTimeMillis();
					String strStatus = "Pass: " + runtotal + " - " + ((curtime-oldtime)/1000) + "s Accept: " + percent_accepted + " - " + ((int)(avg_score[0]*100)) + " " + ((int)(avg_score[1]*100)) + " " + ((int)(avg_score[2]*100));
					System.out.println(strStatus);
					showStatus(strStatus);
					oldtime = curtime;
				
					updateHistogram();
					repaint();	
					calcAvgScore();
					Thread.sleep(400);
				}
				Thread.sleep(10);
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
    
    private void updateHistogram()
	{
    	int red,green,blue;
    	
    	double max_red = Math.pow(FindPeak(img_red), gradient);
    	double max_green = Math.pow(FindPeak(img_green), gradient);
    	double max_blue = Math.pow(FindPeak(img_blue), gradient);
        
        for (int i=0; i<ximlen; i++) 
        {
            for (int j=0; j<yimlen; j++) 
            {
	            red = (int)(( brightness*Math.pow(img_red[i][j], gradient)/max_red)*255);
	            green = (int)(( brightness*Math.pow(img_green[i][j], gradient)/max_green)*255);
	            blue = (int)(( brightness*Math.pow(img_blue[i][j], gradient)/max_blue)*255);

	            if (red > 255) red = 255;
	            if (green > 255) green = 255;
	            if (blue > 255) blue = 255;
	            
                red = red << 16;
				green = green << 8;
                pixels[j*ximlen + i] = 0xff000000 | (red & 0xffff0000) | (green & 0x0000ff00) | (blue & 0xff);
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
    
    
    private long FindPeak(long[][] arr)
    {
    	long max = 0;
    	for (int i = 1; i < ximlen; i++)
    	{
    		for (int j = 1; j < yimlen; j++)
    		{
    			if (arr[i][j] > max)
    			{
    				max = arr[i][j];
    			}
    		}
    	}
    	return max;
    }
	
    private void clearScreenAndReset(double newcx) 
    {
    	img_red = new long[ximlen][yimlen];
    	img_green = new long[ximlen][yimlen];
    	img_blue = new long[ximlen][yimlen];
		
    	lasti = 0;
        runtotal = 0;
        cx = newcx;
        resetScore();
    }
    
    private double getStartX(double xf)
    {
	    if (zoomType.getState())
        {
           	return  4*(xf - 0.5f);
		}
		else
		{
			return  2*(xf - 0.5f)*zoom + xcen;
        }
	 }
	 
	private double getStartY(double yf)
    {
	    if (zoomType.getState())
        {
			return 4*(yf - 0.5f);
		}
		else
		{
            return 2*(yf - 0.5f)*zoom + ycen;
        }
	 }
	 
	 private double[] getNextStartOrbit(int i, int color)
	 {
	 	double[] z = new double[2];
	 	double thread_rnd = Math.random()*0.00001*zoom;

	 	if ( i <= numPixels && i > 6000 && zoom > 1) //do every pixel in order until screen covered
		{
			int xi = i%ximlen;
			double yi = ((double)i / (double)ximlen);
			z[0] = (double)xi / (double)ximlen;
			z[1] = yi / (double)yimlen;
		}
		else
		{
			z[0] = smart_random_points[color][i%points][0] + thread_rnd;
			z[1] = smart_random_points[color][i%points][1] + thread_rnd;
		}
		
		return z;
	}
	
	 
	 
	private int calcOrbit(int i, int depth, double[][] img_cache, int color)
	{
		double a,b,aa,bb,aanew,bbnew,mag2;
        int iter;
         
		double[] z = getNextStartOrbit(i,color);
        a = getStartX(z[0]);
        b = getStartY(z[1]);
			
        aa = a;
        bb = b;
        iter = 0;
        
        if (doHolo)
        {
			z = getNextStartOrbit(i,color);
	        a = getStartX(z[0]);
	        b = getStartY(z[1]);
        }
         
        do
		{
			aanew = calcReal(aa, bb, a , b);
			bbnew = calcImg(aa, bb, a, b);
			aa = aanew;
			bb = bbnew;
			mag2 = aa*aa + bb*bb;
			
			img_cache[iter][0] = aa;
			img_cache[iter][1] = bb;
			
			iter++;
		} 
		while ( (iter < depth) && (mag2 < 4) );
		
		return iter;
	}
	 
	 private double drawOrbit(int i, int iter, double[][] img_cache, Color colorObj, int factor, long[][] img_color, boolean normalZoom)
	 {
	 	int pointx,pointy;
	 	double draw_cutoff = factor*0.2 + 1;
	 	double counter = 0;
	 	double close = 0.001;
		g_g.setColor( colorObj );
		
		for (int j = 0; j < iter; j++)
		{
			if (normalZoom)
			{
				pointx = (int)Math.floor ( ((img_cache[j][0]+xcen)*ximlen)/zoom + ximlen/2 );
				pointy = (int)Math.floor ( ((img_cache[j][1]+ycen)*yimlen)/zoom + yimlen/2 );
			}
			else
			{
				pointx = (int)Math.floor (((img_cache[j][0]+2.0f)/4.0f)*ximlen);
				pointy = (int)Math.floor (((img_cache[j][1]+2.0f)/4.0f)*yimlen);
			}
			
			if (pointx >= -border_buffer && pointy >= -border_buffer && pointx < ximlen + border_buffer && pointy < yimlen + border_buffer)
			{
				counter++;
				if (pointx >= 0 && pointy >= 0 && pointx < ximlen && pointy < yimlen)
				{
					//tag this counter since it was actually drawn to the screen
					counter += 0.001;
					img_color[pointx][pointy]++; 

					//stop accumulating if the point is just repeating itself
					if ( j > 10*factor && doInverse &&
						(  (Math.abs(img_cache[j][0] - img_cache[j-1][0]) < close && Math.abs(img_cache[j][1] - img_cache[j-1][1]) < close ) 
						|| (Math.abs(img_cache[j][0] - img_cache[j-2][0]) < close && Math.abs(img_cache[j][1] - img_cache[j-2][1]) < close ) 
						|| (Math.abs(img_cache[j][0] - img_cache[j-3][0]) < close && Math.abs(img_cache[j][1] - img_cache[j-3][1]) < close )
						)) 
						j = iter;
					
					if (i < numPixels*draw_cutoff && !mDown) 
					{
						g_g.drawLine(pointx,pointy,pointx,pointy);
					}		
				}
			} 
		}
		
		//return the number of points we actually plotted to the screen
		return counter;
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
			border_buffer = Math.min(40000, border_buffer);
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
	
    private void nextpoints()
    {
    	int iter;
		double counter;
		boolean normalZoom = zoomType.getState();
		
		int max_depth = Math.max(depth_green, Math.max(depth_blue, depth_red));
		double[][] img_cache = new double[max_depth][2];
	
		for (int i=lasti; i<lasti+points; i++) 
		{	
			Thread.yield();
			
			//// red  ///////////////////////////////////////////////////////////////////// 
			iter = calcOrbit(i, depth_red, img_cache, 0);
			if ( (iter == depth_red && doInverse) || (iter < depth_red && !doInverse) ) 
			{
				counter = drawOrbit(i, iter, img_cache, Color.red, 3, img_red, normalZoom);
			}
			else
			{
				counter = 0;
			}
			updateScore(counter, 0, i%points);

			//// green  ///////////////////////////////////////////////////////////////////// 
			iter = calcOrbit(i, depth_green, img_cache, 1);
			if ( (iter == depth_green && doInverse) || (iter < depth_green && !doInverse) ) 
			{ 
				counter = drawOrbit(i, iter, img_cache, Color.green, 2, img_green, normalZoom);
			}
			else
			{
				counter = 0;
			}
			updateScore(counter, 1, i%points);

			//// blue  ///////////////////////////////////////////////////////////////////// 
			iter = calcOrbit(i, depth_blue, img_cache, 2);
			if ( (iter == depth_blue && doInverse) || (iter < depth_blue && !doInverse) ) 
			{ 
				counter = drawOrbit(i, iter, img_cache, Color.blue, 1, img_blue, normalZoom);
			}
			else
			{
				counter = 0;
			}
			updateScore(counter, 2, i%points);
		}
		 
		Thread thisThread = Thread.currentThread();
		if (thisThread == runner[0])
		{
			if (lasti <= numPixels*2)
			{
				lasti += points;
			}
			else if (lasti == numPixels*2 + points)
			{
				lasti++;
				for (int k = 1; k < runner.length; k++)
				{
					System.out.println("starting remaining thread " + k);
					if (runner[k] == null)
					{
						runner[k] = new Thread(this);
						runner[k].setPriority(Thread.MIN_PRIORITY);
						runner[k].start();
					}
				}
			}
		}
		runtotal++;
		
		//get all new random points after a while
		if ( (runtotal%200 == 0 && percent_accepted > 40 && zoom > 0.01) || 
		     (runtotal%2000 == 0 && percent_accepted > 20))
			resetScore();
	}
 
 
   public void mouseClicked(MouseEvent e) {
   }

   public void mouseEntered ( MouseEvent e ) {
   }

   public void mouseExited ( MouseEvent e ) {
   }

   public void mousePressed ( MouseEvent e )
   {
   		mDown = true;

        xanchor = e.getX();
		yanchor = e.getY();
		
		int input_modifiers = e.getModifiers(); 
		
		if (e.getButton() == MouseEvent.BUTTON3 && e.isControlDown())
		{
			doInverse = !doInverse;
			
			if (doInverse) gradient = 0.25;
			else gradient = 0.8;
			
			setConstants();
			clearScreenAndReset(cx);
		}
		else if (e.getButton() == MouseEvent.BUTTON3)
		{
			p.setVisible(!p.isVisible());

		}
		else if (e.getButton() == MouseEvent.BUTTON2)
		{
			alg++;
			zoom = 3.7;	   
		    xcen = 0.4;
		    ycen = 0;
		    setConstants();
		    
		    if (alg > 6) alg = 1;
		    
			clearScreenAndReset(0);
		}
		else if (e.getButton() == MouseEvent.BUTTON1 && e.isControlDown())
		{
			doHolo = !doHolo;
			System.out.println("Holographic mode is " + doHolo);
		}
		else if (e.getButton() == MouseEvent.BUTTON1 && e.isShiftDown()) 
		{
       		gradient = 0.25;
       		setConstants();
			cx =  ((double)e.getX()/(double)ximlen-0.5f)*2*zoom + xcen;
	        cy = -((double)e.getY()/(double)yimlen-0.5f)*2*zoom + ycen;
       		clearScreenAndReset(cx);
       		
       		setConstants();
       		
		}
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
        
        //zoom is too small
       	if (dx < 10) return;

        double newxcen = (xanchor + xcurr)/2.0f;
        double newycen = (yanchor + ycurr)/2.0f;
        if ( zoomType.getState() )
        {
        	newxcen = xcen - ((newxcen - ximlen/2)/(ximlen))*zoom;
			newycen = ycen - ((newycen - yimlen/2)/(yimlen))*zoom;
	    }
	    else
	    {
	        newxcen = (newxcen/ximlen-0.5f)*2*zoom + xcen;
	        newycen = (newycen/yimlen-0.5f)*2*zoom + ycen;
		}
		xcen = newxcen;
		ycen = newycen;

        zoom = (double)dx/ximlen*zoom;
        
        clearScreenAndReset(cx);
        setConstants();
        e.consume();
    } 	
	
	public void mouseDragged(MouseEvent e)
	{
		xcurr = e.getX();
		ycurr = e.getY();
		repaint();
		e.consume();
	}
	
	public void mouseMoved(MouseEvent e) 
	{
	}  
     
    private double calcReal(double aa, double bb, double a, double b)
    {
    	if (cx != 0)
    	{
    		a = cx;
    		b = cy;
    	}
    	
    	if (alg == 1)
    	{
    		return aa*aa - Math.sin(bb)*(bb) + Math.atan(a*b);
    	}
    	if (alg == 2)
    	{
    		return (aa*aa*aa - 3*aa*bb*bb) + a;
    	}
    	else if (alg == 3)
    	{
    		return aa*aa - bb*bb + a;
    	}
    	else if (alg == 4)
    	{
    		return Math.sin(aa*aa - bb*bb) + a;
    	}
    	else if (alg == 5)
    	{
    		return -aa*aa + bb*bb + a;
    	}
    	else if (alg == 6)
    	{
    		return Math.tan(aa*aa - bb*bb) - a;
    	}
    	else
    	{
    		return aa*aa - bb*bb + a;
    	}
    }
    
    private double calcImg(double aa, double bb, double a, double b)
    {
    	if (cx != 0)
    	{
    		a = cx;
    		b = cy;
    	}
    	
    	if (alg == 1)
    	{
    		return 2*aa*bb + b - a;
    	}
    	if (alg == 2)
    	{
    		return (3*aa*aa*bb - bb*bb*bb) + b;
    	}
    	else if (alg == 3)
    	{
    		return Math.atan(2*aa*bb) + b;
    	}
    	else if (alg == 4)
    	{
    		return 2*aa*bb + b;
    	}
    	else if (alg == 5)
    	{
    		return 2*aa*bb + b;
    	}
    	else
    	{
    		return 2*aa*bb + b;
    	}
    	
    }
    
    private void getConstants()
	{
		try
		{
			zoom = Double.parseDouble(tzoom.getText());
			xcen = Double.parseDouble(tX.getText());
			ycen = Double.parseDouble(tY.getText());
			cx = Double.parseDouble(tJuliaX.getText());
			cy = Double.parseDouble(tJuliaY.getText());
			gradient = Double.parseDouble(tContrast.getText());
			brightness = Double.parseDouble(tBrightness.getText());
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
			tzoom.setText("" + zoom);
			tX.setText("" + xcen);
			tY.setText("" + ycen);	
			tJuliaX.setText("" + cx);
			tJuliaY.setText("" + cy);	
			tContrast.setText("" + gradient);
			tBrightness.setText("" + brightness);
			tr.setText("" + depth_red);
			tg.setText("" + depth_green);
			tb.setText("" + depth_blue);
		}
		catch(Exception e)
		{
		}
	}
	
	public void keyPressed(KeyEvent e)
	{
	}
	
	public void keyReleased(KeyEvent e)
	{
		getConstants();
		if (e.getSource() != tContrast && e.getSource() != tBrightness)
			clearScreenAndReset(cx);	
	}
	
	public void keyTyped(KeyEvent e)
	{
	}
    
   
    
}







