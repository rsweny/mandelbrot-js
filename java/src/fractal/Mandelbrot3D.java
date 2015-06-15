package fractal;

//**************************************************************************
// 3D Mandelbrot
// Ryan Sweny
// rsweny@alumni.uwaterloo.ca
// 2009-2011
// v1.1
// based on algorithm by Daniel White (http://www.skytopia.com/project/fractal/mandelbulb.html)
// http://www.fractalforums.com/3d-fractal-generation/amazing-fractal/
//**************************************************************************

import java.awt.*;
import java.awt.image.PixelGrabber;
import java.awt.image.MemoryImageSource;
import java.awt.image.BufferedImage;
import java.awt.event.*;
import java.io.*;
import javax.imageio.ImageIO;
import java.util.Vector;

public class Mandelbrot3D extends java.applet.Applet implements Runnable, KeyListener, ActionListener
{
	//detail level
	double rayDetail = 0.003;
	double stepDetail = 0.028;
	double frost = 1.0;

	//ray traced lighting
	double[] LightVector = { 0.12, 0.15, -0.19 };
	double AMBIENT_LIGHT = 8.0;
	int RAY_STEPS = 25;
	double ray_step;
	double primary_light = 28.0;
	double shadow_darkness = 30.0;
	double HORIZON = 20;

	//fog based on path traces
	double fog_factor = 0.04;
	double min_y, max_y;

	//formula variation
	int formula = 0;
	int inverse_azimuth = 1;
	
	//settings that need saving
	int depth = 20;
	int ximlen = 0;
	int yimlen = 0;
	int pal = 1;
	double power = 8;
	double gradient = 0.6;
	double brightness = 1.5;
	double zoom = 2.4;
	double xcen = 0.0;
	double ycen = 0.0;
	double cameraPersp = -0.1;
	double cameraYaw = 0.1;
	double cameraPitch = 0.75;
	double cameraDOF = 0.0;
	double factorDOF = 0.0;
	double focus_depth = 1.0;
	double opacity = 1.4;
	double focus = -0.15; //higher values for further distance in focus.

	//data
	double root_zoom;
	int half_ximlen, half_yimlen;
	double[][] occlusionPositions = null;
	float[][] img_alpha = null;
	float[][] img_red = null;
	float[][] img_green = null;
	float[][] img_blue = null;
	String lastDir = "";
	private int reset = 0;
	Image img;
	MemoryImageSource screenMem;
	int[] pixels;
	Thread[] runner;
	
	//mouse
	int xcurr,ycurr,xanchor,yanchor;
	boolean m_down;
	boolean drawFocus = false;
	
	
	//3D Camera
	double[][] CameraMatrix = new double[3][3];
	//double[][] InverseCameraMatrix = new double[3][3];
	double[][] IrotX = new double[3][3];
	double[][] IrotZ = new double[3][3];
	double[][] rotX = new double[3][3];
	double[][] rotZ = new double[3][3];
	
	//Stats
	long t1 = 0;
	long visiblePixels, allPixels, rayPoints;
	int renderpass = 0;
	double max_alpha = 1;
	
	//UI
	TextField txtZoom = new TextField(8);
	TextField txtGradient = new TextField(5);
	TextField txtBrightness = new TextField(5);
	TextField txtPower = new TextField(5);
	TextField txtDepth = new TextField(5);
	
	Panel panelIndex = new Panel(new FlowLayout(FlowLayout.LEFT));
	Panel panelMain = new Panel(new BorderLayout());
	Panel panelButtons = new Panel(new FlowLayout(FlowLayout.LEFT));
	
	Panel panel3D = new Panel(new GridLayout(5,4));
	TextField txtPerspective = new TextField(5);
	TextField txtZPos = new TextField(5);
	TextField txtYaw = new TextField(5);
	TextField txtPitch = new TextField(5);
	TextField txtDOF = new TextField(5);
	TextField txtOpacity = new TextField(5);
	TextField txtFocus = new TextField(5);
	TextField txtZRes = new TextField(5);
	TextField txtFog = new TextField(5);
	TextField txtLight = new TextField(5);
	
	Button saveFractal = new Button("Save Fractal");
	Button openFractal = new Button("Open Fractal");
	Button exportPNG = new Button("Export PNG");
	
	boolean running = false;
	
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
	
	private void initVars(boolean fresh)
	{
		if (fresh)
		{
			String szheight = getParameter("height");
			yimlen = Integer.parseInt(szheight);
			
			String szwidth = getParameter("width");
			ximlen = Integer.parseInt(szwidth);
		}
		
		half_ximlen = ximlen / 2;
		half_yimlen = yimlen / 2;

		pixels = new int[ ximlen * yimlen ];

		screenMem = new MemoryImageSource(ximlen, yimlen, pixels, 0, ximlen);
		clearScreenAndReset(fresh);
		setConstants();
		getConstants();
	}
	
	public void init()
	{
		Runtime runtime = Runtime.getRuntime();
        	int nrOfProcessors = runtime.availableProcessors();
		runner = new Thread[nrOfProcessors];
		System.out.println(nrOfProcessors + " processors detected.");
		
		initVars(true);
		
		setBackground(Color.black);
		panelIndex.add(new Label("Zoom"));
		panelIndex.add(txtZoom);
		panelIndex.add(new Label("Contrast"));
		panelIndex.add(txtGradient);
		panelIndex.add(new Label("Brightness"));
		panelIndex.add(txtBrightness);
		panelIndex.add(new Label("N"));
		panelIndex.add(txtPower);
		panelIndex.add(new Label("Iter"));
		panelIndex.add(txtDepth);

		txtZoom.addKeyListener(this);
		txtGradient.addKeyListener(this);
		txtBrightness.addKeyListener(this);
		txtPower.addKeyListener(this);
		txtDepth.addKeyListener(this);

		txtPerspective.addKeyListener(this);
		txtZPos.addKeyListener(this);
		txtPitch.addKeyListener(this);
		txtYaw.addKeyListener(this);
		txtDOF.addKeyListener(this);
		txtOpacity.addKeyListener(this);
		txtFocus.addKeyListener(this);
		txtZRes.addKeyListener(this);
		txtFog.addKeyListener(this);
		txtLight.addKeyListener(this);
		
		panel3D.add(new Label("Perspective"));
		panel3D.add(txtPerspective);
		panel3D.add(new Label("Step"));
		panel3D.add(txtZPos);
		panel3D.add(new Label("Yaw"));
		panel3D.add(txtYaw);
		panel3D.add(new Label("Pitch"));
		panel3D.add(txtPitch);
		panel3D.add(new Label("DOF"));
		panel3D.add(txtDOF);
		panel3D.add(new Label("Focus"));
		panel3D.add(txtFocus);
		panel3D.add(new Label("Light Depth"));
		panel3D.add(txtOpacity);
		panel3D.add(new Label("Frost"));
		panel3D.add(txtZRes);
		panel3D.add(new Label("Glow"));
		panel3D.add(txtFog);
		panel3D.add(new Label("Light Intensity"));
		panel3D.add(txtLight);
		
		
		
		saveFractal.addActionListener(this);
		openFractal.addActionListener(this);
		exportPNG.addActionListener(this);
		
		panelButtons.add(saveFractal);
		panelButtons.add(openFractal);
		panelButtons.add(exportPNG);
		
		panelMain.add(panelIndex, BorderLayout.NORTH);
		panelMain.add(panel3D, BorderLayout.CENTER);
		panelMain.add(panelButtons, BorderLayout.SOUTH);

		panelMain.setBackground(Color.WHITE);
		panelIndex.setBackground(new Color(128,128,128));
		
		add(panelMain);
	}

	public void run()
	{
		while (running)
		{
			try
			{
				for(int y = 0; y < yimlen; y++)
				{
					if (reset == 1)
					{
						reset = 0;
						y = 0;
						clearScreenAndReset(true);
						break;
					}
					else if (reset == 2)
					{
						reset = 0;
						updateHistogram();
						repaint();
					}
					
					while (m_down)
					{
						try { Thread.sleep(1000); } catch(Exception e) {}
					}
					renderpass++;
					
					//local variable to keep consistent across threads
					int currentPass = renderpass;
					
					gridpoints(y);
					if (currentPass % (ximlen*runner.length) == 0) updateMinMaxY();
					
					Thread.sleep(10);
					
					if (currentPass < 20 || currentPass % 50 == 0)
					{
						long t2 = System.currentTimeMillis();
						double completeness = Math.round(max_alpha*100)/100.0;
						String strStatus = visiblePixels + " " + rayPoints + " " + ( allPixels / (ximlen*50)) + " Pass: " + currentPass + " max value: " + completeness + " in " + (t2-t1);
						t1 = System.currentTimeMillis();
						System.out.println(strStatus);
						showStatus(strStatus);
						updateHistogram();
						repaint();

						//autosave png and last 5 states.
						if (lastDir != null && lastDir.length() > 0 && (currentPass % 5000 == 0) )
						{
							writePNG(lastDir + "Mandelbulb-" + currentPass + "-" + completeness + "-" + ximlen + "x" + yimlen + ".png");
							writeFractalData(lastDir + "Mandelbulb-" + (currentPass%30000) + "-" + ximlen + "x" + yimlen + ".fractal");
						}
						
						visiblePixels = 0;
						allPixels = 0;
						rayPoints = 0;
					}

				}
			}
			catch (Exception e)
			{
				System.out.println("Thread error: " + e.toString());
			}
		}
		System.out.println("run() exited");
	}

	public void update(Graphics g)
	{
		paint(g);
	}

	public void updateHistogram()
	{
		int red = 0;
		int green = 0;
		int blue = 0;
		
		max_alpha = Math.pow(FindPeak(img_alpha), gradient);
		
		for (int i=0; i<ximlen; i++)
		{
			for (int j=0; j<yimlen; j++)
			{
				double z = Math.pow(img_alpha[i][j], gradient)*brightness/max_alpha;

				red = (int)( (img_red[i][j]*z)/img_alpha[i][j] );
				green = (int)( (img_green[i][j]*z)/img_alpha[i][j] );
				blue = (int)( (img_blue[i][j]*z)/img_alpha[i][j] );
				
				if (drawFocus && occlusionPositions[i][j] > focus)
				{
					green = 20;
					red = 20;
					blue += 50;
				}
				else if (drawFocus && occlusionPositions[i][j] < focus - focus_depth)
				{
					green = 20;
					blue = 20;
					red += 50;
				}
				
				if (red > 255) red = 255;
				if (green > 255) green = 255;
				if (blue > 255) blue = 255;

				red = red << 16;
				green = green << 8;
				pixels[j*ximlen + i] = 0xff000000 |(red & 0x00ff0000) | (green & 0x0000ff00) | (blue & 0xff);
			}
		}
		img = createImage(screenMem);
	}
	
	private void updateMinMaxY()
	{
		min_y = HORIZON;
		max_y = -HORIZON;
		for (int i=0; i<ximlen; i++)
		{
			for (int j=0; j<yimlen; j++)
			{
				if (occlusionPositions[i][j] < min_y)
				{
					min_y = occlusionPositions[i][j]-stepDetail;
				}
				else if (occlusionPositions[i][j] != HORIZON && occlusionPositions[i][j] > max_y)
				{
					max_y = occlusionPositions[i][j]+stepDetail;
				}
			}
		}
		focus_depth = (max_y - min_y) * 0.25; //0.3
		System.out.println(focus_depth + " Y Bounds: " + min_y + " to " + max_y);
	}
	
	public void paint(Graphics g)
	{
		g.drawImage(img, 0, 0, this);
    	if (m_down) 
        {
			g.setColor(Color.white);
			g.drawRect(xanchor,yanchor,xcurr-xanchor,ycurr-yanchor);
        }
	}
	
	public long FindPeak(float[][] arr)
	{
		double max = 0;
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
		return (long)(max+1);
	}

	public void clearScreenAndReset(boolean fresh)
	{
		renderpass = 0;
		max_alpha = 1;
		min_y = -2.0;
		max_y = 2.0;
		
		if (fresh)
		{
			occlusionPositions = new double[ximlen][yimlen];
			img_alpha = new float[ximlen][yimlen];
			img_red = new float[ximlen][yimlen];
			img_green = new float[ximlen][yimlen];
			img_blue = new float[ximlen][yimlen];
		}

		m_down = false;
	}
	
	private void gridpoints(int yRow)
	{
		//make first pass quick
		int found_limit = (int)Math.ceil(frost);
		if (min_y == -2.0) found_limit = 0;
		
		
		double[][] trace_history = new double[depth+1][5];

		double jitter1 = (0.5 - Math.random()) / ximlen;
		double jitter2 = (0.5 - Math.random()) / yimlen;
		
		//step amount for surface detection
		double stepAmount = (stepDetail + Math.random()*stepDetail)*root_zoom;
		
		double z2 = ((double) yRow / yimlen - 0.5 + jitter1) * zoom - ycen;
		double red, green, blue;
		
		for (int x = 0; x < ximlen; x++)
		{
			//convert screen to fractal coordinates
			double x2 = ((double) x / ximlen - 0.5 + jitter2) * zoom - xcen;
			
			//reset occlusion for this pixel
			occlusionPositions[x][yRow] = HORIZON;
			
			//main loop for y values (depth)
			for (double y = min_y; y < max_y; y += stepAmount)
			{
				allPixels++;
		
				//rotate coordinate system from screen to fractal coordinates
				double persp = 1.0 + y * cameraPersp;
				double fractal_x = x2 / persp;
				double fractal_z = z2 / persp;
				double[] point3D = new double[5];
				rotateVector(point3D, CameraMatrix, fractal_x, y, fractal_z);

				
				//do the calculation
				insideFractal(point3D, trace_history);
				int iter = (int)point3D[3];
				double color = point3D[4];
				

				if (iter == depth)
				{
					//plot pixel
					Vector goodPoints = new Vector();
					double rnd1 = Math.random();
					double light_factor = 1.0 + calculateRays(point3D, rnd1, goodPoints);
					plotPixel(x, yRow, y, color, light_factor);
					visiblePixels++;
					
					//plot more pixels near the surface
					int found = 0;
					int counter = 0;
					double newy = y;
					while (found < found_limit && counter < 200*frost)
					{
						newy = y - stepAmount*Math.random()*2.0;
						persp = 1.0 + newy * cameraPersp;
						fractal_x = x2 / persp;
						fractal_z = z2 / persp;
						rotateVector(point3D, CameraMatrix, fractal_x, newy, fractal_z);
						insideFractal(point3D, null);
						iter = (int)point3D[3];
						color = point3D[4];
						if (iter == depth)
						{
							light_factor = 1.0 + calculateRays(point3D, rnd1, goodPoints);
							plotPixel(x, yRow, y, color, light_factor);
							visiblePixels++;
							found++;
						}
						counter++;
					}
					
					//plot additional points found from ray tracing
					if (frost > 0.9)
					{
						for (int i = 0; i < goodPoints.size(); i++)
						{
							double[] traced_point = (double[])goodPoints.elementAt(i);
							color = traced_point[4];
							light_factor = 1.0 + calculateRays(traced_point, rnd1, null);
							double[] screen_point = reversePoint(traced_point);
							
							//if the ray point is not occluded, draw it
							int tempx = (int)screen_point[0];
							int tempy = (int)screen_point[2];
							
							if ( tempx >= 0 && tempy >= 0 && tempx < ximlen && tempy < yimlen && screen_point[1] < occlusionPositions[tempx][tempy])
							{
								//plot the ray pixel
								plotPixel(tempx, tempy, screen_point[1], color, light_factor);
								rayPoints++;
							}
						}
					}
					break;
				}
				else if (min_y != -2.0) //draw glow after first pass
				{
					//increase step accuracy after first pass and near surface
					double rnd = Math.random();
					stepAmount = (stepDetail + rnd*stepDetail) * ((double)(depth/iter) / (double)depth) * root_zoom * 0.5;
					
					int plot_start = 1;
					if (cameraPersp > 0) plot_start = 4;
					if (fog_factor > 0 && rnd > 0.9 && iter > plot_start)
					{
						//double whiteness = 200 - 10*(depth-iter);
						
						//red = 255 - 255/(iter-2);
						//blue = 255/(iter/2.6);
						//green = Math.max(40, 160 - blue);
						
						blue = 250;
						green = 90 + iter*2;
						red = 20 + iter;
						
						for (int c = 0; c < iter; c++)
						{
							//blue = 255/(c+1);
							//green = Math.max(90, 190 - blue);
							//red = 255 - blue;
							plotFogPixel(trace_history[c], fog_factor, red, green, blue);
						}
					}
				}
			}
		}
	}
	
	private void plotPixel(int x, int yRow, double depth, double color, double light_factor)
	{
		if (cameraDOF > 0)
		{
			double blur_factor = focus-depth;
			
			double[] blurred_point = { x, 0, yRow };
			blur(blurred_point, blur_factor);
			
			int tempx = (int)Math.floor(blurred_point[0]);
			int tempy = (int)Math.floor(blurred_point[2]);
			
			
			if (tempx >= 0 && tempy >=0 && tempx < ximlen && tempy < yimlen)
			{
				if (depth > occlusionPositions[tempx][tempy] + stepDetail) return;
				plotShadowPixel(tempx, tempy, depth, color, light_factor);
			}
		}
		else
		{
			plotShadowPixel(x, yRow, depth, color, light_factor);
		}
	}
	
	private void rotateVector(double[] point3D, double[][] rot, double x, double y, double z)
	{
		point3D[0] = rot[0][0] * x + rot[0][1] * y + rot[0][2] * z;
		point3D[1] = rot[1][0] * x + rot[1][1] * y + rot[1][2] * z;
		point3D[2] = rot[2][0] * x + rot[2][1] * y + rot[2][2] * z;
	}
	
	private void blur(double[] orig, double blur_factor)
	{
		if (blur_factor > 0)
		{
			//increase the range that is in focus
			blur_factor -= focus_depth;
			if (blur_factor < 0) return;
		}
		
		double r2 = Math.random()*2*Math.PI;
		double dsin = Math.sin(r2);
		double dcos = Math.cos(r2);
		
		double dr = Math.random() * factorDOF * blur_factor;
		orig[0] = orig[0] + dr*dcos;
		orig[2] = orig[2] + dr*dsin;
		
		/* square blur
		double dr = (Math.random()-0.5) * cameraDOF * blur_factor;
		x = (int)(x + (dr1*800.0));
		if (x < 0) x = 0;
		else if (x > limit - 1) x = limit - 1;
		return x;
		*/
	}
	
	private void insideFractal(double[] data, double[][] trace_history)
	{
		double magnitude, r, theta_power, r_power, phi, phi_cos, phi_sin;
		int iter = 0;
		double x = data[0];
		double y = data[1];
		double z = data[2];
		
		double pixelColor = 0.0;
		
		//box
		//double scale = 2.0;
		//double fixedRadius = 1;
		//double minRadius = 0.5;
		
		do
		{	
			magnitude = x*x + y*y + z*z;
			r = Math.sqrt(magnitude);
			theta_power = Math.atan2(y,x)*power;
			r_power = Math.pow(r,power);
			
			if (formula == 0)
			{
				//2D compatible / sin
				phi = Math.asin(z / r);
				phi_cos = Math.cos(phi*power);
				x = r_power * Math.cos(theta_power) * phi_cos + data[0];
				y = r_power * Math.sin(theta_power) * phi_cos + data[1];
				z = r_power * Math.sin(phi*power)*inverse_azimuth + data[2];
				pixelColor = phi / 3.0;
			}
			else if (formula == 1)
			{
				//wikipedia / original / cos
				phi = Math.atan2(Math.sqrt(x*x + y*y), z);
				phi_sin = Math.sin(phi*power);
				x = r_power * Math.cos(theta_power) * phi_sin + data[0];
				y = r_power * Math.sin(theta_power) * phi_sin + data[1];
				z = r_power * Math.cos(phi*power)*inverse_azimuth + data[2];
				pixelColor = phi / 3.0;
			}
			/*
			else if (formula == 2)
			{
				//mandelbox
				if (x > 1)
					x = 2 - x;
				else if (x < -1)
					x = -2 - x;

				if (y > 1)
					y = 2 - y;
				else if (y < -1)
					y = -2 - y;

				if (z > 1)
					z = 2 - z;
				else if (z < -1)
					z = -2 - z;

					
				double new_magnitude = x*x + y*y + z*z;
				double length = Math.sqrt(new_magnitude);
				
				double fold = 1.0;
				if (length < minRadius)
					fold = (fixedRadius*fixedRadius) / (minRadius*minRadius);
				else if (length < fixedRadius)
					fold = (fixedRadius*fixedRadius) / (length*length);


				x = x*scale*fold + data[0];
				y = y*scale*fold + data[1];
				z = z*scale*fold + data[2];
				
				pixelColor = (double)iter / (double)depth;
			}
			else
			{
				phi = Math.asin(z / r);
				double sign = 1.0;
				double sign2 = 1.0;

				//rotation around Z axis
				if (x < 0) sign2 = -1.0;
				double tempR = Math.sqrt(x * x + y * y);
				
				//z *= (1.0 / tempR);
				x = x*(1.0 / tempR);
				y = y*(1.0 / tempR);
				z = z*(1.0 / tempR);
				
				double temp = x * x - y * y;
				y = 2.0 * x * y;
				x = temp;
				
				//z *= tempR;
				x = x*(tempR);
				y = y*(tempR);
				z = z*(tempR);

				//rotation around X axis
				if (x < 0) sign = -1.0;
				tempR = Math.sqrt(x * x + z * z);
				
				//z *= (1.0 / tempR);
				x = x*(1.0 / tempR);
				y = y*(1.0 / tempR);
				z = z*(1.0 / tempR);
				
				
				temp = x * x - z * z;
				z = 2.0 * x * z * sign2;
				x = temp * sign;
				
				//z *= tempR;
				x = x*(tempR);
				y = y*(tempR);
				z = z*(tempR);

				
				//z = z * r;
				//z += constant;
				x = x*r + data[0];
				y = y*r + data[1];
				z = z*r + data[2];
				
				pixelColor = phi / 3.0;
			}
			*/
			
			//tower
			/*
			phi = Math.acos(z / r);
			phi_cos = Math.cos(phi*power);
			x = r_power * Math.cos(theta_power) * phi_cos + data[0];
			y = r_power * Math.sin(theta_power) * phi_cos + data[1];
			z = r_power * Math.sin(phi*power)*inverse_azimuth + data[2];
			*/
			
			//used for nebula / fog
			if (trace_history != null)
			{
				trace_history[iter][0] = x;
				trace_history[iter][1] = y;
				trace_history[iter][2] = z;
				trace_history[iter][3] = iter;
				trace_history[iter][4] = pixelColor;
			}
			
			iter++;
		}
		while ( iter < depth && r < 8 );

		data[3] = iter;
		data[4] = pixelColor;
	}
	

	private double calculateRays(double[] origPoint, double rnd, Vector goodPoints)
	{
		double light_factor = 1.0;
		//double rndFuzzy = (0.95 + rnd*0.1)*light_depth;
		
		if (rnd > 0.9) rnd *= 2.0;
		double rndFuzzy = Math.max(opacity*rnd, 0.4);
		
		double n0 = ray_step / (15.0*Math.random() + 1);
		double n1 = ray_step / (15.0*Math.random() + 1);
		double n2 = ray_step / (15.0*Math.random() + 1);
		
		//ambient light
		light_factor += calculateRay(origPoint, RAY_STEPS, -n0, -n1, -n2, AMBIENT_LIGHT, rndFuzzy, goodPoints);
		light_factor += calculateRay(origPoint, RAY_STEPS, n0, n1, n2, AMBIENT_LIGHT, rndFuzzy, goodPoints);
		light_factor += calculateRay(origPoint, RAY_STEPS, n0, -n1, -ray_step, AMBIENT_LIGHT, rndFuzzy, goodPoints);
		light_factor += calculateRay(origPoint, RAY_STEPS, n0*LightVector[0]*9.0,  n1*LightVector[1]*9.0, n2*LightVector[2]*9.0, AMBIENT_LIGHT, rndFuzzy, goodPoints);

		//direct light
		light_factor += calculateRay(origPoint, RAY_STEPS*4, ray_step*LightVector[0],  ray_step*LightVector[1], ray_step*LightVector[2], primary_light, rndFuzzy, goodPoints);
		
		return light_factor;
	}
	
	/*
	private double calculateRays(double[] origPoint, double rnd, Vector goodPoints)
	{
		double light_factor = 1.0;
		//double rndFuzzy = (0.95 + rnd*0.1)*light_depth;
		
		if (rnd > 0.9) rnd += 1.0;
		double rndFuzzy = Math.max(opacity*rnd, 0.5);
		
		double n0 = ray_step / (15.0*Math.random() + 1);
		double n1 = ray_step / (15.0*Math.random() + 1);
		double n2 = ray_step / (15.0*Math.random() + 1);
		
		light_factor += calculateRay(origPoint, RAY_STEPS, n0, n1, n2, AMBIENT_LIGHT, rndFuzzy, goodPoints);
		light_factor += calculateRay(origPoint, RAY_STEPS, -n0, -n1, n2, AMBIENT_LIGHT/2, rndFuzzy, goodPoints);
		light_factor += calculateRay(origPoint, RAY_STEPS, n0, -n1, -ray_step, 5.0, rndFuzzy, goodPoints);
		
		//direct light
		light_factor += calculateRay(origPoint, RAY_STEPS*4, ray_step*LightVector[0],  ray_step*LightVector[1], ray_step*LightVector[2], primary_light, rndFuzzy, goodPoints);
		
		return light_factor;
	}
	*/
	
	private double calculateRay(double[] origPoint, double steps, double stepx, double stepy, double stepz, double bright, double rndFuzzy, Vector goodPoints)
	{
		stepx *= rndFuzzy;
		stepy *= rndFuzzy;
		stepz *= rndFuzzy;

		double x = origPoint[0];
		double y = origPoint[1];
		double z = origPoint[2];
		
		for (int i = 1; i < steps; i++)
		{
			origPoint[0] += stepx*i;
			origPoint[1] += stepy*i;
			origPoint[2] += stepz*i;
			
			insideFractal( origPoint, null );
			if ( origPoint[3] == depth )
			{
				//ray hit solid, this is in shadow
				bright = 0.0;
				if (goodPoints != null /*&& i > 1*/)
				{
					//save points that are moderately close to the known surface
					double[] clone_data = new double[5];
					System.arraycopy(origPoint, 0, clone_data, 0, clone_data.length);
					goodPoints.addElement(clone_data);
				}
				break;
			}
			else if (origPoint[3] < 3)
			{
				//ray has left the general area of the solid, stop tracing.
				break;
			}
		}
		origPoint[0] = x;
		origPoint[1] = y;
		origPoint[2] = z;
		
		return bright;
	}
	

	private void plotShadowPixel(int tempx, int tempy, double depth, double colorVal, double light_factor)
	{
		//get color
		int colorIndex = (int)(Math.abs(colorVal)*255.0);
		if (colorIndex > 255) colorIndex = 255;
		
		//apply lighting
		double red = ((double)Pallet.fpalette[pal][colorIndex][0]*light_factor) / shadow_darkness;
		double green = ((double)Pallet.fpalette[pal][colorIndex][1]*light_factor) / shadow_darkness;
		double blue = ((double)Pallet.fpalette[pal][colorIndex][2]*light_factor) / shadow_darkness;
		
		img_red[tempx][tempy] += red;
		img_green[tempx][tempy] += green;
		img_blue[tempx][tempy] += blue;
		img_alpha[tempx][tempy] += 1.0;
		
		//record depth to mask fog
		occlusionPositions[tempx][tempy] = depth;
	}
	
	private void plotFogPixel(double[] origPoint, double factor, double r, double g, double b)
	{
		double[] screenPoint = reversePoint( origPoint );
		int tempx = (int)screenPoint[0];
		int tempy = (int)screenPoint[2];
		double depth = screenPoint[1];
		

		if (tempx >= 0 && tempy >=0 && tempx < ximlen && tempy < yimlen)
		{
			//solid occludes glow
			if (depth > occlusionPositions[tempx][tempy]) return;
		
			if (cameraDOF > 0)
			{
				double blur_factor = focus-depth;
				blur(screenPoint, blur_factor);
				tempx = (int)screenPoint[0];
				tempy = (int)screenPoint[2];
				if (tempx < 0 || tempy < 0 || tempx > ximlen -1 || tempy > yimlen-1) return;
			}
			img_red[tempx][tempy] += (r*factor);
			img_green[tempx][tempy] += (g*factor);
			img_blue[tempx][tempy] += (b*factor);
			img_alpha[tempx][tempy] += factor;
		}
	}
	
	private double[] reversePoint(double[] fracPoint)
	{
		double[] point3D = new double[5];
		point3D[3] = fracPoint[3];
		point3D[4] = fracPoint[4];
		
		rotateVector(point3D, IrotZ, fracPoint[0], fracPoint[1], fracPoint[2]);
		double xx = point3D[0];
		double yy = point3D[1];
		double zz = point3D[2];
		rotateVector(point3D, IrotX, xx, yy, zz);
		
		double persp = 1.0 + point3D[1] * cameraPersp;
		point3D[0] = point3D[0]*persp;
		point3D[2] = point3D[2]*persp;

		point3D[0] = ((point3D[0] + xcen)/zoom)*ximlen + half_ximlen ;
		point3D[2] = ((point3D[2] + ycen)/zoom)*yimlen + half_yimlen ;
		
		return point3D;
	}
	
	private void getConstants()
	{
		try
		{
			power = Double.parseDouble(txtPower.getText());
			depth = Integer.parseInt(txtDepth.getText());
			brightness = Double.parseDouble(txtBrightness.getText());
			gradient = Double.parseDouble(txtGradient.getText());
			setZoom( Double.parseDouble(txtZoom.getText()) );

			cameraPersp = Double.parseDouble(txtPerspective.getText());
			stepDetail = Double.parseDouble(txtZPos.getText());
			cameraPitch = Double.parseDouble(txtPitch.getText());
			cameraYaw = Double.parseDouble(txtYaw.getText());
			cameraDOF = Double.parseDouble(txtDOF.getText());
			opacity = Double.parseDouble(txtOpacity.getText());
			focus = Double.parseDouble(txtFocus.getText());
			frost = Double.parseDouble(txtZRes.getText());
			fog_factor = Double.parseDouble(txtFog.getText());
			primary_light = Double.parseDouble(txtLight.getText());
			
			setCamera();
		}
		catch(Exception e)
		{
			System.out.println("Problem parsing formula: " + e.toString());
		}
	}
	
	private void setConstants()
	{
		try
		{
			txtPower.setText(""+ power);
			txtDepth.setText(""+ depth);
			txtBrightness.setText(""+ brightness);
			txtGradient.setText(""+ gradient);
			txtZoom.setText(""+ zoom);
			
			txtPerspective.setText("" + cameraPersp);
			txtZPos.setText("" + stepDetail);
			txtPitch.setText("" + cameraPitch);
			txtYaw.setText("" + cameraYaw);
			txtDOF.setText("" + cameraDOF);
			txtOpacity.setText("" + opacity);
			txtFocus.setText("" + focus);
			txtZRes.setText("" + frost);
			txtFog.setText("" + fog_factor);
			txtLight.setText("" + primary_light);
		}
		catch(Exception e)
		{
			System.out.println(e.toString());
		}
	}
	
	private void setCamera()
	{
		// 3D camera precalc
		rotX = RotateX(cameraPitch);
		rotZ = RotateZ(cameraYaw);
		CameraMatrix = matrixMult(rotZ, rotX);

		IrotX = RotateX(-cameraPitch);
		IrotZ = RotateZ(-cameraYaw);
		
		factorDOF = cameraDOF * (ximlen/3);
	}
	
	public double[][] RotateX(double angle)
	{
		double[][] rot = new double[3][3];
		double s = Math.sin(angle);
		double c = Math.cos(angle);
		rot[0][0] = 1.0;
		rot[1][1] = c;
		rot[2][2] = c;
		rot[1][2] = -s;
		rot[2][1] = s;
		return rot;
	}
	
	public double[][] RotateY(double angle)
	{
		double[][] rot = new double[3][3];
		double s = Math.sin(angle);
		double c = Math.cos(angle);
		rot[1][1] = 1.0;
		rot[2][2] = c;
		rot[0][0] = c;
		rot[2][0] = -s;
		rot[0][2] = s;
		return rot;
	}
	
	
	public double[][] RotateZ(double angle)
	{
		double[][] rot = new double[3][3];
		double s = Math.sin(angle);
		double c = Math.cos(angle);
		rot[2][2] = 1.0;
		rot[0][0] = c;
		rot[1][1] = c;
		rot[0][1] = -s;
		rot[1][0] = s;
		return rot;
	}
	
	private double[][] matrixMult(double[][] m, double[][] matrix)
	{
		double[][] result = new double[3][3];
		result[0][0] = m[0][0] * matrix[0][0] + m[0][1] * matrix[1][0] + m[0][2] * matrix[2][0];
		result[0][1] = m[0][0] * matrix[0][1] + m[0][1] * matrix[1][1] + m[0][2] * matrix[2][1];
		result[0][2] = m[0][0] * matrix[0][2] + m[0][1] * matrix[1][2] + m[0][2] * matrix[2][2];
		result[1][0] = m[1][0] * matrix[0][0] + m[1][1] * matrix[1][0] + m[1][2] * matrix[2][0];
		result[1][1] = m[1][0] * matrix[0][1] + m[1][1] * matrix[1][1] + m[1][2] * matrix[2][1];
		result[1][2] = m[1][0] * matrix[0][2] + m[1][1] * matrix[1][2] + m[1][2] * matrix[2][2];
		result[2][0] = m[2][0] * matrix[0][0] + m[2][1] * matrix[1][0] + m[2][2] * matrix[2][0];
		result[2][1] = m[2][0] * matrix[0][1] + m[2][1] * matrix[1][1] + m[2][2] * matrix[2][1];
		result[2][2] = m[2][0] * matrix[0][2] + m[2][1] * matrix[1][2] + m[2][2] * matrix[2][2];
		return result;
	}
	
	
	
	private void setZoom(double z)
	{
		zoom = z;
		
		root_zoom = Math.pow(zoom, 0.5);

		//ray tracing
		ray_step = rayDetail*zoom;
	}
	
	public boolean mouseDown(Event evt, int x, int y)
	{
		xanchor = x;
		yanchor = y;

		if (evt.modifiers == Event.ALT_MASK)
		{
			formula  = (formula + 1)%2;
			if (formula == 0) System.out.println("sin mandelbulb");
			else if (formula == 1) System.out.println("cos mandelbulb");
			else System.out.println("symetric mandelbulb");
			reset = 1;
		}
		else if (evt.modifiers == Event.SHIFT_MASK)
		{
			pal  = (pal + 1)%(Pallet.fpalette.length);
			reset = 1;
		}
		else if (evt.modifiers == Event.CTRL_MASK)
		{
			inverse_azimuth = -inverse_azimuth;
			System.out.println("inverse azimuth is " + inverse_azimuth);
			reset = 1;
		}
		else if (evt.modifiers == Event.META_MASK)
		{
			panelMain.setVisible(!panelMain.isVisible());
		}

		return true;
	}

	public boolean mouseDrag(Event evt, int x, int y)
	{
		m_down = true;
		xcurr = x;
		ycurr = y;
		repaint();
		
		int dx = Math.abs(xcurr - xanchor);
		int dy = Math.abs(ycurr - yanchor);
		double newxcen = xanchor + dx/2.0;
		double fx = ((newxcen - (double)half_ximlen)/(double)ximlen)*zoom;
		
		double newycen = yanchor + dy/2.0;
		double fy = ((newycen - (double)half_yimlen)/(double)yimlen)*zoom;
		
		String strStatus = newxcen + " " + newycen + "( " + fx + ","+fy+")";
		showStatus(strStatus);
		
		return true;
	}

	public boolean mouseUp(Event evt, int x, int y)
	{
		xcurr = x;
		ycurr = y;
		m_down = false;

		int dx = Math.abs(xcurr - xanchor);
		int dy = Math.abs(ycurr - yanchor);
		if (dy > dx)  dx = dy;
		
        //make sure zoom isn't too small
       	if (dx > 10)
       	{
			double newxcen = xanchor + dx/2.0;
			newxcen = ((newxcen - half_ximlen)/(double)ximlen)*zoom;
			xcen = xcen - newxcen;

			double newycen = yanchor + dx/2.0;
			newycen = ((newycen - half_yimlen)/(double)yimlen)*zoom;
			ycen = ycen - newycen;
			
			System.out.println("Xcen is " + xcen + " Ycen is " + ycen);
			setZoom( ((double)dx/(double)(ximlen))*zoom);
			reset = 1;
		}

		setConstants();
		getConstants();
		return true;
	}

	
	public void keyPressed(KeyEvent e)
	{
	}
	
	public void keyReleased(KeyEvent e)
	{
		char c = e.getKeyChar();
		int keyCode = (int)c;
		
		//22: ctrl-v
		//8: backspace
		if (keyCode == 22 || keyCode == 8 || c == ' ' || c == '-' || c == '|' || c == '.' || (c >= '0' && c <= '9') )
		{
		
			if (e.getSource() != txtGradient && e.getSource() != txtBrightness && e.getSource() != txtFocus && e.getSource() != txtDOF && e.getSource() != txtFog && e.getSource() != txtLight && e.getSource() != txtZRes && e.getSource() != txtOpacity && e.getSource() != txtZPos) 
				reset = 1;
			else
				reset = 2;
	
			getConstants();
		}
	}
	
	public void keyTyped(KeyEvent e)
	{
		char c = e.getKeyChar();
		int keyCode = (int)c;
		
		//draw focus border for DOF
		if ( c == 'f')
		{
			drawFocus = !drawFocus;
			reset = 2;
			e.consume();
		}
	}
	
	public void actionPerformed(ActionEvent e)
	{
		m_down = true;
		try
		{
			if (e.getSource() == saveFractal)
			{
				String filename = saveFile("Save Fractal", null, "mandelbrot3d_" + ximlen + "x" + yimlen + ".fractal");
				writeFractalData(filename);
			}
			else if (e.getSource() == openFractal)
			{
				String filename = loadFile("Open Fractal", null, ".fractal");
				readFractalData(filename);
			}
			else
			{
				String filename = saveFile("Export PNG", null, ".png");
				writePNG(filename);
			}
		}
		catch(Exception ex)
		{
		}
		m_down = false;
	}
	
	
	/* File Operations */
	public String loadFile(String title, String defDir, String fileType) throws Exception
	{
		Frame parent = new Frame();
		FileDialog fd = new FileDialog(parent, title, FileDialog.LOAD);
		fd.setFile(fileType);
		fd.setDirectory(defDir);
		fd.setLocation(50, 50);
		fd.show();
		lastDir = fd.getDirectory();
		if (lastDir == null) throw new Exception();
		return lastDir + fd.getFile();
	}

	public String saveFile(String title, String defDir, String fileType) throws Exception
	{
		Frame parent = new Frame();
		FileDialog fd = new FileDialog(parent, title, FileDialog.SAVE);
		fd.setFile(fileType);
		fd.setDirectory(defDir);
		fd.setLocation(50, 50);
		fd.show();
		lastDir = fd.getDirectory();
		if (lastDir == null) throw new Exception();
		return lastDir + fd.getFile();
	} 
	
	public void writePNG(String filename)
	{
		try 
		{
			BufferedImage bi = new BufferedImage(ximlen, yimlen, BufferedImage.TYPE_INT_ARGB); 
			bi.setRGB(0, 0, ximlen, yimlen, pixels, 0, ximlen);
			File outputfile = new File(filename);
			ImageIO.write(bi, "png", outputfile);
			System.out.println("Export PNG success: " + outputfile.getAbsolutePath());
			showStatus("Export PNG success: " + outputfile.getAbsolutePath() );
		} 
		catch (Exception e) 
		{
			System.out.println("export: " + e.toString());
		}
	}
	
	public void writeFractalData(String filename)
	{
		try 
		{
			FileOutputStream fos = new FileOutputStream(filename);
			ObjectOutputStream oos = new ObjectOutputStream(fos);

			oos.writeObject(img_alpha);
			oos.writeObject(img_red);
			oos.writeObject(img_green);
			oos.writeObject(img_blue);
			
			//write int pref array
			int[] intPrefs = { depth, ximlen, yimlen, pal, inverse_azimuth, formula };
			oos.writeObject(intPrefs);
			
			//write double pref array
			double[] doublePrefs = { power, gradient, brightness, zoom, xcen, ycen, cameraPersp, stepDetail, cameraYaw, cameraPitch, cameraDOF, opacity, focus, frost, fog_factor, primary_light };
			oos.writeObject(doublePrefs);
			
			oos.flush();
			fos.close();
			System.out.println("write Fractal Data success: " + filename);
			showStatus("write Fractal Data success: " + filename);
		}
		catch (Throwable e) 
		{
			System.out.println("write: " + e.toString());
		} 
	}
	
	public void readFractalData(String filename)
	{
		try 
		{
			FileInputStream fis = new FileInputStream(filename);
			ObjectInputStream ois = new ObjectInputStream(fis);

			img_alpha = (float[][])ois.readObject();
			img_red = (float[][])ois.readObject();
			img_green = (float[][])ois.readObject();
			img_blue = (float[][])ois.readObject();

			//read int pref array
			int[] intPrefs = (int[])ois.readObject();
			depth = intPrefs[0];
			pal = intPrefs[3];
			
			
			//read double pref array
			double[] doublePrefs = (double[])ois.readObject();
			power = doublePrefs[0];
			gradient = doublePrefs[1];
			brightness = doublePrefs[2];
			setZoom( doublePrefs[3] );
			xcen = doublePrefs[4];
			ycen = doublePrefs[5];
			cameraPersp = doublePrefs[6];
			stepDetail = doublePrefs[7];
			cameraYaw = doublePrefs[8];
			cameraPitch = doublePrefs[9];
			cameraDOF = doublePrefs[10];
			opacity = doublePrefs[11];
			focus = doublePrefs[12];
			frost = doublePrefs[13];
			
			try
			{
				inverse_azimuth= intPrefs[4];
				fog_factor = doublePrefs[14];
				primary_light = doublePrefs[15];
				formula = intPrefs[5];
			}
			catch(Exception e)
			{
			}
			
			
			fis.close();
			
			
			if (yimlen != img_alpha[0].length)
			{
				initVars(true);
			}
			else
			{
				ximlen = intPrefs[1];
				yimlen = intPrefs[2];
				initVars(false);
			}
			System.out.println("open success! " + filename);
		}
		catch (Throwable e) 
		{
			System.out.println("read: " + e.toString());
		} 
	}
	
	
	
	
	
		/* Matrix Ops */
	public double determinant(double[][] mat) 
	{
		double result = 0;
		
		if(mat.length == 1) {
			result = mat[0][0];
			return result;
		}

		if(mat.length == 2) {
			result = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
			return result;
		} 

		for(int i = 0; i < mat[0].length; i++) 
		{
			double temp[][] = new double[mat.length - 1][mat[0].length - 1];

			for(int j = 1; j < mat.length; j++)
			{
				for(int k = 0; k < mat[0].length; k++) 
				{
					if(k < i) {
					temp[j - 1][k] = mat[j][k];
					} else if(k > i) {
					temp[j - 1][k - 1] = mat[j][k];
					}
				}
			}

			result += mat[0][i] * Math.pow(-1, (double)i) * determinant(temp);
		}
		return result;
	} 
	
	public double[][] cofactor3x3T(double[][] m) 
	{
		double temp[][] = new double[m[0].length][m[0].length];
		temp[0][0] = m[1][1]*m[2][2] - m[1][2]*m[2][1];
		temp[1][0] = m[0][2]*m[2][1] - m[0][1]*m[2][2];
		temp[2][0] = m[0][1]*m[1][2] - m[0][2]*m[1][1];
		temp[0][1] = m[1][2]*m[2][0] - m[1][0]*m[2][2];
		temp[1][1] = m[0][0]*m[2][2] - m[0][2]*m[2][0];	
		temp[2][1] = m[0][2]*m[1][0] - m[0][0]*m[1][2];
		temp[0][2] = m[1][0]*m[2][1] - m[1][1]*m[2][0];
		temp[1][2] = m[0][1]*m[2][0] - m[0][0]*m[2][1];
		temp[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0];
		return temp;
	} 
	
	public double[][] inverse(double[][] m)
	{
		double d = determinant(m);
		System.out.println("det: " + d);
		if (d == 0) d = 0.00001;  //dubious
		
		System.out.println("Determinant: " + d);
		
		double[][] cofactorT = cofactor3x3T(m);
		d = 1/d;
		constMul(d, cofactorT);
		return cofactorT;
	}
	
	public void constMul(double d, double[][] m)
	{
		m[0][0] *= d;
		m[0][1] *= d;
		m[0][2] *= d;
		m[1][0] *= d;
		m[1][1] *= d;
		m[1][2] *= d;
		m[2][0] *= d;
		m[2][1] *= d;
		m[2][2] *= d;
	}

}


