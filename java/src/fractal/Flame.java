package fractal;

//**************************************************************************
// 3D Flame Fractal Generator
// Ryan Sweny
// rsweny@alumni.uwaterloo.ca
// 2007
// v2.0
// based on Apophysis (http://sourceforge.net/projects/apophysis/)
//**************************************************************************

import java.awt.*;
import java.awt.image.PixelGrabber;
import java.awt.image.MemoryImageSource;
import java.awt.event.*;
import java.util.StringTokenizer;


public class Flame extends java.applet.Applet implements Runnable, KeyListener
{
	private boolean reset = false;
	int CHOOSE_XFORM_GRAIN = 1000;
	double EPS = 1e-12;
	int FUSE = 50;

	XForm[] xp = new XForm[XForm.NVARS];
	Image img;
	MemoryImageSource screenMem;
	int[] pixels;
	Thread runner;

	double zoom = 1;
	double xcen,ycen;
	double gradient = 0.3;
	double brightness = 1.5;

	int ximlen, yimlen;
	double[][] zPositions = null;
	double[][] img_alpha = null;
	double[][] img_red = null;
	double[][] img_green = null;
	double[][] img_blue = null;

	double[] point = new double[4];
	int iter;
	int active = 2;
	int pal = 9;
	int alg = XForm.NVARS-1;
	long renderpass = 0;
	double max_alpha;

	int xcurr,ycurr,xanchor,yanchor;
	boolean m_down;
	
	//3D
	double[][] CameraMatrix = new double[3][3];
	double[] gauss_rnd = new double[4];
	int gauss_N = 0;
	double cameraPersp = 0.9;
	double cameraZpos = 0.0;
	double cameraYaw = 0.3;
	double cameraPitch = 0.86;
	double cameraDOF = 0.01;
	double opacity = 0.9;
	double translucency = 0.1;
	double focus = 0.0;
	double horizon = 1.0;
	double zres = 0.96;
	double newPixels, visiblePixels, occludedPixels;
	
	//UI
	TextField txtZoom = new TextField(5);
	TextField txtGradient = new TextField(5);
	TextField txtBrightness = new TextField(5);
	Panel panelIndex = new Panel(new FlowLayout(FlowLayout.LEFT));
	TextArea txtFormula = new TextArea();
	Panel panelMain = new Panel(new BorderLayout());
	Panel panel3D = new Panel(new GridLayout(2,9));
	TextField txtPerspective = new TextField(5);
	TextField txtZPos = new TextField(5);
	TextField txtYaw = new TextField(5);
	TextField txtPitch = new TextField(5);
	TextField txtDOF = new TextField(5);
	TextField txtOpacity = new TextField(5);
	TextField txtReflect = new TextField(5);
	TextField txtFocus = new TextField(5);
	TextField txtHorizon = new TextField(5);
	TextField txtZRes = new TextField(5);
	
	public void start()
	{
		iter = 80000;
		if (runner == null);
		{
			runner = new Thread(this);
			runner.start();
		}
	}


	public void stop()
	{
		iter = 0;
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
		yimlen = Integer.parseInt(szheight);

		String szwidth = getParameter("width");
		ximlen = Integer.parseInt(szwidth);
		
		pixels = new int[ ximlen * yimlen ];
		point[3] = 0.0;
		screenMem = new MemoryImageSource(ximlen, yimlen, pixels, 0, ximlen);

		
		for (int i = 0; i < XForm.NVARS; i++) xp[i] = new XForm();
		clearScreenAndReset();
		
		loadDefault();
		setConstants();
			
		String formula = getParameter("formula");
		System.out.println(formula);
		if (formula != null && formula.length() > 1) 
		{
			txtFormula.setText(formula);
			getConstants();
			setConstants();
		}
		setCamera();

		setBackground(Color.black);
		panelIndex.add(new Label("Zoom"));
		panelIndex.add(txtZoom);
		panelIndex.add(new Label("Contrast"));
		panelIndex.add(txtGradient);
		panelIndex.add(new Label("Brightness"));
		panelIndex.add(txtBrightness);

		txtZoom.addKeyListener(this);
		txtGradient.addKeyListener(this);
		txtBrightness.addKeyListener(this);
		txtFormula.addKeyListener(this);
		
		txtPerspective.addKeyListener(this);
		txtZPos.addKeyListener(this);
		txtPitch.addKeyListener(this);
		txtYaw.addKeyListener(this);
		txtDOF.addKeyListener(this);
		txtOpacity.addKeyListener(this);
		txtReflect.addKeyListener(this);
		txtFocus.addKeyListener(this);
		txtHorizon.addKeyListener(this);
		txtZRes.addKeyListener(this);
		
		txtFormula.setFont( new Font("Arial", Font.PLAIN, 12) );
		txtFormula.setColumns(80);
		
		panel3D.add(new Label("Perspective"));
		panel3D.add(txtPerspective);
		panel3D.add(new Label("Z"));
		panel3D.add(txtZPos);
		panel3D.add(new Label("Yaw"));
		panel3D.add(txtYaw);
		panel3D.add(new Label("Pitch"));
		panel3D.add(txtPitch);
		panel3D.add(new Label("DOF"));
		panel3D.add(txtDOF);
		panel3D.add(new Label("Focus"));
		panel3D.add(txtFocus);
		panel3D.add(new Label("Opacity"));
		panel3D.add(txtOpacity);
		panel3D.add(new Label("Glow"));
		panel3D.add(txtReflect);
		panel3D.add(new Label("Horizon"));
		panel3D.add(txtHorizon);
		panel3D.add(new Label("Foreground"));
		panel3D.add(txtZRes);
		
		panelMain.add(panelIndex, BorderLayout.NORTH);
		panelMain.add(panel3D, BorderLayout.CENTER);
		panelMain.add(txtFormula, BorderLayout.SOUTH);

		panelMain.setBackground(Color.WHITE);
		panelIndex.setBackground(new Color(128,128,128));
		
		add(panelMain);

		System.gc();
	}

	public void run()
	{
		Thread thisThread = Thread.currentThread();
		thisThread.setPriority(Thread.MIN_PRIORITY);
		while ( runner == thisThread )
		{
			try
			{
				if (reset)
				{
					reset = false;
					clearScreenAndReset();
				}

				if (!m_down) 
				{
					nextpoints();
					
					if (renderpass%2 == 0) 
					{
						String strStatus = visiblePixels + " " + occludedPixels + " Pass: " + renderpass + " max value: " + (Math.round(max_alpha*100)/100.0);
						System.out.println(strStatus);
						showStatus(strStatus);
					
						updateHistogram();
						repaint();
						newPixels = 0;
						visiblePixels = 0;
						occludedPixels = 0;
						
					}
					Thread.sleep(100);
				}
					
				
			}
			catch (Exception e)
			{
				System.out.println(e.toString());
				e.printStackTrace();
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
		
		if (renderpass%4 == 2)
			max_alpha = Math.pow(FindPeak(img_alpha), gradient);
		
		for (int i=0; i<ximlen; i++)
		{
			for (int j=0; j<yimlen; j++)
			{
				double z = Math.pow(img_alpha[i][j], gradient)*brightness/max_alpha;

				red = (int)( (img_red[i][j]*z)/img_alpha[i][j] );
				green = (int)( (img_green[i][j]*z)/img_alpha[i][j] );
				blue = (int)( (img_blue[i][j]*z)/img_alpha[i][j] );
				
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


	public void paint(Graphics g)
	{
		g.drawImage(img, 0, 0, this);
    	if (m_down) 
        {
			g.setColor(Color.white);
			g.drawRect(xanchor,yanchor,xcurr-xanchor,ycurr-yanchor);
        }
	}
	
	public long FindPeak(double[][] arr)
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

	public void clearScreenAndReset()
	{
		renderpass = 0;
		img_alpha = new double[ximlen][yimlen];
		img_red = new double[ximlen][yimlen];
		img_green = new double[ximlen][yimlen];
		img_blue = new double[ximlen][yimlen];
		zPositions = new double[ximlen][yimlen];
		for (int i = 0; i < ximlen; i++)
		{
			for (int j = 0; j < yimlen; j++)
			{
				zPositions[i][j] = -100;
				img_alpha[i][j] = 1;
			}
		}
		
		m_down = false;
	}


	private void nextpoints()
	{
		renderpass++;
		int i,j,k,fn;

		double a, c1, c2, tx, ty, tz, v, nx, ny, nz, r2, t, r, dr, theta;
		double curlx = Math.random();	
		double curly = Math.random();	
		double curlz = Math.random();	
		double curlc2 = curlx*curlx + curly*curly + curlz*curlz;
		
		int[] xform_distrib = new int[CHOOSE_XFORM_GRAIN];


		//----- "set up xform_distrib[], which is an array that converts -----
		//-----    a uniform random variable into one with the distribution -----
		//-----    dictated by the density field -----
		dr = 0.0;
		for (j = 0; j < XForm.NVARS; j++) dr += xp[j].density;
		dr = dr / CHOOSE_XFORM_GRAIN;
		j = 0;
		t = xp[0].density;
		r = 0.0;
		for (k = 0; k < CHOOSE_XFORM_GRAIN; k++)
		{
			while (r >= t)
			{
				j++;
				t += xp[j].density;
			}
			xform_distrib[k] = j;
			r += dr;
		}
		
		//----- iteration loop -----
		for (i = -FUSE; i < iter; i++)
		{
			//----- check for explosion -----
			if ((Math.abs(point[0]) > 1E6) || (Math.abs(point[1]) > 1E6) || (Math.abs(point[2]) > 1E6))
			{
				System.out.println("point out of bounds");
				point[0] = point[1] = point[2] = 0.0;
				break;
			}

			fn = random_distrib(xform_distrib);
			
			
			double rings_dx = xp[fn].c[2][0] * xp[fn].c[2][0] + EPS;
			double fan_dx = Math.PI * xp[fn].c[2][0] * xp[fn].c[2][0] + EPS;
			double fan_dx2 = fan_dx/2;
			
			/* first compute the color coord */
			double s = xp[fn].symmetry;
			point[3] = (point[3] + xp[fn].color) * 0.5 * (1.0 - s) + s * point[3];

			//----- "... then apply the affine part of the function..." -----
			tx = xp[fn].c[0][0] * point[0] + xp[fn].c[1][0] * point[1] + xp[fn].c[2][0];
			ty = xp[fn].c[0][1] * point[0] + xp[fn].c[1][1] * point[1] + xp[fn].c[2][1];
			tz = point[2];
			
			

			
			point[0] = point[1] = point[2] = 0.0;

			//----- "... then add proportional amounts of each -----
			//-----    variation" (S.D.) -----


			//PRE functions that work on tx,ty,tz
			v = xp[fn].var[26];
			if (v != 0.0)
			{
				//PreBlur
				double angle = Math.random()*2*Math.PI;
				double sina = Math.sin(angle);
				double cosa = Math.cos(angle);
				double mag = v * (gauss_rnd[0] + gauss_rnd[1] + gauss_rnd[2] + gauss_rnd[3] - 2);
				gauss_rnd[gauss_N] = Math.random();
				gauss_N = (gauss_N+1) % 4;
				
				tx += mag * cosa;
				ty += mag * sina;
			}

			//NORMAL functions that work on the point
			
			//----- LINEAR -----
			v = xp[fn].var[0];
			if (v > 0.0)
			{
				point[0] += v * tx;
				point[1] += v * ty;
				point[2] += v * tz;
			}

			//----- SINUSOIDAL -----
			v = xp[fn].var[1];
			if (v > 0.0)
			{
				point[0] += v * Math.sin(tx);
				point[1] += v * Math.sin(ty);
			}

			//----- SPHERICAL -----
			v = xp[fn].var[2];
			if (v > 0.0)
			{
				r2 = tx * tx + ty * ty + 1e-6;
				nx = tx / r2;
				ny = ty / r2;
				point[0] += v * nx;
				point[1] += v * ny;
			}

			//----- SWIRL -----
			v = xp[fn].var[3];
			if (v > 0.0)
			{
				r2 = tx * tx + ty * ty;    // "k here is fun" (S.D.)
				c1 = Math.sin(r2);
				c2 = Math.cos(r2);
				nx = c1 * tx - c2 * ty;
				ny = c2 * tx + c1 * ty;
				point[0] += v * nx;
				point[1] += v * ny;
			}

			//----- HORSESHOE -----
			v = xp[fn].var[4];
			if (v > 0.0)
			{
				if (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS)
					a = Math.atan2(tx, ty);      // "times k here is fun" (S.D.)
				else
					a = 0.0;
				c1 = Math.sin(a);
				c2 = Math.cos(a);
				nx = c1 * tx - c2 * ty;
				ny = c2 * tx + c1 * ty;
				point[0] += v * nx;
				point[1] += v * ny;
			}

			//----- POLAR -----
			v = xp[fn].var[5];
			if (v > 0.0)
			{
				if (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS)
					nx = Math.atan2(tx, ty) / Math.PI;
				else
					nx = 0.0;
				ny = Math.sqrt(tx * tx + ty * ty) - 1.0;
				point[0] += v * nx;
				point[1] += v * ny;
			}

			//----- Folded handkerchief -----
			v = xp[fn].var[6];
			if (v > 0.0)
			{
				if (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS)
					a = Math.atan2(tx, ty);
				else
					a = 0.0;

				r2 = Math.sqrt(tx*tx + ty*ty);
				point[0] += v * Math.sin(a+r) * r2;
				point[1] += v * Math.cos(a-r) * r2;
			}

			//----- Heart -----
			v = xp[fn].var[7];
			if (v > 0.0) {
				if (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS)
					a = Math.atan2(tx, ty);
				else
					a = 0.0;
				r2 = Math.sqrt(tx*tx + ty*ty);
				a *= r2;
				point[0] += v * Math.sin(a) * r2;
				point[1] += v * Math.cos(a) * (-r2);
			}

			//----- Disc -----
			v = xp[fn].var[8];
			if (v > 0.0)
			{
				nx = tx * Math.PI;
				ny = ty * Math.PI;
				if (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS)
					a = Math.atan2(nx, ny);
				else
					a = 0.0;
				r2 = Math.sqrt(nx*nx + ny*ny);
				point[0] += v * Math.sin(r2) * a / Math.PI;
				point[1] += v * Math.cos(r2) * a / Math.PI;
			}

			//----- Spiral -----
			v = xp[fn].var[9];
			if (v > 0.0)
			{
				if (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS)
					a = Math.atan2(tx, ty);
				else
					a = 0.0;
				r2 = Math.sqrt(tx*tx + ty*ty) + 1e-6;
				point[0] += v * (Math.cos(a) + Math.sin(r2)) / r2;
				point[1] += v * (Math.sin(a) - Math.cos(r2)) / r2;
			}



			//----- Hyperbolic -----
			v = xp[fn].var[10];
			if (v > 0.0)
			{
				if (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS)
					a = Math.atan2(tx, ty);
				else
					a = 0.0;
				r2 = Math.sqrt(tx*tx + ty*ty) + 1e-6;
				point[0] += v * Math.sin(a) / r2;
				point[1] += v * Math.cos(a) * r2;
			}

			//----- Square / Diamond -----
			v = xp[fn].var[11];
			if (v > 0.0)
			{
				if (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS)
					a = Math.atan2(tx, ty);
				else
					a = 0.0;
				r2 = Math.sqrt(tx*tx + ty*ty);
				point[0] += v * Math.sin(a) * Math.cos(r2);
				point[1] += v * Math.cos(a) * Math.sin(r2);
			}


			//----- Ex -----
			v = xp[fn].var[12];
			if (v > 0.0)
			{
				double n0, n1, m0, m1;
				if (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS)
					a = Math.atan2(tx, ty);
				else
					a = 0.0;
				r2 = Math.sqrt(tx*tx + ty*ty);
				n0 = Math.sin(a+r2);
				n1 = Math.cos(a-r2);
				m0 = n0 * n0 * n0 * r;
				m1 = n1 * n1 * n1 * r;
				point[0] += v * (m0 + m1);
				point[1] += v * (m0 - m1);
			}

			//----- Fan 3D -----
			/*
			v = xp[fn].var[13];
			if (v > 0.0)
			{
				double FAngle = Math.atan2(tx, ty);
				double angle;
				if ( (FAngle + xp[fn].c[2][1]) / fan_dx > 0.5 )
					angle = FAngle - fan_dx2;
				else
					angle = FAngle + fan_dx2;
					
				double sinr = Math.sin(angle);
				double cosr = Math.cos(angle);
				r = v*Math.sqrt(tx*tx + ty*ty);

				point[0] += r*cosr;
				point[1] += r*sinr;
			}
			*/
			
			//----- Julia -----
			v = xp[fn].var[13];
			if (v > 0.0)
			{
				if (tx < -EPS || tx > EPS || ty < -EPS || ty > EPS)
				a = Math.atan2(tx, ty)/2.0;
				else
				a = 0.0;

				if (Math.random() > 0.5) a += Math.PI;

				r2 = Math.pow(tx*tx + ty*ty, 0.25);
				nx = r2 * Math.cos(a);
				ny = r2 * Math.sin(a);
				point[0] += v * nx;
				point[1] += v * ny;
			}

			//----- BENT -----
			v = xp[fn].var[14];
			if (v > 0.0)
			{
				nx = tx;
				ny = ty;
				if (nx < 0.0) nx = nx * 2.0;
				if (ny < 0.0) ny = ny / 2.0;
				point[0] += v * nx;
				point[1] += v * ny;
			}


			v = xp[fn].var[15];
			if (v != 0.0)
			{
				// waves
				double dx = xp[fn].c[2][0];
				double dy = xp[fn].c[2][1];
				nx = tx + xp[fn].c[1][0]*Math.sin(ty/((dx*dx)+EPS));
				ny = ty + xp[fn].c[1][1]*Math.sin(tx/((dy*dy)+EPS));
				point[0] += v * nx;
				point[1] += v * ny;
			}

			v = xp[fn].var[16];
			if (v != 0.0)
			{
				// fisheye
				r = Math.sqrt(tx*tx + ty*ty);
				a = Math.atan2(tx, ty);
				r = 2 * r / (r + 1);
				nx = r * Math.cos(a);
				ny = r * Math.sin(a);
				point[0] += v * nx;
				point[1] += v * ny;
			}

			v = xp[fn].var[17];
			if (v != 0.0)
			{
				// popcorn
				double dx = Math.tan(3*ty);
				double dy = Math.tan(3*tx);

				nx = tx + xp[fn].c[2][0] * Math.sin(dx);
				ny = ty + xp[fn].c[2][1] * Math.sin(dy);
				point[0] += v * nx;
				point[1] += v * ny;
			}
			
			//3D functions
			v = xp[fn].var[18];
			if (v != 0.0)
			{
				//Z Translate
				point[2] += v;
			}
			
			v = xp[fn].var[19];
			if (v != 0.0)
			{
				//Z Scale
  				point[2] = point[2] + point[2]*v;
			}
			
			v = xp[fn].var[20];
			if (v != 0.0)
			{
				//Blur3D
				double angle = Math.random()*2*Math.PI;
				double sina = Math.sin(angle);
				double cosa = Math.cos(angle);
				double mag = v * (gauss_rnd[0] + gauss_rnd[1] + gauss_rnd[2] + gauss_rnd[3] - 2);
				gauss_rnd[gauss_N] = Math.random();
				gauss_N = (gauss_N+1) % 4;
				
				double angle2 = Math.random()*Math.PI;
				double sinb = Math.sin(angle2);
				double cosb = Math.cos(angle2);
				
				point[0] +=  mag * sinb * cosa;
				point[1] +=  mag * sinb * sina;
				point[2] +=  mag * cosb;
			}
			
			v = xp[fn].var[21];
			if (v != 0.0)
			{
				//Z blur
				point[2] = point[2] + v * (gauss_rnd[0] + gauss_rnd[1] + gauss_rnd[2] + gauss_rnd[3] - 2);
				gauss_rnd[gauss_N] = Math.random();
				gauss_N = (gauss_N+1) % 4;
			}
			
			v = xp[fn].var[22];
			if (v != 0.0)
			{
				//Zcone
				point[2] += v * Math.sqrt(tx*tx + ty*ty);
				//point[2] += v * Math.sqrt(point[0]*point[0] + point[1]*point[1]);
			}
			
			v = xp[fn].var[23];
			if (v != 0.0)
			{
				//Bubble
				double mag = (tx*tx + ty*ty)/4 + 1;
				point[2] = point[2] + v * (2 / mag - 1);
				
				mag = v / mag;
				
				point[0] = point[0] + mag * tx;
				point[1] = point[1] + mag * ty;
				
				/*t = 4 / (tx*tx + ty*ty + 4);
				point[0] = point[0] + tx * t * xp[fn].var[19];
				point[1] = point[1] + ty * t * xp[fn].var[19];
				point[2] = point[2] + (2 * t - 1) * xp[fn].var[19];*/
			}
			
			v = xp[fn].var[24];
			if (v != 0.0)
			{
				//Julia 3D
				double r2d = Math.sqrt(tx*tx + ty*ty);
				double mag = v * Math.sqrt(r2d);
				
				point[2] = point[2] + mag * tz / r2d / 2;
				
				int random2 = (int)(Math.random()*2);
				double angle = (Math.atan2(ty, tx)/2 + Math.PI*random2);
				double sina = Math.sin(angle);
				double cosa = Math.cos(angle);
				
				point[0] = point[0] + mag * cosa;
				point[1] = point[1] + mag * sina;
			}
			
			v = xp[fn].var[25];
			if (v != 0.0)
			{
				//Curl 3D
				double curlr2 = tx*tx + ty*ty + tz*tz;
				double mag = v / (curlr2*curlc2 + curlx*2*tx - curly*2*ty + curlz*2*tz + 1);
				
				point[0] += mag * (tx + curlx*curlr2);
				point[1] += mag * (ty - curly*curlr2);
				point[2] += mag * (tz + curlz*curlr2);
			}
			
			
			double[] projectedPoint = (double[])point.clone();
			projectPitchYawDOF(projectedPoint);

			double realx = projectedPoint[0] + xcen;
			double realy = projectedPoint[1] + ycen;
			if (realx > -zoom && realy > -zoom && realx < zoom && realy < zoom && i >= 0)
			{
				//-- Update the histogram --
				int tempx = (int)Math.floor ( (realx/zoom)*(ximlen/2) + (ximlen/2) );
				int tempy = (int)Math.floor ( (realy/zoom)*(ximlen/2) + (yimlen/2) );
				int colorIndex = (int)point[3];
				
				//take into account visibility from camera
				double horizon_dist = projectedPoint[2] + Math.random() / 2;
				if (horizon_dist < -horizon) continue;

				if (projectedPoint[2] > zPositions[tempx][tempy])
				{	
					//this point is visible
					zPositions[tempx][tempy] = projectedPoint[2];
					img_red[tempx][tempy] += Pallet.fpalette[pal][colorIndex][0]*opacity;
					img_green[tempx][tempy] += Pallet.fpalette[pal][colorIndex][1]*opacity;
					img_blue[tempx][tempy] += Pallet.fpalette[pal][colorIndex][2]*opacity;
					img_alpha[tempx][tempy] += opacity;
					visiblePixels++;
				}
				else
				{
					//this point is occluded
					double totalColor = (Pallet.fpalette[pal][colorIndex][0] + Pallet.fpalette[pal][colorIndex][1] + Pallet.fpalette[pal][colorIndex][2]);
					double intensity = totalColor / 765.0;
					
					zPositions[tempx][tempy] = zPositions[tempx][tempy]*zres + projectedPoint[2]*(1.0 - zres);
					img_red[tempx][tempy] += Pallet.fpalette[pal][colorIndex][0]*(1.0-opacity+translucency);
					img_green[tempx][tempy] += Pallet.fpalette[pal][colorIndex][1]*(1.0-opacity+translucency);
					img_blue[tempx][tempy] += Pallet.fpalette[pal][colorIndex][2]*(1.0-opacity+translucency);
					img_alpha[tempx][tempy] += (1.0-opacity+(translucency/2))*intensity;
					occludedPixels++;
				}
			}
		}
	}
	
	private void projectPitchYawDOF(double[] point)
	{
		double x, y, z, zr, dr, dsin, dcos;
		
		z = point[2] - cameraZpos;
		x = CameraMatrix[0][0]*point[0] + CameraMatrix[1][0]*point[1];
		y = CameraMatrix[0][1]*point[0] + CameraMatrix[1][1]*point[1] + CameraMatrix[2][1]*z;
		z = CameraMatrix[0][2]*point[0] + CameraMatrix[1][2]*point[1] + CameraMatrix[2][2]*z;
		zr = 1 - cameraPersp * z;
		
		double r2 = Math.random()*2*Math.PI;
		dsin = Math.sin(r2);
		dcos = Math.cos(r2);
		
		dr = Math.random() * cameraDOF * (z+focus);
		
		point[0] = (x + dr*dcos) / zr;
		point[1] = (y + dr*dsin) / zr;
		point[2] = z + dr*dsin; 
	}
	
	
	/* 
	sym=2 or more means rotational
	sym=1 means identity, ie no symmetry
	sym=0 means pick a random symmetry (maybe none)
	sym=-1 means bilateral (reflection)
	sym=-2 or less means rotational and reflective
	*/
	public int add_symmetry_to_control_point()
	{
		int i, j, k;
		double a;
		int result = 0;

		int sym_distrib[] = {
			-4,
			-3, -3,
			-2, -2, -2, -2,
			-1, -1, -1, -1,
			2, 2, 2, 2,
			3, 3,
			4,
		};

		int sym = random_distrib(sym_distrib);

		System.out.println("sym is " + sym);

		for (i = 0; i < XForm.NVARS; i++)
		if (xp[i].density == 0.0)
		break;
		if (i == XForm.NVARS) return 0;

		if (sym < 0)
		{
			xp[i].density = 1.0;
			xp[i].symmetry = 1.0;
			xp[i].var[0] = 1.0;

			xp[i].color = 1.0;
			xp[i].c[0][0] = -1.0;
			xp[i].c[0][1] = 0.0;
			xp[i].c[1][0] = 0.0;
			xp[i].c[1][1] = 1.0;
			xp[i].c[2][0] = 0.0;
			xp[i].c[2][1] = 0.0;

			i++;
			result++;
			sym = -sym;
		}

		a = 2*Math.PI/sym;

		for (k = 1; (k < sym) && (i < XForm.NVARS); k++)
		{
			xp[i].density = 1.0;
			xp[i].var[0] = 1.0;
			xp[i].symmetry = 1.0;

			xp[i].color = (sym<3) ? 0.0 : ((k-1.0)/(sym-2.0));

			xp[i].c[0][0] = Math.cos(k*a);
			xp[i].c[0][1] = Math.sin(k*a);
			xp[i].c[1][0] = -xp[i].c[0][1];
			xp[i].c[1][1] = xp[i].c[0][0];
			xp[i].c[2][0] = 0.0;
			xp[i].c[2][1] = 0.0;

			i++;
			result++;
		}

		return result;
	}





	public void random_xform(int ivar)
	{
		int i, j, k, var;

		//----- number of xforms to be used in the current image -----
		int xform_distrib[] = {1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6};

		//----- which formula(e) to use in the active xforms -----
		//----- -1 means a random mix; otherwise all xforms will -----
		int var_distrib[] = {24, 24, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};

		//----- weighted choices for a random mix of formulae -----
		int mixed_var_distrib[] = {24, 24, 0, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 14, 14, 15, 16, 17};

		active = random_distrib(xform_distrib);
		var = (ivar < 0) ? random_distrib(var_distrib) : ivar;
		for (i = 0; i < active; i++)
		{
			xp[i].symmetry = 0;
			xp[i].density = 1.0 / (double)active;
			xp[i].color = (int)(((double)i / (double)active)*255);

			for (j = 0; j < 3; j++)
			for (k = 0; k < 2; k++)
			xp[i].c[j][k] = qsrandom11();
			for (j = 0; j < XForm.NVARS; j++)
			xp[i].var[j] = 0.0;
			if (var >= 0) xp[i].var[var] = 1.0;
			else xp[i].var[random_distrib(mixed_var_distrib)] = 1.0;
		}
		for (; i < XForm.NVARS; i++) xp[i].density = 0.0;
	}

	public void predefined_xform(int var)
	{
		if (var == 0)
		{
			var = 1;
			alg++;
		}

		for (int i = 0; i < XForm.NVARS; i++) xp[i].density = 0.0;

		double[] v = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		double[][] d = { {-0.99, -0.05861629} , {0.806310617, -0.688782922} , {-0.368144780, 0.242713706} };
		loadXForm(0, 0.5, 0, 255, v, d);

		double[] v1 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		v1[var] = 1.0;
		double[][] d1 = { {0.191625721, 0.724478896} , {-0.569627979, 0.390423292} , {-0.021881771, -0.208493301} };
		loadXForm(1, 0.5, 0, 0, v1, d1);
	}


	public void loadXForm(int i, double weight, double sym, double color, double[] vars, double[][] coefs)
	{
		xp[i].symmetry = sym;
		xp[i].density = weight;
		xp[i].color = color;
		xp[i].var = vars;
		xp[i].c = coefs;
	}

	public void loadDefault()
	{
		for (int i = 0; i < XForm.NVARS; i++) xp[i].density = 0.0;
		
		String julia = "0.31 0.0 30.0 1.0 0.1 0.1 1.0 0.1 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.2 0.0 0.0 0.19 |5.1 0.4 0.0 0.452548 -0.452548 0.452548 0.452548 0.19 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.1 0.0 0.0 0.0 0.0 0.0 0.9 0.0 0.0 |5.8 0.1 255.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.5 0.0 0.0 |";
		String spiral = "0.25 0.0 0.0 -0.0454 -0.606852 0.620299 -0.212502 0.784553 -0.258353 0.0 0.5 0.0 0.0 0.01 0.0 0.0 0.0 0.0 0.0 0.0 0.49 0.0 0.0 0.0 0.0 0.0 0.0\n 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0|0.25 0.0 128.0 -0.238986 -0.147846 0.147846 -0.238986 -0.07213 -0.801111 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0|0.5 1.0 255.0 0.687624 -0.726067 0.726067 0.687624 0.038586 0.584278 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0|";
		
		if (Math.random() > 0.5)
			parseFormula(julia);
		else
			parseFormula(spiral);
	}

	public void randomize_c()
	{
		int i, j, k;

		for (i = 0; i < 2; i++)
		{
			for (j = 0; j < 3; j++)
			for (k = 0; k < 2; k++)
			xp[i].c[j][k] = qsrandom11();
		}
	}

	public double qsrandom11()
	{
		return Math.random()*2 - 1;
	}

	public int random_distrib(int[] items)
	{
		int r = ((int)Math.floor( Math.random() * items.length ))%items.length;
		return items[r];
	}
	
	
	public void loadFormula(XForm xform, String formula)
	{
		try
		{
			StringTokenizer st = new StringTokenizer(formula, " ");
			
			double[] values = new double[3 + XForm.NVARS + 6];
			int index = 0;
			while (st.hasMoreTokens())
			{
				String token = st.nextToken().trim();
				//System.out.println("token: " + token);
				if (token.length() > 0)
				{
					values[index] = Double.parseDouble( token );
					index++;
				}
			}
				
			xform.density = values[0];
			xform.symmetry = values[1];
			xform.color = values[2];
			
			int counter = 0;	
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					xform.c[i][j] = values[3+counter];
					counter++;
				}
			}
			
			for (int i = 9; i < values.length; i++)
			{
				xform.var[i-9] = values[i];
			}
						
		}
		catch(Exception e)
		{
			System.out.println("Problem parsing xform: " + e.toString());
		}
	}
	
	private void parseFormula(String s)
	{
		for (int i = 0; i < XForm.NVARS; i++) xp[i].density = 0.0;

		StringTokenizer st = new StringTokenizer(s, "|");
		int i = 0;
		while (st.hasMoreTokens())
		{
			System.out.println("Parsing " + i);
			String formula = st.nextToken();
			loadFormula(xp[i], formula);
			i++;
		}
	}
	
	
	private void getConstants()
	{
		try
		{
			brightness = Double.parseDouble(txtBrightness.getText());
			gradient = Double.parseDouble(txtGradient.getText());
			zoom = Double.parseDouble(txtZoom.getText());
			
			cameraPersp = Double.parseDouble(txtPerspective.getText());
			cameraZpos = Double.parseDouble(txtZPos.getText());
			cameraPitch = Double.parseDouble(txtPitch.getText());
			cameraYaw = Double.parseDouble(txtYaw.getText());
			cameraDOF = Double.parseDouble(txtDOF.getText());
			opacity = Double.parseDouble(txtOpacity.getText());
			translucency = Double.parseDouble(txtReflect.getText());	
			focus = Double.parseDouble(txtFocus.getText());
			horizon = Double.parseDouble(txtHorizon.getText());
			zres = Double.parseDouble(txtZRes.getText());
			
			setCamera();
			
			for (int i = 0; i < XForm.NVARS; i++) xp[i] = new XForm();
			parseFormula(txtFormula.getText());
		}
		catch(Exception e)
		{
			System.out.println("Problem parsing formula: " + e.toString());
		}
	}
	
	private void setCamera()
	{
		// 3D camera precalc
		CameraMatrix[0][0] = Math.cos(-cameraYaw);
		CameraMatrix[1][0] = -Math.sin(-cameraYaw);
		CameraMatrix[2][0] = 0;
		CameraMatrix[0][1] = Math.cos(cameraPitch) * Math.sin(-cameraYaw);
		CameraMatrix[1][1] = Math.cos(cameraPitch) * Math.cos(-cameraYaw);
		CameraMatrix[2][1] = -Math.sin(cameraPitch);
		CameraMatrix[0][2] = Math.sin(cameraPitch) * Math.sin(-cameraYaw);
		CameraMatrix[1][2] = Math.sin(cameraPitch) * Math.cos(-cameraYaw);
		CameraMatrix[2][2] = Math.cos(cameraPitch);
	}
	
	private void setConstants()
	{
		try
		{
			txtBrightness.setText(""+ brightness);
			txtGradient.setText(""+ gradient);
			txtZoom.setText(""+ (Math.round(zoom*100))/100.0);
			
			txtPerspective.setText("" + cameraPersp);
			txtZPos.setText("" + cameraZpos);
			txtPitch.setText("" + cameraPitch);
			txtYaw.setText("" + cameraYaw);
			txtDOF.setText("" + cameraDOF);
			txtOpacity.setText("" + opacity);
			txtReflect.setText("" + translucency);
			txtFocus.setText("" + focus);
			txtHorizon.setText("" + horizon);
			txtZRes.setText("" + zres);
			
			StringBuffer sb = new StringBuffer();
			for (int index = 0; index < xp.length; index++)
			{
				if (xp[index].density > 0)
				{
					sb.append( xp[index].density + " " );
					sb.append( xp[index].symmetry + " " );
					sb.append( xp[index].color + " \n" );
					
					for (int i = 0; i < 3; i++)
						for (int j = 0; j < 2; j++)
							sb.append( xp[index].c[i][j] + " " );
							
					sb.append("\n");
							
					for (int i = 0; i < XForm.NVARS; i++)
					{
						sb.append( xp[index].var[i] + " ");
						
						//break at the start of 3D variations
						if (i == 17) sb.append("\n");
					}
							
					sb.append("\n|\n");
				}
			}
			
			txtFormula.setText( sb.toString() );
		}
		catch(Exception e)
		{
			System.out.println(e.toString());
		}
	}
	
	
	public boolean mouseDown(Event evt, int x, int y)
	{
		xanchor = x;
		yanchor = y;
		
		if (evt.modifiers == (Event.META_MASK + Event.CTRL_MASK))
		{
			randomize_c();
			reset = true;
		}
		else if (evt.modifiers == (Event.META_MASK + Event.SHIFT_MASK))
		{
			int num_affected = add_symmetry_to_control_point();
			reset = true;
			System.out.println("affected: " + num_affected);
		}
		else if (evt.modifiers == (Event.ALT_MASK + Event.CTRL_MASK))
		{
			random_xform(-1);
			reset = true;
		}
		else if (evt.modifiers == Event.ALT_MASK)
		{
			alg = (alg+1)%XForm.NVARS;
			predefined_xform(alg);
			reset = true;
		}
		else if (evt.modifiers == Event.SHIFT_MASK)
		{
			pal  = (pal + 1)%(Pallet.fpalette.length);
			reset = true;
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
			double newxcen = (xanchor + xcurr)/2.0f;
			newxcen = xcen - ((newxcen - ximlen/2)/(ximlen))*zoom;
			xcen = newxcen;
	
			double newycen = (yanchor + ycurr)/2.0f;
			newycen = ycen - ((newycen - yimlen/2)/(yimlen))*zoom;
			ycen = newycen;
			
			System.out.println("Xcen is " + xcen + " Ycen is " + ycen);
	
			zoom = (double)dx/ximlen*zoom;
			reset = true;
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
		
			if (e.getSource() != txtGradient && e.getSource() != txtBrightness) 
				reset = true;
	
			getConstants();
		}
	}
	
	public void keyTyped(KeyEvent e)
	{
	}
	
}



class XForm
{
	public static int NVARS = 27;
	double[] var = new double[NVARS];
	public double[][] c = new double[3][2];
	public double density = 0.5;
	public double color = 0;
	public double symmetry = 0;

	public XForm()
	{
	}
}

