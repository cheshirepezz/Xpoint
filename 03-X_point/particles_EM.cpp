#include <iostream> // Declares the standard input/output stream objects
#include <cstdlib>  // Declares a set of general purpose functions
#include <fstream>  // Input/output stream class to operate on files
#include <iomanip>  // Set precision
#include <cmath>    // Declares mathematical library
#include <time.h>   // Random time generator


using namespace std;

// Function prototipes:
void FillIC (double *, double *,double *, double *);
void Ydot (double, double *, double *, double *, double *); 
void RK4 (double , double *, void (*func)(double, double *, double *, double *, double *), double, int);
void Boris (double *, double *, double *, double , int);


// Golobal variable:
const int npart  = 4096; // Number of particles
const int dim    = 3; // Problem dimension
const int neq    = 6; // Dimension of the IC vector
int ntsteps      = 500; // Number of time steps
double k         = 4.0; // Constant
double tbeg      = 0.0; // Begin time
double tend      = 150.0; // End time 
double qtom      = 1.0; // Charge to mass ratio
double vth       = 0.1; // Thermal velocity
double L         = 1000.0; // Lenght domain
double E[dim]; // E field vector
double B[dim]; // B field vector

//****************************************************************************80
//
int main()
//
//****************************************************************************80
//
 /*
  * last update:   05.05.2019
  * code name:     particles_EM.cpp
  * author:        Luca Pezzini
  * e-mail:        luca.pezzini@edu.unito.it
  * website:       https://github.com/lucapezzini
  * license:       GNU LGPL license.
  * 
  * Physics:
  * Study the motion of one or more (non-relativistic) particles in a fixed
  * electromagnetic field.
  * The equation of motion (possibly) involves propagation in all 3 direction.
  * The electric and magnetic field vector are given externally and particles
  * do not interact with each other.
  * The project involves direct integration of the equation of motion using
  * RK-type integrators and/or the generalization of symplectic integration scheme
  * to the case of velocity-dependent force (the Boris algorithm is the progenitor
  * of such schemes).
  * 
  */
{
	//****************************************************************************80
	/*                                   START                                    */
	//****************************************************************************80
	
	cout << "\n";
	cout << " ------------------------------------------------------ \n";
  cout << "+                  PARTICLES_EM.CPP                    +\n";
	cout << "+                                                      +\n";
  cout << "+                    C++ VERSION                       +\n";
  cout << "+           3-D PARTICLES IN EM SIMULATION             +\n";
	cout << "+                                                      +\n";
	cout << " ------------------------------------------------------ \n";
	cout << " \n";
	
	//****************************************************************************80
	/*                                PARAMETER                                   */
	//****************************************************************************80
	
	cout << "                    GLOBAL PARAMETER:                  " << endl;
	cout << "Dimension:                         " <<      dim         << endl;
	cout << "Number of particles:               " <<      npart       << endl;
	cout << "Number of time steps:              " <<      ntsteps     << endl;
	cout << "Time start:                        " <<      tbeg        << endl;
	cout << "Time stop:                         " <<      tend        << endl;
	cout << "Charge to mass ratio:              " <<      qtom        << endl;
	
	double x0[dim]; // Initial position
	double v0[dim]; // Initial velocity
  double x[dim]; // Position
  double v[dim]; // Velocity
	double vd[dim]; // Drift velocity
	double Y[neq]; // Initial condition vector
	double t; // Time parameter
	double h = (tend - tbeg)/double(0.5*ntsteps); // Time step
	
	int pside = 64; // Particles per side
	int pmarker = 0; // Particle marker
	int nframe = 500; // Number of frame 
	double Lstep = 2*L/(double)(pside); // Lenght gap between 2 particles
	double sinalpha; // sin(alpha_random)
	double icvec[npart][neq];
	double Ek;
	double pos[dim];
	double vel[dim];
	
	//****************************************************************************80
	/*                                DATA STORAGE                                */
	//****************************************************************************80
	
	string *fname = new string[nframe];
	char buffer [5];
  sprintf(buffer,"%04d",0);
  string num(buffer);
	fname[0] = "rk4_" + num + ".dat";
	//fname[0] = "boris_" + num + ".dat";

	/* Output File stream
	 * "ios::out" file open for writing mode
	 * c_str() converts a C++ string into a C-style string which is essentially
	 * a null terminated array of bytes. You use it when you want to pass a C++
	 * string into a function that expects a C-style string.
	 */
   ofstream fdata(fname[0].c_str(), ios::out);
	 cout << setiosflags (ios::scientific);
	
	//****************************************************************************80
	/*                              INITIAL CONDITION                             */
	//****************************************************************************80
		
	// def. initial E field comp.:
	E[0] = 0.0;
	E[1] = 0.0;
	E[2] = 0.5;
	
	//****************************************************************************80
	/*                                 RUN FUNCTIONS                              */
	//****************************************************************************80
	
	// drand48() -> returns a pseudo-random number in the range [0.0,1.0)
	srand48(time(NULL)); // Always set a seed value
	
	t = tbeg;

	// Define initial condition for all the particles 
	for (int i = 0; i < pside; i++) // Loop for the y
  {
		for (int j = 0; j < pside; j++) // Loop for the x
		{
			x0[0] = - L + j*Lstep;
			x0[1] = - L + i*Lstep;
			x0[2] = 0.0;
			
			sinalpha = drand48()*2.0 - 1.0; // Returns a pseudo-random number in the range (-1.0, +1.0) 
			
			v0[0] = vth*sqrt(1.0 - sinalpha*sinalpha);
			v0[1] = vth*sinalpha;
			v0[2] = 0.0;
			
			Ek = 0.5*(v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]);
			
			fdata << t << " " << x0[0] << " " << x0[1] << " " << x0[2] << " " << v0[0] << " " << v0[1] << " " << v0[2] << " " << Ek << endl;
			
			for (int k = 0; k < dim; k++)
			{
				icvec[pmarker][k] = x0[k];
				icvec[pmarker][k + 3] = v0[k];
			}
			pmarker += 1;
		}
	}
	
	fdata.close();
	t += h;
	pmarker = 0;
	
	// Time Loop:
	for (int i = 1; i < ntsteps; i++)
	{
		char buffer [5];
    sprintf(buffer,"%04d",i);
    string num(buffer);
    fname[i] = "rk4_" + num + ".dat";
		//fname[i] = "boris_" + num + ".dat";
		ofstream fdata(fname[i].c_str(), ios::out);

		// Particle Loop:
		for (int j = 0; j < npart; j++)
		{
			for (int k = 0; k < dim; k++)
			{
				pos[k] = icvec[pmarker][k];
				vel[k] = icvec[pmarker][k + 3];
			}
			
		  B[0] = pos[1]/L;
		  B[1] = pos[0]/L;
		  B[2] = 0.0;
				
			vd[0] = E[1]*B[2] - E[2]*B[1];
      vd[1] = E[2]*B[0] - E[0]*B[2];
			vd[2] = E[0]*B[1] - E[1]*B[0];
			
			FillIC(Y, pos, vel, vd);
		  RK4(t, Y, Ydot, h, neq);
			//Boris(Y, E, B, h, dim);
			Ek = 0.5*(Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5]);

      fdata << t << " " << Y[0] << " " << Y[1] << " " << Y[2] << " " << Y[3] << " " << Y[4] << " " << Y[5] << " " << Ek << endl;
	    
			for (int k = 0; k < dim; k++)
			{
				icvec[pmarker][k] = Y[k];
			  icvec[pmarker][k + 3] = Y[k + 3];
			}
			pmarker += 1;
		}
		
		pmarker = 0;
		fdata.close();
		t += h;	
	}
	
	//****************************************************************************80
	/*                                   END !                                    */
	//****************************************************************************80

  cout << "\n";
	cout << " ----------------------------------------------------- \n";
	cout << "+                                                     +\n";
  cout << "+               NORMAL END OF EXECUTION !             +\n";
  cout << "+                                                     +\n";
	cout << " ----------------------------------------------------- \n";
	cout << " \n";

  return 0;
}

void FillIC (double *Y, double *x0, double *v0, double *vd)
{
  Y[0] = x0[0];
  Y[1] = x0[1];
  Y[2] = x0[2];
  Y[3] = v0[0] + vd[0];
  Y[4] = v0[1] + vd[1];
  Y[5] = v0[2];
}

void Ydot (double t, double *Y, double *rhs, double *E, double *B) 
{
  double x = Y[0];
  double y = Y[1];
	double z = Y[2];
  double vx = Y[3];
  double vy = Y[4];
  double vz = Y[5];

  rhs[0] = vx;
  rhs[1] = vy;
  rhs[2] = vz;
  rhs[3] = qtom*(E[0] + vy*B[2] - vz*B[1]);
	rhs[4] = qtom*(E[1] - vx*B[2] + vz*B[0]);
	rhs[5] = qtom*(E[2] + vx*B[1] - vy*B[0]);
}

void RK4 (double t, double *Y, void (*Ydot)(double t, double *Y, double *rhs, double *E, double *B), double h, int neq)
{
  double k1[neq], k2[neq], k3[neq], k4[neq];
  double Y1[neq], Y2[neq], Y3[neq];
  
  Ydot(t, Y, k1, E, B);
  for(int n = 0; n < neq; n++)
  {
    Y1[n] = Y[n] + 0.5*h*k1[n];
  }
  Ydot(t + 0.5*h, Y1, k2, E, B);
  for(int n = 0; n < neq; n++)
  {
    Y2[n] = Y[n] + 0.5*h*k2[n];
  }
  Ydot(t + 0.5*h, Y2, k3, E, B);
  for(int n = 0; n < neq; n++)
  {
    Y3[n] = Y[n] + h*k3[n];
  }
  Ydot(t + h, Y3, k4, E, B);
  for(int n = 0; n < neq; n++)
  {
    Y[n] += h*(k1[n] + 2.0*k2[n] + 2.0*k3[n] + k4[n])/6.0;
  }
}

void Boris (double *Y, double *E, double *B, double dt, int neq)
{
	int i;
	double xhalf[neq];
	double xone[neq];
	double vone[neq];
	double vminus[neq];
	double vplus[neq];
	double b[neq];
	double bsquare;
	double cross1[neq];
	double cross2[neq];
	
	// Initialize b: 
	for (i = 0; i < neq; i++)
	{
		b[i] = 0.5*qtom*B[i]*dt;
	}
	
	// Calculate |b|^2: 
	for (i = 0; i < neq; i++)
	{
		bsquare += b[i]*b[i];
	}
	
	// Calculate x(n + 1/2) by calculeting
	// the particles position half a time step further (drift):
	for (i = 0; i < neq; i++)
	{
		xhalf[i] = Y[i] + 0.5*Y[i + 3]*dt;
	}
	
	// Calculate vminus (kick):
	for (i = 0; i < neq; i++)
	{
		vminus[i] = Y[i + 3] + 0.5*E[i]*dt;
	}
	
	// Calculating the 1st wedge product:
  cross1[0] = vminus[1]*b[2] - vminus[2]*b[1];
  cross1[1] = vminus[2]*b[0] - vminus[0]*b[2];
  cross1[2] = vminus[0]*b[1] - vminus[1]*b[0];
	
  // Calculating the 2nd wedge product:
  cross2[0] = ((vminus[1] + cross1[1])*b[2]) - ((vminus[2] + cross1[2])*b[1]);
  cross2[1] = ((vminus[2] + cross1[2])*b[0]) - ((vminus[0] + cross1[0])*b[2]);
  cross2[2] = ((vminus[0] + cross1[0])*b[1]) - ((vminus[1] + cross1[1])*b[0]);
	
  // Calculating v+ (rotate):
  vplus[0] = vminus[0] + 2.0*cross2[0]/(1.0 + bsquare);
  vplus[1] = vminus[1] + 2.0*cross2[1]/(1.0 + bsquare);
  vplus[2] = vminus[2] + 2.0*cross2[2]/(1.0 + bsquare);
	
	// Calculate v(n+1) (kick):
	for (i = 0; i < neq; i++)
	{
		vone[i] = vplus[i] + 0.5*qtom*E[i]*dt;
	}
	
	// Calculate x(n + 1) (drift):
	for (i = 0; i < neq; i++)
	{
		xone[i] = xhalf[i] + 0.5*vone[i]*dt;
	}
	// Save new position and velocity:
  for (i = 0; i < neq; i++)
	{
    Y[i] = xone[i];
		Y[i + 3] = vone[i];
  }
}