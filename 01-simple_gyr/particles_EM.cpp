#include <iostream> // Declares the standard input/output stream objects
#include <cstdlib>  // Declares a set of general purpose functions
#include <fstream>  // Input/output stream class to operate on files
#include <iomanip>  // Set precision
#include <cmath>    // Declares mathematical library
#include <chrono>   // Random time  generator
#include <random>   // Random time  generator

using namespace std;

// Function prototipes:
void FillIC (double *, double *,double *, double *);
void GetParam (double *, double *, double *, double &, double &);
void Analitic (double *, double *, double *, double *, double *, double, double, double);
void Ydot (double, double *, double *, double *, double *); 
void RK4 (double , double *, void (*func)(double, double *, double *, double *, double *), double, int);
void Boris (double *, double *, double *, double , int);

// Golobal variable:
int ntsteps      = 1000; // Number of time steps
int npart        = 1; // Number of particles
const int dim    = 3; // Problem dimension
const int neq    = 6; // Dimension of the IC vector
double k         = 4.0; // Constant
double tbeg      = 0.0; // Begin time
double tend      = k*M_PI; // End time 
double qtom      = 1.0; // Charge to mass ratio
double vth       = 0.001; // Thermal velocity
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
	/*                                DATA STORAGE                                */
	//****************************************************************************80
	
  cout << setiosflags (ios::scientific);
	
	// Create a file .dat for Analitic:
	string fname0   = "analitic_vz.dat";
	ofstream fdata0 (fname0.c_str(), ios::out);
	fdata0 << setiosflags(ios::scientific);
	
	// Create a file .dat for RK4:
	string fname1   = "rk4_vz.dat";
	ofstream fdata1 (fname1.c_str(), ios::out);
	fdata1 << setiosflags(ios::scientific);
	
	// Create a file .dat for Boris:
	string fname2   = "boris_vz.dat";
	ofstream fdata2 (fname2.c_str(), ios::out);
	fdata2 << setiosflags(ios::scientific);
	
	// Create a file .dat for Error:
	string fname3   = "error.dat";
	ofstream fdata3 (fname3.c_str(), ios::out);
	fdata3 << setiosflags(ios::scientific); 
	
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
	
	double x0[dim]; // initial position
	double v0[dim]; // initial velocity
  double x[dim]; // position
  double v[dim]; // velocity
	double vd[dim]; // drift velocity
  double vperp; // perpendicular velocity respect to B
	double Y[neq]; // Initial condition vector
	double t; // Time parameter
	double h = (tend - tbeg)/(ntsteps - 1); // Increment
	double xmoda;
	double vmoda;
	
	//****************************************************************************80
	/*                              INITIAL CONDITION                             */
	//****************************************************************************80
	
	// def. initial position comp.:
	x0[0] = 0.0;
	x0[1] = 0.0;
	x0[2] = 0.0;
	// def. initial velocity comp.:
	v0[0] = 0.0;
	v0[1] = 1.0;
	v0[2] = 0.0;
	// def. initial E field comp.:
	E[0] = 0.0;
	E[1] = 0.0;
	E[2] = 0.0;
	// def. initial B field comp.:
	B[0] = 0.0;
	B[1] = 0.0;
	B[2] = 1.0;
	
	//****************************************************************************80
	/*                                 RUN FUNCTIONS                              */
	//****************************************************************************80
	
	// Fill initial condition inside Y before start integrating
	FillIC (Y, x0, v0, vd);
	// Start integrating from tbeg
	t = tbeg;
	
	double phi;
	// Definien new parameter:
	//GetParam (x0, v0, vd, vperp);
	GetParam (x0, v0, vd, vperp, phi);
	
	// Analitic Integration
	for (int i = 0; i < ntsteps - 1; i++)
	{
		//Analitic (x, v, x0, v0, vd, vperp, t);
		Analitic (x, v, x0, v0, vd, vperp, phi, t);
		fdata0 << t << " " << x[0] << " "<< x[1] << " " << x[2] << endl;
    t += h;
	}
	fdata0 << endl << endl;
	
	xmoda = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	vmoda = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	
	// Restart initial condition
	FillIC (Y, x0, v0, vd);
	// Restart integrating from tbeg
	t = tbeg;
	
	//RK4 Integration
	for (int i = 0; i < ntsteps - 1; i++)
	{
		RK4 (t, Y, Ydot, h, neq);
    fdata1 << t << " " << Y[0] << " " << Y[1] << " "<< Y[2] << endl;;
    fdata1 << endl;
    t += h;
	}
	fdata1 << endl << endl;
	fdata3 << t << " " << (sqrt(Y[0]*Y[0] + Y[1]*Y[1] + Y[2]*Y[2]) - xmoda)/xmoda << " "<< (sqrt(Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5]) - vmoda)/vmoda << endl;
	
	// Restart initial condition
	FillIC (Y, x0, v0, vd);
	// Restart integrating from tbeg
	t = tbeg;
	
	// Boris Integration
	for (int i = 0; i < ntsteps - 1; i++)
	{
		Boris (Y, E, B, h, dim);
		fdata2 << t << " " << Y[0] << " "<< Y[1] << " " << Y[2] << endl;
    t += h;
	}
	fdata2 << endl << endl;
	
	fdata3 << t << " " << (sqrt(Y[0]*Y[0] + Y[1]*Y[1] + Y[2]*Y[2]) - xmoda)/xmoda << " "<< (sqrt(Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5]) - vmoda)/vmoda << endl;
  fdata3 << endl << endl;
	
	fdata0.close();
	fdata1.close();
	fdata2.close();
	fdata3.close();
	
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

void GetParam (double *x0, double *v0, double *vd, double &vperp, double &phi)
{
  double thetaB;
  double thetav;
  double theta;
  double v0abs;
	
	v0abs = sqrt(v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]);
  thetaB = atan(B[2]/(sqrt(B[0]*B[0] + B[1]*B[1])));
  thetav = atan(v0[2]/(sqrt(v0[0]*v0[0] + v0[1]*v0[1])));
  theta = fabs(thetaB - thetav);
  vperp = v0abs*sin(theta);
	
  vd[0] = E[1]*B[2] - E[2]*B[1];
  vd[1] = E[2]*B[0] - E[0]*B[2];
  vd[2] = E[0]*B[1] - E[1]*B[0];
	
  if(v0[0] != 0.0)
  {
    phi = atan(- v0[1]/v0[0]);
  }
  else if(v0[1] != 0.0)
  {
    phi = M_PI/2.0 - atan(- v0[0]/v0[1]);
  }
  else
  {
    cout << "Null vperp!" << endl;
    phi = 0.0;
  }
  phi = - fabs(phi);
}

void Analitic (double *x, double *v, double *x0, double *v0, double *vd, double vperp, double phi, double t)
{
  x[0] = x0[0] + vperp*sin(t + phi) + vd[0]*t - vperp*sin(phi);
  x[1] = x0[1] + vperp*cos(t + phi) + vd[1]*t - vperp*cos(phi);
  x[2] = x0[2] + v0[2]*t;
  v[0] = vperp*cos(t + phi) + vd[0];
  v[1] = -vperp*sin(t + phi) + vd[1];
  v[2] = v0[2];
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