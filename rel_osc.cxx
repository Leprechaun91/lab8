// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx,double* k1, double* k2, double* k3, double* k4);
double interpolaeschion(double* y0, double* yn, double dx, double* k1, double* k2, double* k3, double* k4);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
	const int dim = 2;
	double dx = 0.001,x=0;
	const double L = 100;
	double p0 = 1;
	for (p0 = 0.1; p0 <= 5; p0 += 0.1){
	    double y0[dim] = {p0, 0.0};
	    double yn[dim];
	    double klo = 0;
	    double k1[dim], k2[dim], k3[dim], k4[dim];
	    x = 0;
	    cout << "\\-------------------- p0 = " << p0 << "----------\\" << endl;
	    while(x<=L){
	      x += dx;
	      RKstep(yn, y0, x, dx, k1,k2,k3,k4);
	      if((yn[1] < 0) && (y0[1] > 0)){
		break;
	      }
	      for(int i=0; i<dim; i++) y0[i] = yn[i];
	    
	    }
	    klo = interpolaeschion(y0, yn, dx, k1, k2, k3, k4);
	    out << p0 << "\t" << x - dx + klo * dx << endl;
	}
	out.close();
// 	cout << klo << endl; 
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,const double x, const double dx,double* k1, double* k2, double* k3, double* k4 )
{
	const int dim = 2;
	
	for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1, x);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dx * k1[i];
	f(k2, x+0.5*dx);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dx * k2[i];
	f(k3, x+0.5*dx);

	for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * k3[i];
	f(k4,  x+dx);

	for(int i=0;i<dim; i++)
	 yn[i] = y0[i] + dx*(1/6.0 * k1[i]+1/3.0 * k2[i] + 1/3.0 * k3[i] + 1/6.0 * k4[i]);
	
	
	  
}
//-------------------
// Lorenz model
void f(double* const y0, const double x)
{
	double y[2] = {y0[0], y0[1]};

	y0[0] = y[1];
	y0[1] = -y[0]/sqrt(1+(y[0])*(y[0]));
	}

double interpolaeschion(double* y0, double* yn, double dx, double* k1, double* k2, double* k3, double* k4){
  double theta = 0.5;
  double b[4];
  double thetaL = 0;
  double thetaR = 1;
  b[0] = theta-3*theta*theta*0.5 +2*theta*theta*theta/3.0;
  b[1] = theta*theta - 2*theta*theta*theta/3;
  b[2] = theta*theta - 2*theta*theta*theta/3;
  b[3] = -theta*theta*0.5 + 2*theta*theta*theta/3;
  double ytemp = y0[1] + dx * (b[0] * k1[1] + b[1] * k2[1] + b[2] * k3[1] + b[3]*k4[1]);
//   int imax = 20,i=0;
//     cout << theta <<endl;
  while(abs(ytemp)>(1e-8)){
  if(ytemp > 0)
  thetaL = theta;
  else
  thetaR = theta;
  
  theta = (thetaL + thetaR)/2.0;
  
  b[0] = theta-3*theta*theta*0.5 +2*theta*theta*theta/3.0;
  b[1] = theta*theta - 2*theta*theta*theta/3;
  b[2] = theta*theta - 2*theta*theta*theta/3;
  b[3] = -theta*theta*0.5 + 2*theta*theta*theta/3;
  ytemp = y0[1] + dx * (b[0] * k1[1] + b[1] * k2[1] + b[2] * k3[1] + b[3]*k4[1]);
//   cout << theta << "\t" << ytemp << endl;
 
    
  }
  
  return theta;
  
  
  
  
}