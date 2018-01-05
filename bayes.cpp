/***************************************************************************
 * bayes.c
 *
 * This file contains a collection of routines used in MCMC.
 *
 ***************************************************************************/

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <exception>
#include <vector>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bayes.h"
#include "nrutil.h"

double Rand1(){
  // Return a random number between 0 and 1
  // Most basic random number generator. 
  return( rand()/(1.0*RAND_MAX) );
}


double NormalDev(double mu, double sig){

  double u,v,x,y,q;
  double r;

  do {
    u =  Rand1();
    v = 1.7156*(Rand1()-0.5);
    x = u - 0.449871;
    y = fabs(v) + 0.386595;
    q = x*x + y*(0.19600*y - 0.25472*x);
  } while (q > 0.27597
	   && ( q > 0.27846 || v*v > -4.*log(u)*u*u));

  r = mu + sig*v/u;
  /* std::cout << "NormalDev:  r = " << r 
	    << " mean=" << mu
	    << " sig=" << sig
	    << std::endl;*/

  return ( mu + sig*v/u);

}



double Gaussian(double x, double mu, double sig){
  // Probability of x 
  return (0.398942280401432678/sig * exp(-0.5 * pow( (x-mu)/sig ,2))); 
}

double Gaussian2d(double x, double y, double *mean, double **var){
  // Probability of x,y

  double det = var[1][1] * var[2][2] - pow( var[1][2],2);

  return (0.5/(3.1415192653589793*sqrt(det)) 
	  * exp( - 0.5 * var[1][1] * var[2][2]/det * 
		 ( pow( x-mean[1],2)/var[1][1] + pow(y-mean[2],2)/var[2][2]
		   -2.0 * var[1][2]/(var[1][1]*var[2][2]) * (x-mean[1]) * (y-mean[2]))));

}
