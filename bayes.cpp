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
#include "Struct.h"

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




// Compute the Probability Contours for a 1D distribution
struct Prob1d ProbContours1D(double *xvals, double *probability, double average, unsigned int numbins ){

  /*std::cout << "Hello World! Average = "
	    << average
	    << std::endl;*/

  struct Prob1d x;

  
  double xlo = xvals[1];
  
   double xmed(xlo);
   double xl1(xlo), xl2(xlo), xl3(xlo);
   double xr1(xlo), xr2(xlo), xr3(xlo);

   // Values of probability at the different values of x
   double pmed, pl1, pl2, pl3, pr1, pr2, pr3;
   
   double cumul(0.0);

   // Normalize the histogram
   for (unsigned int i(1); i<=numbins; i++){
     //probability[i] = histogram[i]/(1.0*totprob);
     cumul += probability[i];
     /*std::cout << "i=" << i 
	       << " x[i]= " << xvals[i]
	       << " prob[i]=  " << probability[i]
	       << " cumul= " << cumul
	       << std::endl;*/

     // check for 3sigma
     if ( cumul >= 0.0015 && xl3==xlo ){
       xl3 = xvals[i];
       pl3 = probability[i];
       //std::cout << "3 sig: " << " xl3 = " << xl3 << " cumul = " << cumul << std::endl;
     }
     // check for 2sigma
     if ( cumul >= 0.025 && xl2==xlo ){
       xl2 = xvals[i];
       pl2 = probability[i];
       //std::cout << "2 sig: " << " xl2 = " << xl2 << " cumul = " << cumul << std::endl;
     }
     // check for 1sigma
     if ( cumul >= 0.16 && xl1==xlo ){
       xl1 = xvals[i];
       pl1 = probability[i];
       //std::cout << "1 sig: " << " xl1 = " << xl1 << " cumul = " << cumul << std::endl;
     }

     // check for median value
     if ( cumul >= 0.5 && xmed==xlo ){
       xmed = xvals[i];
       pmed = probability[i];
       //std::cout << "Median: " << " xmed = " << xmed << " cumul = " << cumul << std::endl;
     }

     // check for 1sigma
     if ( cumul >= 0.84 && xr1==xlo ){
       xr1 = xvals[i];
       pr1 = probability[i];
       //std::cout << "1 sig: " << " xr1 = " << xr1 << " cumul = " << cumul << std::endl;
     }
     // check for 2sigma
     if ( cumul >= 0.975 && xr2==xlo ){
       xr2 = xvals[i];
       pr2 = probability[i];
       //std::cout << "2 sig: " << " xr2 = " << xr2 << " cumul = " << cumul << std::endl;
     }
     // check for 3sigma
     if ( cumul >= 0.9985 && xr3==xlo ){
       xr3 = xvals[i];
       pr3 = probability[i];
       //std::cout << "3 sig: " << " xr3 = " << xr3
       //	 << " pr3 = " << pr3
       //	 << " cumul = " << cumul << std::endl;
     }
   } // End for-i-loop

   x.median = xmed;
   x.xl1 = xl1;
   x.xl2 = xl2;
   x.xl3 = xl3;
   x.xr1 = xr1;
   x.xr2 = xr2;
   x.xr3 = xr3;

   x.pmed = pmed;
   x.pl1 = pl1;
   x.pl2 = pl2;
   x.pl3 = pl3;
   x.pr1 = pr1;
   x.pr2 = pr2;
   x.pr3 = pr3;

   return x;
  
}

struct Prob2d ProbContours2D( double *rvals, double *mvals, double **probability, unsigned int numradius, unsigned int nummass){

  struct Prob2d contours;

  double totprob;
  double probs1, probs2, probs3;
  //long int **sigma1 = imatrix(0,numradius,0,nummass);
  
  //  std::cout << "Welcome to ProbContours2D!" << std::endl;

  
   // Find 2D region so probability = 0.68
   double targetprob(0.05);
   double dp(0.0001);
   totprob = 0.0;
   double oldtarget = targetprob;
   double oldtotal = totprob;
   
   while (totprob < 0.68 && targetprob > 0.0){
     oldtarget = targetprob;
     oldtotal = totprob;
     targetprob -= dp;
     totprob = 0.0;
     for (unsigned int i(1); i<=numradius; i++){
       for (unsigned int j(1); j<=nummass; j++){
	 if (probability[i][j] >= targetprob){
	   totprob += probability[i][j];
	   //sigma1[i][j] = 1;
	 }
       }
     }
     //std::cout << "Points with Prob >= " << targetprob
     //	     << " Integrated Prob = " << totprob
     //	       << " OldTarget = " << oldtarget
     //	       << " OldTotal = " << oldtotal
     //	     << std::endl;
   }
   // totprob = 0.68 now
   probs1 = oldtarget + (0.68 - oldtotal) * (targetprob - oldtarget)/(totprob-oldtotal);

   //   std::cout << "1 Sigma Countour at Probability = " << probs1
   //	     << " targetprob = " << targetprob
   //	     << std::endl;
   

   // Find 2D region so probability = 0.95
   //targetprob set from before.
   totprob = 0.0;
   targetprob = probs1;
   dp = 0.0001;
   while (totprob < 0.95 && targetprob > 0.0){
     oldtarget = targetprob;
     oldtotal = totprob;
     targetprob -= dp;
     totprob = 0.0;
     for (unsigned int i(1); i<=numradius; i++){
       for (unsigned int j(1); j<=nummass; j++){
	 if (probability[i][j] >= targetprob){
	   totprob += probability[i][j];
	   //sigma2[i][j] = 1;
	 }
       }
     }
     //std::cout << "Points with Prob >= " << targetprob
     //	     << " Integrated Prob = " << totprob
     //	       << " OldTarget = " << oldtarget
     //	       << " OldTotal = " << oldtotal
     //	     << std::endl;

   }
   // totprob = 0.95 now
 
   probs2 = oldtarget + (0.95 - oldtotal) * (targetprob - oldtarget)/(totprob-oldtotal);

   //std::cout << "2 Sigma Countour at Probability = " << probs2
   //	     << " targetprob = " << targetprob
   //	     << std::endl;
   
   

   // Find 2D region so probability = 0.997
   //targetprob set from before.
   totprob = 0.0;
   targetprob = probs2;
   dp = 0.0001;
   while (totprob < 0.997 && targetprob > 0.0){
     oldtarget = targetprob;
     oldtotal = totprob;
     totprob = 0.0;
     targetprob -= dp;
     for (unsigned int i(1); i<=numradius; i++){
       for (unsigned int j(1); j<=nummass; j++){
	 if (probability[i][j] >= targetprob){
	   totprob += probability[i][j];
	   //sigma3[i][j] = 1;
	 }
       }
     }
     //std::cout << "Points with Prob >= " << targetprob
     //	     << " Integrated Prob = " << totprob
     //	     << std::endl;
     
   }
   // totprob = 0.997 now
  

   probs3 = oldtarget + (0.997 - oldtotal) * (targetprob - oldtarget)/(totprob-oldtotal);

   //   std::cout << "3 Sigma Countour at Probability = " << probs3
   //	     << " targetprob = " << targetprob
   //	     << std::endl;
   
   contours.sigma1 = probs1;
   contours.sigma2 = probs2;
   contours.sigma3 = probs3;

   return contours;


}
