/***************************************************************************************/
/*                                   Metrop.cpp

One Dimensional Hastings-Metropolis Algorithm

*/
/***************************************************************************************/

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <exception>
#include <vector>
#include <string>
#include <string.h>
#include "nrutil.h"
#include "bayes.h"

// MAIN
int main ( int argc, char** argv ) try {  // argc, number of cmd line args; 
                                          // argv, actual character strings of cmd line args

  /*********************************************/
  /* VARIABLE DECLARATIONS AND INITIALIZATIONS */
  /*********************************************/
    
  std::ofstream out;      // output stream; printing information to the output file
  char out_file[256] = "output.txt";

  unsigned int nsteps(100);


  double
    xlo(-9.0), xhi(9.0);
  double
    ylo(-9.0), yhi(9.0);
  int
    numbins(500);




  /*********************************************************/
  /* READING IN PARAMETERS FROM THE COMMAND LINE ARGUMENTS */
  /*********************************************************/
    
    for ( int i(1); i < argc; i++ ) {
        if ( argv[i][0] == '-' ) {  // the '-' flag lets the computer know that we're giving it information from the cmd line
            switch ( argv[i][1] ) {

	    case 'o':  // Name of output file
	                sscanf(argv[i+1], "%s", out_file);
	                break;

	    case 'n': // Number of steps
	      sscanf(argv[i+1], "%u",&nsteps);
	      break;

	                
            } // end switch	
        } // end if
    } // end for

    //    int **histogram = imatrix(1,numradius,1,nummass);
 
    int **histogram = imatrix(1,numbins,1,numbins);
    int **sigma1 = imatrix(1,numbins,1,numbins);
    int **sigma2 = imatrix(1,numbins,1,numbins);
    int **sigma3 = imatrix(1,numbins,1,numbins);

    double **probability = dmatrix(1,numbins,1,numbins);
    double *xprob = dvector(1,numbins);
    double *yprob = dvector(1,numbins);
    double *xvals = dvector(1,numbins);
    double *yvals = dvector(1,numbins);

    double dx, dy;

    dx = (xhi - xlo)/(1.0*numbins);
    dy = (yhi - ylo)/(1.0*numbins);

    for (unsigned int i(1); i<=numbins; i++){

      xvals[i] = xlo + (0.5 + i - 1)*dx;
      yvals[i] = ylo + (0.5 + i - 1)*dy;
      /* std::cout << " bin = " << i
		<< " bin centre = " << xvals[i] 
		<< " left edge = " << xvals[i] - 0.5*dx
		<< std::endl;*/

      //histogram[i] = 0;

    }
    for (unsigned int i(1);i<=numbins; i++)
      for (unsigned int j(1); j<=numbins; j++)
	histogram[i][j] = 0;


    int xbin, ybin;

    // Set up x-bins for storing the histogram

  srand(time(NULL));   // should only be called once

  double x1(0.0), x2(0.0);
  double y1(0.0), y2(0.0);
  double p1, p2;

  double **var = dmatrix(1,2,1,2);
  double *mean = dvector(1,2);
  
  mean[1]=0.0;
  mean[2]=0.0;

  var[1][1] = 2.0;
  var[1][2] = 1.2;
  var[2][1] = 1.2;
  var[2][2] = 2.0;

  double qmux, qmuy, qvar, qsig;

  double xave(0.0), yave(0.0);

  // Create Initial Step 
  x1=-2.0;
  y1=-2.0;

   double r;

   for ( unsigned int i(0); i < nsteps; i++){

     xave += x1;
     yave += y1;

     p1 = Gaussian2d(x1,y1,mean,var);

     // Create Proposal with mean=x1; var=1
     qmux = x1;
     qmuy = y1;
     qvar = 1.0;
     qsig = sqrt(qvar);

     // Draw a new value for x from the Proposal Distribution
     x2 = NormalDev(qmux,qsig);
     y2 = NormalDev(qmuy,qsig);

     p2 = Gaussian2d(x2,y2,mean,var);

     // Draw a random number (uniform between 0 and 1)
     
     r = Rand1();
     
     /*   std::cout << " i = " << i
	       << " x1 = " << x1
	       << " y1 = " << y1
	       << " p1 = " << p1
	       << " x2 = " << x2
	       << " y2 = " << y2
	       << " p2 = " << p2
	       << " r = " << r
	       << " p2/p1 = " << p2/p1
	       << std::endl;*/
    

     // Do we accept the new value?
     if ( p2/p1 > r ){
       x1 = x2;
       y1 = y2;
     }

     // Increment the correct bin of the histogram

     if (x1 < xlo || x1 > xhi)
       std::cout << "x1 = " << x1 << " is out of bounds! Increase them!" << std::endl;
     else{ // Within the x bounds
       if (y1 < ylo || y1 > yhi)
	 std::cout << "y1 = " << y1 << " is out of bounds! Increase them!" << std::endl;
       else{ // Within the y bounds

	 xbin = (x1-xlo)/dx + 1;
	 ybin = (y1-ylo)/dy + 1;
	 /* std::cout << "x = " << x1
	    << " bin = " << bin
	    << " xbin_lo = " << xvals[bin] - 0.5*dx
	    << " xbin_hi = " << xvals[bin] + 0.5*dx
	    << std::endl;*/

	 histogram[xbin][ybin] += 1;

	 /*
	   std::cout << " bin = " << bin
	   << " histogram[bin] = " << histogram[bin]
	   << std::endl;*/
       }
     }

   } // End of MT loop

   xave /= (1.0*nsteps);
   yave /= (1.0*nsteps);

   std::cout << "Average x = " << xave << std::endl;
   std::cout << "Average y = " << yave << std::endl;




   // Integrate the histogram
   double totprob=0.0;
   for (unsigned int i(1); i<=numbins; i++){
     for (unsigned int j(1); j<=numbins; j++){
       totprob += histogram[i][j]*dx*dy;
     }
   }

   // Normalize the probability
   for (unsigned int i(1); i<=numbins; i++){
     for (unsigned int j(1); j<=numbins; j++){
       probability[i][j] = histogram[i][j]/(1.0*totprob);
     }
   }
   // Probability is normalized so integral over x and y yields 1.0

   // totprob = 0.0;
   //Integrate over y to get the xprobability 
   for (unsigned int i(1); i<=numbins; i++){
     xprob[i] = 0.0;
     for (unsigned int j(1); j<=numbins; j++){
       xprob[i] += probability[i][j]*dy;
     }
     //totprob += xprob[i]*dx;
   }
   //   xprob is normalized so integral over x is 1
   //   std::cout << "Total Prob = " << totprob << std::endl;
   //for (unsigned int i(1); i<=numbins; i++)
   // xprob[i] /= totprob;

   
   //totprob = 0.0;
   //Integrate over x to get the yprobability 
   for (unsigned int j(1); j<=numbins; j++){
     yprob[j] = 0.0;
     for (unsigned int i(1); i<=numbins; i++){
       yprob[j] += probability[i][j]*dx;
     }
     //totprob += yprob[j]*dy;
   }
   // yprob is normalized so integral over y = 1
   //   std::cout << "y Total Prob = " << totprob << std::endl;

   //for (unsigned int j(1); j<=numbins; j++)
   //yprob[j] /=totprob;

   // Find 2D region so probability = 0.68
   double targetprob(0.05);
   double dp(0.001);
   totprob = 0.0;
   while (totprob < 0.68 && targetprob > 0.0){
     totprob = 0.0;
     for (unsigned int i(1); i<=numbins; i++){
       for (unsigned int j(1); j<=numbins; j++){
	 if (probability[i][j] >= targetprob){
	   totprob += probability[i][j]*dx*dy;
	   sigma1[i][j] = 1;
	 }
       }
     }
     std::cout << "Points with Prob >= " << targetprob
	     << " Integrated Prob = " << totprob
	     << std::endl;
     targetprob -= dp;
   }
   // totprob = 0.68 now

   // Find 2D region so probability = 0.95
   //targetprob set from before.
   totprob = 0.0;
   while (totprob < 0.95 && targetprob > 0.0){
     totprob = 0.0;
     for (unsigned int i(1); i<=numbins; i++){
       for (unsigned int j(1); j<=numbins; j++){
	 if (probability[i][j] >= targetprob){
	   totprob += probability[i][j]*dx*dy;
	   sigma2[i][j] = 1;
	 }
       }
     }
     std::cout << "Points with Prob >= " << targetprob
	     << " Integrated Prob = " << totprob
	     << std::endl;
     targetprob -= dp;
   }
   // totprob = 0.95 now

   // Find 2D region so probability = 0.997
   //targetprob set from before.
   totprob = 0.0;
   while (totprob < 0.997 && targetprob > 0.0){
     totprob = 0.0;
     for (unsigned int i(1); i<=numbins; i++){
       for (unsigned int j(1); j<=numbins; j++){
	 if (probability[i][j] >= targetprob){
	   totprob += probability[i][j]*dx*dy;
	   sigma3[i][j] = 1;
	 }
       }
     }
     std::cout << "Points with Prob >= " << targetprob
	     << " Integrated Prob = " << totprob
	     << std::endl;
     targetprob -= dp;
   }
   // totprob = 0.997 now


   out.open(out_file);
   for (unsigned int i(1); i<=numbins; i++){
     for (unsigned int j(1); j<=numbins; j++){

       out << xvals[i] << " " 
	   << yvals[j] << " " 
	   << probability[i][j] << " "
	   << Gaussian2d(xvals[i],yvals[j],mean,var) << " "
	   << " " << i << " " << j << std::endl;
     
     }
   }
   out.close();

   out.open("sigma1.txt");
   for (unsigned int i(1); i<=numbins; i++){
     for (unsigned int j(1); j<=numbins; j++){
       if (sigma1[i][j] > 0)
	 out << xvals[i] << " " 
	     << yvals[j] << " "
	     << sigma1[i][j] 
	     << std::endl;
     }
   }
   out.close();

   out.open("sigma2.txt");
   for (unsigned int i(1); i<=numbins; i++){
     for (unsigned int j(1); j<=numbins; j++){
       if (sigma2[i][j] > 0)
	 out << xvals[i] << " " 
	     << yvals[j] << " "
	     << sigma2[i][j] 
	     << std::endl;
     }
   }
   out.close();

   out.open("sigma3.txt");
   for (unsigned int i(1); i<=numbins; i++){
     for (unsigned int j(1); j<=numbins; j++){
       if (sigma3[i][j] > 0)
	 out << xvals[i] << " " 
	     << yvals[j] << " "
	     << sigma3[i][j] 
	     << std::endl;
     }
   }
   out.close();


   
   out.open("xprob.txt");
   for (unsigned int i(1); i<=numbins; i++){
     out << xvals[i] << " "
	 << xprob[i] << std::endl;
   }
   out.close();

   out.open("yprob.txt");
   for (unsigned int i(1); i<=numbins; i++){
     out << yvals[i] << " "
	 << yprob[i] << std::endl;
   }
   out.close();

    return 0;


    //free_ivector(histogram,1,numbins);
    free_dvector(xvals,1,numbins);

 } 
  

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}
