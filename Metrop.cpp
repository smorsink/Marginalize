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
    xlo(-6.0), xhi(10.0);
  int
    numbins(100);




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
 
    long int *histogram = ivector(1,numbins);
    double *probability = dvector(1,numbins);
    double *xvals = dvector(1,numbins);

    double dx;

    dx = (xhi - xlo)/(1.0*numbins);

    for (unsigned int i(1); i<=numbins; i++){

      xvals[i] = xlo + (0.5 + i - 1)*dx;
      /* std::cout << " bin = " << i
		<< " bin centre = " << xvals[i] 
		<< " left edge = " << xvals[i] - 0.5*dx
		<< std::endl;*/

      histogram[i] = 0;

    }

    int bin;

    // Set up x-bins for storing the histogram






  srand(time(NULL));   // should only be called once

  double x1(0.0), x2(0.0);
  double p1, p2;

  double var(2.0);
  double mean(2.0);
  double sig;
  sig = sqrt(var);

  double xave(0.0), yave(0.0);

  double qmu, qvar, qsig;

  // Create Initial Step 
  x1=0.0;

   double r;

   for ( unsigned int i(0); i < nsteps; i++){

     xave += x1;

     p1 = Gaussian(x1,mean,sig);  

     // Create Proposal with mean=x1; var=1
     qmu = x1;
     qvar = 1.0;
     qsig = sqrt(qvar);

     // Draw a new value for x from the Proposal Distribution
     x2 = NormalDev(qmu,qsig);
     p2 = Gaussian(x2,mean,sig);

     // Draw a random number (uniform between 0 and 1)
     
     r = Rand1();
     /*
     std::cout << " i = " << i
	       << " x1 = " << x1
	       << " p1 = " << p1
	       << " x2 = " << x2
	       << " p2 = " << p2
	       << " r = " << r
	       << " p2/p1 = " << p2/p1
	       << std::endl;
     */

     // Do we accept the new value?
     if ( p2/p1 > r )
       x1 = x2;

     // Increment the correct bin of the histogram

     if (x1 < xlo || x1 > xhi)
       std::cout << "x1 = " << x1 << " is out of bounds! Increase them!" << std::endl;
     else{ // Within the bounds

       bin = (x1-xlo)/dx + 1;
       /* std::cout << "x = " << x1
		 << " bin = " << bin
		 << " xbin_lo = " << xvals[bin] - 0.5*dx
		 << " xbin_hi = " << xvals[bin] + 0.5*dx
		 << std::endl;*/

       histogram[bin] += 1;

     }

   } // End of MT loop

   // Compute the average value of x
   xave *= 1.0/(nsteps);
   std::cout << "N steps = " << nsteps << "\t"
	     << " Average x = " << xave << std::endl;

   // Integrate the histogram
   double totprob=0.0;
   for (unsigned int i(1); i<=numbins; i++){
     totprob += histogram[i]*dx;
   }

   //std::cout << "total prob = " << totprob << std::endl;
   double xmed(xlo);
   double xl1(xlo), xl2(xlo), xl3(xlo);
   double xr1(xlo), xr2(xlo), xr3(xlo);
   double cumul(0.0);

   // Normalize the histogram
   for (unsigned int i(1); i<=numbins; i++){
     probability[i] = histogram[i]/(1.0*totprob);
     cumul += probability[i]*dx;
     std::cout << "i=" << i 
	       << " x[i]= " << xvals[i]
	       << " prob[i]=  " << probability[i]
	       << " cumul= " << cumul
	       << std::endl;

     // check for 3sigma
     if ( cumul >= 0.0015 && xl3==xlo ){
       xl3 = xvals[i];
       std::cout << "3 sig: " << " xl3 = " << xl3 << " cumul = " << cumul << std::endl;
     }
     // check for 2sigma
     if ( cumul >= 0.025 && xl2==xlo ){
       xl2 = xvals[i];
       std::cout << "2 sig: " << " xl2 = " << xl2 << " cumul = " << cumul << std::endl;
     }
     // check for 1sigma
     if ( cumul >= 0.16 && xl1==xlo ){
       xl1 = xvals[i];
       std::cout << "1 sig: " << " xl1 = " << xl1 << " cumul = " << cumul << std::endl;
     }

     // check for median value
     if ( cumul >= 0.5 && xmed==xlo ){
       xmed = xvals[i];
       std::cout << "Median: " << " xmed = " << xmed << " cumul = " << cumul << std::endl;
     }

     // check for 1sigma
     if ( cumul >= 0.84 && xr1==xlo ){
       xr1 = xvals[i];
       std::cout << "1 sig: " << " xr1 = " << xr1 << " cumul = " << cumul << std::endl;
     }
     // check for 2sigma
     if ( cumul >= 0.975 && xr2==xlo ){
       xr2 = xvals[i];
       std::cout << "2 sig: " << " xr2 = " << xr2 << " cumul = " << cumul << std::endl;
     }
     // check for 3sigma
     if ( cumul >= 0.9985 && xr3==xlo ){
       xr3 = xvals[i];
       std::cout << "3 sig: " << " xr3 = " << xr3
 << " cumul = " << cumul << std::endl;
     }



   }

      std::cout << "N steps = " << nsteps << "\t"
		<< " dx = " << dx 
	     << " Median x = " << xmed << std::endl;

   // Output histogram
   out.open(out_file);
   for (unsigned int i(1); i<=numbins; i++){

     out << xvals[i] << " " << probability[i] << " " 
	 << Gaussian(xvals[i],mean,sig) << " "
	 << i << std::endl;
   }


    return 0;


    //free_ivector(histogram,1,numbins);
    free_dvector(xvals,1,numbins);

 } 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}
