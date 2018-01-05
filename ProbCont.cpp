/***************************************************************************************/
/*                                   ProbCont.cpp

This code reads in a probability map and creates contours.

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

// MAIN
int main ( int argc, char** argv ) try {  // argc, number of cmd line args; 
                                          // argv, actual character strings of cmd line args

  /*********************************************/
  /* VARIABLE DECLARATIONS AND INITIALIZATIONS */
  /*********************************************/
    
  std::ofstream out1;      // output stream; printing information to the output file
  std::ofstream out2;      // output stream; printing information to the output file
  std::ofstream out3;      // output stream; printing information to the output file

  std::ifstream in;       // probability map


  char out_file1[256] = "sigma1.txt";
  char out_file2[256] = "sigma2.txt";
  char out_file3[256] = "sigma3.txt";


  char in_file1[256] = "No File Specified";
  char in_file2[256] = "No File Specified";


  unsigned int nummass(10),numradius(10);

  unsigned int numsets(1);

  double
    radius_lo(8.0),
    radius_hi(18.0),
    mass_lo(1.0),
    mass_hi(2.2);

  double
    dradius, dmass;

  double level1(0.1),
    level2(0.01),
    level3(0.001);



  /*********************************************************/
  /* READING IN PARAMETERS FROM THE COMMAND LINE ARGUMENTS */
  /*********************************************************/
    
    for ( int i(1); i < argc; i++ ) {
        if ( argv[i][0] == '-' ) {  // the '-' flag lets the computer know that we're giving it information from the cmd line
            switch ( argv[i][1] ) {

	    case 'j': // Number of Mass bins
	      sscanf(argv[i+1], "%u", &nummass);
	      break;

	    case 'i':
	      sscanf(argv[i+1], "%u", &numradius);
	      break;

	    case 'l':
	      sscanf(argv[i+1], "%lf", &level1);
	      break;

	    case 'm':
	      sscanf(argv[i+1], "%lf", &level2);
	      break;

	    case 'n':
	      sscanf(argv[i+1], "%lf", &level3);
	      break;

	    case 'N':
	      sscanf(argv[i+1], "%u", &numsets);
	      break;


	    case 'p':  // Name of 1st input file
	                sscanf(argv[i+1], "%s", in_file1);
	                break;

	    case 'q':  // Name of 1st input file
	                sscanf(argv[i+1], "%s", in_file2);
	                break;

	   
	                
            } // end switch	
        } // end if
    } // end for

    double **probability = dmatrix(1,numradius,1,nummass);
    int **histogram = imatrix(1,numradius,1,nummass);
    int **sigma1 = imatrix(1,numradius,1,nummass);
    int **sigma2 = imatrix(1,numradius,1,nummass);
    int **sigma3 = imatrix(1,numradius,1,nummass);
    double *mvals = dvector(1,nummass);
    double *rvals = dvector(1,numradius);
    
    dradius = (radius_hi - radius_lo)/(numradius);
    dmass = (mass_hi - mass_lo)/(nummass);

    for (unsigned int i(1); i<= numradius; i++){
      rvals[i] = radius_lo + (i-0.5)*dradius;
      //std::cout << "i = " << i << " radius value = " << rvals[i] << std::endl;
    }

    for (unsigned int i(1); i<= nummass; i++){
      mvals[i] = mass_lo + (i-0.5)*dmass;
      //std::cout << "i = " << i << " mass value = " << mvals[i] << std::endl;
    }
    


    //    std::cout << "Opening in_file = " << in_file << std::endl;
    //in.open(in_file);


    std::cout << "Level 1 = " << level1
	      << " Level 2 = " << level2
	      << " Level 3 = " << level3
	      << std::endl;

    // Create Matrices
    
    double totprob(0.0);

    double totprob1(0.0), totprob2(0.0), totprob3(0.0);

    // Read in the datasets

    double lllog, prob, hist, lllog2;

    for(unsigned int k(1); k<=numsets; k++){

      if (k==1){
	std::cout << "Opening in_file = " << in_file1 << std::endl;
	in.open(in_file1);
      }
      if (k==2){
	std::cout << "Opening in_file = " << in_file2 << std::endl;
	in.open(in_file2);
      }
      for (unsigned int i(1); i<numradius; i++){
	for (unsigned int j(1); j<nummass; j++){

	  in >> rvals[i];
	  in >> mvals[j];
	  in >> lllog;
	  in >> prob;
	  in >> lllog2;
	  in >> hist;

	  prob = exp(-lllog2+2306);

	  if (hist !=0)
	  std::cout << "i=" << i << " j=" << j 
		    << " R = " << rvals[i] << " M = " << mvals[j]
		    << " histo = " << hist
		    << " log(prob) = " << lllog2
		    << "   prob = " << prob << "  totprob = " << totprob << std::endl;


	  probability[i][j] += prob;
	  histogram[i][j] += hist;

	  totprob += prob;

	}
      }
      in.close();

    }

    std::cout.precision(10);
    std::cout << "Total Probability = " << totprob << std::endl;
      

    for (unsigned int i(1); i<numradius; i++){
      for (unsigned int j(1); j<nummass; j++){

	probability[i][j] *= 1.0/totprob;

	if ( probability[i][j] >= level1){
	  totprob1 += probability[i][j];
	  sigma1[i][j] = 1;
	}
	else {
	  if ( probability[i][j] >= level2){
	    totprob2 += probability[i][j];
	    sigma2[i][j] = 1;
	      
	  }
	  else{
	    if ( probability[i][j] >= level3){
	      totprob3 += probability[i][j];
	      sigma3[i][j] = 1;
	    }
	  }
	}
      }
    }




      std::cout << "1 sigma = " << totprob1 << std::endl;
      std::cout << "2 sigma = " << totprob2 + totprob1 << std::endl;
      std::cout << "3 sigma = " << totprob3 + totprob2 + totprob1 << std::endl;



      // std::cout << "Opening out_file = " << out_file1 << std::endl;

    out1.open(out_file1);
    out2.open(out_file2);
    out3.open(out_file3);



      for (unsigned int i(1); i<numradius; i++){
	for (unsigned int j(1); j<nummass; j++){

	  if (sigma1[i][j] > 0 ) 
	    out1 << rvals[i] << "\t"
	      << mvals[j] << "\t"
	      << sigma1[i][j] << "\t"
	      << std::endl;

	  if (sigma2[i][j] > 0 ) 
	    out2 << rvals[i] << "\t"
	      << mvals[j] << "\t"
	      << sigma2[i][j] << "\t"
	      << std::endl;

	  if (sigma3[i][j] > 0 ) 
	    out3 << rvals[i] << "\t"
	      << mvals[j] << "\t"
	      << sigma3[i][j] << "\t"
	      << std::endl;


	}
      }

      in.close();
      out1.close();
      out2.close();
      out3.close();




      free_dmatrix(probability,1,numradius,1,nummass);
      free_imatrix(histogram,1,numradius,1,nummass);
      free_imatrix(sigma1,1,numradius,1,nummass);
      free_imatrix(sigma2,1,numradius,1,nummass);
      free_imatrix(sigma3,1,numradius,1,nummass);
      free_dvector(rvals,1, numradius);
      free_dvector(mvals,1, nummass);


    return 0;
 } 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}
