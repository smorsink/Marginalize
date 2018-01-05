/***************************************************************************************/
/*                                   Marginalize.cpp

This code reads in a set of Optimal Models.
Parameters are sorted into bins.

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
  std::ifstream fdim;       // dimensions of optimals
  std::ifstream fopt;       // optimals
  std::ifstream ffit;       // fitness values
  std::ifstream frnk;       // rank of optimal

  char out_file[256] = "output.txt";

  char dim_file[256] = "No File Specified";
  char opt_file[256] = "No File Specified";
  char fit_file[256] = "No File Specified";
  char rnk_file[256] = "No File Specified";
  		
  unsigned int nopt,
    npara,
    nummass(40),numradius(40),numbins(12);

  double
    radius_lo(8.0),
    radius_hi(20.5),
    mass_lo(0.8),
    mass_hi(2.25);

  double
    dradius, dmass;

  double 
    inc_lo(60.0), inc_hi(125.0), dinc;
  double
    theta_lo(50.0),theta_hi(125.0),dtheta;
  double 
    time_lo(0.46),time_hi(0.55),dtime;
  double
    rho_lo(0.2),rho_hi(1.6),drho;
  double
    temp_lo(0.1),temp_hi(0.52),dtemp;
  double
    distance_lo(0.1),distance_hi(0.6),ddistance;

  double
    normalization(2315.0);

  normalization = 2500.0;

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
	      
	    case 'k':
	      sscanf(argv[i+1], "%u", &numbins);
	      break;
	      

	    case 'd':  // Name of input file
	                sscanf(argv[i+1], "%s", dim_file);
	                break;

	    case 'p':  // Name of input file
	                sscanf(argv[i+1], "%s", opt_file);
	                break;

	    case 'f':  // Name of input file
	                sscanf(argv[i+1], "%s", fit_file);
	                break;

	    case 'r':  // Name of input file
	                sscanf(argv[i+1], "%s", rnk_file);
	                break;

	    case 'o':  // Name of output file
	                sscanf(argv[i+1], "%s", out_file);
	                break;
	                
            } // end switch	
        } // end if
    } // end for

    //    double **probability = dmatrix(1,numradius,1,nummass);
    //int **histogram = imatrix(1,numradius,1,nummass);
    //double **minchi = dmatrix(1,numradius,1,nummass);

    double *loglikelihood = dvector(0,numradius*nummass*pow(numbins,6));

    /*    for ( unsigned int i(0); i<=numradius*nummass*pow(numbins,6); i++){
      loglikelihood[i] = 100000.0;
      }*/



    long int *histo = ivector(0,numradius*nummass*pow(numbins,6));

    // Create the Mass and Radius Intervals, and intervals for other parameters
    double *rvals = dvector(0,numradius);    
    double *mvals = dvector(0,nummass);
    double *ivals = dvector(0,numbins);
    double *thetavals = dvector(0,numbins);
    double *timevals = dvector(0,numbins);
    double *rhovals = dvector(0,numbins);
    double *tempvals = dvector(0,numbins);
    double *distancevals = dvector(0,numbins);


    dradius = (radius_hi - radius_lo)/(numradius);
    dmass = (mass_hi - mass_lo)/(nummass);
    dinc = (inc_hi - inc_lo)/numbins;
    dtheta = (theta_hi - theta_lo)/numbins;
    dtime = (time_hi - time_lo)/numbins;
    drho = (rho_hi - rho_lo)/numbins;
    dtemp = (temp_hi - temp_lo)/numbins;
    ddistance = (distance_hi - distance_lo)/numbins;

    for (unsigned int i(1); i<= numradius; i++){
      rvals[i] = radius_lo + (i-0.5)*dradius;
      std::cout << "i = " << i << " radius value = " << rvals[i] << std::endl;
      //for (unsigned int j(1); j<=nummass;j++){
      //minchi[i][j] = 1e8;
      //}
    }
    for (unsigned int i(1); i<= nummass; i++){
      mvals[i] = mass_lo + (i-0.5)*dmass;
      std::cout << "i = " << i << " mass value = " << mvals[i] << std::endl;
    }

    for (unsigned int i(1); i<= numbins; i++){
      ivals[i] = inc_lo + (i-0.5)*dinc;
      thetavals[i] = theta_lo + (i-0.5)*dtheta;
      timevals[i] = time_lo + (i-0.5)*dtime;
      rhovals[i] = rho_lo + (i-0.5)*drho;
      tempvals[i] = temp_lo + (i-0.5)*dtemp;
      distancevals[i] = distance_lo + (i-0.5)*ddistance;

      std::cout << "i = " << i << " time value = " << timevals[i] << std::endl;

    }

      fdim.open(dim_file);
      fopt.open(opt_file);
      ffit.open(fit_file);
      frnk.open(rnk_file);


      
      char line[265]; // line of the data file being read in
      //double get_t, get_e, get_f1, get_f2, get_eem, get_et;

      fdim.getline(line,265);
      sscanf( line, "%u %u", &nopt, &npara);
      fdim.close();

      std::cout << "Number of Optimals = " << nopt << std::endl;

      unsigned int rank;
      //double fitness;

      // Create Matrices
      double **optimals = dmatrix(1,nopt,1,npara);
      double *fitness = dvector(1,nopt);
      double *radius = dvector(1,nopt);
      double *mass = dvector(1,nopt);

    


      unsigned int rbin, mbin;
      unsigned maxrank(1);

      double inc, theta, time, rho, temp, distance;
      int ibin,thetabin,timebin,rhobin,tempbin, distancebin;
      int bin;

      for (unsigned int i(1);i<=nopt;i++){ // loop through the optimals

	  frnk.getline(line,265);  
	  sscanf( line, "%u", &rank );

	  ffit.getline(line,265);
	  sscanf( line, "%lf", &fitness[rank]);

	  //fopt.getline(line,265);
	  
	  fopt >> radius[rank];
	  fopt >> mass[rank];
	  
	  for (unsigned int j(1); j <= npara-2; j++){

	    fopt >> optimals[rank][j];

	  }

	  optimals[rank][3] += 0.5;
	  if (optimals[rank][3] > 1.0) optimals[rank][3]-= 1.0;

	   // Compute the bins for each parameter
	  rbin = 1 + (radius[rank] - radius_lo)/dradius;
	  mbin = 1 + (mass[rank] - mass_lo)/dmass;
	  ibin = 1 + (optimals[rank][1] - inc_lo)/dinc;
	  thetabin = 1 + (optimals[rank][2] - theta_lo)/dtheta;
	  timebin = 1 + (optimals[rank][3] - time_lo)/dtime;
	  rhobin = 1 + (optimals[rank][4] - rho_lo)/drho;
	  tempbin = 1 + (optimals[rank][5] - temp_lo)/dtemp;
	  distancebin = 1 + (optimals[rank][6] - distance_lo)/ddistance;

	  if ( rbin*mbin*ibin*thetabin*rhobin*tempbin*distancebin <= 0.0)
	    std::cout << "ERROR bin=0!!!!! " << std::endl;

	  if ( rbin > numradius || mbin > nummass || ibin > numbins || thetabin > numbins || 
	       timebin > numbins || rhobin > numbins || tempbin > numbins || distancebin > numbins){
	      std::cout << "ERROR bin>numbins!!!!! " << std::endl;
	      std::cout << "rbin = " << rbin
		    << " mbin = " << mbin
		    << " ibin = " << ibin
		    << " thetabin = " << thetabin
		    << " timebin = " << timebin
		    << " rhobin = " << rhobin
		    << " tempbin = " << tempbin
		    << " distancebin = " << distancebin
			<< std::endl;


	  }

	  bin = (rbin-1) + 
	    numradius*((mbin-1) + 
		       nummass*((ibin-1) + 
				numbins*((thetabin-1) + 
					 numbins*((timebin-1) + 
						  numbins*((rhobin-1) 
							   + numbins*((tempbin-1) + 
								      numbins*(distancebin-1)))))));
	 
	  /*  std::cout << "i = " << i
		    << " temperature = " << optimals[rank][5]
		    << " tempbin = " << tempbin
		    << " bin=" << bin
		    << " maxbin = " << numradius*nummass*pow(numbins,6)
		    << std::endl;

	  std::cout << "rbin = " << rbin
		    << " mbin = " << mbin
		    << " ibin = " << ibin
		    << " thetabin = " << thetabin
		    << " timebin = " << timebin
		    << " rhobin = " << rhobin
		    << " tempbin = " << tempbin
		    << " distancebin = " << distancebin
		    << std::endl;*/


	  histo[bin] += 1;

	  if (histo[bin] == 1)
	    loglikelihood[bin] = fitness[rank];
	  else{
	    if ( fitness[rank] <= loglikelihood[bin] )
	      loglikelihood[bin] = fitness[rank];
	  }
	  
	  if (fitness[rank] > 1e5){
	  std::cout << std::endl
		    << " i = " << i 
		    << " rank = " << rank 
		    << " radius = " << radius[rank]
		    << " mass = " << mass[rank]
		    << " log(likelihood) = " << fitness[rank] 
		    << std::endl;

	  std::cout << "rbin = " << rbin
		    << " mbin = " << mbin
		    << " ibin = " << ibin
		    << " thetabin = " << thetabin
		    << " timebin = " << timebin
		    << " rhobin = " << rhobin
		    << " tempbin = " << tempbin
		    << " distancebin = " << distancebin
		    << std::endl;

	  std::cout << "bin = " << bin << " log(Likelihood) = " << loglikelihood[bin] << " Histogram[bin] = " << histo[bin] << std::endl;

	  std::cout << "radius[rank] = " << radius[rank]
		    << " rbin = " << rbin
		    << " rvals[rbin] = " << rvals[rbin]
		    << " mass[rank] = " << mass[rank]
		    << " mbin = " << mbin
		    << " mvals[mbin] = " << mvals[mbin]
		    << std::endl;
	  }

	  //if ( rank > maxrank)
	  // maxrank = rank;

      } // Finished reading in the optimals

      std::cout << "Finished reading in the optimals! " << std::endl;



      free_dmatrix(optimals,1,nopt,1,npara);
      free_dvector(fitness,1,nopt);
      free_dvector(radius,1,nopt);
      free_dvector(mass,1,nopt);



      frnk.close();
      fopt.close();
      ffit.close();

      // Now do Metropolis-Hastings search. 

      // Choose initial values

      double r1, m1, i1, th1, ti1, rh1, te1, d1;
      int bin1;
      double ll1;

      double r2, m2, i2, th2, ti2, rh2, te2, d2;
      int bin2;
      double ll2;

      int yes(0);
      double acceptance(0.0);

      r1 = rvals[19];
      m1 = mvals[18];
      i1 = ivals[12];
      th1 = thetavals[12];
      ti1 = timevals[5];
      rh1 = rhovals[7];
      te1 = tempvals[5];
      d1 = distancevals[7];

 
   double r;
   int nsteps = 1000;

   for ( unsigned int i(1); i < nsteps; i++){

      // Compute loglikelihood for initial case.
	   // Compute the bins for each parameter
	  rbin = 1 + (r1 - radius_lo)/dradius;
	  mbin = 1 + (m1 - mass_lo)/dmass;
	  ibin = 1 + (i1 - inc_lo)/dinc;
	  thetabin = 1 + (th1 - theta_lo)/dtheta;
	  timebin = 1 + (ti1 - time_lo)/dtime;
	  rhobin = 1 + (rh1 - rho_lo)/drho;
	  tempbin = 1 + (te1 - temp_lo)/dtemp;
	  distancebin = 1 + (d1 - distance_lo)/ddistance;

	  bin1 = (rbin-1) + 
	    numradius*((mbin-1) + 
		       nummass*((ibin-1) + 
				numbins*((thetabin-1) + 
					 numbins*((timebin-1) + 
						  numbins*((rhobin-1) 
							   + numbins*((tempbin-1) + 
								      numbins*(distancebin-1)))))));

	  ll1 = loglikelihood[bin1];
	  if (ll1 == 0.0) ll1 = 1e4;

	  std::cout << " i = " << i 
		    << " r1 = " << r1
		    << " m1 = " << m1
		    << " i1 = " << i1
		    << " th1 = " << th1
		    << " ti1 = " << ti1 
		    << " rh1 = " << rh1
		    << " te1 = " << te1
		    << " d1 = " << d1
		    << std::endl;

	  std::cout << " i = " << i 
		    << " bin1 = " << bin1
		    << " ll1 = " << ll1 
		    << std::endl;

	  // Create a proposal
	  r2 = NormalDev(r1,dradius*2);
	  m2 = NormalDev(m1,dmass*2);
	  i2 = NormalDev(i1,dinc*2);
	  th2 = NormalDev(th1,dtheta*2);
	  ti2 = NormalDev(ti1,dtime*0.1);
	  rh2 = NormalDev(rh1,drho);
	  te2 = NormalDev(te1,dtemp);
	  d2 = NormalDev(d1,ddistance);

	  // Ensure that 0<t<1
	  if (ti2 > 1.0)
	    ti2 -= 1.0;
	  if (ti2 < 0.0)
	    ti2 += 1.0;


	   // Compute the bins for each parameter
	  rbin = 1 + (r2 - radius_lo)/dradius;
	  mbin = 1 + (m2 - mass_lo)/dmass;
	  ibin = 1 + (i2 - inc_lo)/dinc;
	  thetabin = 1 + (th2 - theta_lo)/dtheta;
	  timebin = 1 + (ti2 - time_lo)/dtime;
	  rhobin = 1 + (rh2 - rho_lo)/drho;
	  tempbin = 1 + (te2 - temp_lo)/dtemp;
	  distancebin = 1 + (d2 - distance_lo)/ddistance;



	  bin2 = (rbin-1) + 
	    numradius*((mbin-1) + 
		       nummass*((ibin-1) + 
				numbins*((thetabin-1) + 
					 numbins*((timebin-1) + 
						  numbins*((rhobin-1) 
							   + numbins*((tempbin-1) + 
								      numbins*(distancebin-1)))))));

	  ll2 = loglikelihood[bin2];
	  if (ll2 == 0.0) ll2 = 1e4;


	  std::cout << " i = " << i 
		    << " r2 = " << r2
		    << " m2 = " << m2
		    << " i2 = " << i2
		    << " th2 = " << th2
		    << " ti2 = " << ti2 
		    << " rh2 = " << rh2
		    << " te2 = " << te2
		    << " d2 = " << d2
		    << std::endl;

	  std::cout << "rbin = " << rbin
		    << " mbin = " << mbin
		    << " ibin = " << ibin
		    << " thetabin = " << thetabin
		    << " timebin = " << timebin
		    << " rhobin = " << rhobin
		    << " tempbin = " << tempbin
		    << " distancebin = " << distancebin
		    << std::endl;

	  std::cout << " i = " << i 
		    << " bin2 = " << bin2
		    << " ll2 = " << ll2 
		    << std::endl;

	  // Draw a random number
	  r = Rand1();

	  std::cout << "r = " << r
		    << " log(r) = " << log(r)
		    << std::endl;

	  if ( ll2 - ll1 < log(r)){
	    std::cout << "Accept new step!" << std::endl;
	    r1 = r2;
	    m1 = m2;
	    i1 = i2;
	    th1 = th2;
	    ti1 = ti2;
	    rh1 = rh2;
	    te1 = te2;
	    d1 = d2;
	    yes += 1;
	    
	  }
	  else{
	    std::cout << "Reject new step!" << std::endl;
	  }
	    
	  acceptance = yes/(1.0*i);
	  std::cout << "Acceptance ratio = " << acceptance << std::endl;

   }




      

    return 0;
 } 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}
