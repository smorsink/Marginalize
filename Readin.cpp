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
#include "Struct.h"
#include "bayes.h"
#include "fitness.h"


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
  char directory[256] = "WRONG";
  		
  unsigned int nopt(200000),
    npara, rbin, mbin, ibin, thetabin, timebin, rhobin, tempbin, distancebin,
    nummass(40),numradius(40),numbins(12);

  int nzrank;

  unsigned int chainlength(10);

  // Non-zero values for fitness
  double *fit_nz = dvector(1,nopt);
  long int *bin_nz = ivector(1,nopt);

  double
    radius_lo(7.5),
    radius_hi(21.5),
    mass_lo(0.8),
    mass_hi(2.25);

  double 
    inc_lo(60.0), inc_hi(125.0);
  double
    theta_lo(50.0),theta_hi(125.0);
  double 
    time_lo(0.20),time_hi(0.80);
  double
    rho_lo(0.1),rho_hi(1.6);
  double
    temp_lo(0.1),temp_hi(0.52);
  double
    distance_lo(0.1),distance_hi(0.6);

  struct Parameters para, lo, hi, delta;

  double min_loglik(10000000), max_loglik(-1000000);
  double fudge(1.0);
 
  // Initialize lowest values of all parameters.
  lo.radius=radius_lo;
  lo.mass = mass_lo;
  lo.inc = inc_lo;
  lo.theta = theta_lo;
  lo.time = time_lo;
  lo.rho = rho_lo;
  lo.temp = temp_lo;
  lo.distance = distance_lo;

  // Initialize largest values of all parameters.
  hi.radius=radius_hi;
  hi.mass = mass_hi;
  hi.inc = inc_hi;
  hi.theta = theta_hi;
  hi.time = time_hi;
  hi.rho = rho_hi;
  hi.temp = temp_hi;
  hi.distance = distance_hi;


    struct Parameters para1,para2;
    double ll1,ll2;

    //Initialization of starting values for MT loop
    para1.radius = 8.0;
    para1.mass = 0.8;
    para1.inc = 70.1;
    para1.theta = 50.0;
    para1.time = 0.473;
    para1.rho = 0.72;
    para1.temp = 0.437;
    para1.distance = 0.31;


  /*********************************************************/
  /* READING IN PARAMETERS FROM THE COMMAND LINE ARGUMENTS */
  /*********************************************************/
    
    for ( int i(1); i < argc; i++ ) {
        if ( argv[i][0] == '-' ) {  // the '-' flag lets the computer know that we're giving it information from the cmd line
            switch ( argv[i][1] ) {

	    case 'c': // Number of Mass bins
	      sscanf(argv[i+1], "%u", &chainlength);
	      break;

	    case 'j': // Number of Mass bins
	      sscanf(argv[i+1], "%u", &nummass);
	      break;

	    case 'i':
	      sscanf(argv[i+1], "%u", &numradius);
	      break;
	      
	    case 'k':
	      sscanf(argv[i+1], "%u", &numbins);
	      break;
	      
	    case 'b':
	      sscanf(argv[i+1], "%lf", &fudge);
	      break;
	      
	    case 'a': // Name of input file directory
	      sscanf(argv[i+1], "%s", directory);
	      break;

	    case 'o':  // Name of output file
	                sscanf(argv[i+1], "%s", out_file);
	                break;

	    case 'r':
	      sscanf(argv[i+1], "%lf", &para1.radius);
	      break;			

	    case 'm':
	      sscanf(argv[i+1], "%lf", &para1.mass);
	      break;

	    case 'n':
	      sscanf(argv[i+1], "%lf", &para1.inc);
	      break;

	    case 'p':
	      sscanf(argv[i+1], "%lf", &para1.theta);
	      break;

	    case 'q':
	      sscanf(argv[i+1], "%lf", &para1.time);
	      break;

	    case 's':
	      sscanf(argv[i+1], "%lf", &para1.rho);
	      break;

	    case 't':
	      sscanf(argv[i+1], "%lf", &para1.temp);
	      break;
			
	    case 'u':
	      sscanf(argv[i+1], "%lf", &para1.distance);
	      break;

	    

	                
            } // end switch	
        } // end if
    } // end for


    std::cout << "Number of bins = " << numradius*nummass*pow(numbins,6) << std::endl;

    long int *histo = ivector(0,numradius*nummass*pow(numbins,6));
    long int *binrank = ivector(0,numradius*nummass*pow(numbins,6));

	//std::cout << "histo[500] = " << histo[500] << std::endl;

    
    // Create the Mass and Radius Intervals, and intervals for other parameters
    double *rvals = dvector(0,numradius);    
    double *mvals = dvector(0,nummass);
    double *ivals = dvector(0,numbins);
    double *thetavals = dvector(0,numbins);
    double *timevals = dvector(0,numbins);
    double *rhovals = dvector(0,numbins);
    double *tempvals = dvector(0,numbins);
    double *distancevals = dvector(0,numbins);

    // Create histograms for each parameter
    long int *rhist = ivector(0,numradius);
    long int *mhist = ivector(0,nummass);
    long int *ihist = ivector(0,numbins);
    long int *thetahist = ivector(0,numbins);
    long int *timehist = ivector(0,numbins);
    long int *rhohist = ivector(0,numbins);
    long int *temphist = ivector(0,numbins);
    long int *distancehist = ivector(0,numbins);

    delta.radius = (radius_hi - radius_lo)/(numradius);
    delta.mass = (mass_hi - mass_lo)/(nummass);
    delta.inc = (inc_hi - inc_lo)/numbins;
    delta.theta = (theta_hi - theta_lo)/numbins;
    delta.time = (time_hi - time_lo)/numbins;
    delta.rho = (rho_hi - rho_lo)/numbins;
    delta.temp = (temp_hi - temp_lo)/numbins;
    delta.distance = (distance_hi - distance_lo)/numbins;


    for (unsigned int i(1); i<= numradius; i++){
      //rvals[i] = radius_lo + (i-0.5)*dradius;
      rvals[i] = lo.radius + (i-0.5)*delta.radius;
    }
    for (unsigned int i(1); i<= nummass; i++){
      mvals[i] = lo.mass + (i-0.5)*delta.mass;
    }

    for (unsigned int i(1); i<= numbins; i++){
      ivals[i] = lo.inc + (i-0.5)*delta.inc;
      thetavals[i] = lo.theta + (i-0.5)*delta.theta;
      timevals[i] = lo.time + (i-0.5)*delta.time;
      rhovals[i] = lo.rho + (i-0.5)*delta.rho;
      tempvals[i] = lo.temp + (i-0.5)*delta.temp;
      distancevals[i] = lo.distance + (i-0.5)*delta.distance;
    }

    // Number of non-zero fitness values
    int nz(0);

    // Loop through the data sets
    for (unsigned int k(1); k<=7; k++){

      // Set-up the input files.

      sprintf(dim_file, "%s%d/%s", directory,k,"dimensions.txt");
      std::cout << "dim_file = " << dim_file << std::endl;
      sprintf(opt_file, "%s%d/%s", directory,k,"optimals.txt");
      std::cout << "opt_file = " << opt_file << std::endl;
      sprintf(fit_file, "%s%d/%s", directory,k,"fitness.txt");
      std::cout << "fit_file = " << fit_file << std::endl;

      fdim.open(dim_file);
      fopt.open(opt_file);
      ffit.open(fit_file);
     
	std::cout << "Opened Data File: " << k << std::endl;

      char line[265]; // line of the data file being read in

      fdim.getline(line,265);
      sscanf( line, "%u %u", &nopt, &npara);
      fdim.close();

      std::cout << "Number of Optimals = " << nopt << std::endl;    

      double fitness;

      //for (unsigned int i(1);i<=1;i++){ // loop through the optimals
      for (unsigned int i(1);i<=nopt;i++){ // loop through the optimals

	//std::cout <<std::endl << "i = " << i << std::endl;

	ffit.getline(line,265);
	sscanf( line, "%lf", &fitness);
	fitness *= -1.0;
	  
	fopt >> para.radius;
	fopt >> para.mass;
	fopt >> para.inc;
	fopt >> para.theta;
	fopt >> para.time;
	fopt >> para.rho;
	fopt >> para.temp;
	fopt >> para.distance;

	para.time += 0.5;
	if (para.time > 1.0) para.time-= 1.0;

	unsigned int bin = BinValue(para,lo,hi,delta,numradius,nummass,numbins);
	//std::cout << " bin = " << bin << std::endl;

	histo[bin] += 1;
	//std::cout << " histo[" << bin << "] = " << histo[bin] << std::endl;

	if (histo[bin] == 1){
	    // We have a new non-zero fitness bin!
	    nz += 1;
	    bin_nz[nz] = bin;
	    fit_nz[nz] = fitness;
	    binrank[bin] = nz;

	/*std::cout << "New Bin: fitness = " << fitness 
		    << " radius = " << para.radius 
		    << " mass = " << para.mass
		    << " inc = " << para.inc
		    << " theta = " << para.theta
		    << " time = " << para.time
		    << " rho = " << para.rho
		    << " temperature = " << para.temp
		    << " distance = " << para.distance
		    << std::endl; */

	}
	else{
	  nzrank = binrank[bin];
	  if ( fitness <= fit_nz[nzrank] ){
	    fit_nz[nzrank] = fitness;
	  }
	}

	if ( para.radius > 12.0 && para.radius < 12.5){
	  std::cout << "Fitness: fitness = " << fitness 
		    << " radius = " << para.radius 
		    << " mass = " << para.mass
		    << " inc = " << para.inc
		    << " theta = " << para.theta
		    << " time = " << para.time
		    << " rho = " << para.rho
		    << " temperature = " << para.temp
		    << " distance = " << para.distance
		    << std::endl;
	}

	if (fitness < min_loglik){
	  min_loglik = fitness;
	  std::cout << "Lowest Log(likelihood) = " << min_loglik << std::endl;
	  std::cout << " radius = " << para.radius << " mass = " << para.mass << " fitness = " << fitness <<std::endl;
std::cout << "Fitness: radius = " << para.radius 
	    << " mass = " << para.mass
	    << " inc = " << para.inc
	    << " theta = " << para.theta
	    << " time = " << para.time
	    << " rho = " << para.rho
	    << " temperature = " << para.temp
	    << " distance = " << para.distance
	    << std::endl;

	}
	if (fitness > max_loglik){
	  max_loglik = fitness;
	  std::cout << "               Largest Log(likelihood) = " << max_loglik << std::endl;
	  std::cout << "               radius = " << para.radius << " mass = " << para.mass << " fitness = " << fitness <<std::endl;
	  std::cout << " inc = " << para.inc
		    << " theta = " << para.theta
		    << " time = " << para.time
		    << " rho = " << para.rho
		    << " temp = " << para.temp
		    << " distance = " << para.distance 
		    << std::endl;

	}

      } // Finished reading in the optimals

      std::cout << "Finished reading in the optimals! " << std::endl;

      std::cout << "Number of optimals = " << nopt << std::endl;

      std::cout << "number of non-zeros = " << nz << std::endl;


      frnk.close();
      fopt.close();
      ffit.close();

    } // End of the for-k-loop

 
    // Initialize Metropolis-Teller Algorithm


    double fudge2(0.82);

    long int **histogram = imatrix(0,numradius,0,nummass);

    out.open("histtestRM.txt");    
    for (unsigned int i(1);i<=numradius;i++){ //Uncomment this line if you want to run the loop
      for (unsigned int j(1);j<=nummass;j++){
	histogram[i][j] = 0;
	out << rvals[i] << " " 
	    << mvals[j] << " " 
	    << histogram[i][j] << " "
		<< i << " " << j 
	    << std::endl;
      }
    }
    out.close();

      ll1 = Fitness(para1,lo,hi,delta,numradius,nummass,numbins, binrank, fit_nz, histo);
   
    std::cout <<"Step 1: \t"
		  << para1.radius << "\t"
		  << para1.mass << "\t"
		  << para1.inc << "\t"
		  << para1.theta << "\t"
		  << para1.time << "\t"
		  << para1.rho << "\t"
		  << para1.temp << "\t"
		  << para1.distance << "\t"
		  << ll1 << "\t"
		  << std::endl;

    //unsigned int nsteps(2e10);
    double r;
    long int yes(0);
    double acceptance(0.0);
    //double fudge(1.0);

    out.open("trace.txt");


    std::cout << "chainlength = " << chainlength << std::endl;

    for (unsigned int i(1);i<chainlength;i++){

      //std::cout << "i = " << i << std::endl;

      // std::cout 
      //<< "MT Step: fitness1 = " << ll1 << std::endl;

      // Create a proposal
      para2.radius = NormalDev(para1.radius,delta.radius*fudge);
      para2.mass = NormalDev(para1.mass,delta.mass*fudge);
      para2.inc = NormalDev(para1.inc,delta.inc*fudge*fudge2);
      para2.theta = NormalDev(para1.theta,delta.theta*fudge*fudge2);
      para2.time = NormalDev(para1.time,delta.time*fudge*fudge2);
      para2.rho = NormalDev(para1.rho,delta.rho*fudge*fudge2);
      para2.temp = NormalDev(para1.temp,delta.temp*fudge*fudge2);
      para2.distance = NormalDev(para1.distance,delta.distance*fudge*fudge2);

      //std::cout << " Proposal: " << std::endl;
      ll2 = Fitness(para2,lo,hi,delta,numradius,nummass,numbins, binrank, fit_nz, histo);
      //std::cout << " Log(L2) = " << ll2 << std::endl;

      // Find a random number
      r = Rand1();

      if ( ll2 - ll1 > log(r)){ // Accept the new step!

	rbin = 1 + (para2.radius - lo.radius)/delta.radius;
	mbin = 1 + (para2.mass - lo.mass)/delta.mass;
	ibin = 1 + (para2.inc - lo.inc)/delta.inc;
	thetabin = 1 + (para2.theta - lo.theta)/delta.theta;
	timebin = 1 + (para2.time - lo.time)/delta.time;
	rhobin = 1 + (para2.rho - lo.rho)/delta.rho;
	tempbin = 1 + (para2.temp - lo.temp)/delta.temp;
	distancebin = 1 + (para2.distance - lo.distance)/delta.distance;
  
	histogram[rbin][mbin] += 1;
	rhist[rbin] += 1;
	mhist[mbin] += 1;
	ihist[ibin] += 1;
	thetahist[thetabin] += 1;
	timehist[timebin] += 1;
	rhohist[rhobin] += 1;
	temphist[tempbin] += 1;
	distancehist[distancebin] += 1;
	

	para1 = para2;
	ll1 = ll2;
	    
	yes += 1;
	acceptance = yes/(1.0*i);
	  
	

	if ( i%1000 == 0){
	      out << i <<"\t"
		  << para1.radius << "\t"
		  << para1.mass << "\t"
		  << para1.inc << "\t"
		  << para1.theta << "\t"
		  << para1.time << "\t"
		  << para1.rho << "\t"
		  << para1.temp << "\t"
		  << para1.distance << "\t"
		  << ll1 << "\t"
		  << acceptance 
		  << std::endl;




	    }
	    
	  }
	  else{
	    //std::cout << "Reject new step!" << std::endl;
	  }
	    
    }
    out.close();

	// Normalize the probability functions to 1.0
	

    long int totprob = 0;
    double **probRM = dmatrix(0,numradius,0,nummass);
    double *probR = dvector(0,numradius);
    double *probM = dvector(0,numradius);
    double *probi = dvector(0,numbins);
    double *probtheta = dvector(0,numbins);
    double *probtime = dvector(0,numbins);
    double *probrho = dvector(0,numbins);
    double *probtemp = dvector(0,numbins);
    double *probdistance = dvector(0,numbins);

    for (unsigned int i(1);i<=numradius;i++){ 
      for (unsigned int j(1);j<=nummass;j++){
	
	totprob += histogram[i][j];

      }
    }
	std::cout << "Total RM probability = " << totprob << std::endl;
    for (unsigned int i(1);i<=numradius;i++){ 
      for (unsigned int j(1);j<=nummass;j++){
	
	 probRM[i][j] = histogram[i][j]/(1.0*totprob);
	
      }
    }

    for (unsigned int i(1);i<=numradius;i++){ 
	probR[i] = rhist[i]/(1.0*totprob);
    }
    for (unsigned int i(1);i<=numradius;i++){ 
	probM[i] = mhist[i]/(1.0*totprob);
    }
    for (unsigned int i(1);i<=numbins;i++){ 
	probi[i] = ihist[i]/(1.0*totprob);
	probtheta[i] = thetahist[i]/(1.0*totprob);
	probtime[i] = timehist[i]/(1.0*totprob);
	probrho[i] = rhohist[i]/(1.0*totprob);
	probtemp[i] = temphist[i]/(1.0*totprob);
	probdistance[i] = distancehist[i]/(1.0*totprob);
    }


    std::cout << " Open Output files " << std::endl;
    out.open("histRM.txt");    
    for (unsigned int i(1);i<=numradius;i++){ //Uncomment this line if you want to run the loop
      for (unsigned int j(1);j<=nummass;j++){
	out << rvals[i] << " " 
	    << mvals[j] << " " 
	    << probRM[i][j] << " "
		<< i << " " << j 
	    << std::endl;
      }
    }
    out.close();

    out.open("histR.txt");    
    for (unsigned int i(1);i<=numradius;i++){ //Uncomment this line if you want to run the loop
	out << rvals[i] << " " 
	    << probR[i] << " "
	    << std::endl;
    }
    out.close();

    out.open("histM.txt");    
    for (unsigned int i(1);i<=nummass;i++){ //Uncomment this line if you want to run the loop
	out << mvals[i] << " " 
	    << probM[i] << " "
	    << std::endl;
    }
    out.close();


        out.open("histi.txt");    
    for (unsigned int i(1);i<=numbins;i++){ //Uncomment this line if you want to run the loop
	out << ivals[i] << " " 
	    << probi[i] << " "
	    << std::endl;
    }
    out.close();

        out.open("histtheta.txt");    
    for (unsigned int i(1);i<=numbins;i++){ //Uncomment this line if you want to run the loop
	out << thetavals[i] << " " 
	    << probtheta[i] << " "
	    << std::endl;
    }
    out.close();

        out.open("histtime.txt");    
    for (unsigned int i(1);i<=numbins;i++){ //Uncomment this line if you want to run the loop
	out << timevals[i] << " " 
	    << probtime[i] << " "
	    << std::endl;
    }
    out.close();

        out.open("histrho.txt");    
    for (unsigned int i(1);i<=numbins;i++){ //Uncomment this line if you want to run the loop
	out << rhovals[i] << " " 
	    << probrho[i] << " "
	    << std::endl;
    }
    out.close();

        out.open("histtemp.txt");    
    for (unsigned int i(1);i<=numbins;i++){ //Uncomment this line if you want to run the loop
	out << tempvals[i] << " " 
	    << probtemp[i] << " "
	    << std::endl;
    }
    out.close();

        out.open("histdistance.txt");    
    for (unsigned int i(1);i<=numbins;i++){ //Uncomment this line if you want to run the loop
	out << distancevals[i] << " " 
	    << probdistance[i] << " "
	    << std::endl;
    }
    out.close();
    

      // Note this loop takes a long time to execute. In most cases you shouldn't run it.
      //out.open(out_file,std::ios_base::app);
      for (unsigned int i(1);i<1;i++){ // This skips the loop
    //for (unsigned int i(1);i<=numradius;i++){ //Uncomment this line if you want to run the loop
	for (unsigned int j(1);j<=nummass;j++){
	  for (unsigned int k(1);k<=numbins;k++){
	    for (unsigned int l(1);l<=numbins;l++){
	      for (unsigned int m(1);m<=numbins;m++){
		for (unsigned int n(1);n<=numbins;n++){
		  for (unsigned int o(1);o<=numbins;o++){
		    for (unsigned int p(1);p<=numbins;p++){

		    int  bin = + 1 + (i-1) + 
			numradius*((j-1) + 
				   nummass*((k-1) + 
					    numbins*((l-1) + 
						     numbins*((m-1) + 
							      numbins*((n-1) 
								       + numbins*((o-1) + 
										  numbins*(p-1)))))));

		      if (histo[bin] > 0){

			nzrank = binrank[bin];
			double fit = fit_nz[nzrank];
	      
			std::cout << rvals[i] << "\t"
				<< mvals[j] << "\t"
				<< ivals[k] << "\t"
				<< thetavals[l] << "\t"
				<< timevals[m] << "\t"
				<< rhovals[n] << "\t"
				<< tempvals[o] << "\t"
				<< distancevals[p] << "\t"
			  // << loglikelihood[bin] << "\t"
				  << fit << "\t"
				  << bin << "\t"
				  << histo[bin]
				<< std::endl;
		      

		      }
		      else{
			double fit = 1e5;
			std::cout << rvals[i] << "\t"
				  << mvals[j] << "\t"
				  << ivals[k] << "\t"
				  << thetavals[l] << "\t"
				  << timevals[m] << "\t"
				  << rhovals[n] << "\t"
				  << tempvals[o] << "\t"
				  << distancevals[p] << "\t"
			  // << loglikelihood[bin] << "\t"
				  << fit << "\t"
				  << bin << "\t"
				  << histo[bin]
				  << std::endl;

		      }

		    }
		  }
		}
	      }
	    }
	  }
	}
      } // End of for-i-loop

      
      

    return 0;
 } 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}
