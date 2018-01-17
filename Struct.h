/***************************************************************************************/
/*                                    Struct.h

    This code creates the many data storage structures used in other pieces of code.
    
*/
/***************************************************************************************/

#ifndef STRUCT_H
#define STRUCT_H

#include <exception>
#include <vector>
#include <float.h>


struct Parameters {
  double radius;
  double mass;
  double inc;
  double theta;
  double time;
  double rho;
  double temp;
  double distance;
};

struct Prob1d {

  double median;
  double xl1;
  double xl2;
  double xl3;
  double xr1;
  double xr2;
  double xr3;
  double pmed;
  double pl1;
  double pl2;
  double pl3;
  double pr1;
  double pr2;
  double pr3;
  
};

struct Prob2d{
  double sigma1;
  double sigma2;
  double sigma3;
};





#endif // STRUCT_H
