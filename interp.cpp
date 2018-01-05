/***************************************************************************
 * interp.c
 *
 * This file contains a collection of useful interpolation routines.
 *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interp.h"
#include "Struct.h"
#include "nrutil.h"

/* Do a linear interpolation between two points. */

double * interplin(double *xp, double **yp, int numangles, double xb, int *n_nearest_pt){

  int k;
  double *yint;

  yint = dvector(0,3*numangles+1);

  k=*n_nearest_pt;

  for (unsigned int i(0);i<=3*numangles;i++){

    yint[i] =  yp[k][i] + (xb-xp[k])*(yp[k+1][i]-yp[k][i])/(xp[k+1]-xp[k]);

  }

  //return yp[k] + (xb-xp[k])*(yp[k+1]-yp[k])/(xp[k+1]-xp[k]);
 
  return (yint);
}

/* Polynomial interpolation.  This is based on the interpolation routine in
 * numerical recipes.  If P(x) is the polynomial of order order-1 passing 
 * through the points (xp[i],yp[i]), then polint returns P(xb) and err points
 * to an error estimate.  Note xp and yp are arrays indexed from 1 to n>=order.
 */
double polint(double *xp, double *yp, int order, double xb, double *err){
  int i, m, ns=1, j;
  double tmp, diff, den, dnum, cnum, yb;
  double *c, *d;

  c = (double *) malloc((order+1)*sizeof(double));
  d = (double *) malloc((order+1)*sizeof(double));
  diff = fabs(xb-xp[1]);
  for (i=1; i<=order; i++){
    if ( (tmp=fabs(xb-xp[i])) < diff ){
      ns=i;
      diff = tmp;
    }
    c[i] = yp[i];
    d[i] = yp[i];
  }

  yb = yp[ns--];
  for (m=1; m<order; m++){
    for (i=1; i<=order-m; i++){
      cnum = xp[i]-xb;
      dnum = xp[i+m]-xb;
      if ( (den=cnum-dnum) == 0.0 ){
        /*Two values of xp are equal:no polynomial passes through the points.*/
        printf("error in polint: xp[%d]==xp[%d]\n", i, i+m);
	//printf(" xb = %lf ", xb); 
	//for (j=1; j<=order; j++){	
	//	printf(" xp[%d]=%lf",j,xp[j]);
        //}			
        //exit(1);
      }
      tmp = (c[i+1]-d[i])/den;
      c[i] = cnum*tmp;
      d[i] = dnum*tmp;
    }
    *err = 2*ns<order-m ? c[ns+1] : d[ns--];
    yb += *err;
  }
  free(c);
  free(d);
  return yb;
}

/* Do a 1-dimensional interpolation. */
double interp1(double *xp, double *yp, int first, int np, double xb, int *x_nearest_pt){
  const int order = 4; /*Order of interpolation*/
  int lo;

  double err;
  double ans;

  hunt_surf(xp,first,np,xb,x_nearest_pt);

  lo=IMIN(IMAX(*x_nearest_pt-(order-1)/2, 1), np+1-order);

  ans = polint(&xp[lo-1], &yp[lo-1], order, xb, &err);

  //printf("ans = %g err=%g \n",ans,err);

  return ans;

}


// hunt_surf searches the region outside of the star to find the appropriate value
void hunt_surf(double xx[], int first, int n, double x, int *jlo)
{ 
	int jm,jhi,inc,ascnd;

	if ( x >= xx[*jlo] &&  x <= xx[*jlo+1]) {
	  //printf("jlo = %d Perfect! \n",*jlo);
	  //printf("x = %g  ; xx[jlo] = %g ; xx[jlo+1] = %g \n",
	  //	 x, xx[*jlo], xx[*jlo+1]);
	  return;
	}
	else{

	ascnd=(xx[n] > xx[first]);
	if (*jlo <= 0 || *jlo > n) {
		*jlo=first-1;
		jhi=n+1;
	} else {
		inc=1;
		if ((x >= xx[*jlo]) == ascnd) {
			if (*jlo == n) return;
			jhi=(*jlo)+1;
			while ((x >= xx[jhi]) == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n) {
					jhi=n+1;
					break;
				}
			}
		} else {
			if (*jlo == first) {
				*jlo=first-1;
				return;
			}
			jhi=(*jlo);
			*jlo -= 1;
			while ((x < xx[*jlo]) == ascnd) {
				jhi=(*jlo);
				inc += inc;
				*jlo=jhi-inc;
				if (*jlo < 1) {
					*jlo=0;
					break;
				}
			}
		}
	}

	//	printf("x = %g jhi = %d xx[jhi] = %g \n", x, jhi, xx[jhi]);

	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if ((x > xx[jm]) == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
	}
}











/*C*/










/* Polynomial Extrapolation.  This is based on the interpolation routine in
 * numerical recipes.  If P(x) is the polynomial of order order-1 passing 
 * through the points (xp[i],yp[i]), then polint returns P(xb) and err points
 * to an error estimate.  Note xp and yp are arrays indexed from 1 to n>=order.
 */
double extrap(double *xp, double *yp, int order, double xb, double *err){
  int i, m, ns;
  double tmp, den, dnum, cnum, yb;
  double *c, *d;

  c = (double *) malloc((order+1)*sizeof(double));
  d = (double *) malloc((order+1)*sizeof(double));

  ns = order;


  for (i=1; i<=order; i++){
    c[i] = yp[i];
    d[i] = yp[i];
  }

  yb = yp[ns--];
  for (m=1; m<order; m++){
    for (i=1; i<=order-m; i++){
      cnum = xp[i]-xb;
      dnum = xp[i+m]-xb;
      if ( (den=cnum-dnum) == 0.0 ){
        /*Two values of xp are equal:no polynomial passes through the points.*/
        printf("error in polint: xp[%d]==xp[%d]\n", i, i+m);
        exit(1);
      }
      tmp = (c[i+1]-d[i])/den;
      c[i] = cnum*tmp;
      d[i] = dnum*tmp;
    }
    *err = 2*ns<order-m ? c[ns+1] : d[ns--];
    yb += *err;
  }
  free(c);
  free(d);
  return yb;
}
