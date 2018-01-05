#ifndef INTERP_H
#define INTERP_H

#include "Struct.h"

#define IMAX(imaxarg1,imaxarg2) ((imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

#define IMIN(iminarg1,iminarg2) ((iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

double polint(double *xp, double *yp, int order, double xb, double *err);

double polint2(double *xp, double *yp, int order, double xb, double *err);


double printpolint(double *xp, double *yp, int order, double xb, double *err);

double interp1(double *xp, double **yp, int first, int np, double xb, int *x_nearest_pt);

void hunt_surf(double xx[], int first, int n, double x, int *jlo);

double interp2(double *, double *, double **, int, int, double, double, int *,
               int *);

double * interplin(double *, double **yp, int, double, int *);

double interp_b(double *, double **, int, int, double, int *);

/* New version of interp2 Do a 2-dimensional interpolation.  */

double interp2A(double *xp, double *yp, double **zp, 
		double xb, double yb, int xlo, int ylo, int max_order);




double interpB(double xp[], 
              double yp[], 
              double    xb ,
	       int    xlo);


void bcucof(double *y,
	    double *y1,
	    double *y2,
	    double *y12,
	    double d1,
	    double d2,
	    double **c);


/* Copyright (C) 1987,1988 Numerical Recipes Software -- BCUINT */


void bcuint(double *y,
	    double *y1,
	    double *y2,
	    double *y12,
	    double x1l,
	    double x1u,
	    double x2l,
	    double x2u,
	    double x1,
	    double x2,
	    double *ansy,
	    double *ansy1,
	    double *ansy2);


void pot_interp2(  double **pot, 
		   double **pot_s, 
		   double **pot_m, 
		   double **pot_ms,
		   double ss, 
		   double mm, 
		   double *ans, 
		   double *ans_s, 
		   double *ans_m);


double extrap(double *xp, double *yp, int order, double xb, double *err);

double interp_equat(double **pot, double ss);

#endif
