unsigned int BinValue( struct Parameters para, struct Parameters lo, struct Parameters hi, struct Parameters delta, int numradius, int nummass, int numbins);


double Fitness( struct Parameters para, struct Parameters lo, struct Parameters hi, struct Parameters delta, int numradius, int nummass, int numbins,
		long *binrank, double *fit_nz, long *histo);
