/**
 * Initial Data for the EulerFlow.
 * 
 * As we don't have parameters for that in the moment, we use environment
 * variables to distuingish at run time which initial data to load.
 * 
 **/

#include "InitialData.h"
#include "Primitives.h"

#include <stdlib.h>
#include <math.h>
#include <cstring>
using namespace std;

/**
 * Source: MD ADERDG F90 CODE
 **/
void ShuVortex2D(const double* const  x, double* Q, double t=0.0) {
	static const double epsilon = 5.0;
	static const double pi = acos(-1.0);

	double r = sqrt(SQ(x[0]-t-5.)+SQ(x[1]-t-5.));
	double du = epsilon/2./pi*exp(0.5*(1.-r*r))*(5. - x[1] + t);
	double dv = epsilon/2./pi*exp(0.5*(1.-r*r))*(x[0]  - 5.- t);
	double dTemp = -(eos_gamma-1.)*SQ(epsilon)/8./eos_gamma/SQ(pi)*exp(1.-r*r);
	double drho = pow(1.+dTemp, 1./(eos_gamma-1.))-1.;
	double dp   = pow(1.+dTemp, eos_gamma/(eos_gamma-1.))-1.;

	Q[0] = 1. + drho;
	Q[1] = 1. + du;
	Q[2] = 1. + dv;
	Q[3] = 0.0;
	Q[4] = 1. + dp;
}


void gauss(const double* const  x, double* Q) {
    const double GAMMA = 1.4;

    Q[0] = 1.;
    Q[1] = 0.;
    Q[2] = 0.;
    Q[3] = 0.;
#if DIMENSIONS == 2
    Q[4] =
        1. / (GAMMA - 1) +
        exp(-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5)) /
                 (0.05 * 0.05)) *
            1.0e-3;
#else
    Q[4] =
        1. / (GAMMA - 1) +
        exp(-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5) 
                  + (x[2] - 0.5) * (x[2] - 0.5)) / (0.05 * 0.05 * 0.05)) *
            1.0e-3;
#endif  
}

static bool wroteAboutInitialData(false);
#define logInitialData(txt,...) { if(!wroteAboutInitialData) printf(txt); }

void InitialData(const double* const  x, double* Q) {
	const char* id = getenv("EXAHYPE_INITIALDATA");
	logInitialData("Have read '%s'\n", id);
	if(strcmp(id, "ShuVortex")) {
		// ShuVortex
		ShuVortex2D(x, Q);
		logInitialData("Loading ShuVortex Initial Data\n");
	} else {
		// default:
		gauss(x, Q);
		logInitialData("Loading Gaussian Initial Data\n");
	}
	wroteAboutInitialData = true;
}
