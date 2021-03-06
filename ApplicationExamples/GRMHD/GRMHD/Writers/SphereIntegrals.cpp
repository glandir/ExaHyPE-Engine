// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "SphereIntegrals.h"

#include "AbstractGRMHDSolver_ADERDG.h"
#include "kernels/GaussLegendreBasis.h" // kernels::legendre::interpolate
#include <cmath> // std::sin, std::cos, std::sqrt
#include "peano/utils/Loop.h" // dfor; not yet used
using namespace tarch::la;
using namespace std;
typedef Vector<DIMENSIONS, double> dvec;

#include "Fortran/MassAccretionRate.h"


GRMHD::SphereIntegrals::SphereIntegrals( ) :
	exahype::plotters::ADERDG2UserDefined::ADERDG2UserDefined(),
	spherewriter("output/sphere-stuff")
{
  // @TODO Please insert your code here.
}


void GRMHD::SphereIntegrals::plotPatch(const dvec& offsetOfPatch, const dvec& sizeOfPatch, double* u, double timeStamp) {
	// First comment: This is for 3D. This will not work in 2D as the GRMHD AccretionDisk2D itself
	// uses spherical coordinates. We assume cartesian coordinates here in the simulation domain.
	const int nVar  = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
	const int order = GRMHD::AbstractGRMHDSolver_ADERDG::Order;

	// Compute on a sphere:
	const double rSphere = 2.1;
	
	// How much points to compute
	const int ntheta = 10;
	const int nphi = ntheta;

	// Here, we just determine the cartesian positions of all coordinates on the sphere.
	// This is dumb and inefficient, and there should be intelligent checks instead to
	// reduce the loops.
	
	// Comment 2: Michael mentioned there are much more intelligent ways, ie. with gauss
	// points on the sphere and so.
	
	// Comment 3: This should be tested with integration of a constant field (value 1) which
	// should give the value 4*M_PI*rSphere. Also, this should be tested with the EulerFlow
	// standard Gauss test to see whether we get a nice curve.
	
	const double scaling = 1./(2.*M_PI*ntheta) * 1./(2.*M_PI*nphi);
	dvec ip; // integration point
	for(int itheta=0; itheta<ntheta; itheta++) {
		for(int iphi=0; iphi<nphi; iphi++) {
			// compute the position of this integration point
			ip(0) = rSphere * sin(2*M_PI*itheta/ntheta) * cos(2*M_PI*iphi/nphi);
			ip(1) = rSphere * sin(2*M_PI*itheta/ntheta) * sin(2*M_PI*iphi/nphi);
			ip(2) = rSphere * cos(2*M_PI*itheta/ntheta);
			
			// check if it is inside the current cell
			bool isinside = true;
			for(int d=0; d<DIMENSIONS; d++) {
				isinside = isinside && offsetOfPatch(d) < ip(d) && ip(d) < (offsetOfPatch(d)+sizeOfPatch(d));
			}
			
			if(isinside) {
				// for the time being, interpolate all quantities
				double intp[nVar];
				for(int k=0; k<nVar; k++) {
					intp[k] = kernels::legendre::interpolate(offsetOfPatch.data(), sizeOfPatch.data(), ip.data(), nVar, k, order, u);
				}

				// map the interpolated values to something
				double masschange;
				double vx;
				double vy;
				double vz;								
				double ur;
				massaccretionrate_(intp, &masschange, &vx, &vy, &vz);
				// reduce the scalar field on the sphere.

				ur = sin(2*M_PI*itheta/ntheta) * cos(2*M_PI*iphi/nphi)*vx
				   + sin(2*M_PI*itheta/ntheta) * sin(2*M_PI*iphi/nphi)*vy
				   + cos(2*M_PI*itheta/ntheta)*vz
				  - ((2.0/rSphere)/(1.0+2.0/rSphere))/sqrt(1.0/(1.0+2.0/rSphere));

                                masschange = rSphere * rSphere * masschange*ur;
				
				// the "sum" entry in the reductions ASCII file is the integral value.
				spherewriter.addValue(masschange, scaling);
			}
		}
	}

	// For debugging purposes, we should write out this 
	// test at to an ASCII file to ensure whether it works.
}


void GRMHD::SphereIntegrals::startPlotting( double time) {
  spherewriter.startRow(time);
}


void GRMHD::SphereIntegrals::finishPlotting() {
  spherewriter.finishRow();
}
