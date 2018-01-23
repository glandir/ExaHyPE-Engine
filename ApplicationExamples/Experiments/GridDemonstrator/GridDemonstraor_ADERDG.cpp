#include "GridDemonstraor_ADERDG.h"

#include "GridDemonstraor_ADERDG_Variables.h"

#include <math.h>
#include <random>


tarch::logging::Log GridDemonstrator::GridDemonstraor_ADERDG::_log( "GridDemonstrator::GridDemonstraor_ADERDG" );


void GridDemonstrator::GridDemonstraor_ADERDG::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

bool GridDemonstrator::GridDemonstraor_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,const int numberOfObservables,
      const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t, const double dt) const {

	// admissible means the solution is acceptable as being a physical solution.
	bool isAdmissible;
	
	typedef tarch::la::Vector <DIMENSIONS,double> vec;
	
	vec abscenter = tarch::la::abs(center);
	
	// lower left, upper right radius of cell
	double l = tarch::la::norm2(abscenter - dx/2.0);
	double r = tarch::la::norm2(abscenter + dx/2.0);
	// center radius of cell
	double c = tarch::la::norm2(center);

	// works well with maximum-mesh-depth = 3
	// and grid size = [-20,20]
	bool singleNS = true;
	if(singleNS) {
		constexpr double radius_shell = 13; // limit around this shell
		isAdmissible = (l > radius_shell) || (r <= radius_shell);
		//printf("Cell has l=%f,r=%f => isAdmissible=%s\n", l, r, isAdmissible?"true":"false");
	}
	
	// a single NS with some atmosphere thickness
	bool atmoNS = true;
	if(atmoNS) {
		double rnd = ((double) rand() / (RAND_MAX)); // random \in [0,1]
		
		// limit around this shell:
		double radius_shell = 10 + rnd * 4;
		double atmo_thickness = 2;
		isAdmissible = (l > radius_shell+atmo_thickness/2) || (r <= radius_shell-atmo_thickness/2);
		//printf("Cell has l=%f,r=%f => isAdmissible=%s\n", l, r, isAdmissible?"true":"false");
	}
	
	// works well with maximum-mesh-depth = 5
	// and grid size = 20
	bool singleBH = false;
	if(singleBH) {
		constexpr double horizon_size = .5;
		isAdmissible = c > horizon_size;
	}
	
	// two points rotating around the center, roughly representing
	// the inspiral phase of two NS.
	// Status: Does not work in the moment. There is a bug somewhere.
	bool inspiralNS = false;
	constexpr double pi = M_PI;
	constexpr double radius = 6; // size of each NS
	if(inspiralNS) {
		vec pol[2];  // centers in polar coordinates {r,phi}

		double rstart = 15, rziel = 2, tstart = 0, tziel = 10*dt;
		double r = rstart + (rstart - rziel) / (tstart - tziel) * t;
		
		double runden = 5, T = tziel / runden, omega = 2*pi / T, offset = pi;
		double phi = omega * t;

		pol[0] = r, phi;
		pol[1] = r, phi + offset;
		
		bool Admissible[2];
		for(int i=0; i<2; i++) {
			double r=pol[i][0], phi = pol[i][1];
			
			// determine center of star in cartesian coordinates
			vec star;
			star[0] = r * std::cos(phi); // x
			star[1] = r * std::sin(phi); // y

			// determine radius positions of current cell center relative to star center
			
			vec relcenter = tarch::la::abs(center - star);
			// lower left, upper right radius in these coordinates
			l = tarch::la::norm2(relcenter - dx/2.0);
			r = tarch::la::norm2(relcenter + dx/2.0);
			
			Admissible[i] = (l > radius) || (r <= radius);
		}
		
		isAdmissible = (Admissible[0] && Admissible[1]);		
	}

	return isAdmissible;

}

void GridDemonstrator::GridDemonstraor_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 1 + 0
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    Q[0] = 0.0;
  }
}

void GridDemonstrator::GridDemonstraor_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 1 + 0

  // @todo Please implement/augment if required
  stateOut[0] = 0.0;

  fluxOut[0] = 0.0;
}

exahype::solvers::Solver::RefinementControl GridDemonstrator::GridDemonstraor_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void GridDemonstrator::GridDemonstraor_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 1 + 0
  
  // @todo Please implement/augment if required
  lambda[0] = 1.0;
}




//You can either implement this method or modify fusedSource
void GridDemonstrator::GridDemonstraor_ADERDG::algebraicSource(const double* const Q,double* S) {
  // @todo Please implement/augment if required
  S[0] = 0.0;
}


