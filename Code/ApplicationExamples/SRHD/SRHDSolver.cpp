#include "SRHDSolver.h"
#include "GeneratedConstants.h"

#include <memory>

using std::endl;
using std::cout;

extern "C" {
void hastoadjustsolution_(double* time, bool* refine);
void adjustedsolutionvalues_(const double* const x,const double* w,const double* t,const double* dt,double* Q);
void pdeflux_(double* F, const double* const Q);
void pdeeigenvalues_(double* lambda, const double* const Q, const double* nv);
}

void SRHD::SRHDSolver::init(){
  // implement if wanted
}

void SRHD::SRHDSolver::flux(const double* const Q, double** F) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  pdeflux_(F[0], Q);
}



void SRHD::SRHDSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  // normal vector: Allocate for 3 dimensions for convenience
  double nv[3] = {0.};
  nv[normalNonZeroIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}



bool SRHD::SRHDSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {
  return (t < 1e-10);

  // This would be the alternative invocation via Fortran.
  // However this crashes for some reason and is also unneccessary.
  bool refine;
  hastoadjustsolution_(&t, &refine);
  return refine;
}



void SRHD::SRHDSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  adjustedsolutionvalues_(x, &w, &t, &dt, Q);
}

void SRHD::SRHDSolver::source(const double* const Q, double* S){
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
}


exahype::solvers::Solver::RefinementControl SRHD::SRHDSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void SRHD::SRHDSolver::boundaryValues(const double* const x,const double t, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)


  // @todo Please implement
  // fluxOut
  fluxOut[0] = fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = fluxIn[2];
  fluxOut[3] = fluxIn[3];
  fluxOut[4] = fluxIn[4];
  // stateOut
  // @todo Please implement
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];
}


