// This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
// 
// ========================
//   www.exahype.eu
// ========================

#include "MHDSolver.h"

#include "MHDSolver_Variables.h"

#include "InitialDataAdapter.h"
#include "PDE.h"

#include <memory>
#include <cstring>
#include <stdio.h>
#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h" // matrix indexing

tarch::logging::Log MHD::MHDSolver::_log( "MHD::MHDSolver" );


void MHD::MHDSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
  // just pass the pointer to the crazy Fortran glue code. Should be improved.
  // constants = &_constants;
}

void MHD::MHDSolver::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 9 + 0
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    adjustedsolutionvalues_(x, &t, &dt, Q);
  }
}

void MHD::MHDSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 9 + 0

  const int nVar = MHD::AbstractMHDSolver::NumberOfVariables;
  const int order = MHD::AbstractMHDSolver::Order;
  const int basisSize = order + 1;
  const int nDim = DIMENSIONS;

  double Qgp[nVar];
  std::memset(stateOut, 0, nVar * sizeof(double));
  std::memset(fluxOut, 0, nVar * sizeof(double));

  double F[nDim][nVar];

  // Integrate solution in gauss points (Qgp) in time
  for(int i=0; i < basisSize; i++)  { // i == time
     const double weight = kernels::gaussLegendreWeights[order][i];
     const double xi = kernels::gaussLegendreNodes[order][i];
     double ti = t + xi * dt;

     alfenwave_(x, Qgp, &ti);
     pdeflux_(F[0], F[1], (nDim==3) ? F[2] : nullptr, Qgp);
     for(int m=0; m < nVar; m++) {
  //if(m==checkm) printf("fluxOut[%d] += %.20e\n", m, weight * F[normalNonZero][m]);
        stateOut[m] += weight * Qgp[m];
        fluxOut[m] += weight * F[normalNonZero][m];
     }
  }
}

exahype::solvers::Solver::RefinementControl MHD::MHDSolver::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void MHD::MHDSolver::eigenvalues(const double* const Q,const int d,double* const lambda) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 9 + 0
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}





void MHD::MHDSolver::flux(const double* const Q,double** const F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 9 + 0
  // Caveats: Fortran accepts a uniform array of size (nVar*nDim), however C passes an array of pointers.
  // This Fortran interface works only if F is a continous array and F[1]==F[nDim+1] etc!
  
  // Allow non-continous storage:
  const int nVar = MHD::AbstractMHDSolver::NumberOfVariables;
  const int nDim = DIMENSIONS;
  pdeflux_(F[0], F[1], (nDim==3) ? F[2] : nullptr, Q);
}


