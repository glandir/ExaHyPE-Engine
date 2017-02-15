#include "GRMHDSolver_ADERDG.h"

#include "GRMHDSolver_ADERDG_Variables.h"
#include "Fortran/InitialData.h"
#include "Fortran/PDE.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

const double excision_radius = 1.0;
  const int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;


void GRMHD::GRMHDSolver_ADERDG::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

bool GRMHD::GRMHDSolver_ADERDG::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) {
   bool insideExcisionBall = std::sqrt(center[0]*center[0] + center[1]*center[1] + center[2]*center[2]) < excision_radius;
  insideExcisionBall = false;
  return tarch::la::equals(t,0.0) || insideExcisionBall;
}

void GRMHD::GRMHDSolver_ADERDG::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Fortran
  initialdata_(x, &t, Q);
}

void GRMHD::GRMHDSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void GRMHD::GRMHDSolver_ADERDG::flux(const double* const Q,double** F) {
  pdeflux_(F[0], Q);
}


void GRMHD::GRMHDSolver_ADERDG::source(const double* const Q,double* S) {
  pdesource_(S, Q);
}


void GRMHD::GRMHDSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int d,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {

  const int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
  const int order = GRMHD::AbstractGRMHDSolver_ADERDG::Order;
  const int basisSize = order + 1;
  const int nDim = DIMENSIONS;

  double Qgp[nVar];
  std::memset(stateOut, 0, nVar * sizeof(double));
  std::memset(fluxOut, 0, nVar * sizeof(double));

  double F[3 * nVar]; // Fortran needs continous storage!
                      // Use always 3 dimensions here since the kernels works with those internally; see nDim in PDE.f90;

  kernels::idx2 F_idx(nDim, nVar);

  for(int i=0; i < basisSize; i++)  { // i == time
     const double weight = kernels::gaussLegendreWeights[order][i];
     const double xi = kernels::gaussLegendreNodes[order][i];
     double ti = t + xi * dt;

     initialdata_(x, &ti, Qgp);
     pdeflux_(F, Qgp);
     for(int m=0; m < nVar; m++) {
        stateOut[m] += weight * Qgp[m];
        fluxOut[m] += weight * F[F_idx(d, m)];
     }
  }
}


exahype::solvers::Solver::RefinementControl GRMHD::GRMHDSolver_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


bool GRMHD::GRMHDSolver_ADERDG::physicalAdmissibilityDetection(const double* const QMin,const double* const QMax) {
  if (QMin[0] < 0.0) return false;
  if (QMin[4] < 0.0) return false;

  for (int i=0; i<nVar; ++i) {
    if (!std::isfinite(QMin[i])) return false;
    if (!std::isfinite(QMax[i])) return false;
  }

  return true;
}


void GRMHD::GRMHDSolver_ADERDG::ncp(const double* const Q,const double* const gradQ,double* BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
}


void GRMHD::GRMHDSolver_ADERDG::matrixb(const double* const Q,const int d,double* Bn) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);
}

bool GRMHD::GRMHDSolver_ADERDG::isDummyKRequired() const { 
  return false;
}


void GRMHD::GRMHDSolver_ADERDG::dummyK_Value(const double* const x,const double t,const double dt, double* forceVector, double* x0) {
  // ...
}