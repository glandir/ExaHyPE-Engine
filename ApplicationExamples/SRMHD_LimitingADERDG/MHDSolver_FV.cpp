#include "MHDSolver_FV.h"

#include "PDE.h"
#include "InitialData.h"
#include "BoundaryConditions.h"

#include <memory>

void MHD::MHDSolver_FV::init(std::vector<std::string>& cmdlineargs) {
  // All the work is done in MHDSolver_ADERDG::init.
}

bool MHD::MHDSolver_FV::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) {
  return (t < 1e-10);
}

void MHD::MHDSolver_FV::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  if (tarch::la::equals(t, 0.0)) {
    idfunc(x, Q);
  }
}

exahype::solvers::Solver::RefinementControl MHD::MHDSolver_FV::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void MHD::MHDSolver_FV::eigenvalues(const double* const Q,const int normalNonZeroIndex,double* lambda) {
  double nv[3] = {0.};
  nv[normalNonZeroIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void MHD::MHDSolver_FV::flux(const double* const Q,double** F) {
  // Caveats: Fortran accepts a uniform array of size (nVar*nDim), however C passes an array of pointers.
  // This Fortran interface works only if F is a continous array and F[1]==F[nDim+1] etc!
  pdeflux_(F[0], Q);
}


void MHD::MHDSolver_FV::algebraicSource(const double* const Q,double* S) {
  pdesource_(S, Q);
}


void MHD::MHDSolver_FV::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
const double* const stateIn,double* stateOut) {
  // we abuse the BC func for the ADERDG here, thus we reserve arbitrary
  // storage for the fluxes.
  double fluxIn[9], fluxOut[9];
  double nv[3] = {0.};
  nv[normalNonZero] = 1;
  bcfunc(x, &t, &dt, &faceIndex, nv, fluxIn, stateIn, fluxOut, stateOut);

//  const double* normalVelocity = stateIn+1+normalNonZero;
//  const double sign = (faceIndex-2*normalNonZero)==0 ? -1.0 : 1.0;
//  const double outwardDirectedVelocity = sign * (*normalVelocity);
//
//  std::cout << "normalVelocity="<<normalVelocity<<",faceIndex="<<faceIndex<<",outwardDirectedVelocity="<<outwardDirectedVelocity << std::endl;
//
//  if (outwardDirectedVelocity>0) { // outflow; take inside value
//    // do nothing
//  } else { // inflow; take outside values
//    stateOut[0] = 1.0e-4; // density for InitialBlast
//    stateOut[4] = 5.0e-4; // pressure for InitialBlast
//  }
}


//void MHD::MHDSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
//  constexpr int nVar = 9;
//  std::memset(BgradQ, 0, nVar * sizeof(double));
//}
//
//
//void MHD::MHDSolver_FV::coefficientMatrix(const double* const Q,const int normalNonZero,double* Bn) {
//  constexpr int nVar = 9;
//  std::memset(Bn, 0, nVar * nVar * sizeof(double));
//}
