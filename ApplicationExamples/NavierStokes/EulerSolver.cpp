// This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
// 
// ========================
//   www.exahype.eu
// ========================

#include <cmath>

#include "EulerSolver.h"
#include "stableDiffusiveTimeStepSize.h"
#include "diffusiveRiemannSolver.h"

#include "kernels/aderdg/generic/Kernels.h"
#include "EulerSolver_Variables.h"
#include "kernels/KernelUtils.h"


tarch::logging::Log Euler::EulerSolver::_log( "Euler::EulerSolver" );


void Euler::EulerSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants)
{
  // TODO: Fix referenceT, referenceViscosity and sutherlandC
  const auto referenceT = 23.;
  const auto referenceViscosity = 0.0;
  const auto sutherlandC = 0.0;
  ns = NavierStokes::NavierStokes(referenceT, referenceViscosity, sutherlandC);
  // @todo Please implement/augment if required
}

void Euler::EulerSolver::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 5 + 0
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    double p = 1.0;
    Variables vars(Q);
    // Sod shock tube
    if (x[0] < 0.5) {
      vars.rho() = 1.0;
      const double u = 0.0; 
      vars.j(u * vars.rho(), 0.0, 0.0);
      p  = 1.0;
    } else {
      vars.rho() = 0.125;
      const double u = 0;
      vars.j(u * vars.rho(), 0.0, 0.0);
      p  = 0.1;
    }

    // Initial value for E given by equation of state!
    const auto scale = 1.0/(vars.rho() * (ns.GAMMA - 1));
    const auto norm = vars.j(0) * vars.j(0) + vars.j(1) * vars.j(1) + vars.j(2) * vars.j(2);
    vars.E() = scale * (vars.rho() * p + 0.5 * (ns.GAMMA - 1) * norm);

  }

}

void Euler::EulerSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
					const double * const fluxIn,const double* const stateIn, const double* const gradStateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 5 + 0
  // @todo Please implement/augment if required

  // Set wall boundary conditions.

  // First set the state.
  std::copy_n(stateIn, NumberOfVariables, stateOut);
  stateOut[1+normalNonZero]= -stateOut[1+normalNonZero];

  // Compute gradient of state.
  constexpr auto gradSize = NumberOfVariables * DIMENSIONS;
  auto gradStateOut = std::array<double, gradSize>{{0.0}};

  const auto mirroredIdx = 1 + normalNonZero;

  kernels::idx2 idx_gradQ(DIMENSIONS, NumberOfVariables);
  for (int i = 0; i < DIMENSIONS; ++i) {
    for (int j = 0; j < NumberOfVariables; ++j) {
      const auto idx = idx_gradQ(i, j);
      // Mirror gradient where vel. is mirrored.
      // Mirror for all dims or only dim == normalNonZero
      if (j == 1 + normalNonZero) {
	gradStateOut[idx] = gradStateIn[idx];
      } else if(j == 1 || j == 2 || j == 3) {
       	//gradStateOut[idx] = -gradStateIn[idx];
	assert(idx < NumberOfVariables * DIMENSIONS);
	gradStateOut[idx] = -gradStateIn[idx];
      } else {
	gradStateOut[idx] = gradStateIn[idx];
      }
    }
  }


  // Compute appropiate flux.
  double _F[3][NumberOfVariables]={0.0};
  double* F[3] = {_F[0], _F[1], _F[2]};
  // TODO(Lukas): Compute gradient for gradStateOut
  flux(stateOut, gradStateOut.data(), F);

  std::copy_n(F[normalNonZero], NumberOfVariables, fluxOut);
 }

exahype::solvers::Solver::RefinementControl Euler::EulerSolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void Euler::EulerSolver::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 5 + 0
  ns.evaluateEigenvalues(Q, d, lambda);
}

void Euler::EulerSolver::diffusiveEigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 5 + 0
  ns.evaluateDiffusiveEigenvalues(Q, d, lambda);
}

// TODO(Lukas) remove, currently called in boundaryValues!
void Euler::EulerSolver::flux(const double* const Q, double** F) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 5 + 0
  ReadOnlyVariables vars(Q);
  Fluxes f(F);
  assert(false);
}

void Euler::EulerSolver::flux(const double* const Q,const double* const gradQ, double** F) {
  ns.evaluateFlux(Q, gradQ, F);
}

double Euler::EulerSolver::stableTimeStepSize(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx) {
  //return std::min(1., stableDiffusiveTimeStepSize<EulerSolver>(*static_cast<EulerSolver*>(this),luh,dx));
  return kernels::aderdg::generic::c::stableTimeStepSize<EulerSolver>(*static_cast<EulerSolver*>(this),luh,dx);
}

void Euler::EulerSolver::riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,const double dt,const int direction, bool isBoundaryFace, int faceIndex) {
  assertion2(direction>=0,dt,direction);
  assertion2(direction<DIMENSIONS,dt,direction);
  const auto i = 27;
  const auto invDx = tarch::la::Vector<DIMENSIONS, double>({i, i, i});
  //riemannSolverNonlinear<false,EulerSolver>(*static_cast<EulerSolver*>(this),FL,FR,QL,QR,i, dt,direction);
  kernels::aderdg::generic::c::riemannSolverNonlinear<false,EulerSolver>(*static_cast<EulerSolver*>(this),FL,FR,QL,QR,dt,direction);

}
