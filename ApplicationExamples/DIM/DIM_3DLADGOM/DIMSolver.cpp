// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
// ==============================================
// Please do not change the implementations below
// ==============================================
#include "DIMSolver.h"

#include "kernels/limiter/generic/Limiter.h"

DIM::DIMSolver::DIMSolver(
        const double maximumMeshSize,
        const int maximumMeshDepth,
        const int haloCells,
        const int haloBufferCells,
        const int limiterBufferCells,
        const int regularisedFineGridLevels,
        const exahype::solvers::Solver::TimeStepping timeStepping,
        const int DMPObservables,
        const double DMPRelaxationParameter,
        const double DMPDifferenceScaling,
        const int iterationsToCureTroubledCell 
) :
  exahype::solvers::LimitingADERDGSolver::LimitingADERDGSolver(
      "DIMSolver",
    new DIM::DIMSolver_ADERDG(
      maximumMeshSize,maximumMeshDepth,haloCells,haloBufferCells,limiterBufferCells,regularisedFineGridLevels,timeStepping,DMPObservables),
    new DIM::DIMSolver_FV(
      maximumMeshSize, timeStepping),
    DMPRelaxationParameter,
    DMPDifferenceScaling,
    iterationsToCureTroubledCell) {}

void DIM::DIMSolver::projectOnFVLimiterSpace(const double* const luh, double* const lim) const {
  kernels::limiter::generic::c::projectOnFVLimiterSpace<Order+1,NumberOfVariables+NumberOfParameters,GhostLayerWidth>(luh, lim);
}

void DIM::DIMSolver::projectOnDGSpace(const double* const lim, double* const luh) const {
  kernels::limiter::generic::c::projectOnDGSpace<Order+1,NumberOfVariables+NumberOfParameters,GhostLayerWidth>(lim, luh);
}

bool DIM::DIMSolver::discreteMaximumPrincipleAndMinAndMaxSearch(const double* const luh, double* const boundaryMinPerVariables, double* const boundaryMaxPerVariables) {
  return kernels::limiter::generic::c::discreteMaximumPrincipleAndMinAndMaxSearch<AbstractDIMSolver_ADERDG, NumberOfDMPObservables, GhostLayerWidth>(luh, *static_cast<AbstractDIMSolver_ADERDG*>(_solver.get()), _DMPMaximumRelaxationParameter, _DMPDifferenceScaling, boundaryMinPerVariables, boundaryMaxPerVariables);
}

void DIM::DIMSolver::findCellLocalMinAndMax(const double* const luh, double* const localMinPerVariables, double* const localMaxPerVariable) {
  kernels::limiter::generic::c::findCellLocalMinAndMax<AbstractDIMSolver_ADERDG, NumberOfDMPObservables>(luh, *static_cast<AbstractDIMSolver_ADERDG*>(_solver.get()), localMinPerVariables, localMaxPerVariable);
}
void DIM::DIMSolver::findCellLocalLimiterMinAndMax(const double* const lim, double* const localMinPerObservable, double* const localMaxPerObservable) {
  kernels::limiter::generic::c::findCellLocalLimiterMinAndMax<AbstractDIMSolver_ADERDG, NumberOfDMPObservables, GhostLayerWidth>(lim, *static_cast<AbstractDIMSolver_ADERDG*>(_solver.get()), localMinPerObservable,localMaxPerObservable);
}