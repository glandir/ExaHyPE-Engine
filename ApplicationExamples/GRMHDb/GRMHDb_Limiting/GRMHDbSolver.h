#ifndef __GRMHDbSolver_CLASS_HEADER__
#define __GRMHDbSolver_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================

#include <string>

#include "exahype/solvers/LimitingADERDGSolver.h"
#include "GRMHDbSolver_ADERDG.h"
#include "GRMHDbSolver_FV.h"

namespace GRMHDb{
  class GRMHDbSolver;
}

class GRMHDb::GRMHDbSolver: public exahype::solvers::LimitingADERDGSolver {  
  public:
    static constexpr int NumberOfVariables      = GRMHDb::AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
    static constexpr int NumberOfParameters     = GRMHDb::AbstractGRMHDbSolver_ADERDG::NumberOfParameters;
    static constexpr int Order                  = GRMHDb::AbstractGRMHDbSolver_ADERDG::Order;
    static constexpr int NumberOfDMPObservables = GRMHDb::AbstractGRMHDbSolver_ADERDG::NumberOfDMPObservables;
    static constexpr int GhostLayerWidth        = GRMHDb::AbstractGRMHDbSolver_FV::GhostLayerWidth;
  
    GRMHDbSolver(
        const double maximumMeshSize,
        const int maximumMeshDepth,
        const int haloCells,
        const int haloBufferCells,
        const int limiterBufferCells,
        const int regularisedFineGridLevels,
        const exahype::solvers::Solver::TimeStepping timeStepping,
        const int DMPObservables,
        const double DMPRelaxationParameter,
        const double DMPDifferenceScaling
);
    
    
    void projectOnFVLimiterSpace(const double* const luh, double* const lim) const override;
    void projectOnDGSpace(const double* const lim, double* const luh) const override;
    bool discreteMaximumPrincipleAndMinAndMaxSearch(const double* const luh, double* const boundaryMinPerVariables, double* const boundaryMaxPerVariables) override;
    void findCellLocalMinAndMax(const double* const luh, double* const localMinPerVariables, double* const localMaxPerVariable) override;
    void findCellLocalLimiterMinAndMax(const double* const lim, double* const localMinPerObservable, double* const localMaxPerObservable) override;
};

#endif // __GRMHDbSolver_CLASS_HEADER__