// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================


#include <memory>

#include "exahype/parser/ParserView.h"
#include "exahype/profilers/Profiler.h"
#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"



namespace Euler{
  class SecondEulerSolver;
}


class Euler::SecondEulerSolver: public exahype::solvers::ADERDGSolver {
  public:
    SecondEulerSolver(double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler);

    void spaceTimePredictor(double* const lQhbnd,double* const lFhbnd, double** const tempSpaceTimeUnknowns, double** const tempSpaceTimeFluxUnknowns, double* const  tempUnknowns, double* const  tempFluxUnknowns, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt) override; 
    void solutionUpdate(double* const luh, const double* const lduh, const double dt) override;
    void volumeIntegral(double* const lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void surfaceIntegral(double* const lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void riemannSolver(double* const FL, double* const FR, const double* const QL, const double* const QR, double* const tempFaceUnknownsArray, double** const tempStateSizedVectors, double** const tempStateSizedSquareMatrices, const double dt, const int normalNonZeroIndex) override;
    void boundaryConditions(double* const fluxOut, double* const stateOut, const double* const fluxIn, const double* const stateIn, const tarch::la::Vector<DIMENSIONS, double>& cellCentre, const tarch::la::Vector<DIMENSIONS,double>& cellSize, const double t,const double dt, const int faceIndex, const int normalNonZero) override;
    double stableTimeStepSize(const double* const luh, double* const tempEigenvalues, const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void solutionAdjustment(double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) override;
    bool hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t) override;
    exahype::solvers::Solver::RefinementControl refinementCriterion(const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) override;
    void faceUnknownsProlongation(double* const lQhbndFine,double* const lFhbndFine,const double* const lQhbndCoarse,const double* const lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) override;
    void faceUnknownsRestriction(double* const lQhbndCoarse,double* const lFhbndCoarse,const double* const lQhbndFine,const double* const lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) override;
    void volumeUnknownsProlongation(double* const luhFine, const double* const luhCoarse, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) override;
    void volumeUnknownsRestriction(double* const luhCoarse, const double* const luhFine, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) override;
  private:
    void init();
    static void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* const lambda);
    static void flux(const double* const Q, double** const F);
    static void source(const double* const Q, double* const S);
    static void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut);
    static void adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* const Q);
};


