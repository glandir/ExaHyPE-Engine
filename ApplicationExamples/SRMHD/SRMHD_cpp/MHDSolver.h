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



namespace MHDSolver{
  class MHDSolver;
}


class MHDSolver::MHDSolver: public exahype::solvers::ADERDGSolver {
  public:
    MHDSolver(double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler, exahype::parser::ParserView constants);

    void spaceTimePredictor(double* const lQhbnd,double* const lFhbnd, double** tempSpaceTimeUnknowns, double** tempSpaceTimeFluxUnknowns, double* const  tempUnknowns, double* const  tempFluxUnknowns, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt) override; 
    void solutionUpdate(double* const luh, const double* const lduh, const double dt) override;
    void volumeIntegral(double* const lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void surfaceIntegral(double* const lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void riemannSolver(double* const FL, double* const FR, const double* const QL, const double* const QR, double* tempFaceUnknownsArray, double** tempStateSizedVectors, double** tempStateSizedSquareMatrices, const double dt, const int normalNonZeroIndex) override;
    void boundaryConditions(double* const fluxOut, double* stateOut, const double* const fluxIn, const double* const stateIn, const tarch::la::Vector<DIMENSIONS, double>& cellCentre, const tarch::la::Vector<DIMENSIONS,double>& cellSize, const double t,const double dt, const int faceIndex, const int normalNonZero) override;
    double stableTimeStepSize(const double* const luh, double* tempEigenvalues, const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) override;
    bool hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t) override;
    exahype::solvers::Solver::RefinementControl refinementCriterion(const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) override;
    void faceUnknownsProlongation(double* const lQhbndFine,double* const lFhbndFine,const double* const lQhbndCoarse,const double* const lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) override;
    void faceUnknownsRestriction(double* const lQhbndCoarse,double* const lFhbndCoarse,const double* const lQhbndFine,const double* const lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) override;
    void volumeUnknownsProlongation(double* const luhFine, const double* const luhCoarse, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) override;
    void volumeUnknownsRestriction(double* const luhCoarse, const double* const luhFine, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) override;

    static void Prim2Cons(const double* V,double* Q);
    static void Cons2Prim(const double* Q,double* V);
    static void AlfvenWave(const double* x, double* Q, const double t);
  private:
    void init();
    static void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);
    static void flux(const double* const Q, double** F);
    static void source(const double* const Q, double* S);
    static void boundaryValues(const double* const x,const double t, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut);
    static void adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q);
};


