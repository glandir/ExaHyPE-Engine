#ifndef __SRHDSolver_FV_CLASS_HEADER__
#define __SRHDSolver_FV_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "exahype/parser/ParserView.h"
#include "exahype/solvers/FiniteVolumesSolver.h"


namespace SRHD{
  class SRHDSolver_FV;
}

class SRHD::SRHDSolver_FV : public exahype::solvers::FiniteVolumesSolver {
  public:
    SRHDSolver_FV(int cellsPerCoordinateAxis,double maximumMeshSize,exahype::solvers::Solver::TimeStepping timeStepping);
    
    double stableTimeStepSize(const double* const luh,double* tempEigenvalues,const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void solutionAdjustment(double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) override;
    bool hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) override;
    exahype::solvers::Solver::RefinementControl refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) override;
    void solutionUpdate(double* luhNew,const double* luh,double** tempStateSizedArrays,double** tempUnknowns,const tarch::la::Vector<DIMENSIONS,double>& dx,const double dt,double& maxAdmissibleDt) override;
    void ghostLayerFilling(double* luh,const double* luhNeighbour,const tarch::la::Vector<DIMENSIONS,int>& neighbourPosition) override;
    void ghostLayerFillingAtBoundary(double* luh,const double* luhbnd,const tarch::la::Vector<DIMENSIONS,int>& boundaryPosition) override;
    void boundaryLayerExtraction(double* luhbnd,const double* luh,const tarch::la::Vector<DIMENSIONS,int>& boundaryPosition) override;
    void boundaryConditions(double* stateOut,const double* const stateIn,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const int faceIndex,const int normalNonZero) override;
	
    void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const stateIn,double* stateOut);
  private:
  	void init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) final override;
    static void adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q);
    static void eigenvalues(const double* const Q,const int normalNonZeroIndex,double* lambda);
    static void flux(const double* const Q,double** F);
    static void source(const double* const Q,double* S);
};


#endif // __SRHDSolver_FV_CLASS_HEADER__