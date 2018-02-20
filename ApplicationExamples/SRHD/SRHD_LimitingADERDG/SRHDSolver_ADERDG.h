#ifndef __SRHDSolver_ADERDG_CLASS_HEADER__
#define __SRHDSolver_ADERDG_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "exahype/solvers/ADERDGSolver.h"




namespace SRHD{
  class SRHDSolver_ADERDG;
}

class SRHD::SRHDSolver_ADERDG: public exahype::solvers::ADERDGSolver {
  public:
  	// Sorry for being inconsistent here: While AderDGSolver offers the methods getNumberOfVariables() etc.,
    // in static context they cannot be accessed. Thus the toolkit offers you access to the variables here.
    // Thank you,Toolkit!
    static constexpr int nVar      = 5;
    static constexpr int nParams   = 0;
    static constexpr int nDim      = 2;
    static constexpr int order     = 3;
    static constexpr int basisSize = 3 + 1;
  
    SRHDSolver_ADERDG(double maximumMeshSize,int maximumAdaptiveMeshDepth,int DMPObservables,int limiterHelperLayers,exahype::solvers::Solver::TimeStepping timeStepping,std::vector<std::string>& cmdlineargs);

    void spaceTimePredictor(double* lQhbnd,double* lFhbnd,double** tempSpaceTimeUnknowns,double** tempSpaceTimeFluxUnknowns,double* tempUnknowns,double* tempFluxUnknowns,double* tempStateSizedVectors,const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx,const double dt) override; 
    void solutionUpdate(double* luh,const double* const lduh,const double dt) override;
    void volumeIntegral(double* lduh,const double* const lFhi,const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void surfaceIntegral(double* lduh,const double* const lFhbnd,const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,double* tempFaceUnknownsArray,double** tempStateSizedVectors,double** tempStateSizedSquareMatrices,const double dt,const int normalNonZeroIndex) override;
    void boundaryConditions(double* fluxOut,double* stateOut,const double* const fluxIn,const double* const stateIn,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const int faceIndex,const int normalNonZero) override;
    double stableTimeStepSize(const double* const luh,double* tempEigenvalues,const tarch::la::Vector<DIMENSIONS,double>& dx) override;
    void solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) override;
    bool hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) override;
    exahype::solvers::Solver::RefinementControl refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) override;
    void faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1,int>& subfaceIndex) override;
    void faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1,int>& subfaceIndex) override;
    void volumeUnknownsProlongation(double* luhFine,const double* luhCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) override;
    void volumeUnknownsRestriction(double* luhCoarse,const double* luhFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) override;

    void init(const std::vector<std::string>& cmdlineargs,const exahype::Parser::ParserView& constants) final override;
    void eigenvalues(const double* const Q,const int normalNonZeroIndex,double* lambda);
    void flux(const double* const Q,double** F);
    void source(const double* const Q,double* S);
    void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double * const fluxIn,const double* const stateIn,double *fluxOut,double* stateOut);
    void adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q);
    void ncp(const double* const Q,const double* const gradQ,double* BgradQ);
    void matrixb(const double* const Q,const int normalNonZero,double* Bn);

    bool physicalAdmissibilityDetection(const double* const QMin,const double* const QMax) override;
};

#endif // __SRHDSolver_ADERDG_CLASS_HEADER__
