    // This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
// 
// ========================
//   www.exahype.eu
// ========================

#include "GPRDRSolver_ADERDG.h"
#include "GPRDRSolver_FV.h"

#include <algorithm>

#include "GPRDRSolver_ADERDG_Variables.h"
#include "kernels/GaussLegendreQuadrature.h"

#include "PDE.h"
#include "InitialData.h"
#include "Tools.h"

#include "kernels/KernelUtils.h"
#include "peano/utils/Loop.h"

#include "GPRDRSolver_ADERDG_Variables.h"

#include "kernels/KernelUtils.h"
#include "peano/utils/Loop.h"

#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/multicore/Lock.h"

tarch::logging::Log GPRDR::GPRDRSolver_ADERDG::_log( "GPRDR::GPRDRSolver_ADERDG" );


void GPRDR::GPRDRSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  const int order = GPRDR::GPRDRSolver_ADERDG::Order;
  	constexpr int basisSize = AbstractGPRDRSolver_FV::PatchSize;
	constexpr int Ghostlayers = AbstractGPRDRSolver_FV::GhostLayerWidth;
	
	static tarch::multicore::BooleanSemaphore initializationSemaphoreDG;
  
  
    tarch::multicore::Lock lock(initializationSemaphoreDG);	
	
  //inittecplot_(&order,&order);
  if (constants.isValueValidString("reference")) {
    std::string reference = constants.getValueAsString("reference");
	const int length=reference.length();
	logInfo("init(...)","Reference setup:"<<reference);
  	printf("\n******************************************************************");
	printf("\n**************<<<  INIT TECPLOT    >>>****************************");
	printf("\n******************************************************************");
    inittecplot_(&order,&order,&basisSize,&Ghostlayers);
	//inittecplot_(&order,&order);
	printf("\n******************************************************************");
	printf("\n**************<<<  INIT PDE SETUP  >>>****************************");
	printf("\n******************************************************************\n");
    initparameters_(&length,&reference[0]);
	printf("\n******************************************************************");
	printf("\n**************<<<       DONE       >>>****************************");
	printf("\n******************************************************************");
  } else {
    logInfo("init(...)","Not recognized setup.");
	std::abort();
  }	
	
	lock.free();
	
	
}

void GPRDR::GPRDRSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  if (tarch::la::equals(t,0.0)) {
    int md = exahype::solvers::Solver::getMaximumAdaptiveMeshDepth();
    double cms = exahype::solvers::Solver::getCoarsestMeshSize();
    const int order = GPRDR::GPRDRSolver_ADERDG::Order;
    std::fill_n(Q,24,0.0);
    
    //    initialdata_(x, &ti, Qgp,&md,&cms,&order);
    double x_3[3];
    x_3[2]=0;
    std::copy_n(&x[0],DIMENSIONS,&x_3[0]);
    
    initialdata_(x_3, &t, Q);
  }
  for(int i = 0; i< 24 ; i++){
    assert(std::isfinite(Q[i]));
  }
}

void GPRDR::GPRDRSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int direction,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  const int nVar = GPRDR::GPRDRSolver_ADERDG::NumberOfVariables;
  const int order = GPRDR::GPRDRSolver_ADERDG::Order;
  const int basisSize = order + 1;
  const int nDim = DIMENSIONS;
  double Qgp[nVar],*F[nDim], Fs[nDim][nVar];

  double x_3[3];
  x_3[2]=0;
  std::copy_n(&x[0],DIMENSIONS,&x_3[0]);
  
  int md=0;
  double cms=0;
	
  std::memset(stateOut, 0, nVar * sizeof(double));
  std::memset(fluxOut , 0, nVar * sizeof(double));
	
  for(int dd=0; dd<nDim; dd++) F[dd] = Fs[dd];

  for(int i=0; i < basisSize; i++)  { // i == time
    const double weight = kernels::gaussLegendreWeights[order][i];
    const double xi = kernels::gaussLegendreNodes[order][i];
    double ti = t + xi * dt;

    //    initialdata_(x, &ti, Qgp,&md,&cms,&order);
    initialdata_(x_3, &ti, Qgp);
    flux(Qgp, F);
    for(int m=0; m < nVar; m++) {
      stateOut[m] += weight * Qgp[m];
      fluxOut[m] += weight * Fs[direction][m];
    }
  }
  std::copy_n(stateIn,nVar,stateOut);
  std::copy_n(fluxIn,nVar,fluxOut);
}

exahype::solvers::Solver::RefinementControl GPRDR::GPRDRSolver_ADERDG::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,double t,const int level) {
  if (tarch::la::equals(t,0.0)) {
  if(DIMENSIONS == 2){
    if(std::abs(cellCentre[0]) < 1000){
      if(std::abs(cellCentre[1]) < 50){
	return exahype::solvers::Solver::RefinementControl::Refine;
      }
    }
  }else{
    if(std::abs(cellCentre[0]) < 50){
      if(std::abs(cellCentre[1]) < 50){
		  if(std::abs(cellCentre[2]) < 50){
			return exahype::solvers::Solver::RefinementControl::Refine;
		  }
      }
    }	  
	  
  };
  }	 
  
  //return exahype::solvers::Solver::RefinementControl::Keep;
  if ( level > getCoarsestMeshLevel() ) {
    return exahype::solvers::Solver::RefinementControl::Erase;
  }
return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void GPRDR::GPRDRSolver_ADERDG::eigenvalues(const double* const Q,const int direction,double* const lambda) {
  double nv[3] = {0.};
  nv[direction] = 1;
  pdeeigenvalues_(lambda, Q, nv);

  for(int i = 0; i< 24 ; i++){
    assert(std::isfinite(lambda[i]));
  }
}





void GPRDR::GPRDRSolver_ADERDG::flux(const double* const Q,double** const F) {
  constexpr int nVar = GPRDR::GPRDRSolver_ADERDG::NumberOfVariables;
  if(DIMENSIONS == 2){
    double F_3[nVar];
    pdeflux_(F[0], F[1],F_3, Q);
  }else{
    pdeflux_(F[0], F[1],F[2], Q);
  }

  for(int d = 0; d< DIMENSIONS ; d++){
    for(int i = 0; i< 24 ; i++){
      assert(std::isfinite(F[d][i]));
    }
  }
}


//You can either implement this method or modify fusedSource
void GPRDR::GPRDRSolver_ADERDG::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  pdesource_(S, Q);

  for(int i = 0; i< 24 ; i++){
    assert(std::isfinite(S[i]));
  }
}

void GPRDR::GPRDRSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  //std::fill_n(BgradQ,24,0.0);
  pdencp_(BgradQ, Q, gradQ);
  for(int i = 0; i< 24 ; i++){
    assert(std::isfinite(BgradQ[i]));
  }
}

bool GPRDR::GPRDRSolver_ADERDG::isPhysicallyAdmissible(
      const double* const                         solution,
      const double* const                         localDMPObservablesMin,
      const double* const                         localDMPObservablesMax,
      const bool                                  wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>& cellSize,
      const double                                timeStamp) const
{
  		  
  int limvalue;
   
  pdelimitervalue_(&limvalue,&cellCentre[0],&NumberOfDMPObservables, localDMPObservablesMin, localDMPObservablesMax);
  bool ret_value;
  limvalue > 0 ? ret_value=false : ret_value=true;
  return ret_value;
}


void GPRDR::GPRDRSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* const observables,
								       const double* const Q) const {
  VariableShortcuts s;
  observables[0] = Q[17]; //alpha
  observables[1] = Q[20]; //xi
  
  double criticalStress;

  pdecritialstress_(&criticalStress,Q);
  observables[2] = criticalStress; //criticalStress
}




