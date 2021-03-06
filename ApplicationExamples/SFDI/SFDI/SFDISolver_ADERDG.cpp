// This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
// 
// ========================
//   www.exahype.eu
// ========================

#include "SFDISolver_ADERDG.h"
#include "SFDISolver_FV.h"

#include <algorithm>
#include "kernels/GaussLegendreBasis.h"
#include "PDE.h"
#include "InitialData.h"
#include "Tools.h"

#include "SFDISolver_ADERDG_Variables.h"

#include "SFDISolver_FV.h"

#include "kernels/KernelUtils.h"
#include "peano/utils/Loop.h"

#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/multicore/Lock.h"

#include <cstring>

tarch::logging::Log SFDI::SFDISolver_ADERDG::_log( "SFDI::SFDISolver_ADERDG" );


void SFDI::SFDISolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  constexpr int order = SFDI::AbstractSFDISolver_ADERDG::Order;
  constexpr int basisSize = AbstractSFDISolver_FV::PatchSize;
  constexpr int Ghostlayers = AbstractSFDISolver_FV::GhostLayerWidth;
  
  static tarch::multicore::BooleanSemaphore initializationSemaphoreDG;
  tarch::multicore::Lock lock(initializationSemaphoreDG);

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
    printf("\n******************************************************************");
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

void SFDI::SFDISolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  if (tarch::la::equals(t,0.0)) {
    int md = exahype::solvers::Solver::getMaximumAdaptiveMeshDepth();
    double cms = exahype::solvers::Solver::getCoarsestMeshSize();
    const int order = Order;
    std::fill_n(Q,NumberOfVariables,0.0);
    
    //    initialdata_(x, &ti, Qgp,&md,&cms,&order);
    double x_3[3];
    x_3[2]=0;
    std::copy_n(&x[0],DIMENSIONS,&x_3[0]);
    
    initialdata_(x_3, &t, Q);
  }
  for(int i = 0; i< NumberOfVariables ; i++){
    assertion2(std::isfinite(Q[i]),i,Q[i]);
  }
}

void SFDI::SFDISolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int direction,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  const int nVar = NumberOfVariables;
  const int order = Order;
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
    const double weight = kernels::legendre::weights[order][i];
    const double xi = kernels::legendre::nodes[order][i];
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

exahype::solvers::Solver::RefinementControl SFDI::SFDISolver_ADERDG::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,double t,const int level) {
//  return exahype::solvers::Solver::RefinementControl::Keep;
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


void SFDI::SFDISolver_ADERDG::eigenvalues(const double* const Q,const int direction,double* const lambda) {
  double nv[3] = {0.};
  nv[direction] = 1;
  pdeeigenvalues_(lambda, Q, nv);
  
  for(int i = 0; i< NumberOfVariables ; i++){
    assertion2(std::isfinite(lambda[i]),i,lambda[i]);
  }
}





void SFDI::SFDISolver_ADERDG::flux(const double* const Q,double** const F) {
  constexpr int nVar = NumberOfVariables;
  if(DIMENSIONS == 2){
    double F_3[nVar];
    pdeflux_(F[0], F[1],F_3, Q);
  }else{
    pdeflux_(F[0], F[1],F[2], Q);
  }
  
  for(int d = 0; d< DIMENSIONS ; d++){
    for(int i = 0; i< NumberOfVariables ; i++){
      assertion3(std::isfinite(F[d][i]),d,i,F[d][i]);
    }
  }
}


void SFDI::SFDISolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* const observables,
                       const double* const Q) const {
  observables[0] = Q[0];   
  observables[1] = Q[1];
  observables[2] = Q[2];
  observables[3] = Q[3];
  observables[4] = Q[4];
  observables[5] = Q[5];
}

bool SFDI::SFDISolver_ADERDG::isPhysicallyAdmissible(
      const double* const                         solution,
      const double* const                         localDMPObservablesMin,
      const double* const                         localDMPObservablesMax,
      const bool                                  wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>& cellSize,
      const double                                timeStamp) const
{
  int limvalue;
   
  const int nObs = NumberOfDMPObservables;
  pdelimitervalue_(&limvalue,&cellCentre[0],&nObs, localDMPObservablesMin, localDMPObservablesMax);
  bool ret_value;
  limvalue > 0 ? ret_value=false : ret_value=true;
  return ret_value;
}

//You can either implement this method or modify fusedSource
void SFDI::SFDISolver_ADERDG::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  
  pdesource_(S, Q);
  for(int i = 0; i< NumberOfVariables ; i++){
    assertion2(std::isfinite(S[i]),i,S[i]);
  }
  
}

void  SFDI::SFDISolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
  for(int i = 0; i< NumberOfVariables ; i++){
    assertion2(std::isfinite(BgradQ[i]),i,BgradQ[i]);
  }
}


