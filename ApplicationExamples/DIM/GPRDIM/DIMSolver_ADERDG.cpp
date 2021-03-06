// This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
// 
// ========================
//   www.exahype.eu
// ========================

#include "DIMSolver_ADERDG.h"

#include "DIMSolver_ADERDG_Variables.h"

// User defined calls
#include "PDE.h"
#include "InitialData.h"
#include "C2P-DIM.h"
#include "TECPLOTinterface.h"
// Used for the rieman-solver part
#include "peano/utils/Loop.h"
#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreBasis.h"
tarch::logging::Log GPRDIM::DIMSolver_ADERDG::_log( "GPRDIM::DIMSolver_ADERDG" );


void GPRDIM::DIMSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
	const int order = GPRDIM::AbstractDIMSolver_ADERDG::Order;
	inittecplot_(&order,&order);
}

void GPRDIM::DIMSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
	if (tarch::la::equals(t,0.0)) {
		int md = exahype::solvers::Solver::getMaximumAdaptiveMeshDepth();
		double cms = exahype::solvers::Solver::getCoarsestMeshSize();
		const int order = GPRDIM::AbstractDIMSolver_ADERDG::Order;
		initialdata_(x, &t, Q,&md,&cms,&order);
	}
	// Call the dynamic rupture subroutine

	dynamicrupture_(x, &t, Q);
	/*if (t>0.005) {
		//printf("Crack!");
		if(x[0]>0. && x[0]<0.1) {
		Q[17]=0.0;
		}
	}*/
}

void GPRDIM::DIMSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut)
	// Local variables
	const int nVar = GPRDIM::AbstractDIMSolver_ADERDG::NumberOfVariables;
	const int order = GPRDIM::AbstractDIMSolver_ADERDG::Order;
	const int basisSize = order + 1;
	const int nDim = DIMENSIONS;
	double Qgp[nVar],*F[nDim], Fs[nDim][nVar];
	
	int md=0;
	double cms=0;
	
	std::memset(stateOut, 0, nVar * sizeof(double));
	std::memset(fluxOut , 0, nVar * sizeof(double));
	
	for(int dd=0; dd<nDim; dd++) F[dd] = Fs[dd];

	for(int i=0; i < basisSize; i++)  { // i == time
		const double weight = kernels::legendre::weights[order][i];
		const double xi = kernels::legendre::nodes[order][i];
		double ti = t + xi * dt;

		initialdata_(x, &ti, Qgp,&md,&cms,&order);
		flux(Qgp, F);
		for(int m=0; m < nVar; m++) {
			stateOut[m] += weight * Qgp[m];
			fluxOut[m] += weight * Fs[normalNonZero][m];
		}
	}
}

exahype::solvers::Solver::RefinementControl GPRDIM::DIMSolver_ADERDG::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
	const int order = GPRDIM::AbstractDIMSolver_ADERDG::Order;
	int rupture_flag;
	

	ruptureflag_(&rupture_flag,&order,luh,&center[0],&dx[0]);
	if(rupture_flag>0){
		return exahype::solvers::Solver::RefinementControl::Refine;	
	}else{
		return exahype::solvers::Solver::RefinementControl::Keep;
	}
	
	//if ( level > getCoarsestMeshLevel() )
	//	return exahype::solvers::Solver::RefinementControl::Erase;
	//	else return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void GPRDIM::DIMSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* const lambda) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 24 + 0
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void GPRDIM::DIMSolver_ADERDG::flux(const double* const Q,double** const F) {
	const int nVar = GPRDIM::AbstractDIMSolver_ADERDG::NumberOfVariables;
	if(DIMENSIONS == 2){
		double F_3[nVar];
		pdeflux_(F[0], F[1],F_3, Q);
	}else{
		pdeflux_(F[0], F[1],F[2], Q);
	}
}


//You can either implement this method or modify fusedSource
void GPRDIM::DIMSolver_ADERDG::algebraicSource(const double* const Q,double* const S) {
  // @todo Please implement/augment if required
  //S[0] = 0.0;
  //S[1] = 0.0;
  //S[2] = 0.0;
  //S[3] = 0.0;
  //S[4] = 0.0;
  //S[5] = 0.0;
  //S[6] = 0.0;
  //S[7] = 0.0;
  //S[8] = 0.0;
  //S[9] = 0.0;
  //S[10] = 0.0;
  //S[11] = 0.0;
  //S[12] = 0.0;
  //S[13] = 0.0;
  //S[14] = 0.0;
  //S[15] = 0.0;
  //S[16] = 0.0;
  //S[17] = 0.0;
  //S[18] = 0.0;
  //S[19] = 0.0;
  //S[20] = 0.0;
  //S[21] = 0.0;
  //S[22] = 0.0;
  //S[23] = 0.0;
  pdesource_(S, Q);
}

void  GPRDIM::DIMSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
	pdencp_(BgradQ, Q, gradQ);
}





void GPRDIM::DIMSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* const observables, const double* const Q) const {
	assertion(NumberOfDMPObservables==1);
	ReadOnlyVariables vars(Q);
	observables[0]=Q[13];
	observables[1]=Q[13];
}

bool GPRDIM::DIMSolver_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {
		  
	int limvalue, LocNumberOfDMPObservables;
	LocNumberOfDMPObservables=1;
	pdelimitervalue_(&limvalue,&center[0],&LocNumberOfDMPObservables, observablesMin, observablesMax);
	if(limvalue>0){
		return false;
	}else{
		return true;
	};
	/*
	if (tarch::la::equals(t,0.0)) {
	  pdelimitervalue_(&limvalue,&center[0],&NumberOfDMPObservables, observablesMin, observablesMax);
	  //pdegeometriclimitervalue_(&limvalue,&center[0]);
	  if(limvalue>0){
		  return false;
	  }else{
		  return true;
	  };
	}else{
		return !wasTroubledInPreviousTimeStep;
	};
	*/
}


