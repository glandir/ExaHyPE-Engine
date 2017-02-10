#include "MHDSolver_ADERDG.h"
//#include "fortran.h" _ltob

#include "InitialData.h"
#include "BoundaryConditions.h"
#include "PDE.h"

#include <memory>
#include <cstring>
#include <stdio.h>
#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h" // matrix indexing

// storage
InitialDataHandler idfunc;
BoundaryConditionHandler bcfunc;

void MHD::MHDSolver_ADERDG::init(std::vector<std::string>& cmdargs) {
	/**
	 * InitialData with consistent BC can be set by environment:
	 *    export EXAHYPE_INITIALDATA="Blast"
	 * or by command line:
	 *    ./ExaHyPE-MHD ../MHD_Limiting.exahype Blast
	 * You can also set the default_id to trigger another ID.
	 *
	 * Obey that you need the matching grid for your simulation.
	 **/

	static tarch::logging::Log _log("MHDSolver_ADERDG::init");
	// Don't rely on exahype specfile parameters, use Environment
	// variables or command line variables instead.
	const char* env_id = std::getenv("EXAHYPE_INITIALDATA");
	const char* cmd_id = cmdargs.size()>1 ? cmdargs[1].c_str() : nullptr;
	std::string default_id = "Jet";

	std::string id, bc;	
	std::map<std::string, InitialDataHandler> ids;
	std::map<std::string, BoundaryConditionHandler> bcs;
	
	ids["Jet"] = initialjet_;
	bcs["Jet"] = boundaryjet_;

	ids["AlfenWave"] = initialalfenwave_;
	bcs["AlfenWave"] = BoundaryAlfenWave;
	
	ids["Rotor"] = initialrotor_;
	bcs["Rotor"] = boundaryoutflow_;
	
	ids["Blast"] = initialblast_;
	bcs["Blast"] = boundaryoutflow_;
	
	ids["OrsagTang"] = initialorsagtang_;
	bcs["OrsagTang"] = boundaryoutflow_;
	
	ids["ShockTube"] = initialshocktube_;
	bcs["ShockTube"] = boundaryoutflow_;
	
	if(env_id) id = env_id;
	else if(cmd_id) id = cmd_id;
	else id = default_id;
	
	if(id==default_id)
		logError("ID", "Loading default Initial Data: " << id);
	if(!ids.count(id) && !bcs.count(id)) {
		logError("ID", "Cannot understand requested ID: " << id);
		exit(-1);
	} else {
		logInfo("ID", "Will load requested Initial Data: " << id);
	}
	
	// TODO: Lowercase key lookup
	idfunc = ids[id];
	bcfunc = bcs[id];
}

void MHD::MHDSolver_ADERDG::flux(const double* const Q, double** F) {
  // Caveats: Fortran accepts a uniform array of size (nVar*nDim), however C passes an array of pointers.
  // This Fortran interface works only if F is a continous array and F[1]==F[nDim+1] etc!
  pdeflux_(F[0], Q);
}

void MHD::MHDSolver_ADERDG::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  double nv[3] = {0.};
  nv[normalNonZeroIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

bool MHD::MHDSolver_ADERDG::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t, double dt) {
  return (t < 1e-10);
}

void MHD::MHDSolver_ADERDG::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  if (tarch::la::equals(t, 0.0)) {
    idfunc(x, Q);
  }
}

void MHD::MHDSolver_ADERDG::source(const double* const Q, double* S) {
  pdesource_(S, Q);
}

exahype::solvers::Solver::RefinementControl MHD::MHDSolver_ADERDG::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void MHD::MHDSolver_ADERDG::boundaryValues(const double* const x,const double t, const double dt, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut) {
  // pass to Fortran/generic function pointer
  double nv[3] = {0.};
  nv[normalNonZero] = 1;
  bcfunc(x, &t, &dt, &faceIndex, nv, fluxIn, stateIn, fluxOut, stateOut);
}

void MHD::MHDSolver_ADERDG::ncp(const double* const Q, const double* const gradQ, double* BgradQ) {
  constexpr int nVar = 9;
  std::memset(BgradQ, 0, nVar * sizeof(double));
}

void MHD::MHDSolver_ADERDG::matrixb(const double* const Q, const int normalNonZero, double* Bn) {
  constexpr int nVar = 9;
  std::memset(Bn, 0, nVar * nVar * sizeof(double));
}


bool MHD::MHDSolver_ADERDG::physicalAdmissibilityDetection(const double* const QMin,const double* const QMax) {
  if (QMin[0] < 0.0) return false;
  if (QMin[4] < 0.0) return false;

  for (int i=0; i<5; ++i) {
    if (!std::isfinite(QMin[i])) return false;
    if (!std::isfinite(QMax[i])) return false;
  }

  return true;
}

bool MHD::MHDSolver_ADERDG::isDummyKRequired() const {
  return false;
}

void MHD::MHDSolver_ADERDG::dummyK_Value(const double* const x,const double t,const double dt, double* forceVector, double* x0) {
  // do nothing
}


