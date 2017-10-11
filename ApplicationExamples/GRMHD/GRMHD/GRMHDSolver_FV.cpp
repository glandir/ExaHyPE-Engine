#include "GRMHDSolver_FV.h"

#include "GRMHDSolver_FV_Variables.h"

#include "AbstractGRMHDSolver_ADERDG.h"
#include "InitialData/InitialData.h"
#include "Fortran/PDE.h"

#include <stdio.h>
#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

#include "exahype/disableOptimization.h" // bugs when limiting is on. whatevers

#include "BoundaryConditions/BoundaryConditions_FV.h"

const double excision_radius = 1.0;

tarch::logging::Log GRMHD::GRMHDSolver_FV::_log("GRMHDSolver_FV");
constexpr int nVar = GRMHD::AbstractGRMHDSolver_FV::NumberOfVariables;

typedef BoundaryConditions<GRMHD::GRMHDSolver_FV> FV_BC;
FV_BC* fvbc;

// enable nan tracker
#include <fenv.h>

void GRMHD::GRMHDSolver_FV::init(std::vector<std::string>& cmdlineargs,exahype::Parser::ParserView& constants) {
  // try NaN catcher
  // feenableexcept(FE_INVALID | FE_OVERFLOW);  // Enable all floating point exceptions but FE_INEXACT

  // Just copy and pasted from GRMHDSolver_ADERDG
  // Todo: Move this to specfile once we have working constants.
  // Todo: Move this to specfile once we have working constants.
  std::string id_default = "Fortran";
  std::string bc_default = "left:exact,right:exact,top:exact,bottom:exact,front:exact,back:exact";

  // alternatives:
  //std::string id_RNSID = "RNSID";
  //std::string bc_RNSID_octant = "left:refl,right:exact,bottom:refl,top:exact,front:refl,back:exact";

  // try to obtain requested initial data and boundary conditions from the
  // environment variables, as the specfile parameter system is still broken.
  std::string tid = getenv("EXAHYPE_ID") ? getenv("EXAHYPE_ID") : id_default;
  std::string tbc = getenv("EXAHYPE_BC") ? getenv("EXAHYPE_BC") : bc_default;

  if(!prepare_id(tid)) {
	  logError("prepare_id", "Could not setup Initial Data '" << tid << "', probably misspelled.");
	  std::abort();
  }

  fvbc = new FV_BC(this);
  if(!fvbc->setFromSpecFile<StringMapView>(StringMapView(tbc))) {
	logError("boundaryValues", "Some Boundary faces are missing in Specfile. Need: left,right,top,bottom,front,back. Got:" << tbc);
	std::abort();
  }
  
  /*
  prepare_id();
  
  fvbc = new FV_BC(this);

  // face indices: 0 x=xmin 1 x=xmax, 2 y=ymin 3 y=ymax 4 z=zmin 5 z=zmax
  // corresponding 0-left, 1-right, 2-front, 3-back, 4-bottom, 5-top
  fvbc->left = fvbc->parseFromString("exact");
  fvbc->right = fvbc->parseFromString("exact");
  fvbc->front = fvbc->parseFromString("exact");
  fvbc->back = fvbc->parseFromString("exact");
  fvbc->bottom = fvbc->parseFromString("exact");
  fvbc->top = fvbc->parseFromString("exact");
  
  if(!fvbc->allFacesDefined()) {
	logError("boundaryValues", "Some Boundary faces are not defined");
	std::abort();
  }
  */

}


void GRMHD::GRMHDSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* Q) {
  using namespace tarch::la;
  // Do excision only in 3D.
  //bool insideExcisionBall = norm2(x) < excision_radius;
  bool insideExcisionBall = false;
  bool useAdjustSolution = equals(t,0.0) || insideExcisionBall;

  if(useAdjustSolution) id->Interpolate(x, t, Q);
}

void GRMHD::GRMHDSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  double nv[3] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void GRMHD::GRMHDSolver_FV::flux(const double* const Q, double** F) {
  pdeflux_(F[0], F[1], (DIMENSIONS==3)?F[2]:nullptr, Q);
}

// Source is exactly 0
/*
void GRMHD::GRMHDSolver_FV::algebraicSource(const double* const Q, double* S) {
  pdesource_(S, Q);
}
*/

void GRMHD::GRMHDSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateIn,
    double* stateOut) {
	fvbc->apply(FV_BOUNDARY_CALL);
	// Use for the time being: Exact BC
	//id->Interpolate(x, t, stateOut);
}


void GRMHD::GRMHDSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
}

/*
void GRMHD::GRMHDSolver_FV::coefficientMatrix(const double* const Q,const int d,double* Bn) {
  // new FV scheme has no coefficient matrix
  static tarch::logging::Log _log("GRMHDSolver");
  logError("coefficientMatrix()", "Coefficient Matrix invoked");
  exit(-2);

  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);
}
*/
