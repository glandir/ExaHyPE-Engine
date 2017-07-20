#include "GRMHDSolver_ADERDG.h"

#include "GRMHDSolver_ADERDG_Variables.h"
#include "InitialData/InitialData.h"
#include "Fortran/PDE.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

#include "exahype/disableOptimization.h" // we experience compiler bugs sometimes.

#include "BoundaryConditions/BoundaryConditions_ADERDG.h"


const double excision_radius = 1.0;
constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
constexpr int order = GRMHD::AbstractGRMHDSolver_ADERDG::Order;
constexpr int basisSize = order + 1;
constexpr int nDim = DIMENSIONS;

tarch::logging::Log GRMHD::GRMHDSolver_ADERDG::_log("GRMHDSolver_ADERDG");

typedef BoundaryConditions<GRMHD::GRMHDSolver_ADERDG> ADERDG_BC;
ADERDG_BC* abc;

void GRMHD::GRMHDSolver_ADERDG::init(std::vector<std::string>& cmdlineargs,exahype::Parser::ParserView constants) {
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

  abc = new ADERDG_BC(this);
  //if(!abc->setFromSpecFile<exahype::Parser::ParserView>(constants)) {
  if(!abc->setFromSpecFile<StringMapView>(StringMapView(tbc))) {
	logError("boundaryValues", "Some Boundary faces are missing in Specfile. Need: left,right,top,bottom,front,back. Got:" << tbc);
	std::abort();
  }
}

exahype::solvers::ADERDGSolver::AdjustSolutionValue  __attribute__((optimize("O0"))) GRMHD::GRMHDSolver_ADERDG::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
 // disable this in 2D:
 // bool insideExcisionBall = std::sqrt(center[0]*center[0] + center[1]*center[1] + center[2]*center[2]) < excision_radius;
  bool insideExcisionBall = false;
  bool hastoadjust = tarch::la::equals(t,0.0) || insideExcisionBall;
  return hastoadjust ? AdjustSolutionValue::PointWisely : AdjustSolutionValue::No;
}

void __attribute__((optimize("O0"))) GRMHD::GRMHDSolver_ADERDG::adjustPointSolution(const double* const x,const double w,const double t,const double dt,double* Q) {
  id->Interpolate(x, t, Q);
  //printf("Interpoalted at x=[%f,%f,%f], t=%f, Q0=%f\n", x[0],x[1],x[2], t, Q[0]);
  for(int i=0; i<NumberOfVariables; i++) {
	if(!std::isfinite(Q[i])) {
		printf("NAN in i=%d at t=%f, x=[%f,%f,%f], Q[%d]=%f\n", i, t, x[0],x[1],x[2], i, Q[i]);
	}
  }
}

void __attribute__((optimize("O0"))) GRMHD::GRMHDSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void GRMHD::GRMHDSolver_ADERDG::flux(const double* const Q,double** F) {
  pdeflux_(F[0], F[1], (DIMENSIONS==3)?F[2]:nullptr, Q);
}


void __attribute__((optimize("O0"))) GRMHD::GRMHDSolver_ADERDG::algebraicSource(const double* const Q,double* S) {
  pdesource_(S, Q);
}


void GRMHD::GRMHDSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int d,
  const double * const fluxIn,const double* const stateIn, double *fluxOut,double* stateOut) {
	 // for debugging, to make sure BC are set correctly
	double snan = std::numeric_limits<double>::signaling_NaN();
	double weird_number = -1.234567;
	std::memset(stateOut, snan, nVar * sizeof(double));
	std::memset(fluxOut,  snan, nVar * sizeof(double));
	
	//abc->apply(ADERDG_BOUNDARY_CALL);
	abc->exact(ADERDG_BOUNDARY_CALL);
}


exahype::solvers::Solver::RefinementControl GRMHD::GRMHDSolver_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

// only evaluated in Limiting context
void GRMHD::GRMHDSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(
    double* observables,const int numberOfObservables,
    const double* const Q) const {
  assertion(numberOfObservables==2);

  observables[0] = Q[0]; // rho
  observables[1] = Q[4]; // dens
}


bool GRMHD::GRMHDSolver_ADERDG::isPhysicallyAdmissible(
  const double* const solution,
  const double* const observablesMin,const double* const observablesMax,const int numberOfObservables,
  const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,
  const double t, const double dt) const {

  // geometric criterion:
  //  if ((center[0]-0.5)*(center[0]-0.5)+(center[1]-0.5)*(center[1]-0.5)<0.25*dx[0]*dx[0]) return false;

  // Static criterium for startup: When the density makes a large jump,
  // ie. at the star crust
  //if ( QMin[0] != 0.0 && QMax[0]/QMin[0] > 1e3 ) return false;

  if (observablesMin[0] < 0.0) return false;
  if (observablesMin[1] < 0.0) return false;

  // what about this kind of check?
  /*

  for (int i=0; i<nVar; ++i) {
    if (!std::isfinite(QMin[i])) return false;
    if (!std::isfinite(QMax[i])) return false;
  }
  
  */
  return true;
}

void __attribute__((optimize("O0"))) GRMHD::GRMHDSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
}


void GRMHD::GRMHDSolver_ADERDG::coefficientMatrix(const double* const Q,const int d,double* Bn) {
  // new scheme has no coefficient matrix
  static tarch::logging::Log _log("GRMHDSolver");
  logError("coefficientMatrix()", "ADERDG Coefficient Matrix invoked");
  exit(-2);

  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);
}
