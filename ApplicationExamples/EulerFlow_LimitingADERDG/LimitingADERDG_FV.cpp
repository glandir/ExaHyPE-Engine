#include "LimitingADERDG_FV.h"

#include "InitialData.h"
#include "LimitingADERDG_FV_Variables.h"

#include <cstring>

void Euler::LimitingADERDG_FV::init(std::vector<std::string>& cmdlineargs) {
  // This function is called inside the constructur.
  // @todo Please implement/augment if required.
}

void Euler::LimitingADERDG_FV::flux(const double* const Q, double** F) {
  ReadOnlyVariables vars(Q);
  Fluxes f(F);

  tarch::la::Matrix<3,3,double> I;
  I = 1, 0, 0,
      0, 1, 0,
      0, 0, 1;

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  f.rho ( vars.j()                                 );
  f.j   ( irho * outerDot(vars.j(),vars.j()) + p*I );
  f.E   ( irho * (vars.E() + p) * vars.j()         );
}


void Euler::LimitingADERDG_FV::source(const double* const Q, double* S) {
  Variables s(S);
  s.rho()=0;
  s.E()=0;
  s.j(0,0,0);
}


void Euler::LimitingADERDG_FV::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  double u_n = vars.j(normalNonZeroIndex) * irho;
  double c   = std::sqrt(GAMMA * p * irho);

  eigs.rho()=u_n - c;
  eigs.E()  =u_n + c;
  eigs.j(u_n,u_n,u_n);
}


bool Euler::LimitingADERDG_FV::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt) {
  if ( tarch::la::equals(t,0.0) ) {
    return true;
  }
  return false; 
}


void Euler::LimitingADERDG_FV::adjustedSolutionValues(const double* const x, const double w, const double t, const double dt, double* Q) {  
  initialData(x,Q);
}

exahype::solvers::Solver::RefinementControl Euler::LimitingADERDG_FV::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void Euler::LimitingADERDG_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int normalNonZero,
    const double* const stateIn,
    double* stateOut) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];
}
