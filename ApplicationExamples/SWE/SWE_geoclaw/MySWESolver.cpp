#include "MySWESolver.h"
#include "InitialData.h"
#include "MySWESolver_Variables.h"
#include "swe.hpp"

#include "kernels/KernelUtils.h"

using namespace kernels;

double grav;
double epsilon;
int scenario;

enum SolverTypes {
  SOLVERS_ORIGINAL,
  SOLVERS_MINE,
  SOLVERS_SAMOA,
};
int solverType;

tarch::logging::Log SWE::MySWESolver::_log("SWE::MySWESolver");

void SWE::MySWESolver::init(const std::vector<std::string>& cmdlineargs,
                            const exahype::parser::ParserView& constants) {
  if (constants.isValueValidDouble("grav")) {
    grav = constants.getValueAsDouble("grav");
  }
  if (constants.isValueValidDouble("epsilon")) {
    epsilon = constants.getValueAsDouble("epsilon");
  }
  if (constants.isValueValidInt("scenario")) {
    scenario = constants.getValueAsInt("scenario");
  }
  if (constants.isValueValidInt("solverType")) {
    solverType = constants.getValueAsInt("solverType");
  }
}

void SWE::MySWESolver::adjustSolution(const double* const x, const double t,
                                      const double dt, double* const Q) {
  // Dimensions             = 2
  // Number of variables    = 4 + #parameters

  if (tarch::la::equals(t, 0.0)) {
    initialData(x, Q);
  } else {
    if (Q[0] < epsilon) {
      Q[1] = 0;
      Q[2] = 0;
    }
  }
}

void SWE::MySWESolver::eigenvalues(const double* const Q, const int dIndex,
                                   double* const lambda) {
  // Dimensions             = 2
  // Number of variables    = 4 + #parameters

  return swe::eigenvalues(Q, dIndex, lambda, grav, epsilon);
}

void SWE::MySWESolver::boundaryValues(const double* const x, const double t,
                                      const double dt, const int faceIndex,
                                      const int d,
                                      const double* const stateInside,
                                      double* const stateOutside) {
  // Dimensions             = 2
  // Number of variables    = 4 + #parameters

  // OL Code
  //  stateOutside[0] = stateInside[0];
  //  stateOutside[1] = 0.0;
  //  stateOutside[2] = 0.0;
  //  stateOutside[3] = 0.0;

  // normal Code
  stateOutside[0] = stateInside[0];
  stateOutside[1] = stateInside[1];
  stateOutside[2] = stateInside[2];
  stateOutside[3] = stateInside[3];

  // Wall
  stateOutside[d + 1] = -stateInside[d + 1];
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

// to add new PDEs specify them in the specification file, delete this file and
// its header and rerun the toolkit

void SWE::MySWESolver::flux(const double* const Q, double** const F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0

  swe::flux(Q, F, epsilon, grav);
}

void SWE::MySWESolver::nonConservativeProduct(const double* const Q,
                                              const double* const gradQ,
                                              double* const BgradQ) {
  idx2 idx_gradQ(DIMENSIONS, NumberOfVariables);

  BgradQ[0] = 0.0;
  BgradQ[1] = grav * Q[0] * gradQ[idx_gradQ(0, 3)];
  BgradQ[2] = grav * Q[0] * gradQ[idx_gradQ(1, 3)];
  BgradQ[3] = 0.0;
}

double SWE::MySWESolver::riemannSolver(double* fL, double* fR, const double* qL,
                                       const double* qR, const double* gradQL,
                                       const double* gradQR,
                                       const double* cellSize, int direction) {
  if (solverType == SOLVERS_ORIGINAL) {
    return swe::originalRiemannSolver(fL, fR, qL, qR, direction, grav, epsilon);
  } else if (solverType == SOLVERS_MINE) {
    return swe::riemannSolver(fL, fR, qL, qR, direction, grav, epsilon);
  } else if (solverType == SOLVERS_SAMOA) {
    return swe::samoaRiemannSolver(fL, fR, qL, qR, direction, grav, epsilon);
  } else {
    std::abort();
  }
}
