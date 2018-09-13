// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "Plotter.h"
#include <kernels/GaussLegendreQuadrature.h>
#include "NavierStokesSolverDG.h"
#include "NavierStokesSolverDG_Variables.h"
//#include "NavierStokes.h"

NavierStokes::Plotter::Plotter(NavierStokes::NavierStokesSolverDG& solver) {
  // @TODO Please insert your code here.
  order = solver.Order;
}

NavierStokes::Plotter::~Plotter() {
}

void NavierStokes::Plotter::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void NavierStokes::Plotter::finishPlotting() {
  // @TODO Please insert your code here.
}

void NavierStokes::Plotter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  const int writtenUnknowns = 5;
  for (int i=0; i<DIMENSIONS + 2; i++){
    outputQuantities[i] = Q[i];
  }

  /*
  // TODO(Lukas): Make sure we use the correct constants
  // As we only consider air, this should be a given!
  const auto ns = NavierStokes::NavierStokes();
  AbstractEulerSolver::Variables vars(Q);

  const auto pressure = ns.evaluatePressure(vars.E(), vars.rho(), vars.j());
  const auto temperature = ns.evaluateTemperature(vars.rho(), pressure);

  const auto potT = temperature * std::pow((ns.referencePressure/pressure), (ns.gasConstant/ns.c_p));

  // Write potential temperature
  outputQuantities[DIMENSIONS + 2] = potT;
*/
  const auto& weights = kernels::gaussLegendreWeights[order];

  double weight = 1.0;
  for (int i = 0; i < DIMENSIONS; ++i) {
    weight *= weights[pos[i]];
  }
  outputQuantities[DIMENSIONS + 2] = weight;
}