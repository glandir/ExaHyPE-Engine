// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================

#include <sstream>
#include <ostream>
#include "exahype/plotters/Plotter.h"
#include "exahype/profilers/ProfilerFactory.h"
#include "exahype/solvers/Solver.h"
#include "exahype/solvers/SolverCoupling.h"

#include "kernels/KernelCalls.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/GaussLobattoQuadrature.h"
#include "kernels/LimiterProjectionMatrices.h"
#include "kernels/DGMatrices.h"
#include "kernels/DGBasisFunctions.h"
#include "buildinfo.h"

#include "ElasticWaveEquation.h"



void kernels::initSolvers(exahype::Parser& parser, std::vector<std::string>& cmdlineargs) {
  {
  // Create and register solver
  exahype::solvers::RegisteredSolvers.push_back( new ElasticWaveEquation3D::ElasticWaveEquation(parser.getMaximumMeshSize(0), parser.getMaximumMeshDepth(0), 0, 0, parser.getTimeStepping(0), cmdlineargs  ));
  parser.checkSolverConsistency(0);

  
  }

  std::set<int> orders;
  for (const auto p : exahype::solvers::RegisteredSolvers) {
    orders.insert(p->getNodesPerCoordinateAxis()-1);
  }
  kernels::initGaussLegendreNodesAndWeights(orders);
  kernels::initGaussLobattoNodesAndWeights(orders);
  kernels::initLimiterProjectionMatrices(orders);
  kernels::initDGMatrices(orders);
  kernels::initBasisFunctions(orders);
}


void kernels::finalise() {
  std::set<int> orders;
  for (const auto p : exahype::solvers::RegisteredSolvers) {
    orders.insert(p->getNodesPerCoordinateAxis()-1);
  }
  kernels::freeGaussLegendreNodesAndWeights(orders);
  kernels::freeGaussLobattoNodesAndWeights(orders);
  kernels::freeLimiterProjectionMatrices(orders);
  kernels::freeDGMatrices(orders);
  kernels::freeBasisFunctions(orders);

  for (auto solver : exahype::solvers::RegisteredSolvers) {
    delete solver;
  }
  exahype::solvers::RegisteredSolvers.clear();

  for (auto plotter : exahype::plotters::RegisteredPlotters) {
    delete plotter;
  }
  exahype::plotters::RegisteredPlotters.clear();
  for (auto coupling : exahype::solvers::RegisteredSolverCouplings) {
    delete coupling;
  }
  exahype::solvers::RegisteredSolverCouplings.clear();
}



void kernels::toString(std::ostream& ostream) {
/* Generated SolverRegistration code by the toolkit */

  ostream << "inputFileName: Benchmarks/coolmuc2/Elastic3D/single-core/Elastic3D-no-output.exahype\n";
  ostream << "projectName: ElasticWaveEquation3D\n";
  ostream << "Kernel[1].registration: AderdgSolver\n";
  ostream << "Kernel[1].type: ElasticWaveEquation3D::ElasticWaveEquation\n";
  ostream << "Kernel[1].parent: ";
  ElasticWaveEquation3D::AbstractElasticWaveEquation::constantsToString(ostream);
  ostream << "\n";
  ostream << "Kernel[1].hasConstants: null\n";
  ostream << "Kernel[1].variables: v 3 sigma 6 \n";
  ostream << "Kernel[1].kernel: (type: [linear ], terms: [pointsources, flux, materialparameters, ncp, source ], opt: [patchwiseadjust, generic ])\n";
}


