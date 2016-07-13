// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================

#include <sstream>

#include "exahype/plotters/Plotter.h"
#include "exahype/profilers/ProfilerFactory.h"
#include "exahype/solvers/Solver.h"
#include "kernels/KernelCalls.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

#include "SRHDSolver.h"



void kernels::initSolvers(const exahype::Parser& parser) {
  std::string profiler_identifier = parser.getProfilerIdentifier();
  std::string metrics_identifier_list = parser.getMetricsIdentifierList();

  assertion1(metrics_identifier_list.find_first_of("{") == 0,
           metrics_identifier_list);
  assertion1(metrics_identifier_list.find_last_of("}") ==
                 metrics_identifier_list.size() - 1,
             metrics_identifier_list);

  // Split "{metric1,metric2...}" into {"metric1", "metric2", ...}
  std::vector<std::string> metrics_vector;
  std::stringstream ss;
  ss << metrics_identifier_list.substr(1, metrics_identifier_list.size() - 2);
  std::string metric;
  while (std::getline(ss, metric, ',')) {
    metrics_vector.emplace_back(std::move(metric));
  }

  // Create profiler
  auto profiler = exahype::profilers::ProfilerFactory::getInstance().create(
    profiler_identifier, metrics_vector);

  // Create and register solver
  exahype::solvers::RegisteredSolvers.push_back( new SRHD::SRHDSolver(0, parser.getMaximumMeshSize(0), parser.getTimeStepping(0), std::move(profiler)));
  parser.checkSolverConsistency(0);

  
  exahype::plotters::RegisteredPlotters.push_back( new exahype::plotters::Plotter(0,0,parser));


  std::set<int> orders;
  for (const auto p : exahype::solvers::RegisteredSolvers) {
    orders.insert(p->getNodesPerCoordinateAxis()-1);
  }
  kernels::initGaussLegendreNodesAndWeights(orders);
  kernels::initDGMatrices(orders);
}


void kernels::finalise() {
  std::set<int> orders;
  for (const auto p : exahype::solvers::RegisteredSolvers) {
    orders.insert(p->getNodesPerCoordinateAxis()-1);
  }
  kernels::freeGaussLegendreNodesAndWeights(orders);
  kernels::freeDGMatrices(orders);

  for (auto solver : exahype::solvers::RegisteredSolvers) {
    delete solver;
  }
  exahype::solvers::RegisteredSolvers.clear();

  for (auto plotter : exahype::plotters::RegisteredPlotters) {
    delete plotter;
  }
  exahype::plotters::RegisteredPlotters.clear();
}



