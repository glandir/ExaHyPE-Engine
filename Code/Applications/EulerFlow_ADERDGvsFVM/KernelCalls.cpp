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
#include "kernels/GaussLobattoQuadrature.h"
#include "kernels/LimiterProjectionMatrices.h"
#include "kernels/DGMatrices.h"
#include "kernels/DGBasisFunctions.h"

#include "exahype/solvers/SingleSolverCoupling.h"
#include "exahype/solvers/ADERDGAPosterioriSubcellLimiter.h"

#include "ADERDG.h"
#include "ADERDG_Plotter0.h"
#include "FVM.h"
#include "FVM_Plotter0.h"



void kernels::initSolvers(exahype::Parser& parser) {
  {
  std::string profiler_identifier = parser.getProfilerIdentifier();
  std::string metrics_identifier_list = parser.getMetricsIdentifierList();
  std::string profiling_output = parser.getProfilingOutputFilename();

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
    profiler_identifier, metrics_vector, profiling_output);

  // Create and register solver
  exahype::solvers::RegisteredSolvers.push_back( new Euler::ADERDG(parser.getMaximumMeshSize(0), parser.getTimeStepping(0), std::move(profiler)
  ));
//  exahype::solvers::RegisteredSolverCouplings.push_back( new exahype::solvers::SingleSolverCoupling(0) );
  parser.checkSolverConsistency(0);

  
  }
  exahype::plotters::RegisteredPlotters.push_back( new exahype::plotters::Plotter(0,0,parser,new Euler::ADERDG_Plotter0(  *static_cast<Euler::ADERDG*>(exahype::solvers::RegisteredSolvers[0])) ));

  {
  std::string profiler_identifier = parser.getProfilerIdentifier();
  std::string metrics_identifier_list = parser.getMetricsIdentifierList();
  std::string profiling_output = parser.getProfilingOutputFilename();

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
    profiler_identifier, metrics_vector, profiling_output);

  // Create and register solver
  exahype::solvers::RegisteredSolvers.push_back( new Euler::FVM(7, parser.getMaximumMeshSize(1), parser.getTimeStepping(1), std::move(profiler)  ));
//  exahype::solvers::RegisteredSolverCouplings.push_back( new exahype::solvers::SingleSolverCoupling(1) );
  parser.checkSolverConsistency(1);

  
  }
  exahype::plotters::RegisteredPlotters.push_back( new exahype::plotters::Plotter(1,0,parser,new Euler::FVM_Plotter0(  *static_cast<Euler::FVM*>(exahype::solvers::RegisteredSolvers[1])) ));

  // Create and register couplings
  exahype::solvers::RegisteredSolverCouplings.push_back( new exahype::solvers::ADERDGAPosterioriSubcellLimiter(0,1) );


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



