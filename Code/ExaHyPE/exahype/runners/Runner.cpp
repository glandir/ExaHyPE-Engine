/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#include "exahype/runners/Runner.h"

#include "../../../Peano/mpibalancing/HotspotBalancing.h"

#include "exahype/repositories/Repository.h"
#include "exahype/repositories/RepositoryFactory.h"
#include "exahype/mappings/TimeStepSizeComputation.h"
#include "exahype/mappings/Sending.h"

#include "tarch/Assertions.h"

#include "tarch/logging/CommandLineLogger.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "tarch/parallel/FCFSNodePoolStrategy.h"

#include "tarch/multicore/Core.h"
#include "tarch/multicore/MulticoreDefinitions.h"

#include "peano/parallel/JoinDataBufferPool.h"
#include "peano/parallel/JoinDataBufferPool.h"
#include "peano/parallel/SendReceiveBufferPool.h"
#include "peano/parallel/loadbalancing/Oracle.h"
#include "peano/parallel/loadbalancing/OracleForOnePhaseWithGreedyPartitioning.h"

#include "peano/geometry/Hexahedron.h"

#include "peano/utils/UserInterface.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/datatraversal/autotuning/OracleForOnePhaseDummy.h"

#include "sharedmemoryoracles/OracleForOnePhaseWithGrainSizeSampling.h"
#include "sharedmemoryoracles/OracleForOnePhaseWithShrinkingGrainSize.h"

#ifdef Parallel
#include "mpibalancing/GreedyBalancing.h"
#include "mpibalancing/FairNodePoolStrategy.h"
#endif
#include "exahype/plotters/Plotter.h"

#include "exahype/solvers/ADERDGSolver.h"

tarch::logging::Log exahype::runners::Runner::_log("exahype::runners::Runner");

exahype::runners::Runner::Runner(const Parser& parser) : _parser(parser) {}

exahype::runners::Runner::~Runner() {}

void exahype::runners::Runner::initDistributedMemoryConfiguration() {
  #ifdef Parallel
  std::string configuration = _parser.getMPIConfiguration();
  if (_parser.getMPILoadBalancingType()==Parser::MPILoadBalancingType::Static) {
    if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      if (configuration.find( "FCFS" )!=std::string::npos ) {
        tarch::parallel::NodePool::getInstance().setStrategy(
            new tarch::parallel::FCFSNodePoolStrategy()
        );
        logInfo("initDistributedMemoryConfiguration()", "load balancing relies on FCFS answering strategy");
      }
      else if (configuration.find( "fair" )!=std::string::npos ) {
        int ranksPerNode = static_cast<int>(exahype::Parser::getValueFromPropertyString(configuration,"ranks_per_node"));
        if (ranksPerNode<=0) {
          logError( "initDistributedMemoryConfiguration()", "please inform fair balancing how many ranks per node you use through value \"ranks_per_node:XXX\". Read value " << ranksPerNode << " is invalid" );
          ranksPerNode = 1;
        }
        if ( ranksPerNode>=tarch::parallel::Node::getInstance().getNumberOfNodes() ) {
          logWarning( "initDistributedMemoryConfiguration()", "value \"ranks_per_node:XXX\" exceeds total rank count. Reset to 1" );
          ranksPerNode = 1;
        }
        tarch::parallel::NodePool::getInstance().setStrategy(
            new mpibalancing::FairNodePoolStrategy(ranksPerNode)
        );
        logInfo("initDistributedMemoryConfiguration()", "load balancing relies on fair answering strategy with " << ranksPerNode << " rank(s) per node") ;
      }
      else {
        logError("initDistributedMemoryConfiguration()", "no valid load balancing answering strategy specified: use FCFS");
        tarch::parallel::NodePool::getInstance().setStrategy(
            new tarch::parallel::FCFSNodePoolStrategy()
        );
      }
    }

    if ( configuration.find( "greedy" )!=std::string::npos ) {
      logInfo("initDistributedMemoryConfiguration()", "use greedy load balancing without joins (mpibalancing/GreedyBalancing)");
      peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
        new mpibalancing::GreedyBalancing(
          getCoarsestGridLevelOfAllSolvers(),
          getCoarsestGridLevelOfAllSolvers()+1
        )
      );
    }
    else if ( configuration.find( "hotspot" )!=std::string::npos ) {
      logInfo("initDistributedMemoryConfiguration()", "use global hotspot elimination without joins (mpibalancing/StaticBalancing)");
      peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
          new mpibalancing::HotspotBalancing(false,getCoarsestGridLevelOfAllSolvers()+1)
      );
    }
    else {
      logError("initDistributedMemoryConfiguration()", "no valid load balancing configured. Use greedy load balancing without joins");
      peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
          new peano::parallel::loadbalancing::OracleForOnePhaseWithGreedyPartitioning(false)
      );
    }
  }
  else {
    logError("initDistributedMemoryConfiguration()", "only MPI static load balancing supported so far. ");
  }

  // end of static load balancing


  tarch::parallel::NodePool::getInstance().restart();
  tarch::parallel::NodePool::getInstance().waitForAllNodesToBecomeIdle();

  tarch::parallel::Node::getInstance().setDeadlockTimeOut(_parser.getMPITimeOut());
  tarch::parallel::Node::getInstance().setTimeOutWarning(_parser.getMPITimeOut()/2);
  logInfo("initDistributedMemoryConfiguration()", "use MPI time out of " << _parser.getMPITimeOut() << " (warn after half the timeout span)");

  const int bufferSize = _parser.getMPIBufferSize();
  peano::parallel::SendReceiveBufferPool::getInstance().setBufferSize(bufferSize);
  peano::parallel::JoinDataBufferPool::getInstance().setBufferSize(bufferSize);
  logInfo("initDistributedMemoryConfiguration()", "use MPI buffer size of " << bufferSize);

  if ( _parser.getSkipReductionInBatchedTimeSteps() ) {
    logInfo("initDistributedMemoryConfiguration()", "allow ranks to skip reduction" );
    exahype::mappings::Sending::SkipReductionInBatchedTimeSteps = true;
  }
  else {
    logWarning("initDistributedMemoryConfiguration()", "ranks are not allowed to skip any reduction (might harm performance). Use optimisation section to switch feature on" );
    exahype::mappings::Sending::SkipReductionInBatchedTimeSteps = false;
  }
  #endif
}


void exahype::runners::Runner::shutdownDistributedMemoryConfiguration() {
#ifdef Parallel
  tarch::parallel::NodePool::getInstance().terminate();
  exahype::repositories::RepositoryFactory::getInstance().shutdownAllParallelDatatypes();
#endif
}

void exahype::runners::Runner::initSharedMemoryConfiguration() {
#ifdef SharedMemoryParallelisation
  const int numberOfThreads = _parser.getNumberOfThreads();
  tarch::multicore::Core::getInstance().configure(numberOfThreads);

  switch (_parser.getMulticoreOracleType()) {
  case Parser::MulticoreOracleType::Dummy:
    logInfo("initSharedMemoryConfiguration()",
        "use dummy shared memory oracle");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
        new peano::datatraversal::autotuning::OracleForOnePhaseDummy(
            true, false,             // useMultithreading, measureRuntimes
            0,                       // grainSizeOfUserDefinedRegions
            peano::datatraversal::autotuning::OracleForOnePhaseDummy::SplitVertexReadsOnRegularSubtree::Split,
            true, true,              // pipelineDescendProcessing, pipelineAscendProcessing
            tarch::la::aPowI(DIMENSIONS,3*3*3*3/2), 3,  // smallestGrainSizeForAscendDescend
            1,1                      // smallestGrainSizeForEnterLeaveCell, grainSizeForEnterLeaveCell
            )
    );
    break;
  case Parser::MulticoreOracleType::Autotuning:
    logInfo("initSharedMemoryConfiguration()",
        "use autotuning shared memory oracle");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
        new sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize());
    break;
  case Parser::MulticoreOracleType::GrainSizeSampling:
    logInfo("initSharedMemoryConfiguration()",
        "use shared memory oracle sampling");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
        new sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling(
            32,
            false,  // useThreadPipelining,
            true    // logarithmicDistribution
        ));
    break;
  }

  std::ifstream f(_parser.getMulticorePropertiesFile().c_str());
  if (f.good()) {
    peano::datatraversal::autotuning::Oracle::getInstance().loadStatistics(
        _parser.getMulticorePropertiesFile());
  }
  f.close();
#endif
}


void exahype::runners::Runner::initDataCompression() {
  exahype::solvers::ADERDGSolver::CompressionAccuracy = _parser.getDoubleCompressionFactor();

  if (exahype::solvers::ADERDGSolver::CompressionAccuracy==0.0) {
    logInfo( "initDataCompression()", "switched off any data compression");
  }
  else {
    if (!_parser.getFuseAlgorithmicSteps()) {
      logError( "initDataCompression()", "data compression is not supported if you don't use the fused time stepping");
      exahype::solvers::ADERDGSolver::CompressionAccuracy = 0.0;
    }
    else {
      exahype::solvers::ADERDGSolver::SpawnCompressionAsBackgroundThread = _parser.getSpawnDoubleCompressionAsBackgroundTask();
      logInfo( "initDataCompression()", "store all data with accuracy of " << exahype::solvers::ADERDGSolver::CompressionAccuracy << ". Use background threads for data conversion=" << exahype::solvers::ADERDGSolver::SpawnCompressionAsBackgroundThread);
    }
  }
}


void exahype::runners::Runner::shutdownSharedMemoryConfiguration() {
#ifdef SharedMemoryParallelisation
  switch (_parser.getMulticoreOracleType()) {
  case Parser::MulticoreOracleType::Dummy:
    break;
  case Parser::MulticoreOracleType::Autotuning:
  case Parser::MulticoreOracleType::GrainSizeSampling:
#ifdef Parallel
    if (tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getNumberOfNodes()-1) {
      logInfo("shutdownSharedMemoryConfiguration()",
          "wrote statistics into file " << _parser.getMulticorePropertiesFile()
          << ". Dump from all other ranks subpressed to avoid file races"
      );
      peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics(
          _parser.getMulticorePropertiesFile());
    }
#else
    logInfo("shutdownSharedMemoryConfiguration()",
        "wrote statistics into file "
        << _parser.getMulticorePropertiesFile());
    peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics(
        _parser.getMulticorePropertiesFile());
#endif
    break;
  }
#endif
}


int exahype::runners::Runner::getCoarsestGridLevelOfAllSolvers() const {
  double boundingBox = _parser.getBoundingBoxSize()(0);
  double hMax        = exahype::solvers::Solver::getCoarsestMeshSizeOfAllSolvers();

  int    result      = 1;
  double currenthMax = std::numeric_limits<double>::max();
  while (currenthMax>hMax) {
    currenthMax = boundingBox / threePowI(result);
    result++;
  }

  logDebug( "getCoarsestGridLevelOfAllSolvers()", "regular grid depth of " << result << " creates grid with h_max=" << currenthMax );
  return std::max(3,result);
}


exahype::repositories::Repository* exahype::runners::Runner::createRepository() const {
  // Geometry is static as we need it survive the whole simulation time.
  static peano::geometry::Hexahedron geometry(
      _parser.getDomainSize(),
      _parser.getOffset());

  logDebug(
      "createRepository(...)",
      "create computational domain at " << _parser.getOffset() <<
      " of width/size " << _parser.getDomainSize() <<
      ". bounding box has size " << _parser.getBoundingBoxSize() <<
      ". grid regular up to level " << getCoarsestGridLevelOfAllSolvers() << " (level 1 is coarsest available cell in tree)" );
#ifdef Parallel
  const double boundingBoxScaling = static_cast<double>(getCoarsestGridLevelOfAllSolvers()) / (static_cast<double>(getCoarsestGridLevelOfAllSolvers())-2);
  assertion4(boundingBoxScaling>=1.0, boundingBoxScaling, getCoarsestGridLevelOfAllSolvers(), _parser.getDomainSize(), _parser.getBoundingBoxSize() );
  const double boundingBoxShift   = (1.0-boundingBoxScaling)/2.0;
  assertion5(boundingBoxShift<=0.0, boundingBoxScaling, getCoarsestGridLevelOfAllSolvers(), _parser.getDomainSize(), _parser.getBoundingBoxSize(), boundingBoxScaling );

  logInfo(
      "createRepository(...)",
      "increase domain artificially by " << boundingBoxScaling << " and shift bounding box by " << boundingBoxShift << " to simplify load balancing along boundary");
  return exahype::repositories::RepositoryFactory::getInstance().createWithSTDStackImplementation(
      geometry,
      _parser.getBoundingBoxSize()*boundingBoxScaling,
      _parser.getOffset()+boundingBoxShift*_parser.getBoundingBoxSize()
  );
#else
  return exahype::repositories::RepositoryFactory::getInstance().createWithSTDStackImplementation(
      geometry,
      _parser.getBoundingBoxSize(),
      _parser.getOffset()
  );
#endif
}


int exahype::runners::Runner::run() {
  exahype::repositories::Repository* repository = createRepository();

  initDistributedMemoryConfiguration();
  initSharedMemoryConfiguration();
  initDataCompression();

  int result = 0;
  if ( _parser.isValid() ) {
    if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      result = runAsMaster(*repository);
    }
    #ifdef Parallel
    else {
      result = runAsWorker(*repository);
    }
    #endif
  }
  else {
    logError( "run(...)", "do not run code as parser reported errors" );
    result = 1;
  }

  shutdownSharedMemoryConfiguration();
  shutdownDistributedMemoryConfiguration();

  delete repository;

  return result;
}

void exahype::runners::Runner::createGrid(exahype::repositories::Repository& repository) {
#ifdef Parallel
  const bool UseStationaryCriterion = tarch::parallel::Node::getInstance().getNumberOfNodes()==1;
#else
  const bool UseStationaryCriterion = true;
#endif

  int gridSetupIterations = 0;
  repository.switchToMeshRefinement();

  int gridSetupIterationsToRun = 3;
  while (gridSetupIterationsToRun>0) {
    repository.iterate();
    gridSetupIterations++;

    if ( UseStationaryCriterion && repository.getState().isGridStationary() ) {
      gridSetupIterationsToRun--;
    }
    else if ( !repository.getState().isGridBalanced() && tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()>0 ) {
      gridSetupIterationsToRun=3;  // we need at least 3 sweeps to recover from ongoing balancing
    }
    else if ( !repository.getState().isGridBalanced()  ) {
      gridSetupIterationsToRun=1;  // one additional step to get adjacency right
    }
    else {
      gridSetupIterationsToRun--;
    }

    #if defined(TrackGridStatistics) && defined(Asserts)
    logInfo("createGrid()",
        "grid setup iteration #" << gridSetupIterations <<
        ", max-level=" << repository.getState().getMaxLevel() <<
        ", state=" << repository.getState().toString() <<
        ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
    );
    #elif defined(Asserts)
    logInfo("createGrid()",
        "grid setup iteration #" << gridSetupIterations <<
        ", state=" << repository.getState().toString() <<
        ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
    );
    #elif defined(TrackGridStatistics)
    logInfo("createGrid()",
        "grid setup iteration #" << gridSetupIterations <<
        ", max-level=" << repository.getState().getMaxLevel() <<
        ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
    );
    #else
    logInfo("createGrid()",
        "grid setup iteration #" << gridSetupIterations <<
        ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
    );
    #endif

    #if !defined(Parallel)
    logInfo("createGrid(...)", "memoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");
    #endif

    #ifdef Asserts
    if (exahype::solvers::ADERDGSolver::CompressionAccuracy>0.0) {
      DataHeap::getInstance().plotStatistics();
      peano::heap::PlainCharHeap::getInstance().plotStatistics();
    }
    #endif
  }

  logInfo("createGrid(Repository)", "finished grid setup after " << gridSetupIterations << " iterations" );

  if (
    tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()>0
    &&
    tarch::parallel::Node::getInstance().getNumberOfNodes()>1
  ) {
    logWarning( "createGrid(Repository)", "there are still " << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes() << " ranks idle" )
  }

#ifdef Parallel
  // Might be too restrictive for later runs. Remove but keep warning from above
  assertion( tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()==0 );
#endif
}


int exahype::runners::Runner::runAsMaster(exahype::repositories::Repository& repository) {
  peano::utils::UserInterface::writeHeader();

  initSolverTimeStamps();
  repository.getState().switchToPreAMRContext();
  createGrid(repository);
  /*
   * Set ADER-DG corrector time stamp and finite volumes time stamp.
   * Compute ADER-DG corrector time step size implicitly and finite volumes time step size.
   * (Implicitly means here that we set the predictor time step size but after the next newTimeStep(...)
   * the corrector time step size is set as this predictor time step size.)
   *
   * Note that it is important that we run SolutionAdjustmentAnd
   * GlobalTimeStepComputation directly after the grid setup
   * since we receive here the metadata
   * that was sent in the last iteration of the grid setup.
   */
  initSolverTimeStamps();

  repository.getState().switchToInitialConditionAndTimeStepSizeComputationContext();
  repository.switchToInitialConditionAndTimeStepSizeComputation();
  repository.iterate();
  logInfo( "runAsMaster(...)", "initialised all data and computed first time step size" );

  /*
   * Set the time stamps of the solvers to the initial value again.
   *
   * !!! Rationale
   *
   * The time step size computation
   * sets the predictor time stamp to the value
   * of the predictor time stamp plus the admissible
   * time step size on each patch for each solver.
   */
  initSolverTimeStamps();
  /*
   * Compute current first predictor based on current time step size.
   * Set current time step size as old time step size of next iteration.
   * Compute the current time step size of the next iteration.
   */
  repository.getState().switchToPredictionAndTimeStepSizeComputationContext();
  bool plot = exahype::plotters::isAPlotterActive(
      solvers::Solver::getMinSolverTimeStampOfAllSolvers());
  if (plot) {
    #if DIMENSIONS==2
    repository.switchToPredictionAndPlotAndTimeStepSizeComputation2d();
    #else
    repository.switchToPredictionAndPlotAndTimeStepSizeComputation();
    #endif
  }
  else {
    repository.switchToPredictionAndTimeStepSizeComputation();
  }
  repository.iterate();
  /*
   * Reset the time stamps of the finite volumes solvers.
   *
   * !!! Rationale
   * Unlike for the rearranged ADER-DG scheme, we only
   * need one warm up iteration for the finite volumes solvers.
   * But since we have to perform two time step computations
   * to warm up the ADER-DG schemes,
   * the finite volumes solvers think they are already
   * advanced by one time step after the computation
   * of the second time step size.
   */
  initFiniteVolumesSolverTimeStamps();

  /*
   * Finally print the initial time step info.
   */
  printTimeStepInfo(-1);

  const double simulationEndTime = _parser.getSimulationEndTime();

  logDebug("runAsMaster(...)","min solver time stamp: "     << solvers::Solver::getMinSolverTimeStampOfAllSolvers());
  logDebug("runAsMaster(...)","min solver time step size: " << solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers());

  while ((solvers::Solver::getMinSolverTimeStampOfAllSolvers() < simulationEndTime) &&
      tarch::la::greater(solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers(), 0.0)) {
    bool plot = exahype::plotters::isAPlotterActive(
        solvers::Solver::getMinSolverTimeStampOfAllSolvers());

    if (_parser.getFuseAlgorithmicSteps()) {
      repository.getState().setTimeStepSizeWeightForPredictionRerun(
          _parser.getFuseAlgorithmicStepsFactor());

      int numberOfStepsToRun = 1;
      if (plot) {
        numberOfStepsToRun = 0;
      }
      else if (solvers::Solver::allSolversUseTimeSteppingScheme(solvers::Solver::TimeStepping::GlobalFixed)) {
        /**
         * This computation is optimistic. If we were pessimistic, we had to
         * use the max solver time step size. However, this is not necessary
         * here, as we half the time steps anyway.
         */
        if (solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers()>0.0) {
          const double timeIntervalTillNextPlot = std::min(exahype::plotters::getTimeOfNextPlot(),simulationEndTime) - solvers::Solver::getMaxSolverTimeStampOfAllSolvers();
          numberOfStepsToRun = std::floor( timeIntervalTillNextPlot / solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers() * _parser.getTimestepBatchFactor() );
        }
        numberOfStepsToRun = numberOfStepsToRun<1 ? 1 : numberOfStepsToRun;
      }

      runOneTimeStampWithFusedAlgorithmicSteps(
        repository,
        numberOfStepsToRun,
        _parser.getExchangeBoundaryDataInBatchedTimeSteps() && repository.getState().isGridStationary()
      );
      recomputePredictorIfNecessary(repository);
      printTimeStepInfo(numberOfStepsToRun);
    } else {
      runOneTimeStampWithThreeSeparateAlgorithmicSteps(repository, plot);
      printTimeStepInfo(1);
    }

    logDebug("runAsMaster(...)", "state=" << repository.getState().toString());
  }
  if ( tarch::la::equals(solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers(), 0.0)) {
    logWarning("runAsMaster(...)","Minimum solver time step size is zero (up to machine precision).");
  }

  repository.logIterationStatistics(true);
  repository.terminate();

  return 0;
}

void exahype::runners::Runner::initSolverTimeStamps() {
  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    p->initInitialTimeStamp(0.0);
  }
}

void exahype::runners::Runner::initFiniteVolumesSolverTimeStamps() {
  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    if (p->getType()==exahype::solvers::Solver::Type::FiniteVolumes) {
      p->initInitialTimeStamp(0.0);
    }
  }
}

void exahype::runners::Runner::printTimeStepInfo(int numberOfStepsRanSinceLastCall) {
  double currentMinTimeStamp    = std::numeric_limits<double>::max();
  double currentMinTimeStepSize = std::numeric_limits<double>::max();
  double nextMinTimeStepSize    = std::numeric_limits<double>::max();

  static int n = 0;
  if (numberOfStepsRanSinceLastCall==0) {
    n++;
  }
  else if (numberOfStepsRanSinceLastCall>0) {
    n+=numberOfStepsRanSinceLastCall;
  }

  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    currentMinTimeStamp =
        std::min(currentMinTimeStamp, p->getMinTimeStamp());
    currentMinTimeStepSize =
        std::min(currentMinTimeStepSize, p->getMinTimeStepSize());
    nextMinTimeStepSize =
        std::min(nextMinTimeStepSize, p->getNextMinTimeStepSize());
  }

  logInfo("startNewTimeStep(...)",
      "step " << n << "\tt_min          =" << currentMinTimeStamp);

  logInfo("startNewTimeStep(...)",
      "\tdt_min         =" << currentMinTimeStepSize);


  #if !defined(Parallel)
  // memory consumption on rank 0 would not make any sense
  logInfo("startNewTimeStep(...)",
      "\tmemoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");
  #ifdef Asserts
  if (exahype::solvers::ADERDGSolver::CompressionAccuracy>0.0) {
    DataHeap::getInstance().plotStatistics();
    peano::heap::PlainCharHeap::getInstance().plotStatistics();
  }
  #endif

  #endif

  logDebug("startNewTimeStep(...)",
      "\tnext dt_min    =" << nextMinTimeStepSize); // Only interesting for ADER-DG. Prints MAX_DOUBLE for finite volumes.
#if defined(Debug) || defined(Asserts)
  tarch::logging::CommandLineLogger::getInstance().closeOutputStreamAndReopenNewOne();
#endif
}


void exahype::runners::Runner::runOneTimeStampWithFusedAlgorithmicSteps(
    exahype::repositories::Repository& repository, int numberOfStepsToRun, bool exchangeBoundaryData) {
  /*
   * The adapter below performs the following steps:
   *
   * 1. Exchange the fluctuations using the predictor computed in the previous
   *sweep
   *    and the corrector time stemp size.
   * 2. Perform the corrector step using the corrector update and the corrector
   *time step size.
   *    This is a cell-local operation. Thus we immediately obtain the
   *cell-local current solution.
   * 3. Perform the predictor step using the cell-local current solution and the
   *predictor time step size.
   * 4. Compute the cell-local time step sizes
   */
  repository.getState().switchToADERDGTimeStepContext();
  if (numberOfStepsToRun==0) {
    repository.switchToPlotAndADERDGTimeStep();
    repository.iterate();
  } else {
    repository.switchToADERDGTimeStep();
    repository.iterate(numberOfStepsToRun,exchangeBoundaryData);
  }

  // reduction/broadcast barrier
}

void exahype::runners::Runner::recomputePredictorIfNecessary(
    exahype::repositories::Repository& repository) {
  // Must be evaluated before we start a new time step
  //  bool stabilityConditionWasHarmed = setStableTimeStepSizesIfStabilityConditionWasHarmed(factor);
  // Note that it is important to switch the time step sizes, i.e,
  // start a new time step, before we recompute the predictor.

  if (repository.getState().stabilityConditionOfOneSolverWasViolated()) {
    logInfo("startNewTimeStep(...)",
        "\t\t Space-time predictor must be recomputed.");

    repository.getState().switchToPredictionRerunContext();
    repository.switchToPredictionRerun();
    repository.iterate();
  }
}

void exahype::runners::Runner::runOneTimeStampWithThreeSeparateAlgorithmicSteps(
    exahype::repositories::Repository& repository, bool plot) {
  // Only one time step (predictor vs. corrector) is used in this case.
  repository.getState().switchToNeighbourDataMergingContext();
  repository.switchToNeighbourDataMerging();  // Riemann -> face2face
  repository.iterate();

  // no merging and sending
  // TODO(Dominic):
  // Perform the solution update;
  // mark cells as troubled if applicable, send out limiter data and min/max
  // repository.switchPlotAndSolutionUpdateContext();
  if (plot) {
    repository.switchToPlotAndSolutionUpdate();  // Face to cell + Inside cell
  } else {
    repository.switchToSolutionUpdate();  // Face to cell + Inside cell
  }
  repository.iterate();

//  // TODO(Dominic): Spread limiter status and recompute solution
//  if (repository.getState().recomputeADERDGSolution()) {
//    // Spread the limiter status, add FV subcells
//    repository.getState().switchToSpreadLimiterStatusContext();
//    repository.switchToSpreadLimiterStatus();
//    repository.iterate(2);
//
//    TODO(Dominic): Send out limiter status again at the end.
//    // Maybe plot here again the corrected solution
//    repository.getState().switchToRecomputeSolutionContext();
//    repository.switchToRecomputeSolution();
//    repository.iterate(1);
//  }

//  #if DIMENSIONS==2 && !defined(Parallel)
//  if (plot) {
//    repository.switchToPlotAugmentedAMRGrid();
//    repository.iterate();
//  }
//  #endif


//  repository.getState().switchToPreAMRContext(); // This must stay here.
//
//  // TODO(Dominic): Experimental dyn. AMR grid setup that only works
//  // with standard global time stepping at the moment.
//  //
//  // For the refinement criterion try for example a density based one.
//  //
//  // Uncomment the following lines for testing dyn. AMR.
//  createGrid(repository);
//  const int maxAdaptiveGridDepth = exahype::solvers::Solver::getMaxAdaptiveRefinementDepthOfAllSolvers();
//  logDebug("runOneTimeStampWithThreeSeparateAlgorithmicSteps(...)","maxAdaptiveGridDepth="<<maxAdaptiveGridDepth);
//  repository.iterate(maxAdaptiveGridDepth*2); // We need to two iterations for an erasing event.
//  repository.getState().switchToPostAMRContext();
//  repository.switchToPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation();
//  repository.iterate();

  // TODO(Dominic): Drop Limiter status
  repository.getState().switchToTimeStepSizeComputationContext();
  repository.switchToTimeStepSizeComputation();
  repository.iterate();


  repository.getState().switchToPredictionContext();
  repository.switchToPrediction();  // Cell onto faces
  repository.iterate();
}
