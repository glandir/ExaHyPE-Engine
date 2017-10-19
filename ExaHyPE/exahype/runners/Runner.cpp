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

#include <cmath>

#include "../../../Peano/mpibalancing/HotspotBalancing.h"

#include "exahype/repositories/Repository.h"
#include "exahype/repositories/RepositoryFactory.h"
#include "exahype/mappings/TimeStepSizeComputation.h"
#include "exahype/mappings/Sending.h"
#include "exahype/mappings/LoadBalancing.h"


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
#include "mpibalancing/SFCDiffusionNodePoolStrategy.h"
#endif
#include "exahype/plotters/Plotter.h"

#include "exahype/mappings/MeshRefinement.h"
#include "exahype/mappings/LimiterStatusSpreading.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "tarch/multicore/MulticoreDefinitions.h"


tarch::logging::Log exahype::runners::Runner::_log("exahype::runners::Runner");

exahype::runners::Runner::Runner(Parser& parser) :
    _parser(parser),
    _boundingBoxSize(0.0),
    _meshRefinements(0),
    _localRecomputations(0),
    _globalRecomputations(0),
    _predictorReruns(0) {}

exahype::runners::Runner::~Runner() {}

void exahype::runners::Runner::initDistributedMemoryConfiguration() {
  #ifdef Parallel
  const std::string configuration = _parser.getMPIConfiguration();

  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    //
    // Configure answering behaviour of global node pool
    // =================================================
    //
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
    else if (configuration.find( "sfc-diffusion" )!=std::string::npos ) {
      int ranksPerNode = static_cast<int>(exahype::Parser::getValueFromPropertyString(configuration,"ranks_per_node"));
      if (ranksPerNode<=0) {
        logError( "initDistributedMemoryConfiguration()", "please inform SFC balancing how many ranks per node you use through value \"ranks_per_node:XXX\". Read value " << ranksPerNode << " is invalid" );
        ranksPerNode = 1;
      }
      if (tarch::parallel::Node::getInstance().getNumberOfNodes() % ranksPerNode != 0) {
        logError( "initDistributedMemoryConfiguration()", "Value of \"ranks_per_node:XXX\" does not fit to total number of ranks. ExaHyPE requires homogeneous rank distribution" );
        ranksPerNode = 1;
      }
      int primaryRanksPerNode = static_cast<int>(exahype::Parser::getValueFromPropertyString(configuration,"primary_ranks_per_node"));
      if (primaryRanksPerNode<=0) {
        logError( "initDistributedMemoryConfiguration()", "please inform SFC balancing how many primary ranks per node you use through value \"primary_ranks_per_node:XXX\". Read value " << primaryRanksPerNode << " is invalid" );
        primaryRanksPerNode = 1;
      }
      if ( ranksPerNode<primaryRanksPerNode ) {
        logWarning( "initDistributedMemoryConfiguration()", "value \"ranks_per_node:XXX\" is smaller than primary_ranks_per_node. Reset to 1" );
        primaryRanksPerNode = 1;
      }
      tarch::parallel::NodePool::getInstance().setStrategy(
        new mpibalancing::SFCDiffusionNodePoolStrategy(ranksPerNode,primaryRanksPerNode)
      );
      logInfo("initDistributedMemoryConfiguration()", "load balancing relies on fair answering strategy with " << ranksPerNode << " rank(s) per node") ;
    }
    else {
      logError("initDistributedMemoryConfiguration()", "no valid load balancing answering strategy specified");
      _parser.invalidate();
    }
  }

  if ( configuration.find( "greedy-naive" )!=std::string::npos ) {
    exahype::mappings::LoadBalancing::setLoadBalancingAnalysis( exahype::mappings::LoadBalancing::LoadBalancingAnalysis::Greedy );
  }
  else if ( configuration.find( "greedy-regular" )!=std::string::npos ) {
    exahype::mappings::LoadBalancing::setLoadBalancingAnalysis( exahype::mappings::LoadBalancing::LoadBalancingAnalysis::GreedyWithRegularityAnalysis );
  }
  else if ( configuration.find( "hotspot" )!=std::string::npos ) {
    exahype::mappings::LoadBalancing::setLoadBalancingAnalysis( exahype::mappings::LoadBalancing::LoadBalancingAnalysis::Hotspot );
  }
  else {
    logError("initDistributedMemoryConfiguration()", "no valid load balancing analysis type specified");
    _parser.invalidate();
  }

  //
  // Configure answering behaviour of global node pool
  // =================================================
  //
  if (_parser.getMPILoadBalancingType()==Parser::MPILoadBalancingType::Static) {
    switch ( exahype::mappings::LoadBalancing::getLoadBalancingAnalysis() ) {
      case exahype::mappings::LoadBalancing::LoadBalancingAnalysis::Greedy:
        logInfo("initDistributedMemoryConfiguration()", "use greedy load balancing without joins");
        peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
            new peano::parallel::loadbalancing::OracleForOnePhaseWithGreedyPartitioning(false)
        );
        break;
      case exahype::mappings::LoadBalancing::LoadBalancingAnalysis::GreedyWithRegularityAnalysis:
        logInfo("initDistributedMemoryConfiguration()", "use greedy load balancing without joins (mpibalancing/GreedyBalancing)");
        peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
          new mpibalancing::GreedyBalancing(
            getCoarsestGridLevelForLoadBalancing(_boundingBoxSize),
            getFinestUniformGridLevelForLoadBalancing(_boundingBoxSize)+1 /*boundary regularity*/
          )
        );
        break;
      case exahype::mappings::LoadBalancing::LoadBalancingAnalysis::Hotspot:
        logInfo("initDistributedMemoryConfiguration()", "use global hotspot elimination without joins (mpibalancing/StaticBalancing)");
        peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
            new mpibalancing::HotspotBalancing(
                false,getFinestUniformGridLevelForLoadBalancing(_boundingBoxSize)+1 /*boundary regularity*/
          )
        );
        break;
    }
  }
  // Dynamic load balancing
  // ----------------------
  // @todo Missing
  else {
    logError("initDistributedMemoryConfiguration()", "only MPI static load balancing supported so far. ");
  }

  tarch::parallel::NodePool::getInstance().restart();

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

  tarch::parallel::NodePool::getInstance().waitForAllNodesToBecomeIdle();
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
         true, //   bool useMultithreading                  = true,
         0, //   int  grainSizeOfUserDefinedRegions      = 0,
         peano::datatraversal::autotuning::OracleForOnePhaseDummy::SplitVertexReadsOnRegularSubtree::Split,
         #ifdef SharedOMP // Pipelining does not pay off for OpenMP (yet)
         false, //  bool pipelineDescendProcessing          = false,
         false,  //   bool pipelineAscendProcessing           = false,
         #else
         true, //  bool pipelineDescendProcessing          = true,
         true, //   bool pipelineAscendProcessing           = true,
         #endif
         27, //   int  smallestProblemSizeForAscendDescend  = tarch::la::aPowI(DIMENSIONS,3*3*3*3/2),
         3, //   int  grainSizeForAscendDescend          = 3,
         1, //   int  smallestProblemSizeForEnterLeaveCell = tarch::la::aPowI(DIMENSIONS,9/2),
         1, //   int  grainSizeForEnterLeaveCell         = 2,
         1, //   int  smallestProblemSizeForTouchFirstLast = tarch::la::aPowI(DIMENSIONS,3*3*3*3+1),
         1, //   int  grainSizeForTouchFirstLast         = 64,
         1, //   int  smallestProblemSizeForSplitLoadStore = tarch::la::aPowI(DIMENSIONS,3*3*3),
         1  //   int  grainSizeForSplitLoadStore         = 8,
      )
    );
    break;
  case Parser::MulticoreOracleType::AutotuningWithRestartAndLearning:
    logInfo("initSharedMemoryConfiguration()",
        "use learning autotuning shared memory oracle and allow restarts");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
        new sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize(
          tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getNumberOfNodes()-1,
          true
        ));
    peano::datatraversal::autotuning::Oracle::getInstance().loadStatistics(
        _parser.getMulticorePropertiesFile());
    break;
  case Parser::MulticoreOracleType::AutotuningWithoutLearning:
    logInfo("initSharedMemoryConfiguration()",
        "use autotuning shared memory oracle configuration but disable machine learning algorithm");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
        new sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize(false,false));
    peano::datatraversal::autotuning::Oracle::getInstance().loadStatistics(
        _parser.getMulticorePropertiesFile());
    break;
  case Parser::MulticoreOracleType::AutotuningWithLearningButWithoutRestart:
    logInfo("initSharedMemoryConfiguration()",
        "use autotuning shared memory oracle but disable search restarts");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
        new sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize(
          tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getNumberOfNodes()-1,
          false
        ));
    peano::datatraversal::autotuning::Oracle::getInstance().loadStatistics(
        _parser.getMulticorePropertiesFile());
    break;
  case Parser::MulticoreOracleType::GrainSizeSampling:
    logInfo("initSharedMemoryConfiguration()",
        "use shared memory oracle sampling");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
        new sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling(
            64,
            true    // logarithmicDistribution
        ));
    peano::datatraversal::autotuning::Oracle::getInstance().loadStatistics(
        _parser.getMulticorePropertiesFile());
    break;
  }
  #endif
}


void exahype::runners::Runner::initDataCompression() {
  exahype::solvers::Solver::CompressionAccuracy = _parser.getDoubleCompressionFactor();

  if (exahype::solvers::Solver::CompressionAccuracy==0.0) {
    logInfo( "initDataCompression()", "switched off any data compression");
  }
  else {
    exahype::solvers::Solver::SpawnCompressionAsBackgroundThread = _parser.getSpawnDoubleCompressionAsBackgroundTask();
    logInfo( "initDataCompression()", "store all data with accuracy of " << exahype::solvers::Solver::CompressionAccuracy << ". Use background threads for data conversion=" << exahype::solvers::ADERDGSolver::SpawnCompressionAsBackgroundThread);
  }
}


void exahype::runners::Runner::shutdownSharedMemoryConfiguration() {
#ifdef SharedMemoryParallelisation
  switch (_parser.getMulticoreOracleType()) {
  case Parser::MulticoreOracleType::AutotuningWithoutLearning:
    break;
  case Parser::MulticoreOracleType::Dummy:
  case Parser::MulticoreOracleType::AutotuningWithRestartAndLearning:
  case Parser::MulticoreOracleType::AutotuningWithLearningButWithoutRestart:
  case Parser::MulticoreOracleType::GrainSizeSampling:
  #ifdef Parallel
    if (
      tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getNumberOfNodes()-1
      &&
      tarch::multicore::Core::getInstance().getNumberOfThreads()>1
    ) {
      logInfo("shutdownSharedMemoryConfiguration()",
          "wrote statistics into file " << _parser.getMulticorePropertiesFile()
          << ". Dump from all other ranks subpressed to avoid file races"
      );
      peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics(
          _parser.getMulticorePropertiesFile());
    }
  #else
    if ( tarch::multicore::Core::getInstance().getNumberOfThreads()>1 ) {
      logInfo("shutdownSharedMemoryConfiguration()",
        "wrote statistics into file "
        << _parser.getMulticorePropertiesFile());
      peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics(
        _parser.getMulticorePropertiesFile());
    }
  #endif
    break;
  }
#endif
}


int exahype::runners::Runner::getCoarsestGridLevelOfAllSolvers(
    tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize) const {
  double hMax = exahype::solvers::Solver::getCoarsestMaximumMeshSizeOfAllSolvers();

  const int peanoLevel = exahype::solvers::Solver::computeMeshLevel(hMax,boundingBoxSize[0]);

  logDebug( "getCoarsestGridLevelOfAllSolvers()", "regular grid depth of " << peanoLevel << " (1 means a single cell)");
  return peanoLevel;
}

int exahype::runners::Runner::getCoarsestGridLevelForLoadBalancing(
    tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize) const {
  return std::max( 3, getCoarsestGridLevelOfAllSolvers(boundingBoxSize) );
}

int exahype::runners::Runner::getFinestUniformGridLevelOfAllSolvers(
    tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize) const {
  double hMax = exahype::solvers::Solver::getFinestMaximumMeshSizeOfAllSolvers();

  const int peanoLevel = exahype::solvers::Solver::computeMeshLevel(hMax,boundingBoxSize[0]);

  logDebug( "getCoarsestGridLevelOfAllSolvers()", "regular grid depth of " << peanoLevel << " (1 means a single cell)");
  return peanoLevel;
}

int exahype::runners::Runner::getFinestUniformGridLevelForLoadBalancing(
    tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize) const {
  return std::max( 3, getFinestUniformGridLevelOfAllSolvers(boundingBoxSize) );
}

double
exahype::runners::Runner::determineCoarsestMeshSize(
    const tarch::la::Vector<DIMENSIONS, double>& boundingBoxSize) const {
  const double coarsestUserMeshSize = exahype::solvers::Solver::getCoarsestMaximumMeshSizeOfAllSolvers();
  const double maxBoundingBoxExtent = tarch::la::max(boundingBoxSize);

  const int coarsestMeshLevel =
      exahype::solvers::Solver::computeMeshLevel(coarsestUserMeshSize,maxBoundingBoxExtent);
  return maxBoundingBoxExtent / threePowI(coarsestMeshLevel-1);
}

tarch::la::Vector<DIMENSIONS, double>
exahype::runners::Runner::determineScaledDomainSize(
    const tarch::la::Vector<DIMENSIONS, double>& domainSize,
    const double meshSize) const {
  tarch::la::Vector<DIMENSIONS, double> scaledDomainSize =
      domainSize / meshSize;
  for(int i=0; i<DIMENSIONS; i++) {
    scaledDomainSize[i] = std::ceil(scaledDomainSize[i]) * meshSize;
  }
  return scaledDomainSize;
}

/**
 * @return Bounding box size. If we have a non-cubical domain,
 *         then the bounding box still is cubical and all of its entries are
 *         the biggest dimension along one coordinate axis.
 */
tarch::la::Vector<DIMENSIONS, double> exahype::runners::Runner::determineBoundingBoxSize(
    const tarch::la::Vector<DIMENSIONS, double>& domainSize) const {
    double longestH = tarch::la::max(domainSize);
    return tarch::la::Vector<DIMENSIONS, double>(longestH);
}

exahype::repositories::Repository* exahype::runners::Runner::createRepository() {
  // Geometry is static as we need it to survive the whole simulation time.
  _domainOffset = _parser.getOffset();
  _domainSize   = _parser.getDomainSize();
  _boundingBoxSize  = determineBoundingBoxSize(_domainSize);

  const int coarsestUserMeshLevel = getCoarsestGridLevelOfAllSolvers(_boundingBoxSize);
  int boundingBoxMeshLevel = coarsestUserMeshLevel;
  tarch::la::Vector<DIMENSIONS,double> boundingBoxOffset = _domainOffset;

  #ifdef Parallel
  if (exahype::State::VirtuallyExpandBoundingBox) {
    const double coarsestUserMeshSpacing =
        exahype::solvers::Solver::getCoarsestMaximumMeshSizeOfAllSolvers();
    const double maxDomainExtent = tarch::la::max(_domainSize);

    double boundingBoxScaling         = 0;
    double boundingBoxExtent          = 0;
    double boundingBoxMeshSpacing     = std::numeric_limits<double>::max();

    int level = coarsestUserMeshLevel; // level=1 means a single cell
    while (boundingBoxMeshSpacing > coarsestUserMeshSpacing) {
      const double boundingBoxMeshCells = std::pow(3,level-1);
      boundingBoxScaling                = boundingBoxMeshCells / ( boundingBoxMeshCells - 2 );
      boundingBoxExtent                 = boundingBoxScaling * maxDomainExtent;
      boundingBoxMeshSpacing            = boundingBoxExtent/boundingBoxMeshCells;
      level++;
    }
    level--; // decrement result since boundingBox was computed using level-1

    assertion6(boundingBoxScaling>=1.0,boundingBoxScaling,boundingBoxExtent,boundingBoxMeshSpacing,boundingBoxMeshLevel,coarsestUserMeshSpacing,maxDomainExtent);

    boundingBoxMeshLevel = level;
    _boundingBoxSize    *= boundingBoxScaling;
    boundingBoxOffset   -= boundingBoxMeshSpacing;
  }
  #endif

  const double coarsestUserMeshSize = exahype::solvers::Solver::getCoarsestMaximumMeshSizeOfAllSolvers();
  const double coarsestMeshSize     = determineCoarsestMeshSize(_boundingBoxSize);
  tarch::la::Vector<DIMENSIONS,double> scaledDomainSize =
      determineScaledDomainSize(_domainSize,coarsestMeshSize);
  if (!tarch::la::equals(_domainSize,scaledDomainSize)) {
    logInfo("createRepository(...)",
        "scale domain size artificially to " << scaledDomainSize << " from "
        << _domainSize << " since non-cubic domain was specified");
  }
  logInfo("createRepository(...)",
      "coarsest mesh size was chosen as " << coarsestMeshSize << " based on user's maximum mesh size "<<
      coarsestUserMeshSize << " and domain size " << scaledDomainSize);
  if (boundingBoxMeshLevel!=coarsestUserMeshLevel) {
    logInfo("createRepository(...)",
        "We will need to refine the grid " << boundingBoxMeshLevel-coarsestUserMeshLevel << " more time(s) than expected "
            " in order to satisfy user's maximum mesh size criterion while scaling the bounding box");
  }

  logInfo(
      "createRepository(...)",
      "summary: create computational domain at " << _domainOffset <<
      " of width/size " << scaledDomainSize <<
      ". bounding box has offset " << boundingBoxOffset <<
      " and size " << _boundingBoxSize <<
      ". grid regular up to level " << boundingBoxMeshLevel << " (1 means a single cell)");

  _domainSize = scaledDomainSize;

  static peano::geometry::Hexahedron geometry(
      _domainSize,
      _domainOffset);

  return exahype::repositories::RepositoryFactory::getInstance().createWithSTDStackImplementation(
      geometry,
      _boundingBoxSize,
      boundingBoxOffset
  );
}

void exahype::runners::Runner::initHeaps() {
  exahype::DataHeap::getInstance().setName("DataHeap");
  logInfo("initHeaps()","initialised DataHeap="<<exahype::DataHeap::getInstance().toString());
  exahype::solvers::ADERDGSolver::Heap::getInstance().setName("ADERDGCellDescriptionHeap");
  logInfo("initHeaps()","initialised ADERDGSolver::Heap="<<exahype::solvers::ADERDGSolver::Heap::getInstance().toString());
  exahype::solvers::FiniteVolumesSolver::Heap::getInstance().setName("FiniteVolumesCellDescriptionHeap");
  logInfo("initHeaps()","initialised FiniteVolumesSolver::Heap="<<exahype::solvers::FiniteVolumesSolver::Heap::getInstance().toString());
  #ifdef Parallel
  exahype::MetadataHeap::getInstance().setName("MetadataHeap");
  logInfo("initHeaps()","initialised MetadataHeap="<<exahype::MetadataHeap::getInstance().toString());
  #endif
}

void exahype::runners::Runner::initHPCEnvironment() {
  peano::performanceanalysis::Analysis::getInstance().enable(false);
}

int exahype::runners::Runner::run() {
  int result = 0;
  if ( _parser.isValid() ) {
    initHeaps();

    multiscalelinkedcell::HangingVertexBookkeeper::getInstance().disableInheritingOfCoarseGridIndices();

    exahype::State::FuseADERDGPhases                = _parser.getFuseAlgorithmicSteps();
    exahype::State::WeightForPredictionRerun        = _parser.getFuseAlgorithmicStepsFactor();

    exahype::State::EnableMasterWorkerCommunication = _parser.getMPIMasterWorkerCommunication();
    exahype::State::EnableNeighbourCommunication    = _parser.getMPINeighbourCommunication();
    #ifdef Parallel
    exahype::State::VirtuallyExpandBoundingBox =
        _parser.getMPIConfiguration().find( "virtually-expand-domain")!=std::string::npos;
    #endif

    auto* repository = createRepository();
    // must come after repository creation
    initSolvers();
    initDistributedMemoryConfiguration();
    initSharedMemoryConfiguration();
    initDataCompression();
    initHPCEnvironment();

    exahype::mappings::MeshRefinement::IsInitialMeshRefinement=true;
    #ifdef Parallel
    exahype::mappings::MeshRefinement::IsFirstIteration = false;
    exahype::mappings::LimiterStatusSpreading::IsFirstIteration = false;
    #endif

    if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      result = runAsMaster(*repository);
    }
    #ifdef Parallel
    else {
      result = runAsWorker(*repository);
    }
    #endif

    shutdownSharedMemoryConfiguration();
    shutdownDistributedMemoryConfiguration();

    delete repository;
  }
  else {
    logError( "run(...)", "do not run code as parser reported errors" );
    result = 1;
  }

  return result;
}

void exahype::runners::Runner::initSolvers() const {
  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    p->initSolver(0.0,_domainOffset,_domainSize,_boundingBoxSize);
  }
}

void exahype::runners::Runner::plotMeshSetupInfo(
    exahype::repositories::Repository& repository,
    const int meshSetupIterations) const {
  #if defined(TrackGridStatistics) && defined(Asserts)
  logInfo("createGrid()",
      "grid setup iteration #" << meshSetupIterations <<
      ", max-level=" << repository.getState().getMaxLevel() <<
      ", state=" << repository.getState().toString() <<
      ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
  );
  #elif defined(Asserts)
  logInfo("createGrid()",
      "grid setup iteration #" << meshSetupIterations <<
      ", state=" << repository.getState().toString() <<
      ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
  );
  #elif defined(TrackGridStatistics)
  logInfo("createGrid()",
      "grid setup iteration #" << meshSetupIterations <<
      ", max-level=" << repository.getState().getMaxLevel() <<
      ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
  );
  #else
  logInfo("createGrid()",
      "grid setup iteration #" << meshSetupIterations <<
      ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
  );
  #endif

  #ifdef Asserts
  logInfo("createGrid()",
           "grid setup iteration #" << meshSetupIterations <<
           ", run one more iteration=" <<  repository.getState().continueToConstructGrid() ||
                                            exahype::solvers::Solver::oneSolverHasNotAttainedStableState();
   );
  #endif

  #if !defined(Parallel)
  logInfo("createGrid(...)", "memoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");
  #else
  if (tarch::parallel::Node::getInstance().getNumberOfNodes()==1) {
    logInfo("createGrid(...)", "memoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");
  }
  #endif

  #ifdef Asserts
  if (exahype::solvers::ADERDGSolver::CompressionAccuracy>0.0) {
    DataHeap::getInstance().plotStatistics();
    peano::heap::PlainCharHeap::getInstance().plotStatistics();
  }
  #endif
}

bool exahype::runners::Runner::createMesh(exahype::repositories::Repository& repository) {
  bool meshUpdate = false;

  int meshSetupIterations = 0;
  repository.switchToMeshRefinement();

  while (
      repository.getState().continueToConstructGrid() ||
      exahype::solvers::Solver::oneSolverHasNotAttainedStableState()
  ) {
    repository.iterate();
    meshSetupIterations++;

    repository.getState().endedGridConstructionIteration( getFinestUniformGridLevelOfAllSolvers(_boundingBoxSize) );

    plotMeshSetupInfo(repository,meshSetupIterations);

    std::cout << "repository.getState().continueToConstructGrid()=" <<  repository.getState().continueToConstructGrid() << std::endl;
    std::cout << "exahype::solvers::Solver::oneSolverHasNotAttainedStableState()="<< exahype::solvers::Solver::oneSolverHasNotAttainedStableState() << std::endl;

    meshUpdate = true;
  }

  // a few extra iterations for the cell status flag spreading
  logInfo("createGrid()", "more status spreading.");
  int extraIterations =
      std::max (
          4, // two extra iteration to spread the helper and augmentation status, one to allocate memory
          exahype::solvers::LimitingADERDGSolver::getMaxMinimumHelperStatusForTroubledCell());
  while (
      extraIterations > 0
      || repository.getState().continueToConstructGrid()
      || exahype::solvers::Solver::oneSolverHasNotAttainedStableState() // Further mesh refinement is possible
  ) {
    meshUpdate |=
        repository.getState().continueToConstructGrid()
        || exahype::solvers::Solver::oneSolverHasNotAttainedStableState();

    repository.iterate();
    extraIterations--;
    meshSetupIterations++;

    repository.getState().endedGridConstructionIteration( getFinestUniformGridLevelOfAllSolvers(_boundingBoxSize) );

    plotMeshSetupInfo(repository,meshSetupIterations);
  }

  logInfo("createGrid(Repository)", "finished grid setup after " << meshSetupIterations << " iterations" );

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

  return meshUpdate;
}


int exahype::runners::Runner::runAsMaster(exahype::repositories::Repository& repository) {
  peano::utils::UserInterface::writeHeader();

  if (!exahype::solvers::RegisteredSolvers.empty()) {
    initialiseMesh(repository);

    logInfo( "runAsMaster(...)", "initialised all data and computed first time step size" );

    bool plot = exahype::plotters::startPlottingIfAPlotterIsActive(
        solvers::Solver::getMinSolverTimeStampOfAllSolvers());

    repository.getState().setAlgorithmSection(exahype::records::State::TimeStepping);
    repository.getState().switchToPredictionContext();
    if (plot) {
      repository.switchToPredictionAndPlot();
    } else {
      repository.switchToPrediction();
    }
    repository.iterate();
    logInfo("runAsMaster(...)","plotted initial solution (if specified) and computed first predictor");

    printTimeStepInfo(-1,repository);
    validateInitialSolverTimeStepData(exahype::State::fuseADERDGPhases());

    const double simulationEndTime = _parser.getSimulationEndTime();
    logDebug("runAsMaster(...)","min solver time stamp: "     << solvers::Solver::getMinSolverTimeStampOfAllSolvers());
    logDebug("runAsMaster(...)","min solver time step size: " << solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers());
    while ((solvers::Solver::getMinSolverTimeStampOfAllSolvers() < simulationEndTime) &&
        tarch::la::greater(solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers(), 0.0)) {
      bool plot = exahype::plotters::startPlottingIfAPlotterIsActive(
          solvers::Solver::getMinSolverTimeStampOfAllSolvers());

      if (exahype::State::fuseADERDGPhases()) {
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

        runOneTimeStepWithFusedAlgorithmicSteps(
          repository,
          numberOfStepsToRun,
          _parser.getExchangeBoundaryDataInBatchedTimeSteps() && repository.getState().isGridStationary()
        );
        printTimeStepInfo(numberOfStepsToRun,repository);
      } else {
        runOneTimeStepWithThreeSeparateAlgorithmicSteps(repository, plot);
      }

      #if  defined(SharedMemoryParallelisation) && defined(PerformanceAnalysis) && !defined(Parallel)
      if (sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::hasLearnedSinceLastQuery()) {
        static int dumpCounter = -1;
        dumpCounter++;
        peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics( _parser.getMulticorePropertiesFile() + "-dump-" + std::to_string(dumpCounter) );
      }
      #endif

      logDebug("runAsMaster(...)", "state=" << repository.getState().toString());
    }
    if ( tarch::la::equals(solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers(), 0.0)) {
      logWarning("runAsMaster(...)","Minimum solver time step size is zero (up to machine precision).");
    }

    repository.getState().switchToNeighbourDataDroppingContext();
    repository.switchToNeighbourDataMerging();
    repository.iterate(1,false);

    printStatistics();
    repository.logIterationStatistics(false);
  }

  repository.terminate();

  return 0;
}

void exahype::runners::Runner::updateStatistics() {
  _meshRefinements      += (!exahype::solvers::LimitingADERDGSolver::oneSolverRequestedGlobalRecomputation() &&
                           exahype::solvers::Solver::oneSolverRequestedMeshUpdate()) ? 1 : 0;
  _localRecomputations  +=  (exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalRecomputation()) ? 1 : 0;
  _globalRecomputations +=  (exahype::solvers::LimitingADERDGSolver::oneSolverRequestedGlobalRecomputation()) ? 1 : 0;
  _predictorReruns      +=  (exahype::solvers::Solver::stabilityConditionOfOneSolverWasViolated()) ? 1 : 0;
}

void exahype::runners::Runner::printStatistics() {
  logInfo("printStatistics(...)","number of mesh refinements      = "<<_meshRefinements);
  logInfo("printStatistics(...)","number of local recomputations  = "<<_localRecomputations);
  logInfo("printStatistics(...)","number of global recomputations = "<<_globalRecomputations);
  logInfo("printStatistics(...)","number of predictor reruns      = "<<_predictorReruns);
}

void exahype::runners::Runner::validateInitialSolverTimeStepData(const bool fuseADERDGPhases) const {
  #ifdef Asserts
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    assertionEquals(solver->getMinTimeStamp(),0.0);
    assertion1(std::isfinite(solver->getMinTimeStepSize()),solver->getMinTimeStepSize());
    assertion1(solver->getMinTimeStepSize()>0,solver->getMinTimeStepSize());

    switch(solver->getTimeStepping()) {
      case exahype::solvers::Solver::TimeStepping::Global:
        assertionEquals(solver->getMinNextTimeStepSize(),std::numeric_limits<double>::max());
        break;
      case exahype::solvers::Solver::TimeStepping::GlobalFixed:
        break;
    }
    switch (solver->getType()) {
      case exahype::solvers::Solver::Type::ADERDG: {
        auto* aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
        if (!exahype::State::fuseADERDGPhases()) {
          assertionEquals(aderdgSolver->getPreviousMinCorrectorTimeStepSize(),0.0); // TOOD(Dominic): Revision
        }
        assertion1(std::isfinite(aderdgSolver->getMinPredictorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(std::isfinite(aderdgSolver->getMinCorrectorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertionEquals(aderdgSolver->getMinCorrectorTimeStamp(),0.0);
        if (fuseADERDGPhases) {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),aderdgSolver->getMinPredictorTimeStepSize());
        } else {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),0.0);
        }
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(aderdgSolver->getMinNextPredictorTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }
      }
      break;
      case exahype::solvers::Solver::Type::LimitingADERDG: {
        // ADER-DG
        auto* aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
        if (!exahype::State::fuseADERDGPhases()) {
          assertionEquals(aderdgSolver->getPreviousMinCorrectorTimeStepSize(),0.0); // TODO(Dominic): Revision
        }
        assertion1(std::isfinite(aderdgSolver->getMinPredictorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(std::isfinite(aderdgSolver->getMinCorrectorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertionEquals(aderdgSolver->getMinCorrectorTimeStamp(),0.0);
        if (fuseADERDGPhases) {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),aderdgSolver->getMinPredictorTimeStepSize());
        } else {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),0.0);
        }
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(aderdgSolver->getMinNextPredictorTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }
      } break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        auto* finiteVolumesSolver = static_cast<exahype::solvers::FiniteVolumesSolver*>(solver);
        if (!exahype::State::fuseADERDGPhases()) {
          assertionEquals(finiteVolumesSolver->getPreviousMinTimeStepSize(),0.0);
        }
        break;
    }
  }
  #endif
}

void exahype::runners::Runner::initialiseMesh(exahype::repositories::Repository& repository) {
  // We refine here using the previous solution (which is valid)
  logInfo("initialiseMesh(...)","create initial grid");
  repository.getState().setAlgorithmSection(exahype::records::State::TimeStepping);

  repository.getState().switchToUpdateMeshContext();
  createMesh(repository);
  logInfo("initialiseMesh(...)","finalise mesh refinement and compute first time step size");
  repository.getState().switchToTimeStepSizeComputationContext();
  repository.switchToFinaliseMeshRefinementAndTimeStepSizeComputation();
  repository.iterate();
}

void exahype::runners::Runner::updateMeshAndSubdomains(
    exahype::repositories::Repository& repository, const bool fusedTimeStepping) {
  // 1. All solvers drop their MPI messages
  assertion(repository.getState().getAlgorithmSection()==exahype::records::State::AlgorithmSection::TimeStepping);

  if (fusedTimeStepping) {
    repository.getState().switchToNeighbourDataDroppingContext();
    repository.switchToNeighbourDataMerging();
    repository.iterate(1,false);
  }

  // 1. Only the solvers with irregular limiter domain change do the limiter status spreading.
  if (exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLimiterStatusSpreading()) {
    repository.getState().setAlgorithmSection(exahype::records::State::AlgorithmSection::LimiterStatusSpreading);
    logInfo("updateMeshAndSubdomains(...)","pre-spreading of limiter status");
    repository.switchToLimiterStatusSpreading();
    repository.iterate(
        exahype::solvers::LimitingADERDGSolver::getMaxMinimumHelperStatusForTroubledCell(),false);
  }
  if (exahype::solvers::LimitingADERDGSolver::oneSolverRequestedGlobalRecomputation()) {
    assertion(exahype::solvers::Solver::oneSolverRequestedMeshUpdate());
    assertion(exahype::solvers::LimitingADERDGSolver::oneSolverHasNotAttainedStableState());
    logInfo("updateMeshAndSubdomains(...)","global recomputation requested by at least one solver");
    repository.switchToGlobalRollback();
    repository.iterate(1,false);
  }

  // 2. Perform a grid update for those solvers that requested refinement
  if (exahype::solvers::Solver::oneSolverRequestedMeshUpdate()) {
    repository.getState().setAlgorithmSection(exahype::records::State::AlgorithmSection::MeshRefinement);
    logInfo("updateMeshAndSubdomains(...)","perform mesh refinement");
    repository.getState().switchToUpdateMeshContext(); // TODO(Dominic): Adjust context for MPI
    createMesh(repository);
  }

  // 3. Drop the MPI metadata for all solvers that requested grid refinement
  // Further reinitialse solvers that reported an irregular limiter domain change
  if (exahype::solvers::Solver::oneSolverRequestedMeshUpdate() ||
      exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalOrGlobalRecomputation()) {
    repository.getState().setAlgorithmSection(exahype::records::State::AlgorithmSection::MeshRefinementOrLocalOrGlobalRecomputation);

    // TODO(Dominic): Only overwrite the current time step size with the previous values if global recomputation
    // Rename Reinitialisation to LocalReinitialisation. Do only call reinitialiseSolversLocally
    logInfo("updateMeshAndSubdomains(...)","finalise mesh refinement (if applicable)");
    logInfo("updateMeshAndSubdomains(...)","reinitialise cells and send data to neigbours (if applicable)");
    repository.getState().switchToReinitialisationContext();
    repository.switchToFinaliseMeshRefinementAndReinitialisation();
    repository.iterate(1,false);

    // 4. Perform a local recomputation of the solution of the solvers that requested one.
    // Perform a time
    if (!exahype::solvers::Solver::oneSolverRequestedMeshUpdate() &&
        !exahype::solvers::LimitingADERDGSolver::oneSolverRequestedGlobalRecomputation()) {
      repository.getState().setAlgorithmSection(exahype::records::State::AlgorithmSection::LocalRecomputationAllSend);
    }

    // Do not advance the time stamp if global recomputation/mesh refinement
    // Advance time stamp if local recomputation
    logInfo("updateMeshAndSubdomains(...)","recompute solution locally (if applicable) and compute new time step size");
    repository.getState().switchToLocalRecomputationAndTimeStepSizeComputationFusedTimeSteppingContext();
    repository.switchToLocalRecomputationAndTimeStepSizeComputation(); // do not roll forward here if global recomp.; we want to stay at the old time step
    repository.iterate(1,false); // local recomputation: has now recomputed predictor in interface cells
  } // LocalRecomputation is done here

  assertion(!exahype::solvers::LimitingADERDGSolver::oneSolverRequestedGlobalRecomputation() ||
            exahype::solvers::Solver::oneSolverRequestedMeshUpdate());

  if (fusedTimeStepping &&
      exahype::solvers::Solver::oneSolverRequestedMeshUpdate()) {
    repository.getState().setAlgorithmSection(exahype::records::State::AlgorithmSection::MeshRefinementOrGlobalRecomputationAllSend);

    logInfo("updateMeshAndSubdomains(...)","recompute predictor globally and reinitialise fused time stepping");
    repository.getState().switchToPredictionContext();
    repository.switchToPrediction();
    repository.iterate(1,false); // At this stage all solvers that required a mesh update, have
                          // recomputed the predictor
    if (exahype::solvers::LimitingADERDGSolver::oneSolverRequestedGlobalRecomputation()) {
      logInfo("updateMeshAndSubdomains(...)","redo time step since mesh was not prepared");
    }
  } // MeshUpdate and GlobalRecomputation is done here
}

void exahype::runners::Runner::printTimeStepInfo(int numberOfStepsRanSinceLastCall, const exahype::repositories::Repository& repository) {
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
        std::min(nextMinTimeStepSize, p->getMinNextTimeStepSize());

    #if defined(Debug) || defined(Asserts)
    switch(p->getType()) {
      case exahype::solvers::Solver::Type::ADERDG:
        logInfo("startNewTimeStep(...)",
                "\tADER-DG correction: t_min         =" << static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinCorrectorTimeStamp());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG correction: dt_min         =" << static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinCorrectorTimeStepSize());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG prediction: t_min         =" << static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinPredictorTimeStamp());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG prediction: dt_min         =" << static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinPredictorTimeStepSize());
        break;
      case exahype::solvers::Solver::Type::LimitingADERDG:
        logInfo("startNewTimeStep(...)",
                 "\tADER-DG prev correction*:  t_min   =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getPreviousMinCorrectorTimeStamp());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG prev correction*:  dt_min  =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getPreviousMinCorrectorTimeStepSize());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG correction: t_min         =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getMinCorrectorTimeStamp());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG correction: dt_min        =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getMinCorrectorTimeStepSize());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG prediction: t_min         =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getMinPredictorTimeStamp());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG prediction: dt_min        =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getMinPredictorTimeStepSize());
        break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        break;
    }
    #endif
  }

  logInfo("startNewTimeStep(...)",
      "step " << n << "\tt_min          =" << currentMinTimeStamp);

  logInfo("startNewTimeStep(...)",
      "\tdt_min         =" << currentMinTimeStepSize);

  #if !defined(Parallel)
  // memory consumption on rank 0 would not make any sense
  logInfo("startNewTimeStep(...)",
      "\tmemoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");
  #ifdef TrackGridStatistics
  if (exahype::solvers::ADERDGSolver::CompressionAccuracy>0.0) {
    DataHeap::getInstance().plotStatistics();
    peano::heap::PlainCharHeap::getInstance().plotStatistics();

    logInfo(
      "startNewTimeStep(...)",
      "\tpiped-uncompressed-byes=" << exahype::solvers::ADERDGSolver::PipedUncompressedBytes
      << "\tpiped-compressed-byes=" << exahype::solvers::ADERDGSolver::PipedCompressedBytes
      << "\tcompression-rate=" << (exahype::solvers::ADERDGSolver::PipedCompressedBytes/exahype::solvers::ADERDGSolver::PipedUncompressedBytes)
    );
  }
  #endif

  #if defined(TrackGridStatistics)
  logInfo(
    "startNewTimeStep(...)",
    "\tinner cells/inner unrefined cells=" << repository.getState().getNumberOfInnerCells()
    << "/" << repository.getState().getNumberOfInnerLeafCells() );
  logInfo(
    "startNewTimeStep(...)",
    "\tinner max/min mesh width=" << repository.getState().getMaximumMeshWidth()
    << "/" << repository.getState().getMinimumMeshWidth()
    );
  logInfo(
    "startNewTimeStep(...)",
    "\tmax level=" << repository.getState().getMaxLevel()
    );
  #endif

  #endif

  if (solvers::Solver::getMinSolverTimeStampOfAllSolvers()>std::numeric_limits<double>::max()/100.0) {
    logError("runAsMaster(...)","quit simulation as solver seems to explode" );
    exit(-1);
  }

  #if defined(Debug) || defined(Asserts)
  tarch::logging::CommandLineLogger::getInstance().closeOutputStreamAndReopenNewOne();
  #endif
}

void exahype::runners::Runner::runOneTimeStepWithFusedAlgorithmicSteps(
    exahype::repositories::Repository& repository, int numberOfStepsToRun, bool exchangeBoundaryData) {

  if (numberOfStepsToRun==0) {
    logInfo("runOneTimeStepWithFusedAlgorithmicSteps(...)","plot");
  }
  else if (numberOfStepsToRun>1) {
    logInfo("runOneTimeStepWithFusedAlgorithmicSteps(...)","run "<<numberOfStepsToRun<< " iterations within one batch");
  }

  /*
   * The adapter below performs the following steps:
   *
   * 1. Exchange the fluctuations using the predictor computed in the previous
   *sweep
   *    and the corrector time step size.
   * 2. Perform the corrector step using the corrector update and the corrector
   *time step size.
   *    This is a cell-local operation. Thus we immediately obtain the
   *cell-local current solution.
   * 3. Perform the predictor step using the cell-local current solution and the
   *predictor time step size.
   * 4. Compute the cell-local time step sizes
   */
  repository.getState().setAlgorithmSection(exahype::records::State::AlgorithmSection::TimeStepping);
  repository.getState().switchToADERDGTimeStepContext();

  if (numberOfStepsToRun==0) {
    repository.switchToPlotAndFusedTimeStep();
    repository.iterate(1,false);
  } else {
    repository.switchToFusedTimeStep();
    repository.iterate(numberOfStepsToRun,false);
  }

  if (exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalRecomputation()) {
    logInfo("runOneTimeStepWithFusedAlgorithmicSteps(...)","local recomputation requested by at least one solver");
  }
  if (exahype::solvers::LimitingADERDGSolver::oneSolverRequestedGlobalRecomputation()) {
    assertion(exahype::solvers::Solver::oneSolverRequestedMeshUpdate());
    logInfo("runOneTimeStepWithFusedAlgorithmicSteps(...)","global recomputation requested by at least one solver");
  }
  if (exahype::solvers::Solver::oneSolverRequestedMeshUpdate()) {
    logInfo("runOneTimeStepWithFusedAlgorithmicSteps(...)","mesh update requested by at least one solver");
  }

  if (exahype::solvers::Solver::oneSolverRequestedMeshUpdate() ||
      exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalOrGlobalRecomputation()) {
    updateMeshAndSubdomains(repository,true);
  }

  if (exahype::solvers::Solver::stabilityConditionOfOneSolverWasViolated()) {
    logInfo("runOneTimeStepWithFusedAlgorithmicSteps(...)", "\t\t recompute space-time predictor");
    repository.getState().setAlgorithmSection(exahype::records::State::PredictionRerunAllSend);
    repository.getState().switchToPredictionRerunContext();
    repository.switchToPrediction();
    repository.iterate(1,false);
  }

  updateStatistics();
  // ---- reduction/broadcast barrier ----
}

void exahype::runners::Runner::runOneTimeStepWithThreeSeparateAlgorithmicSteps(
    exahype::repositories::Repository& repository, bool plot) {
  // Only one time step (predictor vs. corrector) is used in this case.
  repository.getState().setAlgorithmSection(exahype::records::State::AlgorithmSection::TimeStepping);

  repository.getState().switchToNeighbourDataMergingContext();
  repository.switchToNeighbourDataMerging();  // Riemann -> face2face
  repository.iterate(1,false);

  repository.getState().switchToTimeStepSizeComputationContext();
  repository.switchToSolutionUpdate();  // Face to cell + Inside cell
  repository.iterate(1,false);

  if (exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalRecomputation()) {
    logInfo("runOneTimeStepWithThreeSeparateAlgorithmicSteps(...)","local recomputation requested by at least one solver");
  }
  if (exahype::solvers::LimitingADERDGSolver::oneSolverRequestedGlobalRecomputation()) {
    assertion(exahype::solvers::Solver::oneSolverRequestedMeshUpdate());
    logInfo("runOneTimeStepWithThreeSeparateAlgorithmicSteps(...)","global recomputation requested by at least one solver");
  }
  if (exahype::solvers::Solver::oneSolverRequestedMeshUpdate()) {
    logInfo("runOneTimeStepWithThreeSeparateAlgorithmicSteps(...)","mesh update requested by at least one solver");
  }

  if (exahype::solvers::Solver::oneSolverRequestedMeshUpdate() ||
      exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalOrGlobalRecomputation()) {
    updateMeshAndSubdomains(repository,false);
  }

  printTimeStepInfo(1,repository);

  /*
   * Compute current first predictor based on current time step size.
   * Set current time step size as old time step size of next iteration.
   * Compute the current time step size of the next iteration.
   */
  repository.getState().setAlgorithmSection(exahype::records::State::AlgorithmSection::TimeStepping);
  repository.getState().switchToPredictionContext();
  if (plot) {
    repository.switchToPredictionAndPlot();
  } else {
    repository.switchToPrediction();   // Cell onto faces
  }
  repository.iterate(1,false);

  updateStatistics();
}

void exahype::runners::Runner::validateSolverTimeStepDataForThreeAlgorithmicPhases(const bool fuseADERDGPhases) const {
  #ifdef Asserts
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    assertionEquals(solver->getMinTimeStamp(),0.0);
    assertion1(std::isfinite(solver->getMinTimeStepSize()),solver->getMinTimeStepSize());
    assertion1(solver->getMinTimeStepSize()>0,solver->getMinTimeStepSize());

    switch(solver->getTimeStepping()) {
      case exahype::solvers::Solver::TimeStepping::Global:
        assertionEquals(solver->getMinNextTimeStepSize(),std::numeric_limits<double>::max());
        break;
      case exahype::solvers::Solver::TimeStepping::GlobalFixed:
        break;
    }
    switch (solver->getType()) {
      case exahype::solvers::Solver::Type::ADERDG: {
        auto* aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
        assertion1(std::isfinite(aderdgSolver->getMinPredictorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(std::isfinite(aderdgSolver->getMinCorrectorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(aderdgSolver->getMinCorrectorTimeStamp() > 0.0,aderdgSolver->getMinCorrectorTimeStamp());
        if (!fuseADERDGPhases) {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),aderdgSolver->getMinCorrectorTimeStamp());
        }
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(aderdgSolver->getMinNextPredictorTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }

      }  break;
      case exahype::solvers::Solver::Type::LimitingADERDG: {
        // ADER-DG
        auto* aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
        assertion1(std::isfinite(aderdgSolver->getMinPredictorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(std::isfinite(aderdgSolver->getMinCorrectorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(aderdgSolver->getMinCorrectorTimeStamp() > 0.0,aderdgSolver->getMinCorrectorTimeStamp());
        if (!fuseADERDGPhases) {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),aderdgSolver->getMinCorrectorTimeStamp());
        }
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(aderdgSolver->getMinNextPredictorTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }

        // Finite Volumes
        auto* finiteVolumesSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiter().get();
        assertionEquals(finiteVolumesSolver->getMinTimeStamp(),0.0);
        assertion1(std::isfinite(finiteVolumesSolver->getMinTimeStepSize()),finiteVolumesSolver->getMinTimeStepSize());
        assertion1(finiteVolumesSolver->getMinTimeStepSize()>0,finiteVolumesSolver->getMinTimeStepSize());
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(finiteVolumesSolver->getMinNextTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }
      } break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        break;
    }
  }
  #endif
}
