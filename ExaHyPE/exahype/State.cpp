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

#include "exahype/State.h"
#include "exahype/Cell.h"
#include "exahype/Vertex.h"

#include "peano/grid/Checkpoint.h"

#include "exahype/solvers/Solver.h"
#include "exahype/plotters/Plotter.h"

#include "tarch/parallel/NodePool.h"

#include <limits>

#include "exahype/mappings/MeshRefinement.h"

tarch::logging::Log exahype::State::_log("exahype::State");

int exahype::State::CurrentBatchIteration   = 0;
int exahype::State::NumberOfBatchIterations = 1;

exahype::State::State() : Base() {
  _stateData.setMaxRefinementLevelAllowed(3);
  // I want the code to lb more aggressively, so it should not wait more than
  Base::IterationsInBetweenRebalancing = 2;
}


exahype::State::State(const Base::PersistentState& argument) : Base(argument) {
  // do nothing
}

void exahype::State::setVerticalExchangeOfSolverDataRequired(bool state) {
  _stateData.setVerticalExchangeOfSolverDataRequired(state);
}

bool exahype::State::getVerticalExchangeOfSolverDataRequired() const {
  return _stateData.getVerticalExchangeOfSolverDataRequired();
}

void exahype::State::setAllSolversAttainedStableState(const bool state) {
  _stateData.setAllSolversAttainedStableState(state);
}

bool exahype::State::getAllSolversAttainedStableState() const {
  return _stateData.getAllSolversAttainedStableState();
}

void exahype::State::setMeshRefinementIsInRefiningMode(const bool state) {
  _stateData.setMeshRefinementIsInRefiningMode(state);
}

bool exahype::State::getMeshRefinementIsInRefiningMode() const {
  return _stateData.getMeshRefinementIsInRefiningMode();
}

void exahype::State::setStableIterationsInARow(const int value) {
  _stateData.setStableIterationsInARow(value);
}

int exahype::State::getStableIterationsInARow() const {
  return _stateData.getStableIterationsInARow();
}

void exahype::State::mergeWithMaster(const exahype::State& anotherState) {
  _stateData.setVerticalExchangeOfSolverDataRequired(
      _stateData.getVerticalExchangeOfSolverDataRequired() ||
      anotherState._stateData.getVerticalExchangeOfSolverDataRequired());
  _stateData.setAllSolversAttainedStableState(
      _stateData.getAllSolversAttainedStableState() &&
      anotherState._stateData.getAllSolversAttainedStableState());
}

void exahype::State::writeToCheckpoint(
    peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>& checkpoint) const {
  // do nothing
}

void exahype::State::readFromCheckpoint(
    const peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>& checkpoint) {
  // do nothing
}

void exahype::State::endedGridConstructionIteration(int finestGridLevelPossible) {
  const bool idleNodesLeft =
      tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()>0;
  const bool nodePoolHasGivenOutRankSizeLastQuery =
      tarch::parallel::NodePool::getInstance().hasGivenOutRankSizeLastQuery();

  #ifdef Debug
  std::cout <<  "!getHasChangedVertexOrCellState="            << !_stateData.getHasChangedVertexOrCellState() << std::endl;
  std::cout <<  "!getHasRefined="                             << !_stateData.getHasRefined() << std::endl;
  std::cout <<  "!getHasErased="                              << !_stateData.getHasErased()  << std::endl;
  std::cout <<  "!getHasTriggeredRefinementForNextIteration=" << !_stateData.getHasTriggeredRefinementForNextIteration() << std::endl;
  std::cout <<  "!getHasTriggeredEraseForNextIteration="      << !_stateData.getHasTriggeredEraseForNextIteration() << std::endl;
  #ifdef Parallel
  std::cout <<  "!getCouldNotEraseDueToDecompositionFlag=" << !_stateData.getCouldNotEraseDueToDecompositionFlag() << std::endl;
  #endif
  std::cout << "isGridStationary()=" << isGridStationary() << std::endl;
  std::cout << "getMaxRefinementLevelAllowed()=" << _stateData.getMaxRefinementLevelAllowed() << std::endl;
  #endif

  // No more nodes left. Start to enforce refinement
  if ( !idleNodesLeft
      && _stateData.getMaxRefinementLevelAllowed()>=0
      && !nodePoolHasGivenOutRankSizeLastQuery
	  && isGridBalanced() ) {
    _stateData.setMaxRefinementLevelAllowed(-1);
  }
  // Seems that max permitted level has exceeded max grid level. We may assume
  // that there are more MPI ranks than available trees.
  else if (isGridStationary()
      && _stateData.getMaxRefinementLevelAllowed()>finestGridLevelPossible
      && isGridBalanced()
      && !nodePoolHasGivenOutRankSizeLastQuery) {
    _stateData.setMaxRefinementLevelAllowed( -1 );
  }
  // Refinement is enforced. So we decrease counter. Once we underrun -2, grid
  // construction can terminate as all enforced refined instructions went
  // through.
  else if (_stateData.getMaxRefinementLevelAllowed()<=-1
      && !nodePoolHasGivenOutRankSizeLastQuery
      && isGridStationary()) {
    _stateData.setMaxRefinementLevelAllowed(
        _stateData.getMaxRefinementLevelAllowed()-1 );
  }
  // Nothing has changed in this grid iteration in the grid and we haven't
  // given out new workers. So increase the permitted maximum grid level by
  // one and give another try whether the grid adds more vertices.
  else if (
      (!nodePoolHasGivenOutRankSizeLastQuery)
	  // @todo TW/DEC
	  // we might want to roll this back to isGridStationary()
      && isGridBalanced()
//	  && isGridStationary()
      && (_stateData.getMaxRefinementLevelAllowed()>=0)
  ) {
	static int stationarySweeps = 0;
	stationarySweeps++;
//	if (stationarySweeps>=Base::IterationsInBetweenRebalancing) {
      _stateData.setMaxRefinementLevelAllowed(
        _stateData.getMaxRefinementLevelAllowed()+1);
      stationarySweeps=0;
//	}
      //IterationsInBetweenRebalancing
  }
}

exahype::State::RefinementAnswer exahype::State::mayRefine(bool isCreationalEvent, int level) const
{
#ifdef Parallel
  if (
      _stateData.getMaxRefinementLevelAllowed()<=-2
      &&
      isCreationalEvent
      &&
      !isInvolvedInJoinOrFork() // A Peano assertion was triggered
  ) {
    return RefinementAnswer::EnforceRefinement;
  }
  else if ( _stateData.getMaxRefinementLevelAllowed()<0 ) {
    return RefinementAnswer::Refine;
  }
  else if (
      _stateData.getMaxRefinementLevelAllowed()>level
      &&
      !isCreationalEvent
      &&
      mayForkDueToLoadBalancing()
  ) {
    return RefinementAnswer::Refine;
  }
  else {
    return RefinementAnswer::DontRefineYet;
  }
#else
  return RefinementAnswer::Refine;
#endif
}


bool exahype::State::continueToConstructGrid() {
  const int levelDelta =
      (exahype::solvers::Solver::getMaximumAdaptiveMeshLevelOfAllSolvers() -
      exahype::solvers::Solver::getCoarsestMeshLevelOfAllSolvers());
  static const int iterationsForErasingToConverge = 1;
      //3*levelDelta  +
      //std::max(exahype::solvers::Solver::getMaxRefinementStatus(),3);
  static const int iterationsForRefiningToConverge = 1;
      //2*levelDelta +
      //std::max(exahype::solvers::Solver::getMaxRefinementStatus(),3);

  // convergence analysis
  if ( getAllSolversAttainedStableState() ) {
    setStableIterationsInARow( getStableIterationsInARow()+1 );
    if (  getMeshRefinementIsInRefiningMode() &&
        getStableIterationsInARow() > iterationsForRefiningToConverge ) {
      setMeshRefinementIsInRefiningMode(false);
      setStableIterationsInARow(0);
    }
  } else {
    setStableIterationsInARow(0);
  }
  const bool meshRefinementHasConverged =
      isGridBalanced()                     &&
      !getMeshRefinementIsInRefiningMode() &&
      (mappings::MeshRefinement::IsInitialMeshRefinement ||
          getStableIterationsInARow() > iterationsForErasingToConverge);

  if (!meshRefinementHasConverged) {
    logInfo( "continueToConstructGrid(...)",
        "grid construction not yet finished. grid balanced=" << isGridBalanced() <<
        ", grid stationary=" << isGridStationary() <<
        ", still in refining mode=" << getMeshRefinementIsInRefiningMode() <<
        ", initial refinement=" << mappings::MeshRefinement::IsInitialMeshRefinement <<
        ", stable iterations in a row=" << getStableIterationsInARow() <<
        ", all solvers attained stable state=" << getAllSolversAttainedStableState() <<
        ", max level="<< getMaxLevel()
    );
  }
  return !meshRefinementHasConverged;
}

bool exahype::State::isEvenBatchIteration() {
  return CurrentBatchIteration % 2 == 0;
}

bool exahype::State::isFirstIterationOfBatchOrNoBatch() {
  return NumberOfBatchIterations==1 || CurrentBatchIteration==0;
}

bool exahype::State::isSecondIterationOfBatchOrNoBatch() {
  return NumberOfBatchIterations==1 || CurrentBatchIteration==1;
}

bool exahype::State::isLastIterationOfBatchOrNoBatch() {
  return NumberOfBatchIterations==1 || CurrentBatchIteration==NumberOfBatchIterations-1;
}

bool exahype::State::isSecondToLastIterationOfBatchOrNoBatch()  {
  return NumberOfBatchIterations==1 || CurrentBatchIteration==NumberOfBatchIterations-2;
}

void exahype::State::kickOffIteration(const exahype::records::RepositoryState::Action& action,const int currentBatchIteration,const int numberOfBatchIterations) {
  switch ( action ) {
  case exahype::records::RepositoryState::UseAdapterFinaliseMeshRefinement:
  case exahype::records::RepositoryState::UseAdapterFinaliseMeshRefinementOrLocalRollback:
    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      if ( solver->getMeshUpdateEvent()==exahype::solvers::Solver::MeshUpdateEvent::RefinementRequested ) {
        solver->rollbackToPreviousTimeStep();
      }
      if ( solver->getMeshUpdateEvent() != solvers::Solver::MeshUpdateEvent::None ) {
        solver->resetAdmissibleTimeStepSize(); // If a local re-computation is performed (IrregularLimiterDomainChange), we also
                                               // compute a new adm. time step size as the current minimum might stem from troubled cells
      }
    }
    break;
  case exahype::records::RepositoryState::UseAdapterInitialPrediction:
  case exahype::records::RepositoryState::UseAdapterPrediction:
    if ( currentBatchIteration==0 ) {
      for (auto* solver : exahype::solvers::RegisteredSolvers) {
        solver->kickOffTimeStep(currentBatchIteration==0);
      }
    }
    break;
  case exahype::records::RepositoryState::UseAdapterFusedTimeStep: {
    const bool beginFusedTimeStep = exahype::solvers::Solver::PredictionSweeps==1 || (currentBatchIteration % 2 == 0);
    if ( beginFusedTimeStep ) {
      for (auto* solver : exahype::solvers::RegisteredSolvers) {
        solver->kickOffTimeStep(currentBatchIteration==0);
      }
    }
  } break;
  default:
    break;
  }
}

bool exahype::State::startAndFinishSynchronousExchangeManually(const exahype::records::RepositoryState::Action& action,bool predictionFusedTimeStepCondition) {
  return
      ((action==exahype::records::RepositoryState::UseAdapterInitialPrediction            ||
      action==exahype::records::RepositoryState::UseAdapterPrediction                     ||
      action==exahype::records::RepositoryState::UseAdapterPredictionRerun                ||
      action==exahype::records::RepositoryState::UseAdapterPredictionOrLocalRecomputation ||
      action==exahype::records::RepositoryState::UseAdapterFusedTimeStep) &&
      predictionFusedTimeStepCondition)
      ||
      action==exahype::records::RepositoryState::UseAdapterBroadcast ||
      action==exahype::records::RepositoryState::UseAdapterBroadcastAndDropNeighbourMessages;
}

bool exahype::State::startAndFinishNeighbourExchangeManually(const exahype::records::RepositoryState::Action& action,bool predictionFusedTimeStepCondition) {
  return
      ((action==exahype::records::RepositoryState::UseAdapterInitialPrediction            ||
      action==exahype::records::RepositoryState::UseAdapterPrediction                     ||
      action==exahype::records::RepositoryState::UseAdapterPredictionRerun                ||
      action==exahype::records::RepositoryState::UseAdapterPredictionOrLocalRecomputation ||
      action==exahype::records::RepositoryState::UseAdapterFusedTimeStep) &&
      predictionFusedTimeStepCondition)
      ||
      action==exahype::records::RepositoryState::UseAdapterBroadcast ||
      action==exahype::records::RepositoryState::UseAdapterBroadcastAndDropNeighbourMessages;
}

void exahype::State::kickOffIteration(exahype::records::RepositoryState& repositoryState, exahype::State& solverState,const int currentBatchIteration) {
  CurrentBatchIteration   = currentBatchIteration;
  NumberOfBatchIterations = repositoryState.getNumberOfIterations();

  // the following must come after the global batch iteration variables are set
  if ( tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    kickOffIteration(repositoryState.getAction(),currentBatchIteration,repositoryState.getNumberOfIterations());
  }

  #ifdef Parallel
  const bool manualSychronousExchange = startAndFinishSynchronousExchangeManually(repositoryState.getAction(),CurrentBatchIteration == 0);
  const bool manualNeighbourExchange  = startAndFinishNeighbourExchangeManually  (repositoryState.getAction(),CurrentBatchIteration % 2 == 0);

  if ( manualNeighbourExchange ) {
    logDebug("kickOffIteration(...)","all heaps start to send boundary data (batch iteration="<<currentBatchIteration<<")");
    peano::heap::AbstractHeap::allHeapsStartToSendBoundaryData(solverState.isTraversalInverted());
  }
  if ( manualSychronousExchange ) {
    logDebug("kickOffIteration(...)","all heaps start to send synchronous data (batch iteration="<<currentBatchIteration<<")");
    peano::heap::AbstractHeap::allHeapsStartToSendSynchronousData();

    assertionEquals(tarch::parallel::Node::getGlobalMasterRank(),0);
    const int masterRank = tarch::parallel::Node::getInstance().getGlobalMasterRank();
    if ( tarch::parallel::Node::getInstance().isGlobalMaster() ) { // TODO might be scalability bottleneck for large rank numbers
      for (int workerRank=1; workerRank<tarch::parallel::Node::getInstance().getNumberOfNodes(); workerRank++) {
        if ( !(tarch::parallel::NodePool::getInstance().isIdleNode(workerRank)) ) {
          exahype::State::broadcastGlobalDataToWorker(workerRank,0.0,0);
        }
      }
    } else {
      exahype::State::mergeWithGlobalDataFromMaster(masterRank,0.0,0);
    }
  }

  // kick off on other ranks
  if ( !tarch::parallel::Node::getInstance().isGlobalMaster() ) {
    kickOffIteration(repositoryState.getAction(),currentBatchIteration,repositoryState.getNumberOfIterations());
  }
  #endif
}


void exahype::State::wrapUpIteration(exahype::records::RepositoryState& repositoryState, exahype::State& solverState, const int currentBatchIteration) {
  #ifdef Parallel
  const bool manualSychronousExchange = startAndFinishSynchronousExchangeManually(repositoryState.getAction(),CurrentBatchIteration == NumberOfBatchIterations-1);
  const bool manuaNeighbourExchange   = startAndFinishNeighbourExchangeManually  (repositoryState.getAction(),CurrentBatchIteration == NumberOfBatchIterations-1 || CurrentBatchIteration % 2 != 0);

  // TODO check if the order of the calls (boundary,synchronous) makes a difference

  if ( manuaNeighbourExchange ) {
    bool isTraversalInverted = solverState.isTraversalInverted(); // is not broadcasted anymore
    if ( CurrentBatchIteration % 2 != 0 ) {
      isTraversalInverted = !isTraversalInverted;
    }
    logDebug("wrapUpIteration(...)","all heaps finish to send boundary data (batch iteration="<<currentBatchIteration<<")");
    peano::heap::AbstractHeap::allHeapsFinishedToSendBoundaryData(isTraversalInverted);
  }
  if ( manualSychronousExchange ) {
    logDebug("wrapUpIteration(...)","all heaps finish to send synchronous data (batch iteration="<<currentBatchIteration<<")");
    peano::heap::AbstractHeap::allHeapsFinishedToSendSynchronousData();
  }
  #endif
}

#ifdef Parallel
void exahype::State::broadcastGlobalDataToWorker(
    const int                                   worker,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const int                                   level) {
  for (auto& solver : exahype::solvers::RegisteredSolvers) {
    solver->sendDataToWorker(worker,cellCentre,level);
  }
  for (auto& plotter : exahype::plotters::RegisteredPlotters) {
    plotter->sendDataToWorker(worker,cellCentre,level);
  }
}

void exahype::State::mergeWithGlobalDataFromMaster(
    const int                                   master,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const int                                   level) {
  for (auto& solver : exahype::solvers::RegisteredSolvers) {
    solver->mergeWithMasterData(master,cellCentre,level);
  }
  for (auto& plotter : exahype::plotters::RegisteredPlotters) {
    plotter->mergeWithMasterData(master,cellCentre,level);
  }
}

void exahype::State::reduceGlobalDataToMaster(
    const int                                   master,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const int                                   level) {
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    solver->sendDataToMaster(master,cellCentre,level);
  }
}

void exahype::State::mergeWithGlobalDataFromWorker(
    const int                                   worker,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const int                                   level) {
  for (auto& solver : exahype::solvers::RegisteredSolvers) {
    solver->mergeWithWorkerData(worker,cellCentre,level);
  }
}
#endif
