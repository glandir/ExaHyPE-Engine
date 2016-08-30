#include "exahype/repositories/RepositoryArrayStack.h"

#include "tarch/Assertions.h"
#include "tarch/timing/Watch.h"

#include "tarch/compiler/CompilerSpecificSettings.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"

#ifdef Parallel
#include "peano/parallel/SendReceiveBufferPool.h"
#include "peano/parallel/loadbalancing/Oracle.h"
#endif

#include "peano/datatraversal/autotuning/Oracle.h"

#include "tarch/compiler/CompilerSpecificSettings.h"

#if !defined(CompilerICC)
#include "peano/grid/Grid.cpph"
#endif


tarch::logging::Log exahype::repositories::RepositoryArrayStack::_log( "exahype::repositories::RepositoryArrayStack" );


exahype::repositories::RepositoryArrayStack::RepositoryArrayStack(
  peano::geometry::Geometry&                   geometry,
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
  int                                          maximumSizeOfCellInOutStack,
  int                                          maximumSizeOfVertexInOutStack,
  int                                          maximumSizeOfVertexTemporaryStack
):
  _geometry(geometry),
  _cellStack(maximumSizeOfCellInOutStack),
  _vertexStack(maximumSizeOfVertexInOutStack, maximumSizeOfVertexTemporaryStack),
  _solverState(),
  _gridWithAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialConditionAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndPlotAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGridErasing(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStepAndPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionRerun(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithRiemannSolver(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPrediction(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrection(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectionAndPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositoryArrayStack(...)" );
  
  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );

  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositoryArrayStack(...)" );
}



exahype::repositories::RepositoryArrayStack::RepositoryArrayStack(
  peano::geometry::Geometry&                   geometry,
  int                                          maximumSizeOfCellInOutStack,
  int                                          maximumSizeOfVertexInOutStack,
  int                                          maximumSizeOfVertexTemporaryStack
):
  _geometry(geometry),
  _cellStack(maximumSizeOfCellInOutStack),
  _vertexStack(maximumSizeOfVertexInOutStack,maximumSizeOfVertexTemporaryStack),
  _solverState(),
  _gridWithAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialConditionAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndPlotAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGridErasing(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStepAndPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionRerun(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithRiemannSolver(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPrediction(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrection(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectionAndPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositoryArrayStack(Geometry&)" );
  
  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );

  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositoryArrayStack(Geometry&)" );
}
    
   
exahype::repositories::RepositoryArrayStack::~RepositoryArrayStack() {
  assertion( _repositoryState.getAction() == exahype::records::RepositoryState::Terminate );
}


void exahype::repositories::RepositoryArrayStack::restart(
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
  int                                          domainLevel,
  const tarch::la::Vector<DIMENSIONS,int>&     positionOfCentralElementWithRespectToCoarserRemoteLevel
) {
  logTraceInWith4Arguments( "restart(...)", domainSize, domainOffset, domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel );
  #ifdef Parallel
  assertion( !tarch::parallel::Node::getInstance().isGlobalMaster());
  #endif
  
  logInfo( "restart(...)", "start node for subdomain " << domainOffset << "x" << domainSize << " on level " << domainLevel << " with master " << tarch::parallel::NodePool::getInstance().getMasterRank() );
  
  assertion( _repositoryState.getAction() == exahype::records::RepositoryState::Terminate );

  _vertexStack.clear();
  _cellStack.clear();

  _gridWithAugmentedAMRGrid.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPlotAugmentedAMRGrid.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithInitialConditionAndTimeStepSizeComputation.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionAndPlotAndTimeStepSizeComputation.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionAndTimeStepSizeComputation.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithGridErasing.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithADERDGTimeStep.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithADERDGTimeStepAndPlot.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionRerun.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithRiemannSolver.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPrediction.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrection.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrectionAndPlot.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPlot.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);

 
   _solverState.restart();
 
  logTraceOut( "restart(...)" );
}


void exahype::repositories::RepositoryArrayStack::terminate() {
  logTraceIn( "terminate()" );
  
  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );
  
  #ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
      _repositoryState,
      peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
    );
  }
  peano::parallel::SendReceiveBufferPool::getInstance().terminate();
  #endif
  
  _gridWithAugmentedAMRGrid.terminate();
  _gridWithPlotAugmentedAMRGrid.terminate();
  _gridWithInitialConditionAndTimeStepSizeComputation.terminate();
  _gridWithPredictionAndPlotAndTimeStepSizeComputation.terminate();
  _gridWithPredictionAndTimeStepSizeComputation.terminate();
  _gridWithGridErasing.terminate();
  _gridWithADERDGTimeStep.terminate();
  _gridWithADERDGTimeStepAndPlot.terminate();
  _gridWithPredictionRerun.terminate();
  _gridWithRiemannSolver.terminate();
  _gridWithPrediction.terminate();
  _gridWithCorrection.terminate();
  _gridWithCorrectionAndPlot.terminate();
  _gridWithPlot.terminate();

 
  logTraceOut( "terminate()" );
}


exahype::State& exahype::repositories::RepositoryArrayStack::getState() {
  return _solverState;
}


const exahype::State& exahype::repositories::RepositoryArrayStack::getState() const {
  return _solverState;
}

   
void exahype::repositories::RepositoryArrayStack::iterate(int numberOfIterations, bool exchangeBoundaryVertices) {
  tarch::timing::Watch watch( "exahype::repositories::RepositoryArrayStack", "iterate(bool)", false);
  
  #ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    _repositoryState.setNumberOfIterations(numberOfIterations);
    _repositoryState.setExchangeBoundaryVertices(exchangeBoundaryVertices);
    tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
      _repositoryState,
      peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
    );
  }
  else {
    assertionEquals( numberOfIterations, 1 );
    numberOfIterations = _repositoryState.getNumberOfIterations();
  }

  peano::parallel::SendReceiveBufferPool::getInstance().exchangeBoundaryVertices(_repositoryState.getExchangeBoundaryVertices());

  if ( numberOfIterations > 1 && ( peano::parallel::loadbalancing::Oracle::getInstance().isLoadBalancingActivated() || _solverState.isInvolvedInJoinOrFork() )) {
    logWarning( "iterate()", "iterate invoked for multiple traversals though load balancing is switched on or grid is not balanced globally. Use activateLoadBalancing(false) to deactivate the load balancing before" );
  }

  peano::datatraversal::autotuning::Oracle::getInstance().switchToOracle(_repositoryState.getAction());

  peano::parallel::loadbalancing::Oracle::getInstance().switchToOracle(_repositoryState.getAction());
  
  _solverState.currentlyRunsMultipleIterations(_repositoryState.getNumberOfIterations()>1);
  #else
  peano::datatraversal::autotuning::Oracle::getInstance().switchToOracle(_repositoryState.getAction());
  #endif
  
  for (int i=0; i<numberOfIterations; i++) {
    switch ( _repositoryState.getAction()) {
      case exahype::records::RepositoryState::UseAdapterAugmentedAMRGrid: watch.startTimer(); _gridWithAugmentedAMRGrid.iterate(); watch.stopTimer(); _measureAugmentedAMRGridCPUTime.setValue( watch.getCPUTime() ); _measureAugmentedAMRGridCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid: watch.startTimer(); _gridWithPlotAugmentedAMRGrid.iterate(); watch.stopTimer(); _measurePlotAugmentedAMRGridCPUTime.setValue( watch.getCPUTime() ); _measurePlotAugmentedAMRGridCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterInitialConditionAndTimeStepSizeComputation: watch.startTimer(); _gridWithInitialConditionAndTimeStepSizeComputation.iterate(); watch.stopTimer(); _measureInitialConditionAndTimeStepSizeComputationCPUTime.setValue( watch.getCPUTime() ); _measureInitialConditionAndTimeStepSizeComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionAndPlotAndTimeStepSizeComputation: watch.startTimer(); _gridWithPredictionAndPlotAndTimeStepSizeComputation.iterate(); watch.stopTimer(); _measurePredictionAndPlotAndTimeStepSizeComputationCPUTime.setValue( watch.getCPUTime() ); _measurePredictionAndPlotAndTimeStepSizeComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionAndTimeStepSizeComputation: watch.startTimer(); _gridWithPredictionAndTimeStepSizeComputation.iterate(); watch.stopTimer(); _measurePredictionAndTimeStepSizeComputationCPUTime.setValue( watch.getCPUTime() ); _measurePredictionAndTimeStepSizeComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterGridErasing: watch.startTimer(); _gridWithGridErasing.iterate(); watch.stopTimer(); _measureGridErasingCPUTime.setValue( watch.getCPUTime() ); _measureGridErasingCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterADERDGTimeStep: watch.startTimer(); _gridWithADERDGTimeStep.iterate(); watch.stopTimer(); _measureADERDGTimeStepCPUTime.setValue( watch.getCPUTime() ); _measureADERDGTimeStepCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterADERDGTimeStepAndPlot: watch.startTimer(); _gridWithADERDGTimeStepAndPlot.iterate(); watch.stopTimer(); _measureADERDGTimeStepAndPlotCPUTime.setValue( watch.getCPUTime() ); _measureADERDGTimeStepAndPlotCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionRerun: watch.startTimer(); _gridWithPredictionRerun.iterate(); watch.stopTimer(); _measurePredictionRerunCPUTime.setValue( watch.getCPUTime() ); _measurePredictionRerunCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterRiemannSolver: watch.startTimer(); _gridWithRiemannSolver.iterate(); watch.stopTimer(); _measureRiemannSolverCPUTime.setValue( watch.getCPUTime() ); _measureRiemannSolverCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPrediction: watch.startTimer(); _gridWithPrediction.iterate(); watch.stopTimer(); _measurePredictionCPUTime.setValue( watch.getCPUTime() ); _measurePredictionCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrection: watch.startTimer(); _gridWithCorrection.iterate(); watch.stopTimer(); _measureCorrectionCPUTime.setValue( watch.getCPUTime() ); _measureCorrectionCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrectionAndPlot: watch.startTimer(); _gridWithCorrectionAndPlot.iterate(); watch.stopTimer(); _measureCorrectionAndPlotCPUTime.setValue( watch.getCPUTime() ); _measureCorrectionAndPlotCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPlot: watch.startTimer(); _gridWithPlot.iterate(); watch.stopTimer(); _measurePlotCPUTime.setValue( watch.getCPUTime() ); _measurePlotCalendarTime.setValue( watch.getCalendarTime() ); break;

      case exahype::records::RepositoryState::Terminate:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case exahype::records::RepositoryState::NumberOfAdapters:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case exahype::records::RepositoryState::RunOnAllNodes:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case exahype::records::RepositoryState::ReadCheckpoint:
        assertionMsg( false, "not implemented yet" );
        break;
      case exahype::records::RepositoryState::WriteCheckpoint:
        assertionMsg( false, "not implemented yet" );
        break;
    }
  }
    
  #ifdef Parallel
  if (_solverState.isJoiningWithMaster()) {
    _repositoryState.setAction( exahype::records::RepositoryState::Terminate );
  }
  #endif
}

 void exahype::repositories::RepositoryArrayStack::switchToAugmentedAMRGrid() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterAugmentedAMRGrid); }
 void exahype::repositories::RepositoryArrayStack::switchToPlotAugmentedAMRGrid() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid); }
 void exahype::repositories::RepositoryArrayStack::switchToInitialConditionAndTimeStepSizeComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterInitialConditionAndTimeStepSizeComputation); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictionAndPlotAndTimeStepSizeComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionAndPlotAndTimeStepSizeComputation); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictionAndTimeStepSizeComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionAndTimeStepSizeComputation); }
 void exahype::repositories::RepositoryArrayStack::switchToGridErasing() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterGridErasing); }
 void exahype::repositories::RepositoryArrayStack::switchToADERDGTimeStep() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterADERDGTimeStep); }
 void exahype::repositories::RepositoryArrayStack::switchToADERDGTimeStepAndPlot() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterADERDGTimeStepAndPlot); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictionRerun() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionRerun); }
 void exahype::repositories::RepositoryArrayStack::switchToRiemannSolver() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterRiemannSolver); }
 void exahype::repositories::RepositoryArrayStack::switchToPrediction() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPrediction); }
 void exahype::repositories::RepositoryArrayStack::switchToCorrection() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrection); }
 void exahype::repositories::RepositoryArrayStack::switchToCorrectionAndPlot() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrectionAndPlot); }
 void exahype::repositories::RepositoryArrayStack::switchToPlot() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlot); }



 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterAugmentedAMRGrid() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterAugmentedAMRGrid; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPlotAugmentedAMRGrid() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterInitialConditionAndTimeStepSizeComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterInitialConditionAndTimeStepSizeComputation; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictionAndPlotAndTimeStepSizeComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionAndPlotAndTimeStepSizeComputation; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictionAndTimeStepSizeComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionAndTimeStepSizeComputation; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterGridErasing() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterGridErasing; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterADERDGTimeStep() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterADERDGTimeStep; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterADERDGTimeStepAndPlot() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterADERDGTimeStepAndPlot; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictionRerun() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionRerun; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterRiemannSolver() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterRiemannSolver; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPrediction() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPrediction; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterCorrection() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrection; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterCorrectionAndPlot() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrectionAndPlot; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPlot() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlot; }



peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>* exahype::repositories::RepositoryArrayStack::createEmptyCheckpoint() {
  return new peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>();
} 


void exahype::repositories::RepositoryArrayStack::writeCheckpoint(peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> * const checkpoint) {
  _solverState.writeToCheckpoint( *checkpoint );
  _vertexStack.writeToCheckpoint( *checkpoint );
  _cellStack.writeToCheckpoint( *checkpoint );
} 


void exahype::repositories::RepositoryArrayStack::setMaximumMemoryFootprintForTemporaryRegularGrids(double value) {
  _regularGridContainer.setMaximumMemoryFootprintForTemporaryRegularGrids(value);
}


void exahype::repositories::RepositoryArrayStack::readCheckpoint( peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> const * const checkpoint ) {
  assertionMsg( checkpoint->isValid(), "checkpoint has to be valid if you call this operation" );

  _solverState.readFromCheckpoint( *checkpoint );
  _vertexStack.readFromCheckpoint( *checkpoint );
  _cellStack.readFromCheckpoint( *checkpoint );
}


#ifdef Parallel
void exahype::repositories::RepositoryArrayStack::runGlobalStep() {
  assertion(tarch::parallel::Node::getInstance().isGlobalMaster());

  exahype::records::RepositoryState intermediateStateForWorkingNodes;
  intermediateStateForWorkingNodes.setAction( exahype::records::RepositoryState::RunOnAllNodes );
  
  tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
    intermediateStateForWorkingNodes,
    peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
  );
  tarch::parallel::NodePool::getInstance().activateIdleNodes();
}


exahype::repositories::RepositoryArrayStack::ContinueCommand exahype::repositories::RepositoryArrayStack::continueToIterate() {
  logTraceIn( "continueToIterate()" );

  assertion( !tarch::parallel::Node::getInstance().isGlobalMaster());

  ContinueCommand result;
  if ( _solverState.hasJoinedWithMaster() ) {
    result = Terminate;
  }
  else {
    int masterNode = tarch::parallel::Node::getInstance().getGlobalMasterRank();
    assertion( masterNode != -1 );

    _repositoryState.receive( masterNode, peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag(), true, ReceiveIterationControlMessagesBlocking );

    result = Continue;
    if (_repositoryState.getAction()==exahype::records::RepositoryState::Terminate) {
      result = Terminate;
    } 
    if (_repositoryState.getAction()==exahype::records::RepositoryState::RunOnAllNodes) {
      result = RunGlobalStep;
    } 
  }
   
  logTraceOutWith1Argument( "continueToIterate()", result );
  return result;
}
#endif


void exahype::repositories::RepositoryArrayStack::logIterationStatistics(bool logAllAdapters) const {
  logInfo( "logIterationStatistics()", "|| adapter name \t || iterations \t || total CPU time [t]=s \t || average CPU time [t]=s \t || total user time [t]=s \t || average user time [t]=s  || CPU time properties  || user time properties " );  
   if (logAllAdapters || _measureAugmentedAMRGridCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| AugmentedAMRGrid \t |  " << _measureAugmentedAMRGridCPUTime.getNumberOfMeasurements() << " \t |  " << _measureAugmentedAMRGridCPUTime.getAccumulatedValue() << " \t |  " << _measureAugmentedAMRGridCPUTime.getValue()  << " \t |  " << _measureAugmentedAMRGridCalendarTime.getAccumulatedValue() << " \t |  " << _measureAugmentedAMRGridCalendarTime.getValue() << " \t |  " << _measureAugmentedAMRGridCPUTime.toString() << " \t |  " << _measureAugmentedAMRGridCalendarTime.toString() );
   if (logAllAdapters || _measurePlotAugmentedAMRGridCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PlotAugmentedAMRGrid \t |  " << _measurePlotAugmentedAMRGridCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.getValue()  << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.getValue() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.toString() << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.toString() );
   if (logAllAdapters || _measureInitialConditionAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| InitialConditionAndTimeStepSizeComputation \t |  " << _measureInitialConditionAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCPUTime.getValue()  << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCalendarTime.getValue() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCPUTime.toString() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionAndPlotAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionAndPlotAndTimeStepSizeComputation \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputationCPUTime.getValue()  << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputationCalendarTime.getValue() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputationCPUTime.toString() << " \t |  " << _measurePredictionAndPlotAndTimeStepSizeComputationCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionAndTimeStepSizeComputation \t |  " << _measurePredictionAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionAndTimeStepSizeComputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndTimeStepSizeComputationCPUTime.getValue()  << " \t |  " << _measurePredictionAndTimeStepSizeComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndTimeStepSizeComputationCalendarTime.getValue() << " \t |  " << _measurePredictionAndTimeStepSizeComputationCPUTime.toString() << " \t |  " << _measurePredictionAndTimeStepSizeComputationCalendarTime.toString() );
   if (logAllAdapters || _measureGridErasingCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| GridErasing \t |  " << _measureGridErasingCPUTime.getNumberOfMeasurements() << " \t |  " << _measureGridErasingCPUTime.getAccumulatedValue() << " \t |  " << _measureGridErasingCPUTime.getValue()  << " \t |  " << _measureGridErasingCalendarTime.getAccumulatedValue() << " \t |  " << _measureGridErasingCalendarTime.getValue() << " \t |  " << _measureGridErasingCPUTime.toString() << " \t |  " << _measureGridErasingCalendarTime.toString() );
   if (logAllAdapters || _measureADERDGTimeStepCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| ADERDGTimeStep \t |  " << _measureADERDGTimeStepCPUTime.getNumberOfMeasurements() << " \t |  " << _measureADERDGTimeStepCPUTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepCPUTime.getValue()  << " \t |  " << _measureADERDGTimeStepCalendarTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepCalendarTime.getValue() << " \t |  " << _measureADERDGTimeStepCPUTime.toString() << " \t |  " << _measureADERDGTimeStepCalendarTime.toString() );
   if (logAllAdapters || _measureADERDGTimeStepAndPlotCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| ADERDGTimeStepAndPlot \t |  " << _measureADERDGTimeStepAndPlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measureADERDGTimeStepAndPlotCPUTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepAndPlotCPUTime.getValue()  << " \t |  " << _measureADERDGTimeStepAndPlotCalendarTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepAndPlotCalendarTime.getValue() << " \t |  " << _measureADERDGTimeStepAndPlotCPUTime.toString() << " \t |  " << _measureADERDGTimeStepAndPlotCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionRerunCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionRerun \t |  " << _measurePredictionRerunCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionRerunCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionRerunCPUTime.getValue()  << " \t |  " << _measurePredictionRerunCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionRerunCalendarTime.getValue() << " \t |  " << _measurePredictionRerunCPUTime.toString() << " \t |  " << _measurePredictionRerunCalendarTime.toString() );
   if (logAllAdapters || _measureRiemannSolverCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| RiemannSolver \t |  " << _measureRiemannSolverCPUTime.getNumberOfMeasurements() << " \t |  " << _measureRiemannSolverCPUTime.getAccumulatedValue() << " \t |  " << _measureRiemannSolverCPUTime.getValue()  << " \t |  " << _measureRiemannSolverCalendarTime.getAccumulatedValue() << " \t |  " << _measureRiemannSolverCalendarTime.getValue() << " \t |  " << _measureRiemannSolverCPUTime.toString() << " \t |  " << _measureRiemannSolverCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Prediction \t |  " << _measurePredictionCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionCPUTime.getValue()  << " \t |  " << _measurePredictionCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionCalendarTime.getValue() << " \t |  " << _measurePredictionCPUTime.toString() << " \t |  " << _measurePredictionCalendarTime.toString() );
   if (logAllAdapters || _measureCorrectionCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Correction \t |  " << _measureCorrectionCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectionCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectionCPUTime.getValue()  << " \t |  " << _measureCorrectionCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectionCalendarTime.getValue() << " \t |  " << _measureCorrectionCPUTime.toString() << " \t |  " << _measureCorrectionCalendarTime.toString() );
   if (logAllAdapters || _measureCorrectionAndPlotCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| CorrectionAndPlot \t |  " << _measureCorrectionAndPlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectionAndPlotCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectionAndPlotCPUTime.getValue()  << " \t |  " << _measureCorrectionAndPlotCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectionAndPlotCalendarTime.getValue() << " \t |  " << _measureCorrectionAndPlotCPUTime.toString() << " \t |  " << _measureCorrectionAndPlotCalendarTime.toString() );
   if (logAllAdapters || _measurePlotCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Plot \t |  " << _measurePlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotCPUTime.getValue()  << " \t |  " << _measurePlotCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotCalendarTime.getValue() << " \t |  " << _measurePlotCPUTime.toString() << " \t |  " << _measurePlotCalendarTime.toString() );

}


void exahype::repositories::RepositoryArrayStack::clearIterationStatistics() {
   _measureAugmentedAMRGridCPUTime.erase();
   _measurePlotAugmentedAMRGridCPUTime.erase();
   _measureInitialConditionAndTimeStepSizeComputationCPUTime.erase();
   _measurePredictionAndPlotAndTimeStepSizeComputationCPUTime.erase();
   _measurePredictionAndTimeStepSizeComputationCPUTime.erase();
   _measureGridErasingCPUTime.erase();
   _measureADERDGTimeStepCPUTime.erase();
   _measureADERDGTimeStepAndPlotCPUTime.erase();
   _measurePredictionRerunCPUTime.erase();
   _measureRiemannSolverCPUTime.erase();
   _measurePredictionCPUTime.erase();
   _measureCorrectionCPUTime.erase();
   _measureCorrectionAndPlotCPUTime.erase();
   _measurePlotCPUTime.erase();

   _measureAugmentedAMRGridCalendarTime.erase();
   _measurePlotAugmentedAMRGridCalendarTime.erase();
   _measureInitialConditionAndTimeStepSizeComputationCalendarTime.erase();
   _measurePredictionAndPlotAndTimeStepSizeComputationCalendarTime.erase();
   _measurePredictionAndTimeStepSizeComputationCalendarTime.erase();
   _measureGridErasingCalendarTime.erase();
   _measureADERDGTimeStepCalendarTime.erase();
   _measureADERDGTimeStepAndPlotCalendarTime.erase();
   _measurePredictionRerunCalendarTime.erase();
   _measureRiemannSolverCalendarTime.erase();
   _measurePredictionCalendarTime.erase();
   _measureCorrectionCalendarTime.erase();
   _measureCorrectionAndPlotCalendarTime.erase();
   _measurePlotCalendarTime.erase();

}
