#include "exahype/mappings/RiemannSolver.h"

#include "tarch/multicore/Lock.h"

#include "peano/utils/Globals.h"

#include "exahype/solvers/Solver.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification   exahype::mappings::RiemannSolver::communicationSpecification() {
  return peano::CommunicationSpecification(peano::CommunicationSpecification::SendDataAndStateBeforeFirstTouchVertexFirstTime,peano::CommunicationSpecification::SendDataAndStateAfterLastTouchVertexLastTime,false);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::RiemannSolver::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::OnlyLeaves,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::RiemannSolver::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::RiemannSolver::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::RiemannSolver::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::RiemannSolver::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidCoarseGridRaces);
}


/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification   exahype::mappings::RiemannSolver::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidCoarseGridRaces);
}


tarch::logging::Log                exahype::mappings::RiemannSolver::_log( "exahype::mappings::RiemannSolver" ); 


exahype::mappings::RiemannSolver::RiemannSolver() {
  // do nothing
}


exahype::mappings::RiemannSolver::~RiemannSolver() {
  // do nothing
}


#if defined(SharedMemoryParallelisation)
exahype::mappings::RiemannSolver::RiemannSolver(const RiemannSolver&  masterThread) {
  logTraceIn( "RiemannSolver(RiemannSolver)" );

  _localState.setOldTimeStepSize(masterThread.getState().getOldTimeStepSize());

  logTraceOut( "RiemannSolver(RiemannSolver)" );
}


void exahype::mappings::RiemannSolver::mergeWithWorkerThread(const RiemannSolver& workerThread) {
  // do nothing
}
#endif


void exahype::mappings::RiemannSolver::createHangingVertex(
    exahype::Vertex&     fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
    exahype::Vertex * const   coarseGridVertices,
    const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
    exahype::Cell&       coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  // do nothing
}


void exahype::mappings::RiemannSolver::destroyHangingVertex(
    const exahype::Vertex&   fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  // do nothing
}


void exahype::mappings::RiemannSolver::createInnerVertex(
    exahype::Vertex&               fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  // do nothing
}


void exahype::mappings::RiemannSolver::createBoundaryVertex(
    exahype::Vertex&               fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  // do nothing
}


void exahype::mappings::RiemannSolver::destroyVertex(
    const exahype::Vertex&   fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  // do nothing
}


void exahype::mappings::RiemannSolver::createCell(
    exahype::Cell&                 fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  // do nothing
}


void exahype::mappings::RiemannSolver::destroyCell(
    const exahype::Cell&           fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::RiemannSolver::mergeWithNeighbour(
    exahype::Vertex&  vertex,
    const exahype::Vertex&  neighbour,
    int                                           fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::RiemannSolver::prepareSendToNeighbour(
    exahype::Vertex&  vertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::RiemannSolver::prepareCopyToRemoteNode(
    exahype::Vertex&  localVertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::RiemannSolver::prepareCopyToRemoteNode(
    exahype::Cell&  localCell,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
    int                                           level
) {
  // do nothing
}

void exahype::mappings::RiemannSolver::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex&  localVertex,
    const exahype::Vertex&  masterOrWorkerVertex,
    int                                       fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&  x,
    const tarch::la::Vector<DIMENSIONS,double>&  h,
    int                                       level
) {
  // do nothing
}

void exahype::mappings::RiemannSolver::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell&  localCell,
    const exahype::Cell&  masterOrWorkerCell,
    int                                       fromRank,
    const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
    int                                       level
) {
  // do nothing
}

bool exahype::mappings::RiemannSolver::prepareSendToWorker(
    exahype::Cell&                 fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
    int                                                                  worker
) {
  // do nothing
  return true;
}

void exahype::mappings::RiemannSolver::prepareSendToMaster(
    exahype::Cell&                       localCell,
    exahype::Vertex *                    vertices,
    const peano::grid::VertexEnumerator&       verticesEnumerator,
    const exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
    const exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
  // do nothing
}


void exahype::mappings::RiemannSolver::mergeWithMaster(
    const exahype::Cell&           workerGridCell,
    exahype::Vertex * const        workerGridVertices,
    const peano::grid::VertexEnumerator& workerEnumerator,
    exahype::Cell&                 fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
    int                                                                  worker,
    const exahype::State&          workerState,
    exahype::State&                masterState
) {
  // do nothing
}


void exahype::mappings::RiemannSolver::receiveDataFromMaster(
    exahype::Cell&                        receivedCell,
    exahype::Vertex *                     receivedVertices,
    const peano::grid::VertexEnumerator&        receivedVerticesEnumerator,
    exahype::Vertex * const               receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator&        receivedCoarseGridVerticesEnumerator,
    exahype::Cell&                        receivedCoarseGridCell,
    exahype::Vertex * const               workersCoarseGridVertices,
    const peano::grid::VertexEnumerator&        workersCoarseGridVerticesEnumerator,
    exahype::Cell&                        workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&    fineGridPositionOfCell
) {
  // do nothing
}


void exahype::mappings::RiemannSolver::mergeWithWorker(
    exahype::Cell&           localCell,
    const exahype::Cell&     receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
    int                                          level
) {
  // do nothing
}


void exahype::mappings::RiemannSolver::mergeWithWorker(
    exahype::Vertex&        localVertex,
    const exahype::Vertex&  receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
  // do nothing
}
#endif

void exahype::mappings::RiemannSolver::touchVertexFirstTime(
    exahype::Vertex&               fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  // do nothing
}

void exahype::mappings::RiemannSolver::solveRiemannProblem(
    tarch::la::Vector<TWO_POWER_D,int>& adjacentADERDGCellDescriptionsIndices,
    const int cellIndexL,
    const int cellIndexR,
    const int faceL,
    const int faceR,
    const double * const normal
) {
  // Only continue if this is an internal face. See multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex.
  if (
      adjacentADERDGCellDescriptionsIndices[cellIndexL] < 0
      ||
      adjacentADERDGCellDescriptionsIndices[cellIndexR] < 0
  ) {
    return;
  }

  logDebug("touchVertexLastTime(...)::solveRiemannProblem(...)",
           "Performing Riemann solve. "
           << "faceL:" << faceL
           << " faceR:" << faceR
           << " cell descr. index L:" << adjacentADERDGCellDescriptionsIndices[cellIndexL]
                                                                               << " cell descr. index R:" << adjacentADERDGCellDescriptionsIndices[cellIndexR]);

  std::vector<records::ADERDGCellDescription>&  cellDescriptionsL =
      ADERDGCellDescriptionHeap::getInstance().getData( adjacentADERDGCellDescriptionsIndices[cellIndexL] );
  std::vector<records::ADERDGCellDescription>&  cellDescriptionsR =
      ADERDGCellDescriptionHeap::getInstance().getData( adjacentADERDGCellDescriptionsIndices[cellIndexR] );

  assertion1WithExplanation( cellDescriptionsL.size() == cellDescriptionsR.size(), cellDescriptionsL.size(), "The number of ADERDGCellDescriptions is not the same for both cells!" );

  for (
      int p = 0;
      p < cellDescriptionsL.size();
      p++
  ) {
    exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[ cellDescriptionsL[p].getSolverNumber() ];


    bool riemannSolveNotPerformed = false;
    {
      // Lock the critical multithreading area.
      tarch::multicore::Lock lock(_semaphore);
      riemannSolveNotPerformed = !cellDescriptionsL[p].getRiemannSolvePerformed(faceL)
                                                               &&
                                                               !cellDescriptionsR[p].getRiemannSolvePerformed(faceR);

      if(riemannSolveNotPerformed==true) {
        cellDescriptionsL[p].setRiemannSolvePerformed(faceL,true);
        cellDescriptionsR[p].setRiemannSolvePerformed(faceR,true);
      }
    } // Unlock the critical multithreading area by letting lock go out of scope.

    const int numberOfFaceDof = tarch::la::aPowI(DIMENSIONS-1,solver->getNodesPerCoordinateAxis());

    if (riemannSolveNotPerformed) {
      double * QL = &(DataHeap::getInstance().getData(cellDescriptionsL[p].getExtrapolatedPredictor())[faceL * numberOfFaceDof]._persistentRecords._u);
      double * QR = &(DataHeap::getInstance().getData(cellDescriptionsR[p].getExtrapolatedPredictor())[faceR * numberOfFaceDof]._persistentRecords._u);

      double * FL = &(DataHeap::getInstance().getData(cellDescriptionsL[p].getFluctuation())[faceL * numberOfFaceDof]._persistentRecords._u);
      double * FR = &(DataHeap::getInstance().getData(cellDescriptionsR[p].getFluctuation())[faceR * numberOfFaceDof]._persistentRecords._u);
/*

      solver->riemannSolver(
          FL,
          FR,
          QL,
          QR,
          _localState.getMaxTimeStepSize(),
          normal
      );
*/
    }
  }
}

void exahype::mappings::RiemannSolver::touchVertexLastTime(
    exahype::Vertex&         fineGridVertex,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
    const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  logTraceInWith6Arguments( "touchVertexLastTime(...)", fineGridVertex, fineGridX, fineGridH, coarseGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfVertex );

  // @todo Tobias Weinzierl
  // Delegate to solver-specific code fragments

  //
  assertion1WithExplanation(_localState.getMaxTimeStepSize() < std::numeric_limits<double>::max(),_localState.toString(),"Old time step size was not initialised correctly!");

  tarch::la::Vector<TWO_POWER_D,int>& adjacentADERDGCellDescriptionsIndices = fineGridVertex.getADERDGCellDescriptionsIndex();
  logDebug("touchVertexLastTime(...)","cell descriptions around vertex. "
           << "coarse grid level: " << coarseGridVerticesEnumerator.getLevel()
           << ", fine grid pos.:" << fineGridPositionOfVertex
           << ", adjac. cell descr.:" << adjacentADERDGCellDescriptionsIndices);
  logDebug("touchVertexLastTime(...)","cell descriptions around vertex. "
           << "fine grid x " << fineGridX);

  // todo: Dominic Etienne Charrier: Reverse engineered indices from
  // PatchInitialisation2MultiscaleLinkedCell_1::touchVertexLastTime(...)
  // Not sure what happens with hanging nodes.
  /* Right cell-left cell   pair indices: 0,1; 2,3;   4,5; 6;7
   * Front cell-back cell   pair indices: 0,2; 1,3;   4,6; 5;7
   * Top   cell-bottom cell pair indices: 0,4; 1,5;   2,6; 3;7
   *
   * Note that from the viewpoint of a cell, the face
   * has always the "opposite" index, i.e., we solve a Riemann
   * problem on the left face of the right cell (which
   * is the right face of the left cell).
   */
  // index maps (
  constexpr int cellIndicesLeft   [4] = { 0, 2, 4, 6 };
  constexpr int cellIndicesRight  [4] = { 1, 3, 5, 7 };
  constexpr int cellIndicesFront  [4] = { 0, 1, 4, 5 };
  constexpr int cellIndicesBack   [4] = { 2, 3, 6, 7 };
#if DIMENSIONS==3
  constexpr int cellIndicesBottom [4] = { 0, 1, 2, 3 };
  constexpr int cellIndicesTop    [4] = { 4, 5, 6, 7 };
#endif

  // normal vectors
  const double nx[3]= { 1., 0., 0. };
  const double ny[3]= { 0., 1., 0. };
#if DIMENSIONS==3
  const double nz[3]= { 0., 0., 1. };
#endif

  // Left/right face
  for (int i=0; i<TWO_POWER_D_DIVIDED_BY_TWO; i++) {
    solveRiemannProblem(
        adjacentADERDGCellDescriptionsIndices,
        cellIndicesLeft [i],
        cellIndicesRight[i],
        EXAHYPE_FACE_RIGHT,
        EXAHYPE_FACE_LEFT,
        nx);

    solveRiemannProblem(
        adjacentADERDGCellDescriptionsIndices,
        cellIndicesFront[i],
        cellIndicesBack [i],
        EXAHYPE_FACE_BACK,
        EXAHYPE_FACE_FRONT,
        ny);

#if DIMENSIONS==3
    solveRiemannProblem(
        adjacentADERDGCellDescriptionsIndices,
        cellIndicesBottom[i],
        cellIndicesTop   [i],
        EXAHYPE_FACE_TOP,
        EXAHYPE_FACE_BOTTOM,
        nz);
#endif
  }
  logTraceOutWith1Argument( "touchVertexLastTime(...)", fineGridVertex );
}

void exahype::mappings::RiemannSolver::enterCell(
    exahype::Cell&                 fineGridCell,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  // do nothing
}


void exahype::mappings::RiemannSolver::leaveCell(
    exahype::Cell&           fineGridCell,
    exahype::Vertex * const  fineGridVertices,
    const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
  // do nothing
}


void exahype::mappings::RiemannSolver::beginIteration(
    exahype::State&  solverState
) {
  logTraceInWith1Argument( "beginIteration(State)", solverState );

  _localState = solverState;

  logTraceOutWith1Argument( "beginIteration(State)", solverState);
}


void exahype::mappings::RiemannSolver::endIteration(
    exahype::State&  solverState
) {
  // do nothing
}



void exahype::mappings::RiemannSolver::descend(
    exahype::Cell * const          fineGridCells,
    exahype::Vertex * const        fineGridVertices,
    const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
    exahype::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
    exahype::Cell&                 coarseGridCell
) {
  // do nothing
}


void exahype::mappings::RiemannSolver::ascend(
    exahype::Cell * const    fineGridCells,
    exahype::Vertex * const  fineGridVertices,
    const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
    exahype::Vertex * const  coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&           coarseGridCell
) {
  // do nothing
}

const exahype::State& exahype::mappings::RiemannSolver::getState() const {
  return _localState;
}
