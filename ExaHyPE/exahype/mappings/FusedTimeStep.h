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
 
#ifndef EXAHYPE_MAPPINGS_FusedTimeStep_H_
#define EXAHYPE_MAPPINGS_FusedTimeStep_H_

#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

#include "peano/CommunicationSpecification.h"
#include "peano/MappingSpecification.h"
#include "peano/grid/VertexEnumerator.h"

#include "tarch/multicore/MulticoreDefinitions.h"

#include "exahype/Cell.h"
#include "exahype/State.h"
#include "exahype/Vertex.h"


namespace exahype {
namespace mappings {
class FusedTimeStep;
}
}

/**
 * Update the solution
 *
 * <h2>ADER-DG</h2>
 * We run over all cells of the local spacetree
 *
 * - take the predicted data within each cell
 * - add the contribution from the Riemann solve to the predicted data
 * - mark the cell with the new time stamp
 *
 * All this is done in enterCell().
 *
 * <h2>Finite Volumes</h2>
 * This is where we do the actual time stepping with the finite volume scheme.
 *
 * @developers:
 * TODO(Dominic): I think we should extract the Riemann solve from the finite
 * volume kernels and move it into the RiemannSolver mapping
 * since I will exchange metadata and face data here between MPI ranks.
 *
 * @author Dominic Charrier Tobias Weinzierl
 */
class exahype::mappings::FusedTimeStep {
private:
  /**
   * Logging device for the trace macros.
   */
  static tarch::logging::Log _log;

  /**
   * This semaphore is used for locking the plotters'
   * plotPatch function which is usually not thread-safe.
   */
  static tarch::multicore::BooleanSemaphore SemaphoreForPlotting;

  /**
   * A local copy of the state set
   * in beginIteration(...).
   */
  exahype::State _stateCopy;

  /**
   * A minimum time step size for each solver.
   */
  std::vector<double> _minTimeStepSizes;
  /**
   * The maximum level occupied by cells of a solver.
   */
  std::vector<int> _maxLevels;
  /**
   * Per solver a flag, indicating if has requested
   * a mesh update request or a limiter domain change.
   */
  std::vector<exahype::solvers::Solver::MeshUpdateEvent> _meshUpdateEvents;

  /**
   * Prepare the vectors _minTimeStepSizes, _maxLevels,
   * _meshUpdateRequests, _limiterDomainChanges.
   */
  void initialiseLocalVariables();

  /**
   * Indicates that the background tasks have terminated.
   * No further checks are required in this case.
   *
   * Is initialised with false for the main thread
   * and for the worker threads.
   * As the worker threads; mappings are destroyed but
   * the main thread's mapping continues to
   * exist we reset this value in endIteration(State) to false.
   *
   * We process background jobs in touchVertexFirstTime(...)
   * and set this flag here as well.
   */
  bool _batchIterationCounterUpdated = false;


  int _batchIteration = 0;
  /**
   * Updates a batch iteration counter and
   * ensures that all background jobs of a certain type have terminated.
   * The type of the background jobs is determined based on the batch iteration counter.
   *
   * Signal the solvers that we start a new time step.
   *
   * \note As the state can only be copied in beginIteration(...) and
   * we need to know the batch iteration in touchVertexFirstTime(...),
   * we cannot rely on its batch iteration counter and
   * have to count ourselves.
   * Note that touchVertexFirstTime(...) might be invoked before beginIteration(...)
   * if state broadcasting is turned off.
   *
   */
  void updateBatchIterationCounter(bool initialiseBatchIterationCounter);

  /**
   * \return if the mappings/adapters
   * FusedTimeStep, Prediction, PredictionRerun, and PredictionOrLocalRecomputation
   * are supposed to send out riemann data in this iteration.
   *
   * \see updatePredictionIterationTag(...)
   */
  bool sendOutRiemannDataInThisIteration() const;

  /**
   * \return if the mappings/adapters
   * FusedTimeStep, Prediction, PredictionRerun, and PredictionOrLocalRecomputation
   * are supposed to issue prediction jobs in this iteration.
   *
   * \see updatePredictionIterationTag(...)
   */
  bool issuePredictionJobsInThisIteration() const;

 public:
  /**
   * Run through the whole tree. Run concurrently on the fine grid.
   *
   * Alters the state if we perform a reduction. This
   * is the case if we perform the last iteration of a batch
   * or no batch iteration at all.
   */
  peano::MappingSpecification enterCellSpecification(int level) const;
  /**
   * Run through the whole tree. Run concurrently on the fine grid.
   */
  peano::MappingSpecification leaveCellSpecification(int level) const;
  /**
   * Run through the whole tree. Avoid fine grid races.
   *
   * Alters the state as we have a counter which checks
   * if we have waited for the background jobs to complete.
   */
  peano::MappingSpecification touchVertexFirstTimeSpecification(int level) const;

  /**
   * Nop.
   */
  peano::MappingSpecification touchVertexLastTimeSpecification(int level) const;
  /**
   * Nop.
   */
  peano::MappingSpecification ascendSpecification(int level) const;
  /**
   * Nop.
   */
  peano::MappingSpecification descendSpecification(int level) const;

  /**
   * No data needs to be synchronised between masters and workers.
   */
  peano::CommunicationSpecification communicationSpecification() const;

  /**
   * Delete temporary variables.
   */
  virtual ~FusedTimeStep();

  #if defined(SharedMemoryParallelisation)
  /**
   * Copy the _backgroundJobsHaveTerminated
   * and _batchIteration fields form the masterThread
   */
  FusedTimeStep(const FusedTimeStep& masterThread);
  /**
   * Merge time step data and flags of the worker thread
   * with the master thread.
   */
  void mergeWithWorkerThread(const FusedTimeStep& workerThread);
  #endif

  /**
   * Merge with the neighbours but check beforehand
   * if all backgrounds have terminated.
   */
  void touchVertexFirstTime(
      exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);

  /**
   * If the fine grid cell functions as compute cell for a solver,
   * we update the solution of the solver within the fine grid cell
   * (at a given time).
   * We further ask the solver if we need to adjust its solution
   * values within the fine grid cell (at a given time).
   * If so, we call the corresponding solver routine
   * that adjusts the solvers' solution values within the fine grid cell.
   *
   * <h2>ADER-DG<h2>
   * For ADER-DG solvers, we call the surfaceIntegral(...) routine before
   * we call the FusedTimeStep(...) routine of the solvers.
   *
   * <h2>Finite volumes<h2>
   * For finite volume solvers, we simply call the
   * FusedTimeStep(...) routine.
   */
  void enterCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);


  /**
   * Perform a restriction for solvers who requested it.
   */
  void leaveCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   *Resets the next mesh update request flag to false and
   * the next limiter domain change to Regular
   * using the "setNext..." methods.
   *
   * Updates the prediction iteration tag in the first iteration
   * of a batch or if no batch is run.
   */
  void beginIteration(exahype::State& solverState);

  /**
   * For all solvers, overwrite the current
   * gridUpdateRequested value with the next value.
   *
   * Update the global solver states (next)limiterDomainHasChanged
   * with values from the temporary variables.
   *
   * Finish plotting if a plotter is active.
   *
   * <h2> Background Jobs </h2>
   * Reset the _backgroundJobsHaveTerminated switch.
   *
   * Notify Peano's tarch that we want to start processing
   * background jobs with all available cores.
   *
   * Updates the prediction iteration tag in every iteration.
   */
  void endIteration(exahype::State& solverState);

  //
  // Below every method is nop.
  //
  // ==================================

#ifdef Parallel
  /**
   * TODO(Dominic: Add docu.
   */
  void mergeWithNeighbour(exahype::Vertex& vertex,
                          const exahype::Vertex& neighbour, int fromRank,
                          const tarch::la::Vector<DIMENSIONS, double>& x,
                          const tarch::la::Vector<DIMENSIONS, double>& h,
                          int level);
  /**
   * TODO(Dominic: Add docu.
   */
  void prepareSendToNeighbour(exahype::Vertex& vertex, int toRank,
                              const tarch::la::Vector<DIMENSIONS, double>& x,
                              const tarch::la::Vector<DIMENSIONS, double>& h,
                              int level);

  /**
   * TODO(Dominic: Add docu.
   *
   * \return Returns true in the first or last iteration of a batch or
   * if no batch is run. This return value indicates that a worker will
   * perform a reduction and will be be broadcasted down to
   * the worker.
   *
   * If we do not broadcast it down to the worker,
   * the worker does not know about it. It thus has be broadcasted
   * down to the worker in the first iteration of the batch.
   *
   * The master must be aware of it in the last iteration of the
   * batch. So we have to return true here as well.
   */
  bool prepareSendToWorker(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      int worker);

  /**
   * TODO(Dominic: Add docu.
   */
  void mergeWithMaster(
      const exahype::Cell& workerGridCell,
      exahype::Vertex* const workerGridVertices,
      const peano::grid::VertexEnumerator& workerEnumerator,
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      int worker, const exahype::State& workerState,
      exahype::State& masterState);

  /**
   * TODO(Dominic: Add docu.
   */
  void prepareSendToMaster(
      exahype::Cell& localCell, exahype::Vertex* vertices,
      const peano::grid::VertexEnumerator& verticesEnumerator,
      const exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * TODO(Dominic: Add docu.
   */
  void receiveDataFromMaster(
      exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
      const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
      exahype::Vertex* const receivedCoarseGridVertices,
      const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
      exahype::Cell& receivedCoarseGridCell,
      exahype::Vertex* const workersCoarseGridVertices,
      const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
      exahype::Cell& workersCoarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * TODO(Dominic: Add docu.
   */
  void mergeWithWorker(exahype::Cell& localCell,
                       const exahype::Cell& receivedMasterCell,
                       const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
                       const tarch::la::Vector<DIMENSIONS, double>& cellSize,
                       int level);

  /**
   * Nop.
   */
  void prepareCopyToRemoteNode(exahype::Vertex& localVertex, int toRank,
                               const tarch::la::Vector<DIMENSIONS, double>& x,
                               const tarch::la::Vector<DIMENSIONS, double>& h,
                               int level);
  /**
   * Nop.
   */
  void prepareCopyToRemoteNode(
      exahype::Cell& localCell, int toRank,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);

  /**
   * Nop.
   */
  void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
      int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
      const tarch::la::Vector<DIMENSIONS, double>& h, int level);

  /**
   * Nop.
   */
  void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
      int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);

  /**
   * Nop.
   */
  void mergeWithWorker(exahype::Vertex& localVertex,
                       const exahype::Vertex& receivedMasterVertex,
                       const tarch::la::Vector<DIMENSIONS, double>& x,
                       const tarch::la::Vector<DIMENSIONS, double>& h,
                       int level);
#endif

  /**
   * Nop.
   */
  FusedTimeStep();

  /**
   * Nop.
   */
  void createInnerVertex(
      exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
  /**
   * Nop.
   */
  void createBoundaryVertex(
      exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
  /**
   * Nop.
   */
  void createHangingVertex(
      exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
  /**
   * Nop.
   */
  void destroyHangingVertex(
      const exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
  /**
   * Nop.
   */
  void destroyVertex(
      const exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
  /**
   * Nop.
   */
  void createCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);
  /**
   * Nop.
   */
  void destroyCell(
      const exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Nop.
   */
  void touchVertexLastTime(
      exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);

  /**
   * Nop.
   */
  void descend(
      exahype::Cell* const fineGridCells,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell);
  /**
   * Nop.
   */
  void ascend(exahype::Cell* const fineGridCells,
              exahype::Vertex* const fineGridVertices,
              const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
              exahype::Vertex* const coarseGridVertices,
              const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
              exahype::Cell& coarseGridCell);
};
#endif
