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
 
#ifndef EXAHYPE_MAPPINGS_LocalRecomputation_H_
#define EXAHYPE_MAPPINGS_LocalRecomputation_H_

#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

#include "peano/CommunicationSpecification.h"
#include "peano/MappingSpecification.h"
#include "peano/grid/VertexEnumerator.h"

#include "tarch/multicore/MulticoreDefinitions.h"

#include "exahype/Cell.h"
#include "exahype/State.h"
#include "exahype/Vertex.h"

#include "exahype/solvers/Solver.h"

#include "exahype/solvers/TemporaryVariables.h"

namespace exahype {
namespace mappings {
class LocalRecomputation;
}
}

/**
 * This mapping is one of three mappings (+adapters) which together perform the limiter
 * status spreading and the local recomputation of troubled cells (and their direct) neighbours.
 *
 * |Mapping                 | Event                  | Action                                            |
 * -------------------------------------------------------------------------------------------------------
 * | LimiterStatusSpreading | touchVertexFirstTime   | Merge the face-wise limiter status between local neighbours.|
 * |                        | enterCell              | Determine a unified value of the merged face-wise limiter status values and write it to every face. (Do not update the cell-wise limiter status.)|
 * |                        | prepareSendToNeighbour | Send the unified ace-wise limiter status to neighbouring ranks.|
 * -------------------------------------------------------------------------------------------------------
 * | Reinitialisation       | mergeWithNeighbour     | Receive the neighbour ranks' cells' merged limiter status. Directly update the unified merged limiter status value of cells at the remote boundary |
 * |                        | touchVertexFirstTime   | Merge the limiter status between local neighbours (again).
 * |                        | enterCell              | Determine a unified value of the merged face-wise limiter status values and write it to every face. |
 * |                        |                        | (Do not update the cell-wise limiter status.)                                                       |
 * |                        |                        | Rollback the solution in troubled cells and their next two neighbours.
 * |                        | prepareSendToNeighbour | Based on the unified face-wise limiter status send interface values to the neighbours. |
 * ------------------------------------------------------------------------------------------------------
 * |LocalRecomputation      | mergeWithNeighbour     | Based on the unified face-wise limiter status receive or drop the interface values send by the neighbours.
 * |                        | touchVertexFirstTIme   | Based on the unified face-wise limiter status merge local direct neighbours.
 * |                        | enterCell              | Recompute the solution in the troubled cells (and their direct neighbours).
 * |                        |                        | Set the cell-wise limiter status to the unified face-wise limiter status.
 * |                        |                        | (The normal time marching does only consider the cell-wise limiter status from now on
 * |                        |                        | and overwrites the face-wise values uniformly with Ok or Troubled after
 * |                        |                        | evaluating the discrete maximum principle (DMP) and the physical admissibility detection (PAD).)
 * ------------------------------------------------------------------------
 *
 * This mapping is only used in the algorithmic section LocalRecomputation.
 * It is not used for in the algorithmic section APosterioriRefinement.
 *
 * @author Dominic Charrier
 */
class exahype::mappings::LocalRecomputation {
 private:
  /**
   * Logging device for the trace macros.
   */
  static tarch::logging::Log _log;

  #ifdef Debug // TODO(Dominic): Exclude shared memory etc.
  /*
   *  Counter for the interior face solves for debugging purposes.
   */
  int _interiorFaceMerges;
  /*
   *  Counter for the boundary face solves for debugging purposes.
   */
  int _boundaryFaceMerges;
  #endif

  /**
   * Local copy of the state which
   * is used to determine if a solver
   * is active in the current algorithm section.
   * (See exahype::runners::Runner for locations
   * where the algorithm section is set. The new
   * state is then broadcasted by Peano to all other ranks.)
   */
  exahype::State _localState;

  /**
   * A minimum time step size for each solver.
   */
  std::vector<double> _minTimeStepSizes;

  /**
   * A minimum cell size for each solver.
   */
  std::vector<double> _minCellSizes;

  /**
   * A maximum cell size for each solver.
   */
  std::vector<double> _maxCellSizes;

  /**
   * Prepare a appropriately sized vector _minTimeStepSizes
   * with elements initiliased to MAX_DOUBLE.
   */
  void prepareLocalTimeStepVariables();

  exahype::solvers::PredictionTemporaryVariables     _predictionTemporaryVariables;

  exahype::solvers::MergingTemporaryVariables        _mergingTemporaryVariables;

  /**
   * Initialises the temporary variables.
   *
   * \note We parallelise over the domain
   * (mapping is copied for each thread) and
   * over the solvers registered on a cell.
   *
   * \note We need to initialise the temporary variables
   * in this mapping and not in the solvers since the
   * solvers in exahype::solvers::RegisteredSolvers
   * are not copied for every thread.
   */
  void initialiseTemporaryVariables();

  /**
   * Deletes the temporary variables.
   *
   * \note We need to initialise the temporary variables
   * in this mapping and not in the solvers since the
   * solvers in exahype::solvers::RegisteredSolvers
   * are not copied for every thread.
   */
  void deleteTemporaryVariables();

  /**
   * \return if we perform a local recomputation for this solver.
   */
  static bool performLocalRecomputation(exahype::solvers::Solver* solver);

  #ifdef Parallel
  /**
   * TODO(Dominic): Add docu.
   */
  static void dropNeighbourData(
      const int                                    fromRank,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const int                                    srcCellDescriptionIndex,
      const int                                    destCellDescriptionIndex,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level,
      const exahype::MetadataHeap::HeapEntries& receivedMetadata);

  /**
   * TODO(Dominic): Add docu.
   */
  void mergeNeighourData(
      const int                                    fromRank,
      const tarch::la::Vector<DIMENSIONS,int>&     src,
      const tarch::la::Vector<DIMENSIONS,int>&     dest,
      const int                                    srcCellDescriptionIndex,
      const int                                    destCellDescriptionIndex,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level,
      const exahype::MetadataHeap::HeapEntries& receivedMetadata);
  #endif

 public:
  /**
   * \see exahype::mappings::Merging::touchVertexLastTimeSpecification()
   */
  peano::MappingSpecification touchVertexLastTimeSpecification(int level) const;
  /**
   * Run through the whole tree. Run concurrently on the fine grid.
   */
  peano::MappingSpecification enterCellSpecification(int level) const;

  /**
   * Nop.
   */
  peano::MappingSpecification touchVertexFirstTimeSpecification(int level) const;
  /**
   * Nop.
   */
  peano::MappingSpecification leaveCellSpecification(int level) const;
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
   * TODO(Dominic): Add docu.
   */
  void enterCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Prolongs adjacency indices down from the
   * coarse grid and performs a merge with boundary data
   * if the indices indicate a domain boundary and
   * a limiting ADER-DG solver computes with FV.
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
   * TODO(Dominic): Add docu.
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
   * Initialise temporary variables
   * if they are not initialised yet (or
   * if a new solver was introuced to the grid.
   * This is why we put the initialisation
   * in beginIteration().
   *
   * In debug mode, further resets counters for Riemann solves at
   * interior
   */
  void beginIteration(exahype::State& solverState);

  /**
   * Frees previously allocated temporary variables.
   *
   * In debug mode, prints the output of counters.
   */
  void endIteration(exahype::State& solverState);

  #if defined(SharedMemoryParallelisation)
  /**
   * Copy the local state object over to the worker thread.
   */
  LocalRecomputation(const LocalRecomputation& masterThread);
  #endif

  /**
   * Free previously allocated temporary variables.
   */
  virtual ~LocalRecomputation();

#ifdef Parallel
  /**
   * TODO(Dominic): Add docu.
   */
  void mergeWithNeighbour(exahype::Vertex& vertex,
                          const exahype::Vertex& neighbour, int fromRank,
                          const tarch::la::Vector<DIMENSIONS, double>& x,
                          const tarch::la::Vector<DIMENSIONS, double>& h,
                          int level);
  /**
   * TODO(Dominic): Add docu.
   */
  void prepareSendToNeighbour(exahype::Vertex& vertex, int toRank,
                              const tarch::la::Vector<DIMENSIONS, double>& x,
                              const tarch::la::Vector<DIMENSIONS, double>& h,
                              int level);



  //
  // Below all methods are nop.
  //
  //===================================



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
  bool prepareSendToWorker(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      int worker);
  /**
   * Nop.
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
   * Nop.
   */
  void prepareSendToMaster(
      exahype::Cell& localCell, exahype::Vertex* vertices,
      const peano::grid::VertexEnumerator& verticesEnumerator,
      const exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);
  /**
   * Nop.
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
   * Nop.
   */
  void mergeWithWorker(exahype::Cell& localCell,
                       const exahype::Cell& receivedMasterCell,
                       const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
                       const tarch::la::Vector<DIMENSIONS, double>& cellSize,
                       int level);
  /**
   * Nop.
   */
  void mergeWithWorker(exahype::Vertex& localVertex,
                       const exahype::Vertex& receivedMasterVertex,
                       const tarch::la::Vector<DIMENSIONS, double>& x,
                       const tarch::la::Vector<DIMENSIONS, double>& h,
                       int level);
#endif

#if defined(SharedMemoryParallelisation)
  /**
   * Nop.
   */
  void mergeWithWorkerThread(const LocalRecomputation& workerThread);
#endif
  /**
   * Nop
   */
  LocalRecomputation();

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
  void leaveCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

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
