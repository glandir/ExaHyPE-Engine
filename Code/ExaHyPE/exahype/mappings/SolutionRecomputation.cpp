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
 
#include "exahype/mappings/SolutionRecomputation.h"

#include "tarch/multicore/Loop.h"

#include "peano/utils/Globals.h"
#include "peano/utils/Loop.h"
#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/VertexOperations.h"
#include "multiscalelinkedcell/HangingVertexBookkeeper.h"


peano::CommunicationSpecification
exahype::mappings::SolutionRecomputation::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange,
      true);
}


peano::MappingSpecification
exahype::mappings::SolutionRecomputation::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification
exahype::mappings::SolutionRecomputation::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification
exahype::mappings::SolutionRecomputation::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification
exahype::mappings::SolutionRecomputation::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::SolutionRecomputation::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}


peano::MappingSpecification
exahype::mappings::SolutionRecomputation::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::SolutionRecomputation::_log(
    "exahype::mappings::SolutionRecomputation");

void exahype::mappings::SolutionRecomputation::prepareTemporaryVariables() {
  assertion(_tempStateSizedVectors       ==nullptr);
  assertion(_tempStateSizedSquareMatrices==nullptr);
  assertion(_tempFaceUnknowns            ==nullptr);
  assertion(_tempUnknowns                ==nullptr);

  int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  _tempStateSizedVectors        = new double**[numberOfSolvers];
  _tempStateSizedSquareMatrices = new double**[numberOfSolvers];
  _tempFaceUnknowns             = new double**[numberOfSolvers];
  _tempUnknowns                 = new double**[numberOfSolvers];
//    _tempSpaceTimeFaceUnknownsArray = new double* [numberOfSolvers]; todo

  int solverNumber=0;
  for (auto solver : exahype::solvers::RegisteredSolvers) {
    int numberOfStateSizedVectors  = 0; // TOOD(Dominic): Check if we need number of parameters too
    int numberOfStateSizedMatrices = 0;
    int numberOfFaceUnknowns       = 0;
    int lengthOfFaceUnknowns       = 0;
    int numberOfUnknowns           = 0;
    int lengthOfUnknowns           = 0;
    switch (solver->getType()) {
      case exahype::solvers::Solver::Type::ADER_DG:
        numberOfStateSizedVectors  = 5; // See riemannSolverLinear
        numberOfStateSizedMatrices = 3; // See riemannSolverLinear
        numberOfFaceUnknowns       = 3; // See exahype::solvers::ADERDGSolver::applyBoundaryConditions
        lengthOfFaceUnknowns       =
            static_cast<exahype::solvers::ADERDGSolver*>(solver)->getUnknownsPerFace();
        break;
      case exahype::solvers::Solver::Type::LimitingADERDG:
        // Needs the same temporary variables as the normal ADER-DG scheme.
        numberOfStateSizedVectors  = 5;
        numberOfStateSizedMatrices = 3;
        numberOfFaceUnknowns       = 3;
        lengthOfFaceUnknowns       = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
            _solver->getUnknownsPerFace();
        break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        // TODO(Dominic): We do not consider high-order FV methods yet;
        // numberOfUnknowns is thus set to zero.
        numberOfUnknowns = 0;
        numberOfStateSizedVectors = 5;
        break;
    }

    if (numberOfStateSizedVectors>0) {
      _tempStateSizedVectors[solverNumber] = new double*[numberOfStateSizedVectors];
      for (int i=0; i<numberOfStateSizedVectors; ++i) { // see riemanSolverLinear
        _tempStateSizedVectors[solverNumber][i] = new double[solver->getNumberOfVariables()];
      }
    }
    //
    if (numberOfStateSizedMatrices>0) {
      _tempStateSizedSquareMatrices[solverNumber] = new double*[numberOfStateSizedMatrices];
      for (int i=0; i<3; ++i) { // see riemanSolverLinear
        _tempStateSizedSquareMatrices[solverNumber][i] =
            new double[solver->getNumberOfVariables() * solver->getNumberOfVariables()];
      }
    }
    //
    if (numberOfFaceUnknowns>0) {
      _tempFaceUnknowns[solverNumber] = new double*[3];
      for (int i=0; i<3; ++i) { // see ADERDGSolver::applyBoundaryConditions(...)
        _tempFaceUnknowns[solverNumber][i] = new double[lengthOfFaceUnknowns];
      }
    }
    //
    if (numberOfUnknowns>0) {
      _tempUnknowns[solverNumber] = new double*[3];
      for (int i=0; i<3; ++i) { // see ADERDGSolver::applyBoundaryConditions(...)
        _tempUnknowns[solverNumber][i] = new double[lengthOfUnknowns];
      }
    }

    ++solverNumber;
  }
}

void exahype::mappings::SolutionRecomputation::deleteTemporaryVariables() {
  if (_tempStateSizedVectors!=nullptr) {
    assertion(_tempStateSizedSquareMatrices!=nullptr);

    int solverNumber=0;
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      int numberOfStateSizedVectors  = 0;
      int numberOfStateSizedMatrices = 0;
      int numberOfFaceUnknowns       = 0;
      int numberOfUnknowns           = 0;
      switch (solver->getType()) {
        case exahype::solvers::Solver::Type::ADER_DG:
        case exahype::solvers::Solver::Type::LimitingADERDG:
          numberOfStateSizedVectors  = 5; // See riemannSolverLinear
          numberOfStateSizedMatrices = 3; // See riemannSolverLinear
          numberOfFaceUnknowns       = 3; // See exahype::solvers::ADERDGSolver::applyBoundaryConditions
          break;
        case exahype::solvers::Solver::Type::FiniteVolumes:
          // TODO(Dominic): We do not consider high-order FV methods yet;
          // numberOfUnknowns is thus set to zero.
          numberOfUnknowns = 0;
          numberOfStateSizedVectors = 5;
          break;
      }

      if (numberOfStateSizedVectors>0) {
        for (int i=0; i<numberOfStateSizedVectors; ++i) { // see riemanSolverLinear
          delete[] _tempStateSizedVectors[solverNumber][i];
        }
        delete[] _tempStateSizedVectors[solverNumber];
        _tempStateSizedVectors[solverNumber] = nullptr;
      }
      //
      if (numberOfStateSizedMatrices>0) {
        for (int i=0; i<numberOfStateSizedMatrices; ++i) { // see riemanSolverLinear
          delete[] _tempStateSizedSquareMatrices[solverNumber][i];
        }
        delete[] _tempStateSizedSquareMatrices[solverNumber];
        _tempStateSizedSquareMatrices[solverNumber] = nullptr;
      }
      //
      if (numberOfFaceUnknowns>0) {
        for (int i=0; i<numberOfFaceUnknowns; ++i) { // see riemanSolverLinear
          delete[] _tempFaceUnknowns[solverNumber][i];
        }
        delete[] _tempFaceUnknowns[solverNumber];
        _tempFaceUnknowns[solverNumber] = nullptr;
      }
      //
      if (numberOfUnknowns>0) {
        for (int i=0; i<numberOfUnknowns; ++i) { // see riemanSolverLinear
          delete[] _tempUnknowns[solverNumber][i];
        }
        delete[] _tempUnknowns[solverNumber];
        _tempUnknowns[solverNumber] = nullptr;
      }
      //
      // _tempSpaceTimeFaceUnknownsArray[solverNumber] = nullptr; // todo

      ++solverNumber;
    }

    delete[] _tempStateSizedVectors;
    delete[] _tempStateSizedSquareMatrices;
    delete[] _tempFaceUnknowns;
    delete[] _tempUnknowns;
//    delete[] _tempSpaceTimeFaceUnknownsArray; todo
    _tempStateSizedVectors        = nullptr;
    _tempStateSizedSquareMatrices = nullptr;
    _tempFaceUnknowns             = nullptr;
    _tempUnknowns                 = nullptr;
//    _tempSpaceTimeFaceUnknownsArray  = nullptr; todo
  }
}

exahype::mappings::SolutionRecomputation::SolutionRecomputation()
 #ifdef Debug
 :
 _interiorFaceMerges(0),
 _boundaryFaceMerges(0)
 #endif
{
  // do nothing
}

exahype::mappings::SolutionRecomputation::~SolutionRecomputation() {
  deleteTemporaryVariables();
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::SolutionRecomputation::SolutionRecomputation(
    const SolutionRecomputation& masterThread)
: _localState(masterThread._localState),
  _tempStateSizedVectors(nullptr),
  _tempUnknowns(nullptr) {
  prepareTemporaryVariables();
}

void exahype::mappings::SolutionRecomputation::mergeWithWorkerThread(
    const SolutionRecomputation& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::SolutionRecomputation::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (fineGridCell.isInitialised()) {
    const int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
    // please use a different UserDefined per mapping/event
    peano::datatraversal::autotuning::MethodTrace methodTrace = peano::datatraversal::autotuning::UserDefined7;
    #ifdef SharedMemoryParallelisation
    int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfSolvers, methodTrace);
    #endif
    pfor(i, 0, numberOfSolvers, grainSize)
      auto solver = exahype::solvers::RegisteredSolvers[i];

      if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
        const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),i);
        if (element!=exahype::solvers::Solver::NotFound) {
          auto* limitingADERSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
          limitingADERSolver->recomputeSolution(
              fineGridCell.getCellDescriptionsIndex(),
              element,
              _tempStateSizedVectors[i],
              _tempUnknowns[i],
              fineGridVertices,
              fineGridVerticesEnumerator);

          // It is important that we update the limiter status only after the recomputation since we use
          // the previous and current limiter status in the recomputation.
          limitingADERSolver->updateLimiterStatus(fineGridCell.getCellDescriptionsIndex(),element);
        }
      }
    endpfor
    peano::datatraversal::autotuning::Oracle::getInstance().parallelSectionHasTerminated(methodTrace);
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::SolutionRecomputation::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  _localState = solverState;

  prepareTemporaryVariables();

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::SolutionRecomputation::endIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("endIteration(State)", solverState);

  deleteTemporaryVariables();

  #if defined(Debug) // TODO(Dominic): Use logDebug if it works with filters
  logInfo("endIteration(...)","interiorFaceSolves: " << _interiorFaceMerges);
  logInfo("endIteration(...)","boundaryFaceSolves: " << _boundaryFaceMerges);
  #endif

  logTraceOutWith1Argument("endIteration(State)", solverState);
}

void exahype::mappings::SolutionRecomputation::touchVertexFirstTime(
  exahype::Vertex& fineGridVertex,
  const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
  const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
  exahype::Vertex* const coarseGridVertices,
  const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
  exahype::Cell& coarseGridCell,
  const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
    dfor2(pos1)
      dfor2(pos2)
        if (fineGridVertex.hasToMergeNeighbours(pos1,pos2)) { // Assumes that we have to valid indices
          const peano::datatraversal::autotuning::MethodTrace methodTrace =
              peano::datatraversal::autotuning::UserDefined2;
          #ifdef SharedMemoryParallelisation
          const int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().
              parallelise(solvers::RegisteredSolvers.size(), methodTrace);
          #endif
          pfor(solverNumber, 0, static_cast<int>(solvers::RegisteredSolvers.size()),grainSize)
            auto solver = exahype::solvers::RegisteredSolvers[solverNumber];
            const int cellDescriptionsIndex1 = fineGridVertex.getCellDescriptionsIndex()[pos1Scalar];
            const int cellDescriptionsIndex2 = fineGridVertex.getCellDescriptionsIndex()[pos2Scalar];
            const int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
            const int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
            if (element2>=0 && element1>=0) {
              solver->mergeNeighbours(
                  cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
                  _tempFaceUnknowns[solverNumber],
                  _tempStateSizedVectors[solverNumber],
                  _tempStateSizedSquareMatrices[solverNumber]); // todo uncomment
            }

            #ifdef Debug // TODO(Dominic)
            _interiorFaceMerges++;
            #endif
          endpfor
          peano::datatraversal::autotuning::Oracle::getInstance()
          .parallelSectionHasTerminated(methodTrace);

          fineGridVertex.setMergePerformed(pos1,pos2,true);
        }
        if (fineGridVertex.hasToMergeWithBoundaryData(pos1,pos2)) {
          const peano::datatraversal::autotuning::MethodTrace methodTrace =
              peano::datatraversal::autotuning::UserDefined3;
          const int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().
              parallelise(solvers::RegisteredSolvers.size(), methodTrace);
          pfor(solverNumber, 0, static_cast<int>(solvers::RegisteredSolvers.size()),grainSize)
            auto solver = exahype::solvers::RegisteredSolvers[solverNumber];
            const int cellDescriptionsIndex1 = fineGridVertex.getCellDescriptionsIndex()[pos1Scalar];
            const int cellDescriptionsIndex2 = fineGridVertex.getCellDescriptionsIndex()[pos2Scalar];
            int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
            int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
            assertion4((element1==exahype::solvers::Solver::NotFound &&
                        element2==exahype::solvers::Solver::NotFound)
                       || (element1 >= 0 && element2==exahype::solvers::Solver::NotFound)
                       || (element2 >= 0 && element1==exahype::solvers::Solver::NotFound),
                       cellDescriptionsIndex1,cellDescriptionsIndex2,element1,element2);

            if (element1 >= 0) {
              solver->mergeWithBoundaryData(cellDescriptionsIndex1,element1,pos1,pos2,
                                            _tempFaceUnknowns[solverNumber],
                                            _tempStateSizedVectors[solverNumber],
                                            _tempStateSizedSquareMatrices[solverNumber]);

              #ifdef Debug
              _boundaryFaceMerges++;
              #endif
            }
            if (element2 >= 0){
              solver->mergeWithBoundaryData(cellDescriptionsIndex2,element2,pos2,pos1,
                                            _tempFaceUnknowns[solverNumber],
                                            _tempStateSizedVectors[solverNumber],
                                            _tempStateSizedSquareMatrices[solverNumber]);
              #ifdef Debug
              _boundaryFaceMerges++;
              #endif
            }
          endpfor
          peano::datatraversal::autotuning::Oracle::getInstance()
          .parallelSectionHasTerminated(methodTrace);

          fineGridVertex.setMergePerformed(pos1,pos2,true);
        }
      enddforx
    enddforx
}


//
// Below all methods are nop.
//
//=====================================


void exahype::mappings::SolutionRecomputation::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::SolutionRecomputation::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::SolutionRecomputation::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  // do nothing
  return false;
}

void exahype::mappings::SolutionRecomputation::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::mergeWithMaster(
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
    exahype::State& masterState) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::SolutionRecomputation::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::SolutionRecomputation::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
