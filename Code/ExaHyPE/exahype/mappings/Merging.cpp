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
 
#include "exahype/mappings/Merging.h"

#include "tarch/multicore/Loop.h"

#include "peano/utils/Loop.h"
#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "peano/utils/UserInterface.h"

peano::CommunicationSpecification
exahype::mappings::Merging::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange,
      true);
}

peano::MappingSpecification
exahype::mappings::Merging::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

// Specifications below are all nop.
peano::MappingSpecification
exahype::mappings::Merging::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification
exahype::mappings::Merging::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

peano::MappingSpecification
exahype::mappings::Merging::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}


peano::MappingSpecification
exahype::mappings::Merging::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}


peano::MappingSpecification
exahype::mappings::Merging::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::Merging::_log(
    "exahype::mappings::Merging");

#ifdef Parallel
bool exahype::mappings::Merging::SkipReductionInBatchedTimeSteps = false;
#endif

void exahype::mappings::Merging::prepareTemporaryVariables() {
  if (_tempStateSizedVectors==nullptr) {
    assertion(_tempStateSizedSquareMatrices==nullptr);

    int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
    _tempStateSizedVectors        = new double**[numberOfSolvers];
    _tempStateSizedSquareMatrices = new double**[numberOfSolvers];
    _tempFaceUnknownsArray        = new double* [numberOfSolvers];

    // TOOD(Dominic): Check if we need number of parameters too
    int solverNumber=0;
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      _tempStateSizedVectors[solverNumber] = new double*[5];
      for (int i=0; i<5; ++i) { // see riemanSolverLinear
        _tempStateSizedVectors[solverNumber][i] = new double[solver->getNumberOfVariables()];
      }
      //
      _tempStateSizedSquareMatrices[solverNumber] = new double*[3];
      for (int i=0; i<3; ++i) { // see riemanSolverLinear
        _tempStateSizedSquareMatrices[solverNumber][i] =
            new double[solver->getNumberOfVariables() * solver->getNumberOfVariables()];
      }

      if (solver->getType()==exahype::solvers::Solver::Type::ADER_DG) {
        auto aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
        _tempFaceUnknownsArray[solverNumber] = new double[aderdgSolver->getUnknownsPerFace()];
        // see spaceTimePredictorLinear
      } else {
        _tempFaceUnknownsArray[solverNumber] = nullptr;
      }

      ++solverNumber;
    }
  }
}

void exahype::mappings::Merging::deleteTemporaryVariables() {
  if (_tempStateSizedVectors!=nullptr) {
    assertion(_tempStateSizedSquareMatrices!=nullptr);

    int solverNumber=0;
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      for (int i=0; i<5; ++i) { // see riemanSolverLinear
        delete[] _tempStateSizedVectors[solverNumber][i];
      }
      delete[] _tempStateSizedVectors[solverNumber];
      //
      for (int i=0; i<3; ++i) { // see riemanSolverLinear
        delete[] _tempStateSizedSquareMatrices[solverNumber][i];
      }
      delete[] _tempStateSizedSquareMatrices[solverNumber];


      if (solver->getType()==exahype::solvers::Solver::Type::ADER_DG) {
        delete[] _tempFaceUnknownsArray[solverNumber];
      }
      _tempFaceUnknownsArray[solverNumber] = nullptr;

      ++solverNumber;
    }

    delete[] _tempStateSizedVectors;
    delete[] _tempStateSizedSquareMatrices;
    delete[] _tempFaceUnknownsArray;
    _tempStateSizedVectors        = nullptr;
    _tempStateSizedSquareMatrices = nullptr;
    _tempFaceUnknownsArray        = nullptr;
  }
}

exahype::mappings::Merging::Merging()
  #ifdef Debug
  :
  _interiorFaceSolves(0),
  _boundaryFaceSolves(0)
  #endif
{
  // do nothing
}

exahype::mappings::Merging::~Merging() {
  deleteTemporaryVariables();
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::Merging::Merging(const Merging& masterThread) :
  _localState(masterThread._localState),
  _tempFaceUnknownsArray(nullptr),
  _tempStateSizedVectors(nullptr),
  _tempStateSizedSquareMatrices(nullptr)
  #ifdef Debug
  ,_interiorFaceSolves(0)
  ,_boundaryFaceSolves(0)
  #endif
  {
  prepareTemporaryVariables();
}
#endif

void exahype::mappings::Merging::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  prepareTemporaryVariables();

  _localState = solverState;

  #ifdef Parallel
  exahype::solvers::ADERDGSolver::Heap::getInstance().finishedToSendSynchronousData();
  exahype::solvers::FiniteVolumesSolver::Heap::getInstance().finishedToSendSynchronousData();
  DataHeap::getInstance().finishedToSendSynchronousData();
  MetadataHeap::getInstance().finishedToSendSynchronousData();

  exahype::solvers::ADERDGSolver::Heap::getInstance().startToSendSynchronousData();
  exahype::solvers::FiniteVolumesSolver::Heap::getInstance().startToSendSynchronousData();
  DataHeap::getInstance().startToSendSynchronousData();
  MetadataHeap::getInstance().startToSendSynchronousData();
  #endif

  #ifdef Debug // TODO(Dominic): And not parallel and not shared memory
  _interiorFaceSolves = 0;
  _boundaryFaceSolves = 0;
  #endif

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::Merging::endIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("endIteration(State)", solverState);

  deleteTemporaryVariables();

  #if defined(Debug) // TODO(Dominic): Use logDebug if it works with filters
  logInfo("endIteration(...)","interiorFaceSolves: " << _interiorFaceSolves);
  logInfo("endIteration(...)","boundaryFaceSolves: " << _boundaryFaceSolves);
  #endif

  logTraceOutWith1Argument("endIteration(State)", solverState);
}

void exahype::mappings::Merging::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  logTraceInWith6Arguments("touchVertexFirstTime(...)", fineGridVertex,
                           fineGridX, fineGridH,
                           coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);
  if (_localState.getMergeMode()==exahype::records::State::MergeFaceData ||
      _localState.getMergeMode()==exahype::records::State::BroadcastAndMergeTimeStepDataAndMergeFaceData) {
    dfor2(pos1)
      dfor2(pos2)
        if (fineGridVertex.hasToMergeNeighbours(pos1,pos2)) { // Assumes that we have to valid indices
          const peano::datatraversal::autotuning::MethodTrace methodTrace =
              peano::datatraversal::autotuning::UserDefined4;
          const int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().
              parallelise(solvers::RegisteredSolvers.size(), methodTrace);
          pfor(solverNumber, 0, static_cast<int>(solvers::RegisteredSolvers.size()),grainSize)
            auto solver = exahype::solvers::RegisteredSolvers[solverNumber];
            const int cellDescriptionsIndex1 = fineGridVertex.getCellDescriptionsIndex()[pos1Scalar];
            const int cellDescriptionsIndex2 = fineGridVertex.getCellDescriptionsIndex()[pos2Scalar];
            const int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
            const int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
            assertion4(element2>=0,element2,cellDescriptionsIndex2,solverNumber,fineGridX);
            assertion4(element1>=0,element1,cellDescriptionsIndex1,solverNumber,fineGridX);

            solver->mergeNeighbours(
                cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
                _tempFaceUnknownsArray[solverNumber],
                _tempStateSizedVectors[solverNumber],
                _tempStateSizedSquareMatrices[solverNumber]);

            #ifdef Debug // TODO(Dominic)
            _interiorFaceSolves++;
            #endif
          endpfor
          peano::datatraversal::autotuning::Oracle::getInstance()
          .parallelSectionHasTerminated(methodTrace);

          fineGridVertex.setMergePerformed(pos1,pos2,true);
        }
        if (fineGridVertex.hasToMergeWithBoundaryData(pos1,pos2)) {
          const peano::datatraversal::autotuning::MethodTrace methodTrace =
              peano::datatraversal::autotuning::UserDefined4; // TODO(Dominic): Change UserDefined.
          const int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().
              parallelise(solvers::RegisteredSolvers.size(), methodTrace);
          pfor(solverNumber, 0, static_cast<int>(solvers::RegisteredSolvers.size()),grainSize)
            auto solver = exahype::solvers::RegisteredSolvers[solverNumber];
            const int cellDescriptionsIndex1 = fineGridVertex.getCellDescriptionsIndex()[pos1Scalar];
            const int cellDescriptionsIndex2 = fineGridVertex.getCellDescriptionsIndex()[pos2Scalar];
            int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
            int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
            assertion4((element1 >= 0 && element2==exahype::solvers::Solver::NotFound)
                       || (element2 >= 0 && element1==exahype::solvers::Solver::NotFound),
                       cellDescriptionsIndex1,cellDescriptionsIndex2,element1,element2);

            if (element1 >= 0) {
              solver->mergeWithBoundaryData(cellDescriptionsIndex1,element1,pos1,pos2,
                                            _tempFaceUnknownsArray[solverNumber],
                                            _tempStateSizedVectors[solverNumber],
                                            _tempStateSizedSquareMatrices[solverNumber]);

              #ifdef Debug // TODO(Dominic)
              _boundaryFaceSolves++;
              #endif
            } else { // element2 >= 0
              solver->mergeWithBoundaryData(cellDescriptionsIndex2,element2,pos2,pos1,
                                            _tempFaceUnknownsArray[solverNumber],
                                            _tempStateSizedVectors[solverNumber],
                                            _tempStateSizedSquareMatrices[solverNumber]);
              #ifdef Debug // TODO(Dominic)
              _boundaryFaceSolves++;
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

  logTraceOutWith1Argument("touchVertexFirstTime(...)", fineGridVertex);
}

#ifdef Parallel
///////////////////////////////////////
// NEIGHBOUR
///////////////////////////////////////
// TODO(Dominic): Add to docu: We receive metadata on all vertices.
// We receive only for Cell/Ancestor/Descendants cell descriptions face data.
// EmptyAncestor/EmptyDescendants/InvalidAdjacencyIndices drop face data
// that was sent to them by Cells/Ancestors/Descendants.
void exahype::mappings::Merging::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {

  if (_localState.getMergeMode()==exahype::records::State::MergeFaceData ||
      _localState.getMergeMode()==exahype::records::State::BroadcastAndMergeTimeStepDataAndMergeFaceData) {
    #if !defined(PeriodicBC)
      if (vertex.isBoundary()) return; // TODO(Dominic): Discussion necessary
    #endif

    dfor2(myDest)
      dfor2(mySrc)
        tarch::la::Vector<DIMENSIONS, int> dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest;
        tarch::la::Vector<DIMENSIONS, int> src  = tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;

        int destScalar = TWO_POWER_D - myDestScalar - 1;
        int srcScalar  = TWO_POWER_D - mySrcScalar  - 1;

        if (vertex.hasToReceiveMetadata(src,dest,fromRank)) {
          int receivedMetadataIndex = MetadataHeap::getInstance().
              createData(0,exahype::solvers::RegisteredSolvers.size());
          MetadataHeap::getInstance().receiveData(
              receivedMetadataIndex,
              fromRank, fineGridX, level,
              peano::heap::MessageType::NeighbourCommunication);
          exahype::MetadataHeap::HeapEntries& receivedMetadata = MetadataHeap::getInstance().getData(receivedMetadataIndex);
          assertion(receivedMetadata.size()==solvers::RegisteredSolvers.size());

          if(vertex.hasToMergeWithNeighbourData(src,dest)) {
            mergeWithNeighbourData(
                fromRank,
                vertex.getCellDescriptionsIndex()[srcScalar],
                vertex.getCellDescriptionsIndex()[destScalar],
                src,dest,
                fineGridX,level,
                receivedMetadata);

            vertex.setFaceDataExchangeCountersOfDestination(src,dest,TWO_POWER_D);
            vertex.setMergePerformed(src,dest,true);
          } else {
            dropNeighbourData(
                fromRank,
                vertex.getCellDescriptionsIndex()[srcScalar],
                vertex.getCellDescriptionsIndex()[destScalar],
                src,dest,
                fineGridX,level,
                receivedMetadata);
          }
          // Clean up
          MetadataHeap::getInstance().deleteData(receivedMetadataIndex);
        }
      enddforx
    enddforx
  }
}

void exahype::mappings::Merging::mergeWithNeighbourData(
        const int fromRank,
        const int srcCellDescriptionIndex,
        const int destCellDescriptionIndex,
        const tarch::la::Vector<DIMENSIONS,int>& src,
        const tarch::la::Vector<DIMENSIONS,int>& dest,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int level,
        const exahype::MetadataHeap::HeapEntries& receivedMetadata) {
  int solverNumber = 0;
  for(auto solver : solvers::RegisteredSolvers) {
    if (receivedMetadata[solverNumber].getU()!=exahype::Vertex::InvalidMetadataEntry) {
      const int element = solver->tryGetElement(destCellDescriptionIndex,solverNumber);
      assertion1(element>=0,element);

      logDebug(
          "mergeWithNeighbour(...)", "receive data for solver " << solverNumber << " from " <<
          fromRank << " at vertex x=" << x << ", level=" << level <<
          ", src=" << src << ", dest=" << dest);

      solver->mergeWithNeighbourData(
          fromRank,receivedMetadata[solverNumber].getU(),
          destCellDescriptionIndex,element,src,dest,
          _tempFaceUnknownsArray[solverNumber],
          _tempStateSizedVectors[solverNumber],
          _tempStateSizedSquareMatrices[solverNumber],
          x,level);
    } else {
      logDebug(
            "mergeWithNeighbour(...)", "drop data for solver " << solverNumber << " from " <<
            fromRank << " at vertex x=" << x << ", level=" << level <<
            ", src=" << src << ", dest=" << dest);

      solver->dropNeighbourData(
          fromRank,src,dest,x,level);
    }
    ++solverNumber;
  }
}

void exahype::mappings::Merging::dropNeighbourData(
    const int fromRank,
    const int srcCellDescriptionIndex,
    const int destCellDescriptionIndex,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int level,
    const exahype::MetadataHeap::HeapEntries& receivedMetadata) {
  assertion(receivedMetadata.size()==solvers::RegisteredSolvers.size());

  for(auto solver : solvers::RegisteredSolvers) {
    logDebug(
        "dropNeighbourData(...)", "drop data from " <<
        fromRank << " at vertex x=" << x << ", level=" << level <<
        ", src=" << src << ", dest=" << dest);

    solver->dropNeighbourData(fromRank,src,dest,x,level);
  }
}

///////////////////////////////////////
// MASTER->WORKER
///////////////////////////////////////
bool exahype::mappings::Merging::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  bool workerHasToSendFaceData = false;

  if (_localState.getMergeMode()==exahype::records::State::MergeMode::BroadcastAndMergeTimeStepData ||
      _localState.getMergeMode()==exahype::records::State::MergeMode::BroadcastAndMergeTimeStepDataAndMergeFaceData) {
    for (auto& p : exahype::solvers::RegisteredSolvers) {
      p->sendDataToWorker(
          worker,
          fineGridVerticesEnumerator.getCellCenter(),
          fineGridVerticesEnumerator.getLevel());
    }
  }

  if ((_localState.getMergeMode()==exahype::records::State::MergeMode::MergeFaceData ||
      _localState.getMergeMode()==exahype::records::State::MergeMode::BroadcastAndMergeTimeStepDataAndMergeFaceData)
      && fineGridCell.isInside()) {
    if (fineGridCell.isInitialised()) {
      exahype::Vertex::sendEncodedMetadata(
          worker,fineGridCell.getCellDescriptionsIndex(),
          peano::heap::MessageType::MasterWorkerCommunication,
          fineGridVerticesEnumerator.getCellCenter(),
          fineGridVerticesEnumerator.getLevel());

      int solverNumber=0;
      for (auto solver : exahype::solvers::RegisteredSolvers) {
        int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);

        if (element!=exahype::solvers::Solver::NotFound) {
          workerHasToSendFaceData |=
              solver->sendDataToWorker(
                  worker,
                  fineGridCell.getCellDescriptionsIndex(),element,
                  fineGridVerticesEnumerator.getCellCenter(),
                  fineGridVerticesEnumerator.getLevel());
        } else {
          solver->sendEmptyDataToWorker(
              worker,
              fineGridVerticesEnumerator.getCellCenter(),
              fineGridVerticesEnumerator.getLevel());
        }

        ++solverNumber;
      }
    } else {
      for (auto solver : exahype::solvers::RegisteredSolvers) {
        solver->sendEmptyDataToWorker(
            worker,
            fineGridVerticesEnumerator.getCellCenter(),
            fineGridVerticesEnumerator.getLevel());
      }
    } // else do nothing
  }

  if (
      (!peano::parallel::loadbalancing::Oracle::getInstance().isLoadBalancingActivated()
      &&
      (!workerHasToSendFaceData)
      &&
      SkipReductionInBatchedTimeSteps)
  ) {
    return false;
  }
  else return true;
}

void exahype::mappings::Merging::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  if (_localState.getMergeMode()==exahype::records::State::MergeMode::BroadcastAndMergeTimeStepData ||
      _localState.getMergeMode()==exahype::records::State::MergeMode::BroadcastAndMergeTimeStepDataAndMergeFaceData) {
    for (auto& solver : exahype::solvers::RegisteredSolvers) {
      solver->mergeWithMasterData(
          tarch::parallel::NodePool::getInstance().getMasterRank(),
          receivedVerticesEnumerator.getCellCenter(),
          receivedVerticesEnumerator.getLevel());
    }
  }

  if ((_localState.getMergeMode()==exahype::records::State::MergeMode::MergeFaceData ||
      _localState.getMergeMode()==exahype::records::State::MergeMode::BroadcastAndMergeTimeStepDataAndMergeFaceData)
      && receivedCell.isInside()) {
    if (receivedCell.isInitialised()) {
      int receivedMetadataIndex = MetadataHeap::getInstance().createData(
          0,exahype::solvers::RegisteredSolvers.size());
      MetadataHeap::getInstance().receiveData(
          receivedMetadataIndex,
          tarch::parallel::NodePool::getInstance().getMasterRank(),
          receivedVerticesEnumerator.getCellCenter(),
          receivedVerticesEnumerator.getLevel(),
          peano::heap::MessageType::MasterWorkerCommunication);
      MetadataHeap::HeapEntries& receivedMetadata =
                MetadataHeap::getInstance().getData(receivedMetadataIndex);

      int solverNumber=0;
      for (auto solver : exahype::solvers::RegisteredSolvers) {
        int element = solver->tryGetElement(receivedCell.getCellDescriptionsIndex(),solverNumber);

        if (element!=exahype::solvers::Solver::NotFound &&
            receivedMetadata[solverNumber].getU()!=exahype::Vertex::InvalidMetadataEntry) {
          solver->mergeWithMasterData(
                    tarch::parallel::NodePool::getInstance().getMasterRank(),
                    receivedMetadata[solverNumber].getU(),
                    receivedCell.getCellDescriptionsIndex(),element,
                    receivedVerticesEnumerator.getCellCenter(),
                    receivedVerticesEnumerator.getLevel());
        } else {
          solver->dropMasterData(
              tarch::parallel::NodePool::getInstance().getMasterRank(),
              receivedVerticesEnumerator.getCellCenter(),
              receivedVerticesEnumerator.getLevel());
        }

        ++solverNumber;
      }
      MetadataHeap::getInstance().deleteData(receivedMetadataIndex);

    } else {
      exahype::Vertex::dropMetadata(
          tarch::parallel::NodePool::getInstance().getMasterRank(),
          peano::heap::MessageType::MasterWorkerCommunication,
          receivedVerticesEnumerator.getCellCenter(),
          receivedVerticesEnumerator.getLevel());

      for (auto solver : exahype::solvers::RegisteredSolvers) {
        solver->dropMasterData(
            tarch::parallel::NodePool::getInstance().getMasterRank(),
            receivedVerticesEnumerator.getCellCenter(),
            receivedVerticesEnumerator.getLevel());
      }
    } // else do nothing
  }
}


//
// Below all methods are nop.
//
// ====================================


void exahype::mappings::Merging::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // TODO(Dominic): 1. Copy the received data to the local cell (Maybe just "swap" the indices and delete previous local cell data).
}

void exahype::mappings::Merging::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Merging::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Merging::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Merging::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Merging::mergeWithMaster(
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

void exahype::mappings::Merging::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Merging::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Merging::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

#endif

#if defined(SharedMemoryParallelisation)
void exahype::mappings::Merging::mergeWithWorkerThread(
    const Merging& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::Merging::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Merging::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Merging::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Merging::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Merging::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Merging::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Merging::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Merging::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Merging::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Merging::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Merging::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::Merging::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
