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
 *
 * \author Dominic E. Charrier, Tobias Weinzierl, Jean-Matthieu Gallard, Fabian Güra, Leonhard Rannabauer
 **/
#include "exahype/solvers/ADERDGSolver.h"

#include <limits>
#include <iomanip>

#include <algorithm>

#include "exahype/Cell.h"
#include "exahype/Vertex.h"
#include "exahype/VertexOperations.h"

#include "tarch/la/VectorVectorOperations.h"
#include "tarch/multicore/Lock.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/amr/AdaptiveMeshRefinement.h"

#include "peano/heap/CompressedFloatingPointNumbers.h"
#include "peano/datatraversal/TaskSet.h"

#include "peano/grid/aspects/VertexStateAnalysis.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "kernels/KernelUtils.h"

#include "tarch/multicore/Jobs.h"

namespace {
  constexpr const char* tags[]{"solutionUpdate",
                             "volumeIntegral",
                             "surfaceIntegral",
                             "riemannSolver",
                             "spaceTimePredictor",
                             "stableTimeStepSize",
                             "solutionAdjustment",
                             "faceUnknownsProlongation",
                             "faceUnknownsRestriction",
                             "volumeUnknownsProlongation",
                             "volumeUnknownsRestriction",
                             "boundaryConditions",
                             "deltaDistribution"
                             };
}

tarch::logging::Log exahype::solvers::ADERDGSolver::_log( "exahype::solvers::ADERDGSolver");

// communication status
int exahype::solvers::ADERDGSolver::CellCommunicationStatus                             = 2;
int exahype::solvers::ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication = 1;
// augmentation status
// On-the fly erasing seems to work with those values
int exahype::solvers::ADERDGSolver::MaximumAugmentationStatus                   = 4;
int exahype::solvers::ADERDGSolver::MinimumAugmentationStatusForVirtualRefining = 3;
int exahype::solvers::ADERDGSolver::MinimumAugmentationStatusForRefining        = 3;

tarch::multicore::BooleanSemaphore exahype::solvers::ADERDGSolver::RestrictionSemaphore;

tarch::multicore::BooleanSemaphore exahype::solvers::ADERDGSolver::CoarseGridSemaphore;

void exahype::solvers::ADERDGSolver::addNewCellDescription(
  const int                                     cellDescriptionsIndex,
  const int                                     solverNumber,
  const CellDescription::Type                   cellType,
  const CellDescription::RefinementEvent        refinementEvent,
  const int                                     level,
  const int                                     parentIndex,
  const tarch::la::Vector<DIMENSIONS, double>&  cellSize,
  const tarch::la::Vector<DIMENSIONS, double>&  cellOffset) {
  
  logDebug("addNewCellDescription(...)","Add cell description: index="<<cellDescriptionsIndex<<", type="<<CellDescription::toString(cellType) <<", level="<<level<<", parentIndex="<<parentIndex
            << " for solver=" << solverNumber);

  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion2(parentIndex == -1 || parentIndex != cellDescriptionsIndex, parentIndex, cellDescriptionsIndex);
  assertion2(parentIndex != cellDescriptionsIndex, parentIndex, cellDescriptionsIndex);

  assertion2(static_cast<unsigned int>(solverNumber) < solvers::RegisteredSolvers.size(),
             solverNumber, exahype::solvers::RegisteredSolvers.size());

  CellDescription newCellDescription;
  newCellDescription.setSolverNumber(solverNumber);

  // Default AMR settings
  newCellDescription.setType(cellType);
  newCellDescription.setParentIndex(parentIndex);
  newCellDescription.setLevel(level);
  newCellDescription.setParentCellLevel(-1);
  newCellDescription.setParentOffset(-1);
  newCellDescription.setRefinementEvent(refinementEvent);
  newCellDescription.setRefinementRequest(CellDescription::RefinementRequest::Pending);

  newCellDescription.setHasVirtualChildren(false);
  newCellDescription.setAugmentationStatus(0);
  newCellDescription.setPreviousAugmentationStatus(0);
  if (cellType==CellDescription::Type::Cell) {
    newCellDescription.setPreviousAugmentationStatus(MaximumAugmentationStatus);
  }
  newCellDescription.setFacewiseAugmentationStatus(0); // implicit conversion
  newCellDescription.setCommunicationStatus(0);
  newCellDescription.setFacewiseCommunicationStatus(0); // implicit conversion
  if (cellType==CellDescription::Type::Cell) {
    newCellDescription.setCommunicationStatus(CellCommunicationStatus);
    newCellDescription.setFacewiseCommunicationStatus(CellCommunicationStatus); // implicit conversion
    // TODO(Dominic): Make sure prolongation and restriction considers this.
  }

  std::bitset<DIMENSIONS_TIMES_TWO> neighbourMergePerformed;  // default construction: no bit set
  newCellDescription.setNeighbourMergePerformed(neighbourMergePerformed);

  // Pass geometry information to the cellDescription description
  newCellDescription.setSize(cellSize);
  newCellDescription.setOffset(cellOffset);

  // Initialise MPI helper variables
  #ifdef Parallel
  newCellDescription.setHasToHoldDataForMasterWorkerCommunication(false);
  for (int faceIndex = 0; faceIndex < DIMENSIONS_TIMES_TWO; faceIndex++) {
    newCellDescription.setFaceDataExchangeCounter(faceIndex,TWO_POWER_D);
  }
  #endif

  // Default field data indices
  newCellDescription.setSolution(-1);
  newCellDescription.setPreviousSolution(-1);
  newCellDescription.setUpdate(-1);
  newCellDescription.setExtrapolatedPredictor(-1);
  newCellDescription.setFluctuation(-1);

  // Limiter meta data (oscillations identificator)
  newCellDescription.setLimiterStatus(0); // 0 is CellDescription::LimiterStatus::Ok
  newCellDescription.setPreviousLimiterStatus(0);
  newCellDescription.setFacewiseLimiterStatus(0);  // implicit conversion
  newCellDescription.setSolutionMin(-1);
  newCellDescription.setSolutionMax(-1);
  newCellDescription.setIterationsToCureTroubledCell(0);

  // Compression
  newCellDescription.setCompressionState(CellDescription::CompressionState::Uncompressed);
  newCellDescription.setSolutionAverages(-1);
  newCellDescription.setPreviousSolutionAverages(-1);
  newCellDescription.setUpdateAverages(-1);
  newCellDescription.setExtrapolatedPredictorAverages(-1);
  newCellDescription.setFluctuationAverages(-1);

  newCellDescription.setSolutionCompressed(-1);
  newCellDescription.setPreviousSolutionAverages(-1);
  newCellDescription.setUpdateCompressed(-1);
  newCellDescription.setExtrapolatedPredictorCompressed(-1);
  newCellDescription.setFluctuationCompressed(-1);

  newCellDescription.setBytesPerDoFInExtrapolatedPredictor(-1);
  newCellDescription.setBytesPerDoFInFluctuation(-1);
  newCellDescription.setBytesPerDoFInPreviousSolution(-1);
  newCellDescription.setBytesPerDoFInSolution(-1);
  newCellDescription.setBytesPerDoFInUpdate(-1);

  tarch::multicore::Lock lock(exahype::HeapSemaphore);
  ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex).push_back(newCellDescription);
  lock.free();
}

/**
 * Returns the ADERDGCellDescription heap vector
 * at address \p cellDescriptionsIndex.
 */
exahype::solvers::ADERDGSolver::Heap::HeapEntries& exahype::solvers::ADERDGSolver::getCellDescriptions(
    const int cellDescriptionsIndex) {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);

  return Heap::getInstance().getData(cellDescriptionsIndex);
}

/**
 * Returns the ADERDGCellDescription with index \p element
 * in the heap vector at address \p cellDescriptionsIndex.
 */
exahype::solvers::ADERDGSolver::CellDescription& exahype::solvers::ADERDGSolver::getCellDescription(
    const int cellDescriptionsIndex,
    const int element) {
  assertion2(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex,element);
  assertion2(element>=0,cellDescriptionsIndex,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),cellDescriptionsIndex,element);

  return Heap::getInstance().getData(cellDescriptionsIndex)[element];
}

/**
 * Returns if a ADERDGCellDescription type holds face data.
 */
bool exahype::solvers::ADERDGSolver::holdsFaceData(const CellDescription& cellDescription) {
  assertion1(cellDescription.getType()!=CellDescription::Type::Cell ||
            cellDescription.getCommunicationStatus()==CellCommunicationStatus,cellDescription.toString());
  return
      cellDescription.getType()!=CellDescription::Type::Ancestor &&
      (
        cellDescription.getCommunicationStatus()>=MinimumCommunicationStatusForNeighbourCommunication 
        #ifdef Parallel
        || cellDescription.getHasToHoldDataForMasterWorkerCommunication()
        #endif
      );
}

void exahype::solvers::ADERDGSolver::ensureNoUnnecessaryMemoryIsAllocated(
    CellDescription& cellDescription) const {

  if (
      cellDescription.getType()!=CellDescription::Type::Cell &&
      DataHeap::getInstance().isValidIndex(cellDescription.getSolution())
  ) {
    tarch::multicore::Lock lock(exahype::HeapSemaphore);

    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolution()));
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolution()));

    if (cellDescription.getSolution()>=0) {
      DataHeap::getInstance().deleteData(cellDescription.getSolution());
      assertion(cellDescription.getSolutionCompressed()==-1);
    }
    else {
      assertion(CompressionAccuracy>0.0);
      assertion(cellDescription.getSolution()==-1);
      CompressedDataHeap::getInstance().deleteData(cellDescription.getSolutionCompressed());
    }

    DataHeap::getInstance().deleteData(cellDescription.getSolutionAverages());
    DataHeap::getInstance().deleteData(cellDescription.getPreviousSolutionAverages());

    cellDescription.setPreviousSolution(-1);
    cellDescription.setSolution(-1);

    cellDescription.setPreviousSolutionAverages(-1);
    cellDescription.setSolutionAverages(-1);

    cellDescription.setPreviousSolutionCompressed(-1);
    cellDescription.setSolutionCompressed(-1);

    lock.free();
  }

  // deallocate update and boundary arrays
  if (
      !holdsFaceData(cellDescription) &&
      DataHeap::getInstance().isValidIndex(cellDescription.getUpdate())
  ) {
    // update
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getUpdate()));
    if ( cellDescription.getUpdate()>=0 ) {
      assertion(cellDescription.getUpdateCompressed()==-1);

      DataHeap::getInstance().deleteData(cellDescription.getUpdate());
      cellDescription.setUpdate(-1);
    }
    else {
      assertion(CompressionAccuracy>0.0);
      assertion(cellDescription.getUpdate()==-1);

      CompressedDataHeap::getInstance().deleteData(cellDescription.getUpdateCompressed());
      cellDescription.setUpdateCompressed(-1);
    }
    DataHeap::getInstance().deleteData(cellDescription.getUpdateAverages());
    cellDescription.setUpdateAverages(-1);

    // extrapolated predictor
    tarch::multicore::Lock lock(exahype::HeapSemaphore);
    if ( cellDescription.getExtrapolatedPredictor()>=0 ) {
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictor()));
      assertion(cellDescription.getExtrapolatedPredictorCompressed()==-1);

      DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictor());
      cellDescription.setExtrapolatedPredictor(-1);
    }
    else {
      assertion(CompressionAccuracy>0.0);
      assertion(cellDescription.getExtrapolatedPredictor()==-1);

      CompressedDataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictorCompressed());
      cellDescription.setExtrapolatedPredictorCompressed(-1);
    }
    DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictorAverages());
    cellDescription.setExtrapolatedPredictorAverages(-1);

    // fluctuations
    if ( cellDescription.getFluctuation()>=0 ) {
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()));
      assertion(cellDescription.getFluctuationCompressed()==-1);

      DataHeap::getInstance().deleteData(cellDescription.getFluctuation());
      cellDescription.setFluctuation(-1);
    }
    else {
      assertion(CompressionAccuracy>0.0);
      assertion(cellDescription.getFluctuation()==-1);

      CompressedDataHeap::getInstance().deleteData(cellDescription.getFluctuationCompressed());
      cellDescription.setFluctuationCompressed(-1);
    }
    DataHeap::getInstance().deleteData(cellDescription.getFluctuationAverages());
    cellDescription.setFluctuationAverages(-1);

    if ( getDMPObservables()>0 ) {
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionMin()));
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionMax()));
      DataHeap::getInstance().deleteData(cellDescription.getSolutionMin());
      DataHeap::getInstance().deleteData(cellDescription.getSolutionMax());

      cellDescription.setSolutionMin(-1);
      cellDescription.setSolutionMax(-1);
    }

    lock.free();
  }
}

void exahype::solvers::ADERDGSolver::ensureNecessaryMemoryIsAllocated(
    CellDescription& cellDescription) const {
  // allocate solution
  if (
      cellDescription.getType()==CellDescription::Type::Cell &&
      !DataHeap::getInstance().isValidIndex(cellDescription.getSolution())
  ) {
    assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolution()));

    tarch::multicore::Lock lock(exahype::HeapSemaphore);
    // Allocate volume DoF for limiter
    const int dataPerNode     = getNumberOfVariables()+getNumberOfParameters();
    const int dataPerCell     = getDataPerCell(); // Only the solution and previousSolution store material parameters
    cellDescription.setPreviousSolution( DataHeap::getInstance().createData( dataPerCell, dataPerCell ) );
    cellDescription.setSolution( DataHeap::getInstance().createData( dataPerCell, dataPerCell ) );
    
    cellDescription.setSolutionCompressed(-1);
    cellDescription.setPreviousSolutionCompressed(-1);

    cellDescription.setPreviousSolutionAverages( DataHeap::getInstance().createData( dataPerNode, dataPerNode ) );
    cellDescription.setSolutionAverages(         DataHeap::getInstance().createData( dataPerNode, dataPerNode ) );

    cellDescription.setCompressionState(CellDescription::Uncompressed);

    lock.free();
  }

  // allocate update and boundary arrays
  if (
      holdsFaceData(cellDescription) &&
      !DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictor())
  ) {
    assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()));

    tarch::multicore::Lock lock(exahype::HeapSemaphore);

    // allocate update dof
    cellDescription.setUpdate( DataHeap::getInstance().createData( getUpdateSize(), getUpdateSize() ) );
    cellDescription.setUpdateAverages( DataHeap::getInstance().createData( getNumberOfVariables(), getNumberOfVariables() ) );
    cellDescription.setUpdateCompressed(-1);

    // extrapolated predictor
    const int dataPerBnd = getBndTotalSize();
    cellDescription.setExtrapolatedPredictor( DataHeap::getInstance().createData(dataPerBnd, dataPerBnd) );
    cellDescription.setExtrapolatedPredictorCompressed(-1);
    const int boundaryData     = (getNumberOfParameters()+getNumberOfVariables()) * DIMENSIONS_TIMES_TWO; //TODO JMG / Dominic adapt for padding with optimized kernels //TODO Tobias: Does it make sense to pad these arrays.
    cellDescription.setExtrapolatedPredictorAverages( DataHeap::getInstance().createData( boundaryData,  boundaryData  ) );

    // fluctuations
    const int dofPerBnd  = getBndFluxTotalSize();
    cellDescription.setFluctuation( DataHeap::getInstance().createData(dofPerBnd,  dofPerBnd) );
    cellDescription.setFluctuationCompressed(-1);
    const int boundaryUnknowns = getNumberOfVariables() * DIMENSIONS_TIMES_TWO;     //TODO JMG / Dominic adapt for padding with optimized kernels //TODO Tobias: Does it make sense to pad these arrays.
    cellDescription.setFluctuationAverages( DataHeap::getInstance().createData( boundaryUnknowns, boundaryUnknowns ) );

    // Allocate volume DoF for limiter (we need for every of the 2*DIMENSIONS faces an array of min values
    // and array of max values of the neighbour at this face).
    const int numberOfObservables = getDMPObservables();
    if ( numberOfObservables>0 ) {
      cellDescription.setSolutionMin(DataHeap::getInstance().createData(
          numberOfObservables * DIMENSIONS_TIMES_TWO, numberOfObservables * DIMENSIONS_TIMES_TWO ));
      cellDescription.setSolutionMax(DataHeap::getInstance().createData(
          numberOfObservables * DIMENSIONS_TIMES_TWO, numberOfObservables * DIMENSIONS_TIMES_TWO ));

      for (int i=0; i<numberOfObservables * DIMENSIONS_TIMES_TWO; i++) {
        DataHeap::getInstance().getData( cellDescription.getSolutionMin() )[i] = std::numeric_limits<double>::max();
        DataHeap::getInstance().getData( cellDescription.getSolutionMax() )[i] = -std::numeric_limits<double>::max();
      }
    }

    lock.free();
  }
}

void exahype::solvers::ADERDGSolver::eraseCellDescriptions(
    const int cellDescriptionsIndex) {
  assertion(Heap::getInstance().isValidIndex(cellDescriptionsIndex));
  for (auto& p : Heap::getInstance().getData(cellDescriptionsIndex)) {
    auto *solver = exahype::solvers::RegisteredSolvers[p.getSolverNumber()];

    ADERDGSolver* aderdgSolver = nullptr;
    if (solver->getType()==Solver::Type::ADERDG) {
      aderdgSolver = static_cast<ADERDGSolver*>(solver);
    }
    else if (solver->getType()==Solver::Type::LimitingADERDG) {
      aderdgSolver =
          static_cast<LimitingADERDGSolver*>(solver)->getSolver().get();
    }
    assertion(aderdgSolver!=nullptr);

    p.setType(CellDescription::Type::Erased);
    aderdgSolver->ensureNoUnnecessaryMemoryIsAllocated(p);
  }

  Heap::getInstance().getData(cellDescriptionsIndex).clear();
}

exahype::solvers::ADERDGSolver::ADERDGSolver(
    const std::string& identifier, int numberOfVariables,
    int numberOfParameters, int DOFPerCoordinateAxis,
    double maximumMeshSize, int maximumAdaptiveMeshDepth,
    int DMPObservables,
    int limiterHelperLayers,
    exahype::solvers::Solver::TimeStepping timeStepping,
    std::unique_ptr<profilers::Profiler> profiler)
    : Solver(identifier, Solver::Type::ADERDG, numberOfVariables,
             numberOfParameters, DOFPerCoordinateAxis,
             maximumMeshSize, maximumAdaptiveMeshDepth,
             timeStepping, std::move(profiler)),
     _previousMinCorrectorTimeStamp( std::numeric_limits<double>::max() ),
     _previousMinCorrectorTimeStepSize( std::numeric_limits<double>::max() ),
     _minCorrectorTimeStamp( std::numeric_limits<double>::max() ),
     _minCorrectorTimeStepSize( std::numeric_limits<double>::max() ),
     _minPredictorTimeStamp( std::numeric_limits<double>::max() ),
     _minPredictorTimeStepSize( std::numeric_limits<double>::max() ),
     _minNextTimeStepSize( std::numeric_limits<double>::max() ),
     _stabilityConditionWasViolated( false ),
     _DMPObservables(DMPObservables),
     _minimumLimiterStatusForPassiveFVPatch(2),
     _minimumLimiterStatusForActiveFVPatch(limiterHelperLayers+_minimumLimiterStatusForPassiveFVPatch),
     _minimumLimiterStatusForTroubledCell (2*limiterHelperLayers+_minimumLimiterStatusForPassiveFVPatch) {

  // register tags with profiler
  for (const char* tag : tags) {
    _profiler->registerTag(tag);
  }

  #ifdef Parallel
  _invalidExtrapolatedPredictor.resize(getBndFaceSize());
  _invalidFluctuations.resize(getBndFluxSize());
  std::fill_n(_invalidExtrapolatedPredictor.data(),_invalidExtrapolatedPredictor.size(),-1);
  std::fill_n(_invalidFluctuations.data(),_invalidFluctuations.size(),-1);

  _receivedExtrapolatedPredictor.resize(getBndFaceSize());
  _receivedFluctuations.resize(getBndFluxSize());

  _receivedUpdate.reserve(getUpdateSize());
  #endif
}

int exahype::solvers::ADERDGSolver::getUnknownsPerFace() const {
  return _numberOfVariables * power(_nodesPerCoordinateAxis, DIMENSIONS - 1);
}

int exahype::solvers::ADERDGSolver::getUnknownsPerCellBoundary() const {
  return DIMENSIONS_TIMES_TWO * getUnknownsPerFace();
}

int exahype::solvers::ADERDGSolver::getUnknownsPerCell() const {
  return _numberOfVariables * power(_nodesPerCoordinateAxis, DIMENSIONS + 0);
}

int exahype::solvers::ADERDGSolver::getFluxUnknownsPerCell() const {
  return (DIMENSIONS + 1) * getUnknownsPerCell(); // +1 for sources
}

int exahype::solvers::ADERDGSolver::getSpaceTimeUnknownsPerCell() const {
  return _numberOfVariables * power(_nodesPerCoordinateAxis, DIMENSIONS + 1);
}

int exahype::solvers::ADERDGSolver::getSpaceTimeFluxUnknownsPerCell() const {
  return (DIMENSIONS + 1) * getSpaceTimeUnknownsPerCell();  // +1 for sources
}

int exahype::solvers::ADERDGSolver::getDataPerFace() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS - 1);
}

int exahype::solvers::ADERDGSolver::getDataPerCellBoundary() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS - 1) * DIMENSIONS_TIMES_TWO;
}

int exahype::solvers::ADERDGSolver::getDataPerCell() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS + 0);
}

int exahype::solvers::ADERDGSolver::getSpaceTimeDataPerCell() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS + 1);
}

int exahype::solvers::ADERDGSolver::getDMPObservables() const {
  return _DMPObservables;
}

int exahype::solvers::ADERDGSolver::getMinimumLimiterStatusForActiveFVPatch() const {
  return _minimumLimiterStatusForActiveFVPatch;
}

int exahype::solvers::ADERDGSolver::getMinimumLimiterStatusForTroubledCell() const {
  return _minimumLimiterStatusForTroubledCell;
}

void exahype::solvers::ADERDGSolver::synchroniseTimeStepping(
    CellDescription& p) const {
  switch (_timeStepping) {
    case TimeStepping::Global:
      p.setPreviousCorrectorTimeStamp(_previousMinCorrectorTimeStamp);
      p.setPreviousCorrectorTimeStepSize(_previousMinCorrectorTimeStepSize);

      p.setCorrectorTimeStamp(_minCorrectorTimeStamp);
      p.setCorrectorTimeStepSize(_minCorrectorTimeStepSize);

      p.setPredictorTimeStamp(_minPredictorTimeStamp);
      p.setPredictorTimeStepSize(_minPredictorTimeStepSize);
      break;
    case TimeStepping::GlobalFixed:
      p.setPreviousCorrectorTimeStamp(_previousMinCorrectorTimeStamp);
      p.setPreviousCorrectorTimeStepSize(_previousMinCorrectorTimeStepSize);

      p.setCorrectorTimeStamp(_minCorrectorTimeStamp);
      p.setCorrectorTimeStepSize(_minCorrectorTimeStepSize);

      p.setPredictorTimeStamp(_minPredictorTimeStamp);
      p.setPredictorTimeStepSize(_minPredictorTimeStepSize);
      break;
  }
}

void exahype::solvers::ADERDGSolver::synchroniseTimeStepping(
      const int cellDescriptionsIndex,
      const int element) const {
  synchroniseTimeStepping(Heap::getInstance().getData(cellDescriptionsIndex)[element]);
}

void exahype::solvers::ADERDGSolver::startNewTimeStep() {
  switch (_timeStepping) {
    case TimeStepping::Global:
      // n-1
      _previousMinCorrectorTimeStamp    = _minCorrectorTimeStamp;
      _previousMinCorrectorTimeStepSize = _minCorrectorTimeStepSize;
      // n
      _minCorrectorTimeStepSize = _minNextTimeStepSize;
      _minCorrectorTimeStamp    = _minCorrectorTimeStamp+_minCorrectorTimeStepSize;

      _minPredictorTimeStepSize = _minCorrectorTimeStepSize;
      _minPredictorTimeStamp    = _minCorrectorTimeStamp;

      _minNextTimeStepSize = std::numeric_limits<double>::max();
      break;
    case TimeStepping::GlobalFixed:
      // n-1
      _previousMinCorrectorTimeStamp    = _minCorrectorTimeStamp;
      _previousMinCorrectorTimeStepSize = _minCorrectorTimeStepSize;
      // n
      _minCorrectorTimeStepSize = _minNextTimeStepSize;
      _minCorrectorTimeStamp    = _minCorrectorTimeStamp+_minNextTimeStepSize;

      _minPredictorTimeStepSize = _minCorrectorTimeStepSize;
      _minPredictorTimeStamp    = _minCorrectorTimeStamp;
      break;
  }

  _maxLevel     = _nextMaxLevel;
  _nextMaxLevel = -std::numeric_limits<int>::max(); // "-", min
}

void exahype::solvers::ADERDGSolver::startNewTimeStepFused(
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch) {
  // n-1
  if ( isFirstIterationOfBatch ) {
    _previousMinCorrectorTimeStamp    = _minCorrectorTimeStamp;
    _previousMinCorrectorTimeStepSize = _minCorrectorTimeStepSize;
  }
  // n
  _minCorrectorTimeStamp    = _minPredictorTimeStamp;
  _minCorrectorTimeStepSize = _minPredictorTimeStepSize;
  // n+1
  _minPredictorTimeStamp    = _minPredictorTimeStamp + _minPredictorTimeStepSize;
  if ( isLastIterationOfBatch ) {
    // TODO(Dominic): Add to docu. We minimise the time step size over all batch iterations
    // Otherwise, we freeze (do not overwrite) the minPredictorTimeStepSize
    switch (_timeStepping) {
      case TimeStepping::Global:
        _minPredictorTimeStepSize     = _minNextTimeStepSize;
        _minNextTimeStepSize = std::numeric_limits<double>::max();
        break;
      case TimeStepping::GlobalFixed:
        _minPredictorTimeStepSize = _minNextTimeStepSize;
        break;
    }

    _maxLevel     = _nextMaxLevel;
    _nextMaxLevel = -std::numeric_limits<int>::max(); // "-", min
  }
}

void exahype::solvers::ADERDGSolver::updateTimeStepSizesFused() {
  switch (_timeStepping) {
  case TimeStepping::Global:
    _minCorrectorTimeStepSize = _minNextTimeStepSize;
    _minPredictorTimeStepSize = _minNextTimeStepSize;

    _minPredictorTimeStamp    =  _minCorrectorTimeStamp+_minNextTimeStepSize;

    _minNextTimeStepSize = std::numeric_limits<double>::max();
    break;
  case TimeStepping::GlobalFixed:
    _minCorrectorTimeStepSize = _minNextTimeStepSize;
    _minPredictorTimeStepSize = _minNextTimeStepSize;

    _minPredictorTimeStamp =  _minCorrectorTimeStamp+_minNextTimeStepSize;
    break;
  }

  _stabilityConditionWasViolated = false;

  _maxLevel     = _nextMaxLevel;
  _nextMaxLevel = -std::numeric_limits<int>::max(); // "-", min
}

void exahype::solvers::ADERDGSolver::updateTimeStepSizes() {
  switch (_timeStepping) {
    case TimeStepping::Global:
      _minCorrectorTimeStepSize = _minNextTimeStepSize;
      _minPredictorTimeStepSize = _minNextTimeStepSize;

      _minPredictorTimeStamp    =  _minCorrectorTimeStamp;

      _minNextTimeStepSize = std::numeric_limits<double>::max();
      break;
    case TimeStepping::GlobalFixed:
      _minCorrectorTimeStepSize = _minNextTimeStepSize;
      _minPredictorTimeStepSize = _minNextTimeStepSize;

      _minPredictorTimeStamp =  _minCorrectorTimeStamp;
      break;
  }

  _maxLevel     = _nextMaxLevel;
  _nextMaxLevel = -std::numeric_limits<int>::max(); // "-", min
}

void exahype::solvers::ADERDGSolver::zeroTimeStepSizes() {
  _previousMinCorrectorTimeStepSize = 0;
  _minCorrectorTimeStepSize         = 0;
  _minPredictorTimeStepSize         = 0;

  _minPredictorTimeStamp = _minCorrectorTimeStamp;
}

void exahype::solvers::ADERDGSolver::rollbackToPreviousTimeStep() {
  switch (_timeStepping) {
    case TimeStepping::Global:
      _minNextTimeStepSize                     = std::numeric_limits<double>::max();

      _minPredictorTimeStamp                    = _previousMinCorrectorTimeStamp;
      _minPredictorTimeStepSize                 = _previousMinCorrectorTimeStepSize;

      _minCorrectorTimeStamp                    = _previousMinCorrectorTimeStamp;
      _minCorrectorTimeStepSize                 = _previousMinCorrectorTimeStepSize;

      _previousMinCorrectorTimeStamp            = std::numeric_limits<double>::max();
      _previousMinCorrectorTimeStepSize         = std::numeric_limits<double>::max();
      break;
    case TimeStepping::GlobalFixed:
      _minPredictorTimeStamp                    = _previousMinCorrectorTimeStamp;
      _minPredictorTimeStepSize                 = _previousMinCorrectorTimeStepSize;

      _minCorrectorTimeStamp                    = _previousMinCorrectorTimeStamp;
      _minCorrectorTimeStepSize                 = _previousMinCorrectorTimeStepSize;

      _previousMinCorrectorTimeStamp            = std::numeric_limits<double>::max();
      _previousMinCorrectorTimeStepSize         = std::numeric_limits<double>::max();
      break;
  }

  _maxLevel     = _nextMaxLevel;
  _nextMaxLevel = -std::numeric_limits<int>::max(); // "-", min
}

void exahype::solvers::ADERDGSolver::rollbackToPreviousTimeStepFused() {
  switch (_timeStepping) {
    case TimeStepping::Global:
      _minNextTimeStepSize                      = std::numeric_limits<double>::max();

//      _minPredictorTimeStamp                    = _minCorrectorTimeStamp;
      _minPredictorTimeStamp                    = _previousMinCorrectorTimeStamp+_previousMinCorrectorTimeStepSize;
      _minPredictorTimeStepSize                 = _minCorrectorTimeStepSize;

      _minCorrectorTimeStamp                    = _previousMinCorrectorTimeStamp;
      _minCorrectorTimeStepSize                 = _previousMinCorrectorTimeStepSize;

      _previousMinCorrectorTimeStamp            = std::numeric_limits<double>::max();
      _previousMinCorrectorTimeStepSize         = std::numeric_limits<double>::max();
      break;
    case TimeStepping::GlobalFixed:
//      _minPredictorTimeStamp                    = _minCorrectorTimeStamp;
      _minPredictorTimeStamp                    = _previousMinCorrectorTimeStamp+_previousMinCorrectorTimeStepSize;
      _minPredictorTimeStepSize                 = _minCorrectorTimeStepSize;

      _minCorrectorTimeStamp                    = _previousMinCorrectorTimeStamp;
      _minCorrectorTimeStepSize                 = _previousMinCorrectorTimeStepSize;

      _previousMinCorrectorTimeStamp            = std::numeric_limits<double>::max();
      _previousMinCorrectorTimeStepSize         = std::numeric_limits<double>::max();
      break;
  }

  _maxLevel     = _nextMaxLevel;
  _nextMaxLevel = -std::numeric_limits<int>::max(); // "-", min
}

void exahype::solvers::ADERDGSolver::updateMinNextPredictorTimeStepSize(
    const double& minNextPredictorTimeStepSize) {
  switch (_timeStepping) {
    case TimeStepping::Global:
      _minNextTimeStepSize =
          std::min(_minNextTimeStepSize, minNextPredictorTimeStepSize);
      break;
    case TimeStepping::GlobalFixed: // TODO(Dominic): Problematic in MPI where we merge with the worker first
      _minNextTimeStepSize =
          _minPredictorTimeStamp == _minCorrectorTimeStamp
              ? std::min(_minNextTimeStepSize,
                         minNextPredictorTimeStepSize)
              : _minNextTimeStepSize;
      break;
  }
}

double exahype::solvers::ADERDGSolver::getMinNextPredictorTimeStepSize() const {
  return _minNextTimeStepSize;
}

double exahype::solvers::ADERDGSolver::getMinCorrectorTimeStamp() const {
  return _minCorrectorTimeStamp;
}

double exahype::solvers::ADERDGSolver::getMinPredictorTimeStamp() const {
  return _minPredictorTimeStamp;
}

double exahype::solvers::ADERDGSolver::getMinCorrectorTimeStepSize() const {
  return _minCorrectorTimeStepSize;
}

double exahype::solvers::ADERDGSolver::getMinPredictorTimeStepSize() const {
  return _minPredictorTimeStepSize;
}

void exahype::solvers::ADERDGSolver::setMinPredictorTimeStepSize(const double value) {
  _minPredictorTimeStepSize = value;
}

double exahype::solvers::ADERDGSolver::getPreviousMinCorrectorTimeStepSize() const {
  return _previousMinCorrectorTimeStepSize;
}

double exahype::solvers::ADERDGSolver::getPreviousMinCorrectorTimeStamp() const {
  return _previousMinCorrectorTimeStamp;
}

double exahype::solvers::ADERDGSolver::getMinTimeStamp() const {
  return getMinCorrectorTimeStamp();
}

double exahype::solvers::ADERDGSolver::getMinTimeStepSize() const {
  return getMinCorrectorTimeStepSize();
}

double exahype::solvers::ADERDGSolver::getMinNextTimeStepSize() const {
  return getMinNextPredictorTimeStepSize();
}

void exahype::solvers::ADERDGSolver::updateMinNextTimeStepSize( double value ) {
  updateMinNextPredictorTimeStepSize(value);
}

void exahype::solvers::ADERDGSolver::initSolver(
    const double timeStamp,
    const tarch::la::Vector<DIMENSIONS,double>& domainOffset,
    const tarch::la::Vector<DIMENSIONS,double>& domainSize,
    const tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize,
    const std::vector<std::string>& cmdlineargs,
    const exahype::parser::ParserView& parserView
) {
  _domainOffset=domainOffset;
  _domainSize=domainSize;
  std::pair<double,int> coarsestMeshInfo =
      exahype::solvers::Solver::computeCoarsestMeshSizeAndLevel(_maximumMeshSize,boundingBoxSize[0]);
  _coarsestMeshSize  = coarsestMeshInfo.first;
  _coarsestMeshLevel = coarsestMeshInfo.second;

  _previousMinCorrectorTimeStepSize = 0.0;
  _minCorrectorTimeStepSize = 0.0;
  _minPredictorTimeStepSize = 0.0;

  _previousMinCorrectorTimeStamp = timeStamp;
  _minCorrectorTimeStamp         = timeStamp;
  _minPredictorTimeStamp         = timeStamp;

  _meshUpdateRequest = true;

  init(cmdlineargs,parserView); // call user define initalisiation
}

bool exahype::solvers::ADERDGSolver::isPerformingPrediction(
    const exahype::State::AlgorithmSection& section) const {
  bool isPerformingPrediction = false;

  switch (section) {
    case exahype::State::AlgorithmSection::TimeStepping:
      isPerformingPrediction = true;
      break;
    case exahype::State::AlgorithmSection::PredictionRerunAllSend:
      isPerformingPrediction = !getMeshUpdateRequest() &&
                               getStabilityConditionWasViolated();
      break;
    case exahype::State::AlgorithmSection::PredictionOrLocalRecomputationAllSend:
      isPerformingPrediction = getMeshUpdateRequest();
      break;
    default:
      break;
  }

  return isPerformingPrediction;
}

bool exahype::solvers::ADERDGSolver::isMergingMetadata(
    const exahype::State::AlgorithmSection& section) const {
  bool isMergingMetadata = false;

  switch (section) {
    case exahype::State::AlgorithmSection::MeshRefinement:
      isMergingMetadata = getMeshUpdateRequest();
      break;
    default:
      break;
  }

  return isMergingMetadata;
}

void exahype::solvers::ADERDGSolver::setStabilityConditionWasViolated(bool state) {
  _stabilityConditionWasViolated = state;
}

bool exahype::solvers::ADERDGSolver::getStabilityConditionWasViolated() const {
  return _stabilityConditionWasViolated;
}

bool exahype::solvers::ADERDGSolver::isValidCellDescriptionIndex(
    const int cellDescriptionsIndex) {
  bool result = cellDescriptionsIndex>=0;
  assertion1(!result || Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  return result;
}

int exahype::solvers::ADERDGSolver::tryGetElement(
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  if ( isValidCellDescriptionIndex(cellDescriptionsIndex) ) {
    int element=0;
    for (auto& p : Heap::getInstance().getData(cellDescriptionsIndex)) {
      if (p.getSolverNumber()==solverNumber) {
        return element;
      }
      ++element;
    }
  }
  return NotFound;
}

///////////////////////////////////
// CELL-LOCAL MESH REFINEMENT
///////////////////////////////////
bool exahype::solvers::ADERDGSolver::progressMeshRefinementInEnterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const bool initialGrid,
    const int solverNumber) {
  bool newComputeCell = false;

  // Fine grid cell based uniform mesh refinement.
  const int fineGridCellElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (
      fineGridCellElement==exahype::solvers::Solver::NotFound &&
      tarch::la::allSmallerEquals(fineGridVerticesEnumerator.getCellSize(),getMaximumMeshSize()) &&
      tarch::la::oneGreater(coarseGridVerticesEnumerator.getCellSize(),getMaximumMeshSize())
  ) {
    logDebug("progressMeshRefinementInEnterCell(...)","Add new uniform grid cell at centre="<<fineGridVerticesEnumerator.getCellCenter() <<", level="<<fineGridVerticesEnumerator.getLevel()
        << " for solver=" << solverNumber);

    addNewCell(
        fineGridCell,fineGridVerticesEnumerator,
        multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex,
        solverNumber);
    newComputeCell = true;
  }
  else if ( fineGridCellElement!=exahype::solvers::Solver::NotFound ) {
    CellDescription& fineGridCellDescription =
        getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridCellElement);
  
    assertion5(tarch::la::equals(fineGridVerticesEnumerator.getCellCenter(),fineGridCellDescription.getOffset()+0.5*fineGridCellDescription.getSize()),fineGridVerticesEnumerator.getCellCenter(),fineGridCellDescription.getOffset()+0.5*fineGridCellDescription.getSize(),fineGridVerticesEnumerator.getLevel(),fineGridCellDescription.getLevel(),tarch::parallel::Node::getInstance().getRank());
    assertionEquals3(fineGridVerticesEnumerator.getLevel(),fineGridCellDescription.getLevel(),fineGridVerticesEnumerator.getCellCenter(),fineGridCellDescription.getOffset()+0.5*fineGridCellDescription.getSize(),tarch::parallel::Node::getInstance().getRank());
    // ensure that the fine grid cell descriptions's parent index is pointing to the
    // coarse grid cell's cell descriptions index; this is important to re-establish
    // the parent-child relations on a new worker after a fork.
    // and to ensure
    ensureConsistencyOfParentInformation(fineGridCellDescription,coarseGridCell.getCellDescriptionsIndex());

    #if defined(Asserts) || defined(Debug)
    const int coarseGridCellElement =
        tryGetElement(coarseGridCell.getCellDescriptionsIndex(),solverNumber);
    assertion5(
        coarseGridCellElement==exahype::solvers::Solver::NotFound ||
        fineGridCellDescription.getParentIndex()==coarseGridCell.getCellDescriptionsIndex(),
        fineGridCellDescription.toString(),
        getCellDescription(coarseGridCell.getCellDescriptionsIndex(),coarseGridCellElement).toString(),
        getCellDescription(fineGridCellDescription.getParentIndex(),0).toString(),
        fineGridCell.toString(),
        coarseGridCell.toString()); // see mergeCellDescriptionsWithRemoteData.
    #endif
    #ifdef Parallel // TODO(Dominic): Still needed?
    fineGridCellDescription.setAdjacentToRemoteRank(
        exahype::Cell::isAtRemoteBoundary(fineGridVertices,fineGridVerticesEnumerator));
    #endif

    updateCommunicationStatus(fineGridCellDescription);
    ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
    ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);

    updateAugmentationStatus(fineGridCellDescription);

    progressCollectiveRefinementOperationsInEnterCell(fineGridCellDescription);

    if ( 
         tarch::parallel::Node::getInstance().getRank()==4 ||
         tarch::parallel::Node::getInstance().getRank()==3
      ) {
       logDebug("progressMeshRefinementInEnterCell(...)","[rank=3,index="<<fineGridCell.getCellDescriptionsIndex()<<"] touch"<< fineGridCellDescription.toString() );
    }

    decideOnRefinement(fineGridCellDescription);
    decideOnVirtualRefinement(fineGridCellDescription);
  }

  // Coarse grid cell based adaptive mesh refinement operations.
  // Add new cells to the grid and veto erasing or erasing virtual children
  // requests if there are cells on the fine level.
  const int coarseGridCellElement =
      tryGetElement(coarseGridCell.getCellDescriptionsIndex(),solverNumber);
  if (coarseGridCellElement!=exahype::solvers::Solver::NotFound) {
    CellDescription& coarseGridCellDescription = getCellDescription(
        coarseGridCell.getCellDescriptionsIndex(),coarseGridCellElement);

    alterErasingRequestsIfNecessary(
        coarseGridCellDescription,
        fineGridCell.getCellDescriptionsIndex());

    addNewDescendantIfVirtualRefiningRequested(
            fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
            coarseGridCellDescription,coarseGridCell.getCellDescriptionsIndex());
    newComputeCell |=
        addNewCellIfRefinementRequested(
            fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
            coarseGridCellDescription,coarseGridCell.getCellDescriptionsIndex(),
            initialGrid);
  }

  return newComputeCell;
}

void exahype::solvers::ADERDGSolver::markForRefinement(CellDescription& cellDescription) {
  assertion1(cellDescription.getType()==CellDescription::Type::Cell,cellDescription.toString());
  assertion1(cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::None,cellDescription.toString());
  assertion(cellDescription.getRefinementRequest()==CellDescription::RefinementRequest::Pending);

  double* solution = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
  exahype::solvers::Solver::RefinementControl refinementControl =
      refinementCriterion(
          solution,cellDescription.getOffset()+0.5*cellDescription.getSize(),
          cellDescription.getSize(),
          cellDescription.getCorrectorTimeStamp(), // is done after the update
          cellDescription.getLevel());

  switch (refinementControl) {
  case exahype::solvers::Solver::RefinementControl::Keep:
    cellDescription.setRefinementRequest(CellDescription::RefinementRequest::Keep);
    break;
  case exahype::solvers::Solver::RefinementControl::Erase:
    cellDescription.setRefinementRequest(CellDescription::RefinementRequest::Erase);
    break;
  case exahype::solvers::Solver::RefinementControl::Refine:
    cellDescription.setRefinementRequest(CellDescription::RefinementRequest::Refine);
    break;
  }
}

void exahype::solvers::ADERDGSolver::decideOnRefinement(
    CellDescription& fineGridCellDescription) {
  // TODO(Dominic): We will balance in the decideOnRefinement (vetoErasingRequests)
  // routines based on the augmentation status flag

  if (
     fineGridCellDescription.getType()==CellDescription::Type::Cell
     &&
     (fineGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::None ||
     fineGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::VirtualRefiningRequested)
     &&
     fineGridCellDescription.getLevel()<getMaximumAdaptiveMeshLevel()
     &&
     fineGridCellDescription.getRefinementRequest()==CellDescription::RefinementRequest::Refine
  ) {
    fineGridCellDescription.setRefinementEvent(CellDescription::RefiningRequested);
  }
  else if (
      fineGridCellDescription.getType()==CellDescription::Type::Ancestor  &&
      fineGridCellDescription.getRefinementEvent()==CellDescription::None &&
      fineGridCellDescription.getRefinementRequest()==CellDescription::RefinementRequest::Pending
      // this means the the former Cell now Ancestor was not yet refined during the current
      // mesh refinement iterations.
  ) {
    fineGridCellDescription.setRefinementEvent(CellDescription::RefinementEvent::ErasingChildrenRequested);
  }

  if ( fineGridCellDescription.getRefinementRequest()!=CellDescription::RefinementRequest::Erase ) {
    const int coarseGridCellElement = tryGetElement(
        fineGridCellDescription.getParentIndex(),fineGridCellDescription.getSolverNumber());
    if ( coarseGridCellElement!=exahype::solvers::Solver::NotFound ) {
      auto& coarseGridCellDescription = getCellDescription(
          fineGridCellDescription.getParentIndex(),coarseGridCellElement);

      tarch::multicore::Lock lock(CoarseGridSemaphore);
      if (
          coarseGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ErasingChildrenRequested ||
          coarseGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ChangeChildrenToVirtualChildrenRequested
      ) {
        coarseGridCellDescription.setRefinementEvent(CellDescription::None);
        assertion1(coarseGridCellDescription.getType()==CellDescription::Type::Ancestor,
                   coarseGridCellDescription.toString());
      }
      lock.free();
    }
  }
}

void exahype::solvers::ADERDGSolver::decideOnVirtualRefinement(
    CellDescription& fineGridCellDescription) {
  // TODO(Dominic): We will balance in the decideOnRefinement (vetoErasingRequests)
  // routines based on the augmentation status flag

  // 1. Check if we can request augmenting or erasing virtual children request.
  bool idleCellOrDescendant =
      (fineGridCellDescription.getType()==CellDescription::Type::Cell ||
      fineGridCellDescription.getType()==CellDescription::Type::Descendant) &&
      fineGridCellDescription.getRefinementEvent()==CellDescription::None;

  if (
      idleCellOrDescendant &&
      fineGridCellDescription.getHasVirtualChildren() &&
      fineGridCellDescription.getAugmentationStatus()<MinimumAugmentationStatusForVirtualRefining
  ) {
    fineGridCellDescription.setRefinementEvent(CellDescription::RefinementEvent::ErasingVirtualChildrenRequested);
  }
  else if (
      idleCellOrDescendant &&
      !fineGridCellDescription.getHasVirtualChildren() &&
      fineGridCellDescription.getAugmentationStatus()>=MinimumAugmentationStatusForVirtualRefining
  ) {
    fineGridCellDescription.setRefinementEvent(CellDescription::VirtualRefiningRequested);
  }

  // 2. Check if we must veto the erasing virtual children request of the parent.
  if (
      fineGridCellDescription.getHasVirtualChildren() ||
      fineGridCellDescription.getAugmentationStatus()>0 // TODO(Dominic): Still necessary?
  ) {
    const int coarseGridCellElement = tryGetElement(fineGridCellDescription.getParentIndex(),
                                              fineGridCellDescription.getSolverNumber());
    if (coarseGridCellElement!=exahype::solvers::Solver::NotFound) {
      auto& coarseGridCellDescription = getCellDescription(fineGridCellDescription.getParentIndex(),
                                                           coarseGridCellElement);
      tarch::multicore::Lock lock(CoarseGridSemaphore);
      if (coarseGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ErasingVirtualChildrenRequested) {
        coarseGridCellDescription.setRefinementEvent(CellDescription::None);

        assertion1(fineGridCellDescription.getType()==CellDescription::Type::Descendant,
                   fineGridCellDescription.toString());
        assertion1(coarseGridCellDescription.getType()==CellDescription::Type::Cell ||
                   coarseGridCellDescription.getType()==CellDescription::Type::Descendant,
                   coarseGridCellDescription.toString());
      }
      lock.free();
    }
  }
}

void exahype::solvers::ADERDGSolver::alterErasingRequestsIfNecessary(
    CellDescription& coarseGridCellDescription,
    const int fineGridCellDescriptionsIndex) const {
  const int fineGridCellElement = tryGetElement(
      fineGridCellDescriptionsIndex,coarseGridCellDescription.getSolverNumber());
  if (fineGridCellElement!=exahype::solvers::Solver::NotFound) {
    CellDescription& fineGridCellDescription =
       getCellDescription(fineGridCellDescriptionsIndex,fineGridCellElement);
    if (
        fineGridCellDescription.getHasVirtualChildren()
        || fineGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::VirtualRefining
        || fineGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::VirtualRefiningRequested
        #ifdef Parallel
        || fineGridCellDescription.getHasToHoldDataForMasterWorkerCommunication()
        #endif
    ) {
      tarch::multicore::Lock lock(CoarseGridSemaphore);
      switch (coarseGridCellDescription.getRefinementEvent()) {
        case CellDescription::RefinementEvent::ErasingVirtualChildrenRequested: {
          assertion1(coarseGridCellDescription.getType()==CellDescription::Type::Cell ||
                     coarseGridCellDescription.getType()==CellDescription::Type::Descendant,
                     coarseGridCellDescription.toString());

          coarseGridCellDescription.setRefinementEvent(CellDescription::None);
        }  break;
        case CellDescription::RefinementEvent::ErasingChildrenRequested: {
          assertion1(coarseGridCellDescription.getType()==CellDescription::Type::Ancestor,
              coarseGridCellDescription.toString());

          coarseGridCellDescription.setRefinementEvent(
              CellDescription::RefinementEvent::ChangeChildrenToVirtualChildrenRequested);
        } break;
        default:
          break;
      }
      lock.free();
    }
  }
}

void exahype::solvers::ADERDGSolver::addNewCell(
    exahype::Cell& fineGridCell,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    const int coarseGridCellDescriptionsIndex,
    const int solverNumber) {
  logDebug("addNewCell(...)","Add new grid cell with center "<<fineGridVerticesEnumerator.getCellCenter() <<
              " at level "<<fineGridVerticesEnumerator.getLevel());

  fineGridCell.addNewCellDescription(
              solverNumber,
              CellDescription::Type::Cell,
              CellDescription::None,
              fineGridVerticesEnumerator.getLevel(),
              coarseGridCellDescriptionsIndex,
              fineGridVerticesEnumerator.getCellSize(),
              fineGridVerticesEnumerator.getVertexPosition());

  const int fineGridCellElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  CellDescription& fineGridCellDescription =
      getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridCellElement);
  ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
}

void exahype::solvers::ADERDGSolver::addNewDescendantIfVirtualRefiningRequested(
     exahype::Cell& fineGridCell,
     exahype::Vertex* const fineGridVertices,
     const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
     CellDescription& coarseGridCellDescription,
     const int coarseGridCellDescriptionsIndex) {
  const int fineGridElement = tryGetElement(
      fineGridCell.getCellDescriptionsIndex(),coarseGridCellDescription.getSolverNumber());

  // read and modify coarse grid
  tarch::multicore::Lock lock(CoarseGridSemaphore);
  const bool virtualRefiningRequested =
      coarseGridCellDescription.getRefinementEvent()==CellDescription::VirtualRefiningRequested;
  const bool virtualRefining =
      coarseGridCellDescription.getRefinementEvent()==CellDescription::VirtualRefining;

  if ( virtualRefining || virtualRefiningRequested ) {
    assertion1(coarseGridCellDescription.getType()==CellDescription::Type::Cell ||
               coarseGridCellDescription.getType()==CellDescription::Type::Descendant,
               coarseGridCellDescription.toString());
    coarseGridCellDescription.setRefinementEvent(CellDescription::None);
    if ( fineGridElement==exahype::solvers::Solver::NotFound ) {
      coarseGridCellDescription.setRefinementEvent(CellDescription::VirtualRefining);
    } else if ( virtualRefiningRequested ) {
      coarseGridCellDescription.setRefinementEvent(CellDescription::None);

      #if defined(Debug) || defined(Asserts)
      CellDescription& fineGridCellDescription =
          getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridElement);
      #endif
      assertion1(fineGridCellDescription.getType()==CellDescription::Type::Descendant,
                 fineGridCellDescription.toString());
    }
  }
  lock.free();

  // work on fine grid
  if (
      (virtualRefiningRequested || virtualRefining) &&
      fineGridElement==exahype::solvers::Solver::NotFound
  ) {
    fineGridCell.addNewCellDescription( // (EmptyDescendant),None
        coarseGridCellDescription.getSolverNumber(),
        CellDescription::Type::Descendant,
        CellDescription::None, // This should be removed from the signature and defaulted to none
        fineGridVerticesEnumerator.getLevel(),
        coarseGridCellDescriptionsIndex,
        fineGridVerticesEnumerator.getCellSize(),
        fineGridVerticesEnumerator.getVertexPosition());

    const int fineGridElement = tryGetElement(
        fineGridCell.getCellDescriptionsIndex(),coarseGridCellDescription.getSolverNumber());
    assertion(fineGridElement!=exahype::solvers::Solver::NotFound);
    CellDescription& fineGridCellDescription =
        getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridElement);
  }
}

// TODO(Dominic): Cannot be called multiple times. Need extra state?
bool exahype::solvers::ADERDGSolver::addNewCellIfRefinementRequested(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    CellDescription& coarseGridCellDescription,
    const int coarseGridCellDescriptionsIndex,
    const bool initialGrid) {
  // read and modify coarse grid
  tarch::multicore::Lock lock(CoarseGridSemaphore);
  bool refiningOrRefiningRequested =
      coarseGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::RefiningRequested ||
      coarseGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::Refining;
  if ( refiningOrRefiningRequested ) {
    assertion1(coarseGridCellDescription.getType()==CellDescription::Type::Cell,
               coarseGridCellDescription.toString());
    coarseGridCellDescription.setRefinementEvent(CellDescription::RefinementEvent::Refining);
  }
  lock.free();

  // work on fine grid
  if ( refiningOrRefiningRequested ) {
    const int fineGridCellElement = tryGetElement(
        fineGridCell.getCellDescriptionsIndex(),coarseGridCellDescription.getSolverNumber());

    if ( fineGridCellElement==exahype::solvers::Solver::NotFound ) {
      addNewCell(fineGridCell,fineGridVerticesEnumerator,
                 coarseGridCellDescriptionsIndex,
                 coarseGridCellDescription.getSolverNumber());
      CellDescription& fineGridCellDescription =
          getCellDescriptions(fineGridCell.getCellDescriptionsIndex()).back();
      fineGridCellDescription.setRefinementEvent(CellDescription::Prolongating);
      fineGridCellDescription.setRefinementRequest(CellDescription::RefinementRequest::Pending);
    } else {
      CellDescription& fineGridCellDescription = getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridCellElement);
      #ifdef Parallel
      assertion4(fineGridCellDescription.getType()==CellDescription::Type::Descendant ||
                 fineGridCellDescription.getType()==CellDescription::Type::Cell,
                 fineGridCellDescription.toString(),coarseGridCellDescription.toString(),
                 coarseGridCellDescriptionsIndex,
                 tarch::parallel::Node::getInstance().getRank());
      #else  
      assertion2(fineGridCellDescription.getType()==CellDescription::Type::Descendant,
                 fineGridCellDescription.toString(),coarseGridCellDescription.toString());
      #endif 
      assertion2(fineGridCellDescription.getParentIndex()==coarseGridCellDescriptionsIndex,
                 fineGridCellDescription.toString(),coarseGridCellDescriptionsIndex);

      fineGridCellDescription.setType(CellDescription::Type::Cell);
      fineGridCellDescription.setRefinementEvent(CellDescription::RefinementEvent::Prolongating);
      fineGridCellDescription.setRefinementRequest(CellDescription::RefinementRequest::Pending);
      fineGridCellDescription.setCommunicationStatus(CellCommunicationStatus);
      fineGridCellDescription.setFacewiseCommunicationStatus(0); // implicit conversion
      ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
    }
    return true;
  }
  return false;
}

void exahype::solvers::ADERDGSolver::prolongateVolumeData(
    CellDescription&       fineGridCellDescription,
    const bool initialGrid) {
  const int coarseGridElement =
      tryGetElement(fineGridCellDescription.getParentIndex(),fineGridCellDescription.getSolverNumber());
  assertion1(coarseGridElement!=exahype::solvers::Solver::NotFound,fineGridCellDescription.toString());
  CellDescription& coarseGridCellDescription =
      getCellDescription(fineGridCellDescription.getParentIndex(),coarseGridElement);

  tarch::la::Vector<DIMENSIONS,int> subcellIndex =
      exahype::amr::computeSubcellIndex(
          fineGridCellDescription.getOffset(),
          fineGridCellDescription.getSize(),coarseGridCellDescription.getOffset());

  const int levelFine = fineGridCellDescription.getLevel();
  const int levelCoarse = coarseGridCellDescription.getLevel();
  assertion(levelCoarse < levelFine);

  // current solution
  double* solutionFine   = DataHeap::getInstance().getData(
      fineGridCellDescription.getSolution()).data();
  double* solutionCoarse = DataHeap::getInstance().getData(
      coarseGridCellDescription.getSolution()).data();
  volumeUnknownsProlongation(
      solutionFine,solutionCoarse,
      levelCoarse,levelFine,
      subcellIndex);

  // previous solution
  assertion(DataHeap::getInstance().isValidIndex(fineGridCellDescription.getPreviousSolution()));
  double* previousSolutionFine   = DataHeap::getInstance().getData(
      fineGridCellDescription.getPreviousSolution()).data();
  double* previousSolutionCoarse = DataHeap::getInstance().getData(
      coarseGridCellDescription.getPreviousSolution()).data();
  volumeUnknownsProlongation(
      previousSolutionFine,previousSolutionCoarse,
      levelCoarse,levelFine,
      subcellIndex);

  fineGridCellDescription.setCorrectorTimeStamp(coarseGridCellDescription.getCorrectorTimeStamp());
  fineGridCellDescription.setPredictorTimeStamp(coarseGridCellDescription.getPredictorTimeStamp());
  fineGridCellDescription.setCorrectorTimeStepSize(coarseGridCellDescription.getCorrectorTimeStepSize());
  fineGridCellDescription.setPredictorTimeStepSize(coarseGridCellDescription.getPredictorTimeStepSize());

  // TODO Dominic: This is a little inconsistent since I orignially tried to hide
  // the limiting from the pure ADER-DG scheme
  fineGridCellDescription.setPreviousLimiterStatus(0);
  fineGridCellDescription.setLimiterStatus(0);

  // TODO Dominic:
  // During the inital mesh build where we only refine
  // according to the PAD, we don't want to have a too broad refined area.
  // We thus do not flag children cells with troubled
  if (
      !initialGrid
      &&
      coarseGridCellDescription.getLimiterStatus()>=_minimumLimiterStatusForTroubledCell
  ) {
    fineGridCellDescription.setLimiterStatus(_minimumLimiterStatusForTroubledCell);
    fineGridCellDescription.setIterationsToCureTroubledCell(coarseGridCellDescription.getIterationsToCureTroubledCell());
  }
  fineGridCellDescription.setFacewiseLimiterStatus(0);
}

bool exahype::solvers::ADERDGSolver::attainedStableState(
        exahype::Cell& fineGridCell,
        exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
        const int solverNumber) const {
  const int element = tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if ( element!=exahype::solvers::Solver::NotFound ) {
    CellDescription& cellDescription = getCellDescription(fineGridCell.getCellDescriptionsIndex(),element);

    // compute flagging gradients in inside cells
    bool flaggingHasConverged = true;
    if (
        (cellDescription.getType()==CellDescription::Type::Cell ||
        cellDescription.getType()==CellDescription::Type::Ancestor)
        &&
        !peano::grid::aspects::VertexStateAnalysis::isOneVertexBoundary(fineGridVertices,fineGridVerticesEnumerator) ) {
      for (int d=0; d<DIMENSIONS; d++) {
        flaggingHasConverged &=
            std::abs(cellDescription.getFacewiseAugmentationStatus(2*d+1)  - cellDescription.getFacewiseAugmentationStatus(2*d+0)) <= 2;
        flaggingHasConverged &=
            std::abs(cellDescription.getFacewiseCommunicationStatus(2*d+1) - cellDescription.getFacewiseCommunicationStatus(2*d+0)) <= 2;
        flaggingHasConverged &=
            std::abs(cellDescription.getFacewiseLimiterStatus(2*d+1)       - cellDescription.getFacewiseLimiterStatus(2*d+0)) <= 2;
      }
    }

    return 
        cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::None
        &&
        (cellDescription.getType()!=CellDescription::Cell ||
        cellDescription.getRefinementRequest()!=CellDescription::RefinementRequest::Pending)
        &&
        flaggingHasConverged;
  } else {
    return true;
  }
}

bool exahype::solvers::ADERDGSolver::progressMeshRefinementInLeaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber) {
  bool newComputeCell = false;

  const int fineGridCellElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (
       fineGridCellElement!=exahype::solvers::Solver::NotFound) {
    CellDescription& fineGridCellDescription = getCellDescription(
            fineGridCell.getCellDescriptionsIndex(),fineGridCellElement);

    // skip remainder if the refinement criterion has not been evaluated yet for a Cell
    // Reading the refinement request might result into data race but this is accepted at this point
    // as we only read and not write
    if (
      fineGridCellDescription.getType()              != CellDescription::Type::Cell ||
      fineGridCellDescription.getRefinementRequest() != CellDescription::RefinementRequest::Pending 
    ) { 
      // start or finish collective operations
      newComputeCell |= progressCollectiveRefinementOperationsInLeaveCell(fineGridCellDescription);

      const int coarseGridElement =
          tryGetElement(coarseGridCell.getCellDescriptionsIndex(),solverNumber);
      if ( coarseGridElement!=exahype::solvers::Solver::NotFound ) {
        assertion3(fineGridCellDescription.getParentIndex()==coarseGridCell.getCellDescriptionsIndex(),
                   fineGridCellDescription.toString(),fineGridCell.toString(),
                   coarseGridCell.toString()); // see mergeCellDescriptionsWithRemoteData.

        CellDescription& coarseGridCellDescription = getCellDescription(
            fineGridCellDescription.getParentIndex(),coarseGridElement);
        assertion1(fineGridCellDescription.getSolverNumber()==
            coarseGridCellDescription.getSolverNumber(),
                       fineGridCellDescription.toString());

        eraseCellDescriptionIfNecessary(
                fineGridCell.getCellDescriptionsIndex(),
                fineGridCellElement,
                coarseGridCellDescription);

        // copy and restrict the limiter status
        if (
          fineGridCellDescription.getType()==CellDescription::Ancestor ||
          fineGridCellDescription.getType()==CellDescription::Type::Cell
        ) {
          restrictToNextParent(fineGridCellDescription,coarseGridElement);
        }
      }
    }
  }
  return newComputeCell;
}

exahype::solvers::Solver::RefinementControl
exahype::solvers::ADERDGSolver::eraseOrRefineAdjacentVertices(
    const int cellDescriptionsIndex,
    const int solverNumber,
    const tarch::la::Vector<DIMENSIONS, double>& cellOffset,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize,
    const bool checkThoroughly) const {
  if ( tarch::la::oneGreater(cellSize,_maximumMeshSize) ) {
     return RefinementControl::Refine;
  } else {
    const int isValidIndex =
        cellDescriptionsIndex > 0 &&
        (!checkThoroughly ||
        Heap::getInstance().getInstance().isValidIndex(cellDescriptionsIndex));

    if ( isValidIndex ) {
      const int element = tryGetElement(cellDescriptionsIndex,solverNumber);
      if (element!=NotFound) {
        CellDescription& cellDescription = getCellDescription(
            cellDescriptionsIndex,element);

        if ( !checkThoroughly || tarch::la::equals( cellDescription.getOffset(), cellOffset ) )  {
          bool refineAdjacentVertices =
              cellDescription.getType()==CellDescription::Type::Ancestor ||
              cellDescription.getHasVirtualChildren() ||
              cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ChangeChildrenToVirtualChildrenRequested ||
              cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ChangeChildrenToVirtualChildren ||
              cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::VirtualRefiningRequested ||
              cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::RefiningRequested ||
              cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::Refining ||
              cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::VirtualRefining;

          #ifdef Asserts
          assertion1(
              cellDescription.getRefinementEvent()!=CellDescription::RefinementEvent::RefiningRequested ||
              cellDescription.getType()==CellDescription::Type::Cell,
              cellDescription.toString());
          assertion1(
              cellDescription.getRefinementEvent()!=CellDescription::RefinementEvent::VirtualRefiningRequested ||
              cellDescription.getType()==CellDescription::Type::Cell ||
              cellDescription.getType()==CellDescription::Type::Descendant,
              cellDescription.toString());
          #endif

          bool eraseAdjacentVertices =
              (cellDescription.getType()==CellDescription::Type::Cell ||
                  cellDescription.getType()==CellDescription::Type::Descendant)
                  &&
                  cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::None
                  &&
                  !cellDescription.getHasVirtualChildren()
                  &&
                  cellDescription.getAugmentationStatus()==0 // TODO(Dominic): Probably can tune here. This is chosen to large
                  &&
                  cellDescription.getLimiterStatus()==0;

          if (refineAdjacentVertices) {
            return RefinementControl::Refine;
          } else if (eraseAdjacentVertices) {
            return RefinementControl::Erase;
          } else {
            return RefinementControl::Keep;
          }
        } else {
          return RefinementControl::Keep; // ?
        }
      } else {
        return RefinementControl::Erase;
      }
    } else {
      return RefinementControl::Erase;
    }
  }
}

void exahype::solvers::ADERDGSolver::prepareVolumeDataRestriction(
    CellDescription& cellDescription) const {
  double* solution =
      DataHeap::getInstance().getData(cellDescription.getSolution()).data();
  std::fill_n(solution,getDataPerCell(),0.0);
  double* previousSolution =
      DataHeap::getInstance().getData(cellDescription.getPreviousSolution()).data();
  std::fill_n(previousSolution,getDataPerCell(),0.0);
}

void exahype::solvers::ADERDGSolver::progressCollectiveRefinementOperationsInEnterCell(
     CellDescription& fineGridCellDescription) {
  switch (fineGridCellDescription.getRefinementEvent()) {
    case CellDescription::Refining:
      assertion1(fineGridCellDescription.getType()==CellDescription::Type::Cell,
                 fineGridCellDescription.toString());
      fineGridCellDescription.setType(CellDescription::Type::Ancestor);
      fineGridCellDescription.setAugmentationStatus(MaximumAugmentationStatus);
      fineGridCellDescription.setFacewiseAugmentationStatus(MaximumAugmentationStatus); // implicit conversion
      fineGridCellDescription.setCommunicationStatus(0);
      fineGridCellDescription.setFacewiseCommunicationStatus(0); // implicit conversion
      fineGridCellDescription.setHasVirtualChildren(false); // since we might replace descendants with cells
      ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);
      fineGridCellDescription.setRefinementEvent(CellDescription::None);
      break;
    case CellDescription::VirtualRefining:
      fineGridCellDescription.setHasVirtualChildren(true);
      fineGridCellDescription.setRefinementEvent(CellDescription::None);
      break;
    default:
      break;
  }
}

bool exahype::solvers::ADERDGSolver::progressCollectiveRefinementOperationsInLeaveCell(
     CellDescription& fineGridCellDescription) {
  bool newComputeCell = false;
  switch (fineGridCellDescription.getRefinementEvent()) {
    case CellDescription::RefinementEvent::ErasingChildrenRequested:
      fineGridCellDescription.setType(CellDescription::Type::Cell);
      fineGridCellDescription.setAugmentationStatus(0);
      fineGridCellDescription.setFacewiseAugmentationStatus(0); // implicit conversion
      fineGridCellDescription.setCommunicationStatus(CellCommunicationStatus);
      fineGridCellDescription.setFacewiseCommunicationStatus(CellCommunicationStatus); // implicit conversion
      ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
      prepareVolumeDataRestriction(fineGridCellDescription);
      fineGridCellDescription.setRefinementEvent(CellDescription::RefinementEvent::ErasingChildren);
      break;
    case CellDescription::RefinementEvent::ErasingChildren:
      fineGridCellDescription.setRefinementEvent(CellDescription::RefinementEvent::None);
      fineGridCellDescription.setRefinementRequest(CellDescription::RefinementRequest::Pending);
      newComputeCell = true;
      break;
    case CellDescription::RefinementEvent::ChangeChildrenToVirtualChildrenRequested:
      fineGridCellDescription.setType(CellDescription::Type::Cell);
      fineGridCellDescription.setAugmentationStatus(0);
      fineGridCellDescription.setFacewiseAugmentationStatus(0); // implicit conversion
      fineGridCellDescription.setCommunicationStatus(CellCommunicationStatus);
      fineGridCellDescription.setFacewiseCommunicationStatus(CellCommunicationStatus); // implicit conversion
      ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
      prepareVolumeDataRestriction(fineGridCellDescription);
      fineGridCellDescription.setRefinementEvent(CellDescription::RefinementEvent::ChangeChildrenToVirtualChildren);
      break;
    case CellDescription::ChangeChildrenToVirtualChildren:
      fineGridCellDescription.setHasVirtualChildren(true);
      fineGridCellDescription.setRefinementEvent(CellDescription::RefinementEvent::None);
      fineGridCellDescription.setRefinementRequest(CellDescription::RefinementRequest::Pending);
      newComputeCell = true;
      break;
    case CellDescription::RefinementEvent::ErasingVirtualChildrenRequested:
      fineGridCellDescription.setRefinementEvent(CellDescription::RefinementEvent::ErasingVirtualChildren);
      break;
    case CellDescription::RefinementEvent::ErasingVirtualChildren:
      fineGridCellDescription.setHasVirtualChildren(false);
      fineGridCellDescription.setRefinementEvent(CellDescription::RefinementEvent::None);
      break;
    default:
      break;
  }
  return newComputeCell;
}

void exahype::solvers::ADERDGSolver::eraseCellDescriptionIfNecessary(
    const int cellDescriptionsIndex,
    const int fineGridCellElement,
    CellDescription& coarseGridCellDescription) {
  tarch::multicore::Lock lock(CoarseGridSemaphore);
  const bool erasingChildren =
      coarseGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ErasingChildren;
  const bool changeChildredToDescendants =
      coarseGridCellDescription.getRefinementEvent()==CellDescription::ChangeChildrenToVirtualChildren;
  const bool deaugmentingChildren =
      coarseGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ErasingVirtualChildren;
  lock.free();

  if ( erasingChildren ) {
    CellDescription& fineGridCellDescription = getCellDescription(
          cellDescriptionsIndex,fineGridCellElement);

    tarch::la::Vector<DIMENSIONS,int> subcellIndex =
        exahype::amr::computeSubcellIndex(
            fineGridCellDescription.getOffset(),
            fineGridCellDescription.getSize(),coarseGridCellDescription.getOffset());

    // restrict values.
    tarch::multicore::Lock lock(RestrictionSemaphore);
    restrictVolumeData(
        coarseGridCellDescription,
        fineGridCellDescription,
        subcellIndex);
    // TODO(Dominic): Reconsider for anarchic time stepping.
    // coarseGridCellDescription.setCorrectorTimeStamp(fineGridCellDescription.getCorrectorTimeStamp());
    // coarseGridCellDescription.setPredictorTimeStamp(fineGridCellDescription.getPredictorTimeStamp());
    lock.free();

    // erase cell description // or change to descendant
    fineGridCellDescription.setType(CellDescription::Erased);
    fineGridCellDescription.setCommunicationStatus(0);
    fineGridCellDescription.setFacewiseCommunicationStatus(0); // implicit conversion
    ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);

    getCellDescriptions(cellDescriptionsIndex).erase(
        getCellDescriptions(cellDescriptionsIndex).begin()+fineGridCellElement);
  }
  else if ( changeChildredToDescendants ) {
    CellDescription& fineGridCellDescription = getCellDescription(
        cellDescriptionsIndex,fineGridCellElement);

    tarch::la::Vector<DIMENSIONS,int> subcellIndex =
        exahype::amr::computeSubcellIndex(
            fineGridCellDescription.getOffset(),
            fineGridCellDescription.getSize(),coarseGridCellDescription.getOffset());

    // restrict values.
    tarch::multicore::Lock lock(RestrictionSemaphore);
    restrictVolumeData(
        coarseGridCellDescription,
        fineGridCellDescription,
        subcellIndex);
    // TODO(Dominic): Reconsider for anarchic time stepping.
    // coarseGridCellDescription.setCorrectorTimeStamp(fineGridCellDescription.getCorrectorTimeStamp());
    // coarseGridCellDescription.setPredictorTimeStamp(fineGridCellDescription.getPredictorTimeStamp());
    lock.free();

    // erase cell description // or change to descendant
    fineGridCellDescription.setType(CellDescription::Type::Descendant);
    fineGridCellDescription.setCommunicationStatus(0);
    fineGridCellDescription.setFacewiseCommunicationStatus(0); // implicit conversion
    ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);
  }
  else if ( deaugmentingChildren ) {
    CellDescription& fineGridCellDescription = getCellDescription(
          cellDescriptionsIndex,fineGridCellElement);

    fineGridCellDescription.setType(CellDescription::Erased);
    fineGridCellDescription.setCommunicationStatus(0);
    fineGridCellDescription.setFacewiseCommunicationStatus(0); // implicit conversion
    ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);

    getCellDescriptions(cellDescriptionsIndex).erase(
        getCellDescriptions(cellDescriptionsIndex).begin()+fineGridCellElement);
  }
}

void exahype::solvers::ADERDGSolver::restrictVolumeData(
    CellDescription&       coarseGridCellDescription,
    const CellDescription& fineGridCellDescription,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {
//  assertion1(coarseGridCellDescription.getLimiterStatus()==CellDescription::LimiterStatus::Ok,
//      coarseGridCellDescription.toString()); // TODO(Dominic): Does not always apply see veto
  assertion1(fineGridCellDescription.getLimiterStatus()==0,
        fineGridCellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(
      fineGridCellDescription.getSolution()),fineGridCellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(
      coarseGridCellDescription.getSolution()),coarseGridCellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(
      fineGridCellDescription.getPreviousSolution()),fineGridCellDescription.toString());
  assertion1(DataHeap::getInstance().isValidIndex(
      coarseGridCellDescription.getPreviousSolution()),coarseGridCellDescription.toString());

  const int levelFine  = fineGridCellDescription.getLevel();
  const int levelCoarse = coarseGridCellDescription.getLevel();
  assertion(levelCoarse < levelFine);

  // restrict current solution
  double* solutionFine   = DataHeap::getInstance().getData(
      fineGridCellDescription.getSolution()).data();
  double* solutionCoarse = DataHeap::getInstance().getData(
      coarseGridCellDescription.getSolution()).data();
  volumeUnknownsRestriction(
      solutionCoarse,solutionFine,
      levelCoarse,levelFine,
      subcellIndex);

  // restrict next solution
  double* previousSolutionFine   = DataHeap::getInstance().getData(
      fineGridCellDescription.getPreviousSolution()).data();
  double* previousSolutionCoarse = DataHeap::getInstance().getData(
      coarseGridCellDescription.getPreviousSolution()).data();
  volumeUnknownsRestriction(
      previousSolutionCoarse,previousSolutionFine,
      levelCoarse,levelFine,
      subcellIndex);

  // Reset the min and max
  const int numberOfObservables = getDMPObservables();
  if (numberOfObservables>0) {
    double* solutionMin = DataHeap::getInstance().getData(
        coarseGridCellDescription.getSolutionMin()).data();
    std::fill_n(solutionMin,DIMENSIONS_TIMES_TWO*numberOfObservables,
        std::numeric_limits<double>::max());
    double* solutionMax = DataHeap::getInstance().getData(
        coarseGridCellDescription.getSolutionMax()).data();
    std::fill_n(solutionMax,DIMENSIONS_TIMES_TWO*numberOfObservables,
        -std::numeric_limits<double>::max()); // Be aware of "-"
  }

  // TODO(Dominic): What to do with the time step data for anarchic time stepping?
  // Tobias proposed some waiting procedure. Until they all have reached
  // a certain time level.
//  coarseGridCellDescription.setCorrectorTimeStamp(fineGridCellDescription.getCorrectorTimeStamp());
//  coarseGridCellDescription.setPredictorTimeStamp(fineGridCellDescription.getPredictorTimeStamp());
//  coarseGridCellDescription.setCorrectorTimeStepSize(fineGridCellDescription.getCorrectorTimeStepSize());
//  coarseGridCellDescription.setPredictorTimeStepSize(fineGridCellDescription.getPredictorTimeStepSize());
}

void exahype::solvers::ADERDGSolver::ensureConsistencyOfParentInformation(
    CellDescription& fineGridCellDescription,
    const int coarseGridCellDescriptionsIndex) {

  const int coarseGridElement = tryGetElement(coarseGridCellDescriptionsIndex,fineGridCellDescription.getSolverNumber());
  if ( coarseGridElement!=exahype::solvers::Solver::NotFound ) {
    assertion2(
        fineGridCellDescription.getParentIndex()==coarseGridCellDescriptionsIndex ||
        fineGridCellDescription.getParentIndex()==multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex,
        coarseGridCellDescriptionsIndex,
        fineGridCellDescription.toString());
    fineGridCellDescription.setParentIndex(coarseGridCellDescriptionsIndex);

    if ( fineGridCellDescription.getType()==CellDescription::Type::Descendant ) {
      CellDescription& coarseGridCellDescription = getCellDescription(coarseGridCellDescriptionsIndex,coarseGridElement);
      if ( coarseGridCellDescription.getType()==CellDescription::Type::Cell ) {
        fineGridCellDescription.setParentCellLevel(coarseGridCellDescription.getLevel());
        fineGridCellDescription.setParentOffset(coarseGridCellDescription.getOffset());
      } else {
        assertion1(coarseGridCellDescription.getType()==CellDescription::Type::Descendant,coarseGridCellDescription.toString());
        fineGridCellDescription.setParentCellLevel(coarseGridCellDescription.getParentCellLevel());
        fineGridCellDescription.setParentOffset(coarseGridCellDescription.getParentOffset());
      }
    }
  } else {
    fineGridCellDescription.setParentIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
  }
}

void exahype::solvers::ADERDGSolver::finaliseStateUpdates(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) {
  const int element =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (element!=exahype::solvers::Solver::NotFound) {
    CellDescription& cellDescription = getCellDescription(fineGridCell.getCellDescriptionsIndex(),element);

    cellDescription.setPreviousAugmentationStatus(cellDescription.getAugmentationStatus());
    cellDescription.setRefinementRequest(CellDescription::RefinementRequest::Pending);
  }
}

////////////////////////////////////////
// CELL-LOCAL
////////////////////////////////////////
void exahype::solvers::ADERDGSolver::validateCellDescriptionData(
  const CellDescription& cellDescription,
  const bool validateTimeStepData,
  const bool afterCompression,
  const std::string& methodTraceOfCaller) const {
  #if defined(Debug) || defined(Asserts)
  if (validateTimeStepData) {
    assertion2(std::isfinite(cellDescription.getPredictorTimeStepSize()),
        cellDescription.toString(),toString());
    assertion3(cellDescription.getPredictorTimeStepSize()<
        std::numeric_limits<double>::max(),
        cellDescription.toString(),toString(),tarch::parallel::Node::getInstance().getRank());
    assertion2(cellDescription.getPredictorTimeStepSize()>0,
        cellDescription.toString(),toString());

    assertion2(std::isfinite(cellDescription.getPredictorTimeStamp()),
        cellDescription.toString(),toString());
    assertion2(cellDescription.getPredictorTimeStamp()<
        std::numeric_limits<double>::max(),
        cellDescription.toString(),toString());
    assertion2(cellDescription.getPredictorTimeStamp()>=0,
        cellDescription.toString(),toString());
  }

  assertion1(cellDescription.getRefinementEvent()==CellDescription::None,cellDescription.toString());
  assertion1(getType()==exahype::solvers::Solver::Type::ADERDG,cellDescription.toString());

  if (afterCompression) {
    // TODO(Dominic)
  } else {
    assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getSolution()),cellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getUpdate()),cellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictor()),cellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()),cellDescription.toString());

    double* luh = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
    double* lduh = DataHeap::getInstance().getData(cellDescription.getUpdate()).data();

    double* lQhbnd = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data();
    double* lFhbnd = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data();

    int dataPerCell             = getDataPerCell();
    int updateSize              = getUpdateSize();

    int dataPerCellBoundary     = getBndTotalSize();
    int unknownsPerCellBoundary = getBndFluxTotalSize();

    for (int i=0; i<dataPerCell; i++) {
      assertion4(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(luh[i]),
          cellDescription.toString(),toString(),methodTraceOfCaller,i);
    }

    for (int i=0; i<updateSize; i++) {
     assertion4(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(lduh[i]),
         cellDescription.toString(),toString(),methodTraceOfCaller,i);
    }

    for (int i=0; i<dataPerCellBoundary; i++) {
      assertion4(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(lQhbnd[i]),
          cellDescription.toString(),toString(),methodTraceOfCaller,i);
    }

    for (int i=0; i<unknownsPerCellBoundary; i++) {
      assertion4(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(lFhbnd[i]),
          cellDescription.toString(),toString(),methodTraceOfCaller,i);
    }
  }
  #endif
}

bool exahype::solvers::ADERDGSolver::evaluateRefinementCriterionAfterSolutionUpdate(
      const int cellDescriptionsIndex,
      const int element) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Type::Cell &&
      cellDescription.getLevel()<getMaximumAdaptiveMeshLevel()) {
    const double* solution = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
    RefinementControl refinementControl = refinementCriterion(
                      solution,cellDescription.getOffset()+0.5*cellDescription.getSize(),
                      cellDescription.getSize(),
                      cellDescription.getCorrectorTimeStamp(), // must be called after advancing in time
                      cellDescription.getLevel());

    return refinementControl==RefinementControl::Refine;
  }

  return false;
}

exahype::solvers::Solver::UpdateResult exahype::solvers::ADERDGSolver::fusedTimeStepBody(
    const int cellDescriptionsIndex,
    const int element,
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch,
    const bool isSkeletonCell,
    const bool mustBeDoneImmediately ) {
  auto& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  // solver->synchroniseTimeStepping(cellDescription); // assumes this was done in neighbour merge
  updateSolution(cellDescription,isFirstIterationOfBatch);

  // This is important to memorise before calling startNewTimeStepFused; TODO(Dominic): Add to docu and/or make cleaner
  UpdateResult result;
  const double predictorTimeStamp    = cellDescription.getPredictorTimeStamp();
  const double predictorTimeStepSize = cellDescription.getPredictorTimeStepSize();
  result._timeStepSize        = startNewTimeStepFused(
      cellDescriptionsIndex,element,isFirstIterationOfBatch,isLastIterationOfBatch);
  result._refinementRequested = evaluateRefinementCriterionAfterSolutionUpdate(cellDescriptionsIndex,element);

  if (
      !SpawnPredictionAsBackgroundJob ||
      mustBeDoneImmediately
  ) {   // TODO(Dominic): Add to docu. This will spawn or do a compression job right afterwards and must thus come last. This order is more natural anyway
    performPredictionAndVolumeIntegralBody(
          cellDescriptionsIndex, element,
          predictorTimeStamp,predictorTimeStepSize,
          false, isSkeletonCell );
  } else {
    PredictionJob predictionJob(
        *this, cellDescriptionsIndex, element, predictorTimeStamp,predictorTimeStepSize,
        false/*is uncompressed*/, isSkeletonCell );
    Solver::submitPredictionJob(predictionJob,isSkeletonCell);
  }
  return result;
}

exahype::solvers::Solver::UpdateResult exahype::solvers::ADERDGSolver::fusedTimeStep(
    const int cellDescriptionsIndex,
    const int element,
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch,
    const bool isAtRemoteBoundary) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  if (cellDescription.getType()==CellDescription::Type::Cell) {
    const bool isAMRSkeletonCell     = ADERDGSolver::belongsToAMRSkeleton(cellDescription,isAtRemoteBoundary);
    const bool isSkeletonCell        = isAMRSkeletonCell || isAtRemoteBoundary;
    const bool mustBeDoneImmediately = isSkeletonCell && PredictionSweeps==1;

    if (
        !SpawnPredictionAsBackgroundJob ||
        isFirstIterationOfBatch ||
        isLastIterationOfBatch
    ) {
      return
          fusedTimeStepBody(
              cellDescriptionsIndex,element,
              isFirstIterationOfBatch,isLastIterationOfBatch,isSkeletonCell, mustBeDoneImmediately );
    } else {
      FusedTimeStepJob fusedTimeStepJob( *this, cellDescriptionsIndex, element, isSkeletonCell );
      Solver::submitPredictionJob(fusedTimeStepJob,isSkeletonCell);
      return UpdateResult();
    }
  } else {
    return UpdateResult();
  }
}

exahype::solvers::Solver::UpdateResult exahype::solvers::ADERDGSolver::update(
      const int cellDescriptionsIndex,
      const int element,
      const bool isAtRemoteBoundary){
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  if (cellDescription.getType()==CellDescription::Type::Cell) {
    uncompress(cellDescription);

    UpdateResult result;
    updateSolution(cellDescriptionsIndex,element,true);
    result._timeStepSize         = startNewTimeStep(cellDescriptionsIndex,element);
    result._refinementRequested |= evaluateRefinementCriterionAfterSolutionUpdate(
                                   cellDescriptionsIndex,element);

    compress(cellDescriptionsIndex,element,isAtRemoteBoundary);
 
    return result;
  } else {
    return UpdateResult();
  }
}

void exahype::solvers::ADERDGSolver::compress(
      const int cellDescriptionsIndex,
      const int element,
      const bool isAtRemoteBoundary) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  if (cellDescription.getType()==CellDescription::Type::Cell) {
    const bool isSkeletonCell = belongsToAMRSkeleton(cellDescription,isAtRemoteBoundary);
    compress(cellDescription,isSkeletonCell);
  }
}


void exahype::solvers::ADERDGSolver::performPredictionAndVolumeIntegral(
    exahype::solvers::Solver* solver,
    const int cellDescriptionsIndex,
    const int element,
    const bool isAtRemoteBoundary) {
  exahype::solvers::ADERDGSolver* aderdgSolver = nullptr;

  switch (solver->getType()) {
    case exahype::solvers::Solver::Type::ADERDG:
      aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
      break;
    case exahype::solvers::Solver::Type::LimitingADERDG:
      aderdgSolver =
          static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
      break;
    case exahype::solvers::Solver::Type::FiniteVolumes:
      break;
  }

  if (aderdgSolver!=nullptr) {
    CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
    if (cellDescription.getType()==CellDescription::Type::Cell) {
      aderdgSolver->synchroniseTimeStepping(cellDescription);
      aderdgSolver->performPredictionAndVolumeIntegral(
          cellDescriptionsIndex,element,
          cellDescription.getPredictorTimeStamp(),
          cellDescription.getPredictorTimeStepSize(),
          true,isAtRemoteBoundary);
    }
  }
}

bool exahype::solvers::ADERDGSolver::belongsToAMRSkeleton(const CellDescription& cellDescription, const bool isAtRemoteBoundary) {
  return cellDescription.getHasVirtualChildren();
}

void exahype::solvers::ADERDGSolver::performPredictionAndVolumeIntegralBody(
    const int    cellDescriptionsIndex,
    const int    element,
    const double predictorTimeStamp,
    const double predictorTimeStepSize,
    const bool   uncompressBefore,
    const bool   isSkeletonCell ) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (uncompressBefore) {
    uncompress(cellDescription);
  }

  validateCellDescriptionData(cellDescription,true,false,"exahype::solvers::ADERDGSolver::performPredictionAndVolumeIntegralBody [pre]");

  double* luh  = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
  double* lduh = DataHeap::getInstance().getData(cellDescription.getUpdate()).data();
  double* lQhbnd = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data();
  double* lFhbnd = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data();

  #if defined(Debug) || defined(Asserts)
  for (int i=0; i<getUnknownsPerCell(); i++) { // cellDescription.getCorrectorTimeStepSize==0.0 is an initial condition
    assertion3(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0) || std::isfinite(luh[i]),cellDescription.toString(),"performPredictionAndVolumeIntegral(...)",i);
  }
  #endif

  fusedSpaceTimePredictorVolumeIntegral(
      lduh,lQhbnd,lFhbnd,
      luh,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      predictorTimeStamp,
      predictorTimeStepSize);

  compress(cellDescription,isSkeletonCell);

  validateCellDescriptionData(cellDescription,true,true,"exahype::solvers::ADERDGSolver::performPredictionAndVolumeIntegralBody [post]");
}

void exahype::solvers::ADERDGSolver::performPredictionAndVolumeIntegral(
    const int    cellDescriptionsIndex,
    const int    element,
    const double predictorTimeStamp,
    const double predictorTimeStepSize,
    const bool   uncompressBefore,
    const bool   isAtRemoteBoundary) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Type::Cell) {
    const bool isAMRSkeletonCell     = ADERDGSolver::belongsToAMRSkeleton(cellDescription,isAtRemoteBoundary);
    const bool isSkeletonCell        = isAMRSkeletonCell || isAtRemoteBoundary;
    const bool mustBeDoneImmediately = isSkeletonCell && PredictionSweeps==1;

    if (
        !SpawnPredictionAsBackgroundJob ||
        mustBeDoneImmediately
    ) {
      performPredictionAndVolumeIntegralBody(
          cellDescriptionsIndex,element,
          predictorTimeStamp,predictorTimeStepSize,
          uncompressBefore,isSkeletonCell);
    }
    else {
      PredictionJob predictionJob( *this,cellDescriptionsIndex,element,predictorTimeStamp,predictorTimeStepSize,
          uncompressBefore,isSkeletonCell );
      Solver::submitPredictionJob(predictionJob,isSkeletonCell);
    }
  }
}

double exahype::solvers::ADERDGSolver::computeTimeStepSize(CellDescription& cellDescription) {
  if (cellDescription.getType()==CellDescription::Type::Cell) {
    assertion1(cellDescription.getRefinementEvent()==CellDescription::None,cellDescription.toString());
    const double* luh = exahype::DataHeap::getInstance().getData(cellDescription.getSolution()).data();

    validateCellDescriptionData(cellDescription,false,false,"computeTimeStepSizes(...)");
    double admissibleTimeStepSize = stableTimeStepSize(luh,cellDescription.getSize());
    assertion2(admissibleTimeStepSize>0,admissibleTimeStepSize,cellDescription.toString());
    assertion3(admissibleTimeStepSize<std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),admissibleTimeStepSize,cellDescription.toString());
    assertion2(std::isfinite(admissibleTimeStepSize),admissibleTimeStepSize,cellDescription.toString());

    return admissibleTimeStepSize;
  }

  return std::numeric_limits<double>::max();
}


double exahype::solvers::ADERDGSolver::startNewTimeStepFused(
    const int cellDescriptionsIndex,
    const int element,
    const bool firstBatchIteration,
    const bool lastBatchIteration) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  if (cellDescription.getType()==CellDescription::Type::Cell) {
    double admissibleTimeStepSize = computeTimeStepSize(cellDescription);

    // Direct update of the cell description time steps
    // Note that these local quantities might
    // be overwritten again by the synchronisation
    // happening in the next time step.

    // n-1
    if (firstBatchIteration) {
      cellDescription.setPreviousCorrectorTimeStamp(cellDescription.getCorrectorTimeStamp());
      cellDescription.setPreviousCorrectorTimeStepSize(cellDescription.getCorrectorTimeStepSize());
    }

    // n
    cellDescription.setCorrectorTimeStamp(cellDescription.getPredictorTimeStamp());
    cellDescription.setCorrectorTimeStepSize(cellDescription.getPredictorTimeStepSize());

    // n+1
    cellDescription.setPredictorTimeStamp(
        cellDescription.getPredictorTimeStamp() + cellDescription.getPredictorTimeStepSize());
    if (lastBatchIteration) { // freeze the predictor time step size until last iteration
      cellDescription.setPredictorTimeStepSize(admissibleTimeStepSize);
    }

    return admissibleTimeStepSize; // still reduce the admissibile time step size over multiple iterations
  }

  return std::numeric_limits<double>::max();
}

double exahype::solvers::ADERDGSolver::startNewTimeStep(const int cellDescriptionsIndex,const int element) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  if (cellDescription.getType()==CellDescription::Type::Cell) {
    double admissibleTimeStepSize = computeTimeStepSize(cellDescription);

    // Direct update of the cell description time steps
    // Note that these local quantities might
    // be overwritten again by the synchronisation
    // happening in the next time step.

    // n-1
    cellDescription.setPreviousCorrectorTimeStamp(cellDescription.getCorrectorTimeStamp());
    cellDescription.setPreviousCorrectorTimeStepSize(cellDescription.getCorrectorTimeStepSize());

    // n
    cellDescription.setCorrectorTimeStamp(
        cellDescription.getCorrectorTimeStamp()+cellDescription.getCorrectorTimeStepSize());
    cellDescription.setCorrectorTimeStepSize(admissibleTimeStepSize);

    cellDescription.setPredictorTimeStamp(cellDescription.getCorrectorTimeStamp()); // just copy the stuff over
    cellDescription.setPredictorTimeStepSize(cellDescription.getCorrectorTimeStepSize());

    return admissibleTimeStepSize;
  }

  return std::numeric_limits<double>::max();
}

double exahype::solvers::ADERDGSolver::updateTimeStepSizes(
      const int cellDescriptionsIndex,
      const int solverElement) {
  CellDescription& cellDescription = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  if (cellDescription.getType()==CellDescription::Type::Cell) {
    const double admissibleTimeStepSize = computeTimeStepSize(cellDescription);

    cellDescription.setCorrectorTimeStepSize( admissibleTimeStepSize );
    cellDescription.setPredictorTimeStepSize( admissibleTimeStepSize );
    cellDescription.setPredictorTimeStamp   ( cellDescription.getCorrectorTimeStamp() ); // Potentially not necessary

    return admissibleTimeStepSize;
  }
  return std::numeric_limits<double>::max();
}

double exahype::solvers::ADERDGSolver::updateTimeStepSizesFused(
      const int cellDescriptionsIndex,
      const int solverElement) {
  CellDescription& cellDescription = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  if (cellDescription.getType()==CellDescription::Type::Cell) {
    const double admissibleTimeStepSize = computeTimeStepSize(cellDescription);

    cellDescription.setCorrectorTimeStepSize( admissibleTimeStepSize );
    cellDescription.setPredictorTimeStepSize( admissibleTimeStepSize );
    cellDescription.setPredictorTimeStamp   ( cellDescription.getCorrectorTimeStamp()+admissibleTimeStepSize );

    return admissibleTimeStepSize;
  }
  return std::numeric_limits<double>::max();
}

void exahype::solvers::ADERDGSolver::zeroTimeStepSizes(
    const int cellDescriptionsIndex,
    const int element) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Type::Cell) {
    cellDescription.setPreviousCorrectorTimeStepSize(0.0);
    cellDescription.setCorrectorTimeStepSize(0.0);
    cellDescription.setPredictorTimeStepSize(0.0);

    cellDescription.setPredictorTimeStamp(
        cellDescription.getCorrectorTimeStamp());
  }
}

void exahype::solvers::ADERDGSolver::rollbackToPreviousTimeStep(
    const int cellDescriptionsIndex,
    const int element) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  // n+1
  cellDescription.setPredictorTimeStamp   (cellDescription.getPreviousCorrectorTimeStamp());
  cellDescription.setPredictorTimeStepSize(cellDescription.getPreviousCorrectorTimeStepSize());

  // n
  cellDescription.setCorrectorTimeStamp   (cellDescription.getPreviousCorrectorTimeStamp());
  cellDescription.setCorrectorTimeStepSize(cellDescription.getPreviousCorrectorTimeStepSize());

  // n-1
  cellDescription.setPreviousCorrectorTimeStamp(std::numeric_limits<double>::max());
  cellDescription.setPreviousCorrectorTimeStepSize(std::numeric_limits<double>::max()); // TODO(Dominic): get rid of the last time level.
}

void exahype::solvers::ADERDGSolver::rollbackToPreviousTimeStepFused(
    const int cellDescriptionsIndex,
    const int element) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  // n+1
  cellDescription.setPredictorTimeStamp   (
      cellDescription.getPreviousCorrectorTimeStamp()+cellDescription.getPreviousCorrectorTimeStepSize());
  cellDescription.setPredictorTimeStepSize(cellDescription.getCorrectorTimeStepSize());

  // n
  cellDescription.setCorrectorTimeStamp(cellDescription.getPreviousCorrectorTimeStamp());
  cellDescription.setCorrectorTimeStepSize(cellDescription.getPreviousCorrectorTimeStepSize());

  // n-1
  cellDescription.setPreviousCorrectorTimeStamp(std::numeric_limits<double>::max());
  cellDescription.setPreviousCorrectorTimeStepSize(std::numeric_limits<double>::max()); // TODO(Dominic): get rid of the last time level.
}

void exahype::solvers::ADERDGSolver::adjustSolutionDuringMeshRefinementBody(
    const int cellDescriptionsIndex,
    const int element,
    const bool isInitialMeshRefinement) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  assertion1(cellDescription.getType()==CellDescription::Type::Cell,cellDescription.toString());
  assertion1(cellDescription.getRefinementRequest()==CellDescription::RefinementRequest::Pending,cellDescription.toString());

  zeroTimeStepSizes(cellDescriptionsIndex,element); // TODO(Dominic): Still necessary?
  synchroniseTimeStepping(cellDescription);

  if (cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::Prolongating) {
    prolongateVolumeData(cellDescription,isInitialMeshRefinement);
    cellDescription.setRefinementEvent(CellDescription::RefinementEvent::None);
  }
  adjustSolution(cellDescription);
  markForRefinement(cellDescription);
}

void exahype::solvers::ADERDGSolver::adjustSolution(CellDescription& cellDescription) {
  assertion1(cellDescription.getType()==CellDescription::Type::Cell,cellDescription.toString());
  assertion1(cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::None,cellDescription.toString());

  cellDescription.setRefinementRequest(CellDescription::RefinementRequest::Pending);

  double* solution = exahype::DataHeap::getInstance().getData(cellDescription.getSolution()).data();
  adjustSolution(
      solution,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      cellDescription.getCorrectorTimeStamp(),
      cellDescription.getCorrectorTimeStepSize());

  #if defined(Debug) || defined(Asserts)
  for (int i=0; i<getUnknownsPerCell(); i++) {
    assertion3(std::isfinite(solution[i]),cellDescription.toString(),"adjustSolution(...)",i);
  }
  #endif
}

void exahype::solvers::ADERDGSolver::updateSolution(
    CellDescription& cellDescription,
    const bool backupPreviousSolution) {
  if (
    cellDescription.getType()==CellDescription::Type::Cell &&
    cellDescription.getRefinementEvent()==CellDescription::None
  ) {
    double* newSolution = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
    if (backupPreviousSolution) {
      double* solution  = DataHeap::getInstance().getData(cellDescription.getPreviousSolution()).data();
      std::copy(newSolution,newSolution+getUnknownsPerCell(),solution); // Copy (current solution) in old solution field.

      #if defined(Debug) || defined(Asserts)
      for (int i=0; i<getDataPerCell(); i++) { // cellDescription.getCorrectorTimeStepSize()==0.0 is an initial condition
        assertion3(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0)  || std::isfinite(solution[i]),cellDescription.toString(),"updateSolution(...)",i);
      } 
      #endif
    }

    double* update       = exahype::DataHeap::getInstance().getData(cellDescription.getUpdate()).data();
    #if defined(Debug) || defined(Asserts)
    for (int i=0; i<getUnknownsPerCell(); i++) {
      assertion3(tarch::la::equals(cellDescription.getCorrectorTimeStepSize(),0.0)  || std::isfinite(update[i]),cellDescription.toString(),"updateSolution",i);
    } 
    #endif

    assertion1(cellDescription.getCorrectorTimeStamp()<std::numeric_limits<double>::max(),cellDescription.toString());
    assertion1(cellDescription.getCorrectorTimeStepSize()<std::numeric_limits<double>::max(),cellDescription.toString());
    solutionUpdate(newSolution,update,cellDescription.getCorrectorTimeStepSize());

    adjustSolution(
        newSolution,
        cellDescription.getOffset()+0.5*cellDescription.getSize(),
        cellDescription.getSize(),
        cellDescription.getCorrectorTimeStamp()+cellDescription.getCorrectorTimeStepSize(),
        cellDescription.getCorrectorTimeStepSize());

    #if defined(Debug) || defined(Asserts)
    for (int i=0; i<getDataPerCell(); i++) {
      assertion3(std::isfinite(newSolution[i]),cellDescription.toString(),"updateSolution(...)",i);
    }
    #endif
  }
  assertion(cellDescription.getRefinementEvent()==CellDescription::None);

  // update helper status // TODO(Dominic): Check if we can work with the reduced values in the neighbour exchange
  updateCommunicationStatus(cellDescription);
  // marking for augmentation
  updateAugmentationStatus(cellDescription);
}

void exahype::solvers::ADERDGSolver::updateSolution(
    const int cellDescriptionsIndex,
    const int element,
    const bool backupPreviousSolution) {
  // reset helper variables
  CellDescription& cellDescription  = getCellDescription(cellDescriptionsIndex,element);
  updateSolution(cellDescription,backupPreviousSolution);
}

void exahype::solvers::ADERDGSolver::swapSolutionAndPreviousSolution(CellDescription& cellDescription) const {
  assertion(cellDescription.getType()==CellDescription::Type::Cell);
  assertion(cellDescription.getRefinementEvent()==CellDescription::None);

  // Simply swap the heap indices
  const int previousSolution = cellDescription.getPreviousSolution();
  cellDescription.setPreviousSolution(cellDescription.getSolution());
  cellDescription.setSolution(previousSolution);
}

void exahype::solvers::ADERDGSolver::prepareFaceDataOfAncestor(CellDescription& cellDescription) {
  logDebug("prepareFaceDataOfAncestor(...)","cell="<<cellDescription.getOffset()+0.5*cellDescription.getSize() <<
          ", level=" << cellDescription.getLevel());

  assertion1(cellDescription.getType()==CellDescription::Type::Ancestor,cellDescription.toString());
  std::fill_n(DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).begin(),
              getBndTotalSize(), 0.0);
  std::fill_n(DataHeap::getInstance().getData(cellDescription.getFluctuation()).begin(),
              getBndFluxTotalSize(), 0.0);

  #if defined(Debug) || defined(Asserts)
  double* Q = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data();
  double* F = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data();
  #endif

  #if defined(Debug) || defined(Asserts)
  for(int i=0; i<getBndTotalSize(); ++i) {
    assertion2(tarch::la::equals(Q[i],0.0),i,Q[i]);
  }

  for(int i=0; i<getBndFluxTotalSize(); ++i) {
    assertion2(tarch::la::equals(F[i],0.0),i,F[i]);
  }
  #endif

  // observables min and max
  const int numberOfObservables = getDMPObservables();
  if (numberOfObservables>0) {
    std::fill_n(DataHeap::getInstance().getData(cellDescription.getSolutionMin()).begin(),
        numberOfObservables*DIMENSIONS_TIMES_TWO, std::numeric_limits<double>::max());
    std::fill_n(DataHeap::getInstance().getData(cellDescription.getSolutionMax()).begin(),
        numberOfObservables*DIMENSIONS_TIMES_TWO, -std::numeric_limits<double>::max()); // !!! Be aware of the "-"
  }
}

void exahype::solvers::ADERDGSolver::prolongateFaceDataToDescendant(
    CellDescription& cellDescription,
    SubcellPosition& subcellPosition) {
  assertion2(Heap::getInstance().isValidIndex(
      subcellPosition.parentCellDescriptionsIndex),
      subcellPosition.parentCellDescriptionsIndex,cellDescription.toString());

  CellDescription& parentCellDescription = getCellDescription(
      subcellPosition.parentCellDescriptionsIndex,subcellPosition.parentElement);

  assertion(parentCellDescription.getSolverNumber() == cellDescription.getSolverNumber());
  assertion(parentCellDescription.getType() == CellDescription::Type::Cell ||
            parentCellDescription.getType() == CellDescription::Type::Descendant);

  const int levelFine   = cellDescription.getLevel();
  const int levelCoarse = parentCellDescription.getLevel();
  assertion(levelCoarse < levelFine);
  const int levelDelta = levelFine - levelCoarse;

  DataHeap::HeapEntries& update = DataHeap::getInstance().getData(cellDescription.getUpdate());
  std::fill(update.begin(),update.end(),0.0);

  for (int d = 0; d < DIMENSIONS; ++d) {
    // Check if cell is at "left" or "right" d face of parent
    if (subcellPosition.subcellIndex[d]==0 ||
        subcellPosition.subcellIndex[d]==tarch::la::aPowI(levelDelta,3)-1) {
      const int faceIndex = 2*d + ((subcellPosition.subcellIndex[d]==0) ? 0 : 1); // Do not remove brackets.

      const int numberOfFaceDof = getBndFaceSize();
      const int numberOfFluxDof = getBndFluxSize();
      
      logDebug("prolongateFaceDataToDescendant(...)","cell=" << cellDescription.getOffset()+0.5*cellDescription.getSize() <<
               ",level=" << cellDescription.getLevel() << ",d=" << d <<
               ",face=" << faceIndex << ",subcellIndex" << subcellPosition.subcellIndex.toString() << " to " <<
               " cell="<<parentCellDescription.getOffset()+0.5*parentCellDescription.getSize()<<
               " level="<<parentCellDescription.getLevel());

      // extrapolated predictor and flux interpolation
      // extrapolated predictor
      assertion1(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictor()),cellDescription.toString());
      assertion1(DataHeap::getInstance().isValidIndex(parentCellDescription.getExtrapolatedPredictor()),parentCellDescription.toString());
      double* lQhbndFine = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
          (faceIndex * numberOfFaceDof);
      const double* lQhbndCoarse = DataHeap::getInstance().getData(parentCellDescription.getExtrapolatedPredictor()).data() +
          (faceIndex * numberOfFaceDof);
      // flux
      double* lFhbndFine = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
          (faceIndex * numberOfFluxDof);
      const double* lFhbndCoarse = DataHeap::getInstance().getData(parentCellDescription.getFluctuation()).data() +
          (faceIndex * numberOfFluxDof);

      faceUnknownsProlongation(lQhbndFine,lFhbndFine,lQhbndCoarse,
                               lFhbndCoarse, levelCoarse, levelFine,
                               exahype::amr::getSubfaceIndex(subcellPosition.subcellIndex,d));

      // time step data TODO(LTS), still need veto
      cellDescription.setPredictorTimeStamp(parentCellDescription.getPredictorTimeStamp());
      cellDescription.setPredictorTimeStepSize(parentCellDescription.getPredictorTimeStepSize());

      prolongateObservablesMinAndMax(cellDescription,parentCellDescription,faceIndex);
    }
  }
}

void exahype::solvers::ADERDGSolver::prolongateObservablesMinAndMax(
    const CellDescription& cellDescription,
    const CellDescription& cellDescriptionParent,
    const int faceIndex) const {
  const int numberOfObservables = getDMPObservables();
  if (numberOfObservables>0) {
    // fine
    double* minFine = DataHeap::getInstance().getData(cellDescription.getSolutionMin()).data() +
        (faceIndex * numberOfObservables);
    double* maxFine = DataHeap::getInstance().getData(cellDescription.getSolutionMax()).data() +
        (faceIndex * numberOfObservables);
    // coarse
    const double* minCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getSolutionMin()).data() +
        (faceIndex * numberOfObservables);
    const double* maxCoarse = DataHeap::getInstance().getData(cellDescriptionParent.getSolutionMax()).data() +
        (faceIndex * numberOfObservables);

    std::copy_n( minCoarse,numberOfObservables, minFine );
    std::copy_n( maxCoarse,numberOfObservables, maxFine );
  }
}

void exahype::solvers::ADERDGSolver::prolongateFaceData(
    const int cellDescriptionsIndex,
    const int element) {
  CellDescription& cellDescription =
      exahype::solvers::ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  if (
      cellDescription.getType()==CellDescription::Type::Descendant &&
      cellDescription.getCommunicationStatus()>=MinimumCommunicationStatusForNeighbourCommunication &&
      isValidCellDescriptionIndex(cellDescription.getParentIndex()) // might be at master-worker boundary
  ) {
    exahype::solvers::Solver::SubcellPosition subcellPosition =
        exahype::amr::computeSubcellPositionOfDescendant<CellDescription,Heap,false>(
            cellDescription);

    prolongateFaceDataToDescendant(cellDescription,subcellPosition);
  }
  assertion2(
      cellDescription.getType()!=CellDescription::Type::Descendant ||
      isValidCellDescriptionIndex(cellDescription.getParentIndex()),
      cellDescription.toString(),
      tarch::parallel::Node::getInstance().getRank());
}

void exahype::solvers::ADERDGSolver::restriction( // TODO(Dominic): Does it still make sense?
    const int fineGridCellDescriptionsIndex,
    const int fineGridElement) {
  CellDescription& cellDescription = getCellDescription(fineGridCellDescriptionsIndex,fineGridElement);

  const int parentElement = tryGetElement(
      cellDescription.getParentIndex(),cellDescription.getSolverNumber());
  if ( parentElement!=exahype::solvers::Solver::NotFound ) {
    // restrict some flags to direct parent
    //restrictToNextParent(cellDescription,parentElement);

    if (
        cellDescription.getType()==CellDescription::Type::Descendant &&
        cellDescription.getCommunicationStatus()>=MinimumCommunicationStatusForNeighbourCommunication
    ) {
      exahype::solvers::Solver::SubcellPosition subcellPosition =
          exahype::amr::computeSubcellPositionOfDescendant<CellDescription,Heap,true>(cellDescription);
      assertion1(subcellPosition.parentElement!=exahype::solvers::Solver::NotFound,cellDescription.toString());

      // restrict update and minMax
      restrictToTopMostParent(cellDescription,
          subcellPosition.parentCellDescriptionsIndex,
          subcellPosition.parentElement);
    }
  }
  // TODO(Dominic): Merge again; Have veto mechanism per face; set at interface with Ancestor
}


// TODO(Dominic): Only use for parallel routines; refactor code
void exahype::solvers::ADERDGSolver::restrictToNextParent(
      const CellDescription& cellDescription,
      const int parentElement) const {
  CellDescription& parentCellDescription =
      ADERDGSolver::getCellDescription(cellDescription.getParentIndex(),parentElement);

  tarch::multicore::Lock lock(CoarseGridSemaphore);
  if (
      cellDescription.getLimiterStatus()>=ADERDGSolver::_minimumLimiterStatusForTroubledCell &&
      parentCellDescription.getType()==CellDescription::Type::Ancestor
  ) {
    parentCellDescription.setLimiterStatus(ADERDGSolver::_minimumLimiterStatusForTroubledCell);
    parentCellDescription.setFacewiseLimiterStatus(0);
  }
  lock.free();
}

void exahype::solvers::ADERDGSolver::restrictToTopMostParent( // TODO must be merged with faceIntegral
                  const CellDescription& cellDescription,
                  const int parentCellDescriptionsIndex,
                  const int parentElement) {
  CellDescription& parentCellDescription =
      getCellDescription(parentCellDescriptionsIndex,parentElement);
  assertion(parentCellDescription.getSolverNumber()==cellDescription.getSolverNumber());
  assertion1(cellDescription.getType()==CellDescription::Type::Descendant &&
             parentCellDescription.getCommunicationStatus()>=MinimumCommunicationStatusForNeighbourCommunication,
             cellDescription.toString());
  #ifdef Parallel
  assertion1(parentCellDescription.getType()==CellDescription::Type::Cell ||
      (parentCellDescription.getType()==CellDescription::Type::Descendant &&
      parentCellDescription.getHasToHoldDataForMasterWorkerCommunication()),
      parentCellDescription.toString());
  #else
  assertion1(parentCellDescription.getType()==CellDescription::Type::Cell,
             parentCellDescription.toString());
  #endif

  DataHeap::HeapEntries& updateFine   = DataHeap::getInstance().getData(cellDescription.getUpdate());
  DataHeap::HeapEntries& updateCoarse = DataHeap::getInstance().getData(parentCellDescription.getUpdate());
  
  tarch::multicore::Lock lock(RestrictionSemaphore);
  for (int i = 0; i < getUpdateSize(); ++i) {
      updateCoarse[i] += updateFine[i];
  }
  lock.free();
  std::fill(updateFine.begin(),updateFine.end(),0.0);

  // For restricting the observables min and max, we can go level by level
  // or directly up to the top-most parent.
  const int levelDelta = cellDescription.getLevel() - parentCellDescription.getLevel();
  const tarch::la::Vector<DIMENSIONS,int> subcellIndex =
      exahype::amr::computeSubcellIndex(
          cellDescription.getOffset(),cellDescription.getSize(),
            parentCellDescription.getOffset());

  logDebug("restriction(...)","cell=" << cellDescription.getOffset()+0.5*cellDescription.getSize() <<
           ",level=" << cellDescription.getLevel() << 
           ",subcellIndex" << subcellIndex.toString() << " to " <<
           " cell="<<parentCellDescription.getOffset()+0.5*parentCellDescription.getSize()<<
           " level="<<parentCellDescription.getLevel());

  for (int d = 0; d < DIMENSIONS; d++) {
    if ( subcellIndex[d]==0 ||
         subcellIndex[d]==tarch::la::aPowI(levelDelta,3)-1 ) {
      const int faceIndex = 2*d + ((subcellIndex[d]==0) ? 0 : 1); // Do not remove brackets.

      restrictObservablesMinAndMax(parentCellDescription,cellDescription,faceIndex);
    }
  }
}

void exahype::solvers::ADERDGSolver::restrictObservablesMinAndMax(
    const CellDescription& cellDescription,
    const CellDescription& parentCellDescription,
    const int faceIndex) const {
  const int numberOfObservables = getDMPObservables();
  if (numberOfObservables>0) {
    // fine
    double* minFine = DataHeap::getInstance().getData(cellDescription.getSolutionMin()).data() +
        (faceIndex * numberOfObservables);
    double* maxFine = DataHeap::getInstance().getData(cellDescription.getSolutionMax()).data() +
        (faceIndex * numberOfObservables);
    // coarse
    double* minCoarse = DataHeap::getInstance().getData(parentCellDescription.getSolutionMin()).data() +
        (faceIndex * numberOfObservables);
    double* maxCoarse = DataHeap::getInstance().getData(parentCellDescription.getSolutionMax()).data() +
        (faceIndex * numberOfObservables);

    tarch::multicore::Lock lock(RestrictionSemaphore);
    for (int i=0; i<numberOfObservables; i++) {
      *(minCoarse+i) = std::min( *(minFine+i), *(minCoarse+i) );
      *(maxCoarse+i) = std::max( *(maxFine+i), *(maxCoarse+i) );
    }
    lock.free();
  }
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void exahype::solvers::ADERDGSolver::mergeWithLimiterStatus(
    CellDescription& cellDescription,
    const int faceIndex,
    const int otherLimiterStatus) const {
  const int croppedOtherLimiterStatus =
      std::min(
          otherLimiterStatus,
          _minimumLimiterStatusForTroubledCell );

  cellDescription.setFacewiseLimiterStatus( faceIndex, croppedOtherLimiterStatus );
}

/**
 * Iterate over the merged limiter statuses per face and
 * determine a unique value.
 */
int
exahype::solvers::ADERDGSolver::determineLimiterStatus(
    const CellDescription& cellDescription,
    const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed) {
  int max = 0;
  for (unsigned int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
    if ( neighbourMergePerformed[i] ) {
      max = std::max( max, cellDescription.getFacewiseLimiterStatus(i)-1 );
    }
  }
  return max;
}

void
exahype::solvers::ADERDGSolver::updateCommunicationStatus(
    exahype::solvers::ADERDGSolver::CellDescription& cellDescription) const {
  cellDescription.setCommunicationStatus(determineCommunicationStatus(cellDescription));
  assertion1(
      cellDescription.getType()!=CellDescription::Type::Cell ||
      cellDescription.getCommunicationStatus()==CellCommunicationStatus,
      cellDescription.toString());
}

int
exahype::solvers::ADERDGSolver::determineCommunicationStatus(
    exahype::solvers::ADERDGSolver::CellDescription& cellDescription) const {
  if ( cellDescription.getType()==CellDescription::Type::Cell ) {
    return CellCommunicationStatus;
  } else {
    int max = 0;
    for (unsigned int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
      if ( cellDescription.getNeighbourMergePerformed(i) ) {
        max = std::max( max, cellDescription.getFacewiseCommunicationStatus(i)-1 );
      }
    }
    return max;
  }
}

void exahype::solvers::ADERDGSolver::mergeWithCommunicationStatus(
    CellDescription& cellDescription,
    const int faceIndex,
    const int otherCommunicationStatus) const {
  assertion3(cellDescription.getCommunicationStatus()<=CellCommunicationStatus,
             cellDescription.getCommunicationStatus(),otherCommunicationStatus,
             cellDescription.getCommunicationStatus());
  cellDescription.setFacewiseCommunicationStatus( faceIndex, otherCommunicationStatus );
}

void
exahype::solvers::ADERDGSolver::updateAugmentationStatus(
    exahype::solvers::ADERDGSolver::CellDescription& cellDescription) const {
  cellDescription.setAugmentationStatus(determineAugmentationStatus(cellDescription));
  assertion1(
      cellDescription.getType()!=CellDescription::Type::Ancestor ||
      cellDescription.getAugmentationStatus()==MaximumAugmentationStatus,
      cellDescription.toString());
}

int
exahype::solvers::ADERDGSolver::determineAugmentationStatus(
    exahype::solvers::ADERDGSolver::CellDescription& cellDescription) const {
  if (cellDescription.getType()==CellDescription::Type::Ancestor) {
    return MaximumAugmentationStatus;
  } else {
    int max = 0;
    for (unsigned int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
      if ( cellDescription.getNeighbourMergePerformed(i) ) {
        max = std::max( max, cellDescription.getFacewiseAugmentationStatus(i)-1 );
      }
    }
    return max;
  }
}

void exahype::solvers::ADERDGSolver::mergeWithAugmentationStatus(
    CellDescription& cellDescription,
    const int faceIndex,
    const int otherAugmentationStatus) const {
  assertion3(
      cellDescription.getAugmentationStatus()<=MaximumAugmentationStatus,
      cellDescription.getAugmentationStatus(),otherAugmentationStatus,
      cellDescription.getAugmentationStatus());
  cellDescription.setFacewiseAugmentationStatus( faceIndex, otherAugmentationStatus );
}

// merge metadata
void exahype::solvers::ADERDGSolver::mergeNeighboursMetadata(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) const {

  CellDescription& cellDescription1  = getCellDescription(cellDescriptionsIndex1,element1);
  CellDescription& cellDescription2 = getCellDescription(cellDescriptionsIndex2,element2);

  const int direction    = tarch::la::equalsReturnIndex(pos1,pos2);
  const int orientation1 = (1 + pos2(direction) - pos1(direction))/2;
  const int orientation2 = 1-orientation1;

  const int faceIndex1 = 2*direction+orientation1;
  const int faceIndex2 = 2*direction+orientation2;

  mergeWithCommunicationStatus(cellDescription1,faceIndex1,cellDescription2.getCommunicationStatus());
  mergeWithAugmentationStatus(cellDescription1,faceIndex1,cellDescription2.getAugmentationStatus());
  mergeWithLimiterStatus(cellDescription1,faceIndex1,cellDescription2.getLimiterStatus());

  mergeWithCommunicationStatus(cellDescription2,faceIndex2,cellDescription1.getCommunicationStatus());
  mergeWithAugmentationStatus(cellDescription2,faceIndex2,cellDescription1.getAugmentationStatus());
  mergeWithLimiterStatus(cellDescription2,faceIndex2,cellDescription1.getLimiterStatus());
}

// merge compute data
void exahype::solvers::ADERDGSolver::mergeNeighbours(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));
  const int direction    = tarch::la::equalsReturnIndex(pos1, pos2);
  const int orientation1 = (1 + pos2(direction) - pos1(direction))/2;

  const int indexOfRightFaceOfLeftCell = 2*direction+1;
  const int indexOfLeftFaceOfRightCell = 2*direction+0;

  CellDescription& cellDescriptionLeft  =
      (orientation1==1) ?
          getCellDescription(cellDescriptionsIndex1,element1)
          :
          getCellDescription(cellDescriptionsIndex2,element2);
  CellDescription& cellDescriptionRight =
      (orientation1==1)
          ?
          getCellDescription(cellDescriptionsIndex2,element2)
          :
          getCellDescription(cellDescriptionsIndex1,element1);

  // synchronise time stepping if necessary
  synchroniseTimeStepping(cellDescriptionLeft);
  synchroniseTimeStepping(cellDescriptionRight);

  if ( CompressionAccuracy > 0.0 ) {
    peano::datatraversal::TaskSet uncompression(
      [&] () -> bool {
        uncompress(cellDescriptionLeft);
        return false;
      },
      [&] () -> bool {
        uncompress(cellDescriptionRight);
        return false;
      },
      peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
      peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
      true
    );
  }

  solveRiemannProblemAtInterface(
      cellDescriptionLeft,cellDescriptionRight,indexOfRightFaceOfLeftCell,indexOfLeftFaceOfRightCell);
}



void exahype::solvers::ADERDGSolver::solveRiemannProblemAtInterface(
    CellDescription& pLeft,
    CellDescription& pRight,
    const int faceIndexLeft,
    const int faceIndexRight) {
  if (
      (pLeft.getCommunicationStatus()==CellCommunicationStatus &&
      pLeft.getFacewiseCommunicationStatus(faceIndexLeft) >= MinimumCommunicationStatusForNeighbourCommunication &&
      pLeft.getFacewiseAugmentationStatus(faceIndexLeft)  <  MaximumAugmentationStatus) // excludes Ancestors
      ||
      (pRight.getCommunicationStatus()==CellCommunicationStatus &&
      pRight.getFacewiseCommunicationStatus(faceIndexRight) >= MinimumCommunicationStatusForNeighbourCommunication &&
      pRight.getFacewiseAugmentationStatus(faceIndexRight)  <  MaximumAugmentationStatus) // excludes Ancestors
  ) {
    assertion1(DataHeap::getInstance().isValidIndex(pLeft.getExtrapolatedPredictor()),pLeft.toString());
    assertion1(DataHeap::getInstance().isValidIndex(pLeft.getFluctuation()),pLeft.toString());
    assertion1(DataHeap::getInstance().isValidIndex(pRight.getExtrapolatedPredictor()),pRight.toString());
    assertion1(DataHeap::getInstance().isValidIndex(pRight.getFluctuation()),pRight.toString());
    assertion1(pLeft.getRefinementEvent()==CellDescription::None,pLeft.toString());
    assertion1(pRight.getRefinementEvent()==CellDescription::None,pRight.toString());
    assertionEquals4(pLeft.getNeighbourMergePerformed(faceIndexLeft),pRight.getNeighbourMergePerformed(faceIndexRight),faceIndexLeft,faceIndexRight,pLeft.toString(),pRight.toString());
    assertion4(std::abs(faceIndexLeft-faceIndexRight)==1,faceIndexLeft,faceIndexRight,pLeft.toString(),pRight.toString());

    const int dataPerFace = getBndFaceSize();
    const int dofPerFace  = getBndFluxSize();

    double* QL = DataHeap::getInstance() .getData(pLeft.getExtrapolatedPredictor()).data() + /// !!! Be aware of the dataPerFace, Left, Right
        (faceIndexLeft * dataPerFace);
    double* QR = DataHeap::getInstance().getData(pRight.getExtrapolatedPredictor()).data() +
        (faceIndexRight * dataPerFace);

    double* FL = DataHeap::getInstance().getData(pLeft.getFluctuation()).data() + /// !!! Be aware of the dofPerFace, Left, Right
        (faceIndexLeft * dofPerFace);
    double* FR = DataHeap::getInstance().getData(pRight.getFluctuation()).data() +
        (faceIndexRight * dofPerFace);

    // todo Time step must be interpolated in local time stepping case
    // both time step sizes are the same, so the min has no effect here.
    assertion1(faceIndexLeft>=0,faceIndexLeft);
    assertion1(faceIndexRight>=0,faceIndexRight);
    assertion1(faceIndexLeft<DIMENSIONS_TIMES_TWO,faceIndexLeft);
    assertion1(faceIndexRight<DIMENSIONS_TIMES_TWO,faceIndexRight);
    assertion1(faceIndexRight%2==0,faceIndexRight);
    const int normalDirection = (faceIndexRight - (faceIndexRight %2))/2;
    assertion3(normalDirection==(faceIndexLeft - (faceIndexLeft %2))/2,normalDirection,faceIndexLeft,faceIndexRight);
    assertion3(normalDirection<DIMENSIONS,normalDirection,faceIndexLeft,faceIndexRight);

    assertion3(std::isfinite(pLeft.getCorrectorTimeStepSize()),pLeft.toString(),faceIndexLeft,normalDirection);
    assertion3(std::isfinite(pRight.getCorrectorTimeStepSize()),pRight.toString(),faceIndexRight,normalDirection);
    assertion3(pLeft.getCorrectorTimeStepSize()>=0.0,pLeft.toString(),faceIndexLeft,normalDirection);
    assertion3(pRight.getCorrectorTimeStepSize()>=0.0,pRight.toString(),faceIndexRight,normalDirection);

    #if defined(Debug) || defined(Asserts)
    for(int i=0; i<dataPerFace; ++i) {
      assertion5(tarch::la::equals(pLeft.getCorrectorTimeStepSize(),0.0) || std::isfinite(QL[i]),pLeft.toString(),faceIndexLeft,normalDirection,i,QL[i]);
      assertion5(tarch::la::equals(pRight.getCorrectorTimeStepSize(),0.0) || std::isfinite(QR[i]),pRight.toString(),faceIndexRight,normalDirection,i,QR[i]);
    }
    for(int i=0; i<dofPerFace; ++i) {
      assertion5(tarch::la::equals(pLeft.getCorrectorTimeStepSize(),0.0) || std::isfinite(FL[i]),pLeft.toString(),faceIndexLeft,normalDirection,i,FL[i]);
      assertion5(tarch::la::equals(pRight.getCorrectorTimeStepSize(),0.0) || std::isfinite(FR[i]),pRight.toString(),faceIndexRight,normalDirection,i,FR[i]);
    }
    #endif
    
    riemannSolver( // TODO(Dominic): Merge Riemann solver directly with the face integral and push the result on update
                   // does not make sense to overwrite the flux when performing local time stepping; coarse grid flux must be constant, or not?
        FL,FR,QL,QR,
        std::min(pLeft.getCorrectorTimeStepSize(),
            pRight.getCorrectorTimeStepSize()),
        normalDirection, false, -1);
    
    #if defined(Debug) || defined(Asserts)
    for(int i=0; i<dofPerFace; ++i) {
      assertion8(tarch::la::equals(pLeft.getCorrectorTimeStepSize(),0.0) || (std::isfinite(FL[i]) && std::isfinite(FR[i])),
                 pLeft.toString(),faceIndexLeft,pRight.toString(),faceIndexRight,normalDirection,i,FL[i],FR[i]);
    }  
    #endif

    // directly perform the face integral afterwards
    int levelDeltaLeft  = 0;
    int levelDeltaRight = 0;
    tarch::la::Vector<DIMENSIONS-1,int> subfaceIndexLeft(0);
    tarch::la::Vector<DIMENSIONS-1,int> subfaceIndexRight(0);

    const int orientationLeft  = faceIndexLeft % 2;
    const int orientationRight = 1-orientationLeft;
    const int direction        = (faceIndexLeft-orientationLeft)/2;
    if ( pLeft.getType()==CellDescription::Type::Descendant ) {
      levelDeltaLeft = pLeft.getLevel() - pLeft.getParentCellLevel();

      const tarch::la::Vector<DIMENSIONS,int> subcellIndex =
          exahype::amr::computeSubcellIndex(
              pLeft.getOffset(),pLeft.getSize(),
                pLeft.getParentOffset());

      subfaceIndexLeft = exahype::amr::getSubfaceIndex(subcellIndex,direction);
    }
    if (  pRight.getType()==CellDescription::Type::Descendant ) {
      levelDeltaRight = pRight.getLevel() - pRight.getParentCellLevel();

      const tarch::la::Vector<DIMENSIONS,int> subcellIndex =
          exahype::amr::computeSubcellIndex(
              pRight.getOffset(),pRight.getSize(),
              pRight.getParentOffset());

      subfaceIndexRight = exahype::amr::getSubfaceIndex(subcellIndex,direction);
    }
    DataHeap::HeapEntries& updateLeft  = DataHeap::getInstance().getData(pLeft.getUpdate());
    DataHeap::HeapEntries& updateRight = DataHeap::getInstance().getData(pRight.getUpdate());

    faceIntegral(updateLeft.data(),FL,direction,orientationLeft,subfaceIndexLeft,levelDeltaLeft,pLeft.getSize());
    faceIntegral(updateRight.data(),FR,direction,orientationRight,subfaceIndexRight,levelDeltaRight,pRight.getSize());
  }
}

void exahype::solvers::ADERDGSolver::mergeWithBoundaryData(
    const int                                 cellDescriptionsIndex,
    const int                                 element,
    const tarch::la::Vector<DIMENSIONS, int>& posCell,
    const tarch::la::Vector<DIMENSIONS, int>& posBoundary) {
  if (tarch::la::countEqualEntries(posCell,posBoundary)!=(DIMENSIONS-1)) {
    return; // We only consider faces; no corners.
  }

  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  synchroniseTimeStepping(cellDescription);

  if (cellDescription.getType()==CellDescription::Type::Cell) {
    const int direction   = tarch::la::equalsReturnIndex(posCell, posBoundary);
    const int orientation = (1 + posBoundary(direction) - posCell(direction))/2;
    const int faceIndex   = 2*direction+orientation;

    uncompress(cellDescription);

    applyBoundaryConditions(cellDescription,faceIndex);

    mergeWithAugmentationStatus(cellDescription,faceIndex,BoundaryStatus);
    mergeWithCommunicationStatus(cellDescription,faceIndex,BoundaryStatus);
    mergeWithLimiterStatus(cellDescription,faceIndex,BoundaryStatus);
  }
}

void exahype::solvers::ADERDGSolver::applyBoundaryConditions(CellDescription& p,const int faceIndex) {
  assertion1(p.getType()==CellDescription::Type::Cell,p.toString());
  assertion1(p.getRefinementEvent()==CellDescription::None,p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getExtrapolatedPredictor()),p.toString());
  assertion1(DataHeap::getInstance().isValidIndex(p.getFluctuation()),p.toString());

  const int dataPerFace = getBndFaceSize();
  const int dofPerFace  = getBndFluxSize();
  double* QIn = DataHeap::getInstance().getData(p.getExtrapolatedPredictor()).data() +
      (faceIndex * dataPerFace);
  double* FIn = DataHeap::getInstance().getData(p.getFluctuation()).data() +
      (faceIndex * dofPerFace);

  double* update = DataHeap::getInstance().getData(p.getUpdate()).data();

  const int orientation = faceIndex % 2;
  const int direction   = (faceIndex - orientation)/2;
  #if defined(Debug) || defined(Asserts)
  assertion2(direction<DIMENSIONS,faceIndex,direction);
  for(int i=0; i<dataPerFace; ++i) {
    assertion5(std::isfinite(QIn[i]), p.toString(),faceIndex, direction, i, QIn[i]);
  }
  for(int i=0; i<dofPerFace; ++i) {
    assertion5(std::isfinite(FIn[i]), p.toString(),faceIndex, direction, i, FIn[i]);
  } 
  #endif

  // TODO(Dominic): Hand in space-time volume data. Time integrate it afterwards
  boundaryConditions(
      update,
      FIn,QIn,
      p.getOffset() + 0.5*p.getSize(),
      p.getSize(),
      p.getCorrectorTimeStamp(),
      p.getCorrectorTimeStepSize(),
      direction,orientation);

  #if defined(Debug) || defined(Asserts)
  assertion4(std::isfinite(p.getCorrectorTimeStamp()),p.toString(),faceIndex,direction,p.getCorrectorTimeStamp());
  assertion4(std::isfinite(p.getCorrectorTimeStepSize()),p.toString(),faceIndex,direction,p.getCorrectorTimeStepSize());
  assertion4(p.getCorrectorTimeStepSize()>=0.0, p.toString(),faceIndex, direction,p.getCorrectorTimeStepSize());
  for(int i=0; i<dofPerFace; ++i) {
    assertion5(tarch::la::equals(p.getCorrectorTimeStepSize(),0.0) || std::isfinite(FIn[i]),p.toString(),faceIndex,direction,i,FIn[i]);
  }
  #endif
}

#ifdef Parallel
const int exahype::solvers::ADERDGSolver::DataMessagesPerNeighbourCommunication    = 2;
const int exahype::solvers::ADERDGSolver::DataMessagesPerForkOrJoinCommunication   = 2;
const int exahype::solvers::ADERDGSolver::DataMessagesPerMasterWorkerCommunication = 2;

/**
 * After the forking the master's cell descriptions
 * are not accessed by enterCell(...) on the master
 * anymore. However we still use them as buffer for saving data.
 */
bool exahype::solvers::ADERDGSolver::sendCellDescriptions(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const bool                                    fromWorkerSide,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  if ( isValidCellDescriptionIndex(cellDescriptionsIndex) ) {
    logDebug("sendCellDescriptions(...)","send "<< Heap::getInstance().getData(cellDescriptionsIndex).size()<<
        " cell descriptions to rank "<<toRank<<" (x="<< x.toString() << ",level="<< level << ")");
    bool oneSolverRequiresVerticalCommunication = false;
    for (auto& cellDescription : Heap::getInstance().getData(cellDescriptionsIndex)) {
      if ( fromWorkerSide ) {
        prepareWorkerCellDescriptionAtMasterWorkerBoundary(cellDescription);
      }
    }
    Heap::getInstance().sendData(cellDescriptionsIndex,toRank,x,level,messageType);
    return oneSolverRequiresVerticalCommunication;
  }
  else {
    logDebug("sendCellDescriptions(...)","send "
        " empty cell descriptions to rank "<<toRank<<" (x="<< x.toString() << ",level="<< level << ")");
    
    sendEmptyCellDescriptions(toRank,messageType,x,level);
    return false;
  }
}

void exahype::solvers::ADERDGSolver::sendEmptyCellDescriptions(
    const int                                     toRank,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  logDebug("sendEmptyCellDescriptions(...)","send empty message to " <<
          "rank "<<toRank <<
          " at (center="<< x.toString() <<
          ",level="<< level << ")");

  Heap::HeapEntries emptyMessage(0);
  Heap::getInstance().sendData(emptyMessage,
      toRank,x,level,messageType);
}

void exahype::solvers::ADERDGSolver::receiveCellDescriptions(
    const int                                    fromRank,
    exahype::Cell&                               localCell,
    const peano::heap::MessageType&              messageType,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  Heap::getInstance().receiveData(
      localCell.getCellDescriptionsIndex(),fromRank,x,level,messageType);

  logDebug("mergeCellDescriptionsWithRemoteData(...)","received " <<
          Heap::getInstance().getData(localCell.getCellDescriptionsIndex()).size() <<
          " cell descriptions for cell (centre="<< x.toString() << ", level="<< level << ")");

  for (auto& cellDescription : Heap::getInstance().getData(localCell.getCellDescriptionsIndex())) {
    resetIndicesAndFlagsOfReceivedCellDescription(
        cellDescription,multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
  }
}

void exahype::solvers::ADERDGSolver::resetIndicesAndFlagsOfReceivedCellDescription(
    CellDescription& cellDescription,const int parentIndex) {
  cellDescription.setParentIndex(parentIndex);

  cellDescription.setAdjacentToRemoteRank(false);
  cellDescription.setFaceDataExchangeCounter(0);

  // Default field data indices
  cellDescription.setSolution(-1);
  cellDescription.setPreviousSolution(-1);
  cellDescription.setUpdate(-1);
  cellDescription.setExtrapolatedPredictor(-1);
  cellDescription.setFluctuation(-1);

  // Limiter meta data (oscillations identificator)
  cellDescription.setSolutionMin(-1);
  cellDescription.setSolutionMax(-1);

  // compression
  cellDescription.setCompressionState(CellDescription::CompressionState::Uncompressed);

  cellDescription.setExtrapolatedPredictorCompressed(-1);
  cellDescription.setFluctuationCompressed(-1);
  cellDescription.setSolutionCompressed(-1);
  cellDescription.setPreviousSolutionCompressed(-1);
  cellDescription.setUpdateCompressed(-1);

  cellDescription.setSolutionAverages(-1);
  cellDescription.setSolutionAverages(-1);
  cellDescription.setUpdateAverages(-1);
  cellDescription.setExtrapolatedPredictorAverages(-1);
  cellDescription.setFluctuationAverages(-1);
  cellDescription.setBytesPerDoFInPreviousSolution(-1);
  cellDescription.setBytesPerDoFInSolution(-1);
  cellDescription.setBytesPerDoFInUpdate(-1);
  cellDescription.setBytesPerDoFInExtrapolatedPredictor(-1);
  cellDescription.setBytesPerDoFInFluctuation(-1);

  // reset the facewise flags
  cellDescription.setFacewiseAugmentationStatus(0);
  cellDescription.setFacewiseCommunicationStatus(0);
  cellDescription.setFacewiseLimiterStatus(0);

  // limiter flagging
  cellDescription.setIterationsToCureTroubledCell(-1);
}

void exahype::solvers::ADERDGSolver::ensureOnlyNecessaryMemoryIsAllocated(CellDescription& cellDescription) {
  auto* solver = RegisteredSolvers[cellDescription.getSolverNumber()];
  switch (solver->getType()) {
  case exahype::solvers::Solver::Type::ADERDG:
    static_cast<ADERDGSolver*>(solver)->ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
    static_cast<ADERDGSolver*>(solver)->ensureNecessaryMemoryIsAllocated(cellDescription);
    break;
  case exahype::solvers::Solver::Type::LimitingADERDG:
    static_cast<LimitingADERDGSolver*>(solver)->
    getSolver()->ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
    static_cast<LimitingADERDGSolver*>(solver)->
        getSolver()->ensureNecessaryMemoryIsAllocated(cellDescription);
    break;
  case exahype::solvers::Solver::Type::FiniteVolumes:
    assertionMsg(false,"Solver type not supported!");
    break;
  }
}

/**
 * Drop cell descriptions received from \p fromRank.
 */
void exahype::solvers::ADERDGSolver::dropCellDescriptions(
    const int                                     fromRank,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  Heap::getInstance().receiveData(fromRank,x,level,messageType);
}

////////////////////////////////////
// MASTER <=> WORKER
////////////////////////////////////
void
exahype::solvers::ADERDGSolver::appendMasterWorkerCommunicationMetadata(
    MetadataHeap::HeapEntries& metadata,
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  const int element = tryGetElement(cellDescriptionsIndex,solverNumber);

  if (element!=exahype::solvers::Solver::NotFound)  {
    CellDescription& cellDescription =
        getCellDescription(cellDescriptionsIndex,element);
    metadata.push_back(static_cast<int>(cellDescription.getType()));
    metadata.push_back(cellDescription.getAugmentationStatus()); // TODO(Dominic): Add to docu: Might be merged multiple times!
    metadata.push_back(cellDescription.getCommunicationStatus());
    metadata.push_back(cellDescription.getLimiterStatus());
    metadata.push_back(
        (cellDescription.getHasToHoldDataForMasterWorkerCommunication()) ? 1 : 0 );
  } else {
    for (int i = 0; i < exahype::MasterWorkerCommunicationMetadataPerSolver; ++i) {
      metadata.push_back(exahype::InvalidMetadataEntry); // implicit conversion
    }
  }
}

void exahype::solvers::ADERDGSolver::prepareWorkerCellDescriptionAtMasterWorkerBoundary(
    CellDescription& cellDescription) {
  if ( 
     cellDescription.getType()==CellDescription::Type::Cell ||
     cellDescription.getType()==CellDescription::Type::Descendant
  ) {
    cellDescription.setHasToHoldDataForMasterWorkerCommunication(cellDescription.getHasVirtualChildren());
  }
}

void exahype::solvers::ADERDGSolver::deduceChildCellErasingEvents(CellDescription& cellDescription) const {
  const int coarseGridElement = tryGetElement(cellDescription.getParentIndex(),cellDescription.getSolverNumber());
  if ( coarseGridElement!=exahype::solvers::Solver::NotFound ) {
    CellDescription& coarseGridCellDescription =
        getCellDescription(cellDescription.getParentIndex(),coarseGridElement);

    tarch::multicore::Lock lock(CoarseGridSemaphore);

    switch (coarseGridCellDescription.getRefinementEvent()) {
    case CellDescription::RefinementEvent::ErasingChildrenRequested: {
      assertion1(coarseGridCellDescription.getType()==CellDescription::Type::Ancestor,
          coarseGridCellDescription.toString());
      cellDescription.setRefinementEvent(CellDescription::RefinementEvent::Erasing);
    } break;
    case CellDescription::RefinementEvent::ErasingVirtualChildrenRequested: {
      assertion1(coarseGridCellDescription.getType()==CellDescription::Type::Cell ||
          coarseGridCellDescription.getType()==CellDescription::Type::Descendant,coarseGridCellDescription.toString());
      cellDescription.setRefinementEvent(CellDescription::RefinementEvent::ErasingVirtualCell);
    }  break;
    case CellDescription::RefinementEvent::ChangeChildrenToVirtualChildrenRequested: {
      assertion1(coarseGridCellDescription.getType()==CellDescription::Type::Ancestor,coarseGridCellDescription.toString());
      cellDescription.setRefinementEvent(CellDescription::RefinementEvent::ChangeToVirtualCell);
    } break;
    default:
      break;
    }
    lock.free();
  }
}

void exahype::solvers::ADERDGSolver::progressMeshRefinementInPrepareSendToWorker(
    const int workerRank,
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const bool initialGrid,
    const int solverNumber) {
  // coarse grid based operations
  const int coarseGridCellDescriptionsIndex = coarseGridCell.getCellDescriptionsIndex();
  const int coarseGridCellElement = tryGetElement(coarseGridCellDescriptionsIndex,solverNumber);
  if (coarseGridCellElement!=exahype::solvers::Solver::NotFound) {
    CellDescription& coarseGridCellDescription = getCellDescription(
        coarseGridCell.getCellDescriptionsIndex(),coarseGridCellElement);

    addNewDescendantIfVirtualRefiningRequested(
        fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
        coarseGridCellDescription,coarseGridCell.getCellDescriptionsIndex());

    bool addedNewCell =
        addNewCellIfRefinementRequested(
            fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
            coarseGridCellDescription,coarseGridCell.getCellDescriptionsIndex(),
            initialGrid);
    
    const int cellDescriptionsIndex = fineGridCell.getCellDescriptionsIndex();
    const int element = tryGetElement(cellDescriptionsIndex,solverNumber);
    
    if ( element!=exahype::solvers::Solver::NotFound ) {
      CellDescription& fineGridCellDescription = getCellDescription(cellDescriptionsIndex,element);
      if ( 
        fineGridCellDescription.getType()==CellDescription::Type::Descendant &&
        fineGridCellDescription.getHasToHoldDataForMasterWorkerCommunication()
      ) {
        exahype::solvers::Solver::SubcellPosition subcellPosition =
            exahype::amr::computeSubcellPositionOfDescendant<CellDescription,Heap,true>(fineGridCellDescription);
        CellDescription& topMostParentCellDescription = 
            getCellDescription(subcellPosition.parentCellDescriptionsIndex,subcellPosition.parentElement);
        if ( topMostParentCellDescription.getType()==CellDescription::Type::Cell ) {
           logDebug( "progressMeshRefinementInPrepareSendToWorker(...)"," try to refine parent " << topMostParentCellDescription.toString());
           topMostParentCellDescription.setRefinementRequest(CellDescription::RefinementRequest::Refine);
        }
      }
    }

    if ( addedNewCell ) {
      CellDescription& fineGridCellDescription = getCellDescription(cellDescriptionsIndex,element);
      prolongateVolumeData(fineGridCellDescription,initialGrid);
      assertion1( fineGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::Prolongating,
        fineGridCellDescription.toString());
    } 
  }
}

void exahype::solvers::ADERDGSolver::sendDataToWorkerIfProlongating(
    const int                                     workerRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  CellDescription& fineGridCellDescription = getCellDescription(cellDescriptionsIndex,element);

  // send out the data
  if ( fineGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::Prolongating ) {
    logDebug( "sendDataToWorkerIfProlongating(...)","send prolongated solution to rank "<<workerRank<< " at x="<<x.toString()<< ",level="<<level << " cell="<<fineGridCellDescription.toString());

    sendDataToWorkerOrMasterDueToForkOrJoin(workerRank,cellDescriptionsIndex,element,
        peano::heap::MessageType::MasterWorkerCommunication,x,level);
  }
}

void exahype::solvers::ADERDGSolver::receiveDataFromMasterIfProlongating(
  const int masterRank,
	const int receivedCellDescriptionsIndex,
  const int receivedElement,
  const tarch::la::Vector<DIMENSIONS,double>& x,
  const int level) const {
  CellDescription& receivedCellDescription = getCellDescription(receivedCellDescriptionsIndex,receivedElement);

  if ( receivedCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::Prolongating ) {
    logDebug( "receiveDataFromMasterIfProlongating(...)","receiving prolongated solution from rank "<<masterRank<< " at x="<<x.toString()<< ",level="<<level);

    mergeWithWorkerOrMasterDataDueToForkOrJoin(
      masterRank,receivedCellDescriptionsIndex,receivedElement,
      peano::heap::MessageType::MasterWorkerCommunication,x,level);
  }
}

void exahype::solvers::ADERDGSolver::ensureSameNumberOfMasterAndWorkerCellDescriptions(
    exahype::Cell& localCell,
    const exahype::Cell& receivedMasterCell) {
  for (CellDescription& pReceived : Heap::getInstance().getData(receivedMasterCell.getCellDescriptionsIndex())) {
    bool found = false;
    for (CellDescription& pLocal : Heap::getInstance().getData(localCell.getCellDescriptionsIndex())) {
      found |= pReceived.getSolverNumber()==pLocal.getSolverNumber();
    }
    if ( !found ) {
      Heap::getInstance().getData(localCell.getCellDescriptionsIndex()).push_back(pReceived); // this copies
    }
  }
}

void exahype::solvers::ADERDGSolver::progressMeshRefinementInMergeWithWorker(
    const int localCellDescriptionsIndex,    const int localElement,
    const int receivedCellDescriptionsIndex, const int receivedElement,
    const bool initialGrid) {
  CellDescription& localCellDescription    = getCellDescription(localCellDescriptionsIndex,localElement);
  CellDescription& receivedCellDescription = getCellDescription(receivedCellDescriptionsIndex,receivedElement); // need because of data

  // finalise prolongation operation started on master
  if ( receivedCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::Prolongating ) {
    logDebug( "progressMeshRefinementInMergeWithWorker(...)","merging prolongated solution");

    assertion( localCellDescription.getType()==CellDescription::Type::Cell ||
               localCellDescription.getType()==CellDescription::Type::Descendant);
    assertion(receivedCellDescription.getType()==CellDescription::Type::Cell);
    assertion(DataHeap::getInstance().isValidIndex(receivedCellDescription.getSolution()));
    assertion(DataHeap::getInstance().isValidIndex(receivedCellDescription.getPreviousSolution()));

    // we know we have received data in this case
    localCellDescription.setType(CellDescription::Type::Cell);
    localCellDescription.setRefinementEvent(CellDescription::RefinementEvent::Prolongating);
    localCellDescription.setRefinementRequest(CellDescription::RefinementRequest::Pending);
    localCellDescription.setCommunicationStatus(CellCommunicationStatus);
    localCellDescription.setFacewiseCommunicationStatus(0); // implicit conversion

    ensureNecessaryMemoryIsAllocated(localCellDescription); // TODO could simply copy index
    std::copy(
        DataHeap::getInstance().getData(receivedCellDescription.getSolution()).begin(),
        DataHeap::getInstance().getData(receivedCellDescription.getSolution()).end(),
        DataHeap::getInstance().getData(localCellDescription.getSolution()).begin());
    std::copy(
        DataHeap::getInstance().getData(receivedCellDescription.getPreviousSolution()).begin(),
        DataHeap::getInstance().getData(receivedCellDescription.getPreviousSolution()).end(),
        DataHeap::getInstance().getData(localCellDescription.getPreviousSolution()).begin());

    // adjust solution
    localCellDescription.setRefinementEvent(CellDescription::RefinementEvent::None);
    Solver::adjustSolutionDuringMeshRefinement(
        localCellDescriptionsIndex,localElement,initialGrid);
  }
}

void exahype::solvers::ADERDGSolver::progressMeshRefinementInPrepareSendToMaster(
    const int masterRank,
    const int cellDescriptionsIndex, const int element,
    const tarch::la::Vector<DIMENSIONS,double>& x,
    const int level) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  // send out data
  if (
      cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::Erasing ||
      cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ChangeToVirtualCell
  ) {
    sendDataToWorkerOrMasterDueToForkOrJoin(masterRank,cellDescriptionsIndex,element,
        peano::heap::MessageType::MasterWorkerCommunication,x,level); // assumes blocking/copy
  }

  // erase or change type of cell descriptions
  if (
      cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::Erasing ||
      cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ErasingVirtualCell
  ) {
    cellDescription.setType(CellDescription::Type::Erased);
    ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
    getCellDescriptions(cellDescriptionsIndex).erase(
        getCellDescriptions(cellDescriptionsIndex).begin()+element);
  } else if ( cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ChangeToVirtualCell ) {
    assertion(cellDescription.getType()==CellDescription::Type::Cell)
            cellDescription.setType(CellDescription::Type::Descendant);
    ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
    ensureNecessaryMemoryIsAllocated(cellDescription);
  }
}

bool exahype::solvers::ADERDGSolver::progressMeshRefinementInMergeWithMaster(
    const int worker,
    const int localCellDescriptionsIndex,
    const int localElement,
    const int coarseGridCellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  CellDescription& cellDescription = getCellDescription(localCellDescriptionsIndex,localElement);
  cellDescription.setParentIndex(coarseGridCellDescriptionsIndex);

  // receive the data
  if (
      cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::Erasing ||
      cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ChangeToVirtualCell
  ) {
    mergeWithWorkerOrMasterDataDueToForkOrJoin(worker,localCellDescriptionsIndex,localElement,
        peano::heap::MessageType::MasterWorkerCommunication,x,level); // assumes blocking/copy
  }

  // work with the data
  if (
      cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::Erasing ||
      cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ChangeToVirtualCell ||
      cellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ErasingVirtualCell
  ) {
    // use the received data
    const int coarseGridElement = tryGetElement(
        cellDescription.getParentIndex(),cellDescription.getSolverNumber());
    assertion1(coarseGridElement!=exahype::solvers::Solver::NotFound,cellDescription.toString());
    CellDescription& coarseGridCellDescription =
        getCellDescription(cellDescription.getParentIndex(),coarseGridElement); // TODO(Dominic): Have helper function for that
    assertion2( coarseGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ErasingChildren ||
        coarseGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ChangeChildrenToVirtualChildren ||
        coarseGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ErasingVirtualChildren,
        cellDescription.toString(),coarseGridCellDescription.toString());

    eraseCellDescriptionIfNecessary(localCellDescriptionsIndex,localElement,coarseGridCellDescription);
  }
  

  progressCollectiveRefinementOperationsInLeaveCell(cellDescription);
  // ignore return value as responsibiliy is still on fine grid.

  // check if any cell description requires vertical communication
  bool solverRequiresVerticalCommunication = false;
 
  ensureConsistencyOfParentInformation(fineGridCellDescription,coarseGridCell.getCellDescriptionsIndex());
 
  if (
      cellDescription.getType()==CellDescription::Ancestor ||
      cellDescription.getType()==CellDescription::Type::Cell
  ) {
    // copy and restrict the limiter status
    const int coarseGridElement = tryGetElement(
        cellDescription.getParentIndex(),cellDescription.getSolverNumber());
    if (coarseGridElement!=exahype::solvers::Solver::NotFound) {
      restrictToNextParent(cellDescription,coarseGridElement);
    }
  }

  // block erasing request of coarse grid cell description if deployed cell
  // does not want to be erased
  decideOnRefinement(cellDescription);

  return solverRequiresVerticalCommunication;
}

///////////////////////////////////
// FORK OR JOIN
///////////////////////////////////

void exahype::solvers::ADERDGSolver::sendDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  
  assertion5(tarch::la::equals(x,cellDescription.getOffset()+0.5*cellDescription.getSize()),x,cellDescription.getOffset()+0.5*cellDescription.getSize(),level,cellDescription.getLevel(),tarch::parallel::Node::getInstance().getRank());
  assertion2(cellDescription.getLevel()==level,cellDescription.getLevel(),level);

  if ( cellDescription.getType()==CellDescription::Type::Cell ) {
    logDebug("sendDataToWorkerOrMasterDueToForkOrJoin(...)",""
            "solution of solver " << cellDescription.getSolverNumber() << " sent to rank "<<toRank<<
                 ", cell: "<< x << ", level: " << level);

    assertion2(DataHeap::getInstance().isValidIndex(cellDescription.getSolution()),
        cellDescriptionsIndex,cellDescription.toString());
    assertion2(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolution()),
            cellDescriptionsIndex,cellDescription.toString());
    DataHeap::getInstance().sendData(
        DataHeap::getInstance().getData(cellDescription.getSolution()).data(),
        getDataPerCell(), toRank, x, level,messageType);
    DataHeap::getInstance().sendData(
        DataHeap::getInstance().getData(cellDescription.getPreviousSolution()).data(),
        getDataPerCell(), toRank, x, level,messageType);
  }
}

void exahype::solvers::ADERDGSolver::mergeWithWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  assertion5(tarch::la::equals(x,cellDescription.getOffset()+0.5*cellDescription.getSize()),x,cellDescription.getOffset()+0.5*cellDescription.getSize(),level,cellDescription.getLevel(),tarch::parallel::Node::getInstance().getRank());
  assertion2(cellDescription.getLevel()==level,cellDescription.getLevel(),level);

  // allocate memory
  if ( messageType==peano::heap::MessageType::ForkOrJoinCommunication ) {
    ensureNecessaryMemoryIsAllocated(cellDescription);
    ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
  } else if ( cellDescription.getType()==CellDescription::Type::Cell ) {
    ensureNecessaryMemoryIsAllocated(cellDescription);
    ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
  }

  // receive data
  if ( cellDescription.getType()==CellDescription::Type::Cell ) {
    logDebug("mergeWithRemoteDataDueToForkOrJoin(...)","[solution] receive from rank "<<fromRank<<
             ", cell: "<< x << ", level: " << level);

    DataHeap::getInstance().getData(cellDescription.getSolution()).clear();
    DataHeap::getInstance().getData(cellDescription.getPreviousSolution()).clear();
    DataHeap::getInstance().receiveData(cellDescription.getSolution(),
        fromRank,x,level,messageType);
    DataHeap::getInstance().receiveData(cellDescription.getPreviousSolution(),
        fromRank,x,level,messageType);
  }
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void
exahype::solvers::ADERDGSolver::appendNeighbourCommunicationMetadata(
    exahype::MetadataHeap::HeapEntries& metadata,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  const int element = tryGetElement(cellDescriptionsIndex,solverNumber);

  if (element!=exahype::solvers::Solver::NotFound)  {
    CellDescription& cellDescription =
        getCellDescription(cellDescriptionsIndex,element);

    metadata.push_back(static_cast<int>(cellDescription.getType()));
    metadata.push_back(cellDescription.getAugmentationStatus()); // TODO(Dominic): Add to docu: Might be merged multiple times!
    metadata.push_back(cellDescription.getCommunicationStatus());
    metadata.push_back(cellDescription.getLimiterStatus());
  } else {
    for (int i = 0; i < exahype::NeighbourCommunicationMetadataPerSolver; ++i) {
      metadata.push_back(exahype::InvalidMetadataEntry); // implicit conversion
    }
  }
}

void exahype::solvers::ADERDGSolver::mergeWithNeighbourMetadata(
    const exahype::MetadataHeap::HeapEntries& neighbourMetadata,
    const tarch::la::Vector<DIMENSIONS, int>& src,
    const tarch::la::Vector<DIMENSIONS, int>& dest,
    const int                                 cellDescriptionsIndex,
    const int                                 element) const {
  if (tarch::la::countEqualEntries(src,dest)==DIMENSIONS-1) { // only consider faces
    const int direction   = tarch::la::equalsReturnIndex(src, dest);
    const int orientation = (1 + src(direction) - dest(direction))/2;
    const int faceIndex   = 2*direction+orientation;

    const int neighbourAugmentationStatus =
        neighbourMetadata[exahype::NeighbourCommunicationMetadataAugmentationStatus];
    const int neighbourCommunicationStatus       =
        neighbourMetadata[exahype::NeighbourCommunicationMetadataCommunicationStatus      ];
    const int neighbourLimiterStatus      =
        neighbourMetadata[exahype::NeighbourCommunicationMetadataLimiterStatus   ];

    CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

    mergeWithAugmentationStatus (cellDescription,faceIndex,neighbourAugmentationStatus);
    mergeWithCommunicationStatus(cellDescription,faceIndex,neighbourCommunicationStatus);
    mergeWithLimiterStatus      (cellDescription,faceIndex,neighbourLimiterStatus);
  }
}

void exahype::solvers::ADERDGSolver::sendDataToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  assertion( tarch::la::countEqualEntries(src,dest)==(DIMENSIONS-1) );

  const int direction    = tarch::la::equalsReturnIndex(src, dest);
  const int orientation  = (1 + dest(direction) - src(direction))/2;
  const int faceIndex    = 2*direction+orientation;

  CellDescription& cellDescription = Heap::getInstance().getData(cellDescriptionsIndex)[element];
  if (
      cellDescription.getCommunicationStatus()>=MinimumCommunicationStatusForNeighbourCommunication &&
      cellDescription.getAugmentationStatus() < MaximumAugmentationStatus // excludes Ancestors
  ) {
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictor()));
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()));

    const int dofPerFace  = getBndFluxSize();
    const int dataPerFace = getBndFaceSize();

    const double* lQhbnd = DataHeap::getInstance().getData(
        cellDescription.getExtrapolatedPredictor()).data() +
        (faceIndex * dataPerFace);
    const double* lFhbnd = DataHeap::getInstance().getData(
        cellDescription.getFluctuation()).data() +
        (faceIndex * dofPerFace);

/*
    #ifdef Asserts
    tarch::la::Vector<DIMENSIONS,double> faceBarycentre =
            exahype::Cell::computeFaceBarycentre(
                cellDescription.getOffset(),cellDescription.getSize(),direction,orientation);
    logInfo("sendDataToNeighbour(...)", "send "<<DataMessagesPerNeighbourCommunication<<" msgs to rank " <<
            toRank << " vertex="<<x.toString()<<" face=" << faceBarycentre.toString());
    #endif
*/

    // Send order: lQhbnd,lFhbnd,observablesMin,observablesMax
    // Receive order: observablesMax,observablesMin,lFhbnd,lQhbnd
    DataHeap::getInstance().sendData(
        lQhbnd, dataPerFace, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
    DataHeap::getInstance().sendData(
        lFhbnd, dofPerFace, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);

    // TODO(Dominic): If anarchic time stepping send the time step over too.
  } else {
    sendEmptyDataToNeighbour(toRank,x,level);
  }
}

void exahype::solvers::ADERDGSolver::sendEmptyDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  // Send order: lQhbnd,lFhbnd,observablesMin,observablesMax
  // Receive order: observablesMax,observablesMin,lFhbnd,lQhbnd
  // TODO(WORKAROUND)
  #if defined(UsePeanosSymmetricBoundaryExchanger)
  const int dofPerFace  = getBndFluxSize();
  const int dataPerFace = getBndFaceSize();
  DataHeap::getInstance().sendData(
      _invalidExtrapolatedPredictor.data(), dataPerFace, toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
  DataHeap::getInstance().sendData(
      _invalidFluctuations.data(), dofPerFace, toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
  #else
  for(int sends=0; sends<DataMessagesPerNeighbourCommunication; ++sends)
    DataHeap::getInstance().sendData(
        exahype::EmptyDataHeapMessage, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
  #endif
}

// TODO(Dominic): Add to docu: We only perform a Riemann solve if a Cell is involved.
void exahype::solvers::ADERDGSolver::mergeWithNeighbourData(
    const int                                    fromRank,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  assertionEquals(tarch::la::countEqualEntries(src,dest),DIMENSIONS-1); // only faces

  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  synchroniseTimeStepping(cellDescription);

  const int direction    = tarch::la::equalsReturnIndex(src, dest);
  const int orientation  = (1 + src(direction) - dest(direction))/2;
  const int faceIndex    = 2*direction+orientation;
  if(
      (cellDescription.getCommunicationStatus()                ==CellCommunicationStatus &&
      cellDescription.getFacewiseCommunicationStatus(faceIndex)>=MinimumCommunicationStatusForNeighbourCommunication &&
      cellDescription.getFacewiseAugmentationStatus(faceIndex) < MaximumAugmentationStatus)
      ||
      (cellDescription.getFacewiseCommunicationStatus(faceIndex)==CellCommunicationStatus &&
      cellDescription.getCommunicationStatus()                  >=MinimumCommunicationStatusForNeighbourCommunication &&
      cellDescription.getAugmentationStatus()                   < MaximumAugmentationStatus)
  ){
    assertion3(cellDescription.getNeighbourMergePerformed(faceIndex),
        faceIndex,cellDescriptionsIndex,cellDescription.toString());
/*
    #ifdef Asserts
    tarch::la::Vector<DIMENSIONS,double> faceBarycentre =
        exahype::Cell::computeFaceBarycentre(
            cellDescription.getOffset(),cellDescription.getSize(),direction,orientation);
    logInfo("mergeWithNeighbourData(...)", "receive "<<DataMessagesPerNeighbourCommunication<<" msgs from rank " <<
        fromRank << " vertex="<<x.toString()<<" face=" << faceBarycentre.toString());
    #endif
*/
    // Send order: lQhbnd,lFhbnd
    // Receive order: lFhbnd,lQhbnd
    // TODO(Dominic): If anarchic time stepping, receive the time step too.
    const int dofPerFace  = getBndFluxSize();
    const int dataPerFace = getBndFaceSize();
    DataHeap::getInstance().receiveData(
        const_cast<double*>(_receivedFluctuations.data()),dofPerFace, // TODO const-correct peano
        fromRank, x, level,peano::heap::MessageType::NeighbourCommunication);
    DataHeap::getInstance().receiveData(                              // TODO const-correct peano
        const_cast<double*>(_receivedExtrapolatedPredictor.data()),dataPerFace,
        fromRank, x, level, peano::heap::MessageType::NeighbourCommunication);

    solveRiemannProblemAtInterface(
        cellDescription,
        faceIndex,
        _receivedExtrapolatedPredictor.data(),
        _receivedFluctuations.data(),
        fromRank);
  } else  {
    dropNeighbourData(fromRank,src,dest,x,level);
  }
}

void exahype::solvers::ADERDGSolver::solveRiemannProblemAtInterface(
    records::ADERDGCellDescription& cellDescription,
    const int faceIndex,
    const double* const lQhbnd,
    const double* lFhbnd,
    const int fromRank) {
  assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictor()));
  assertion(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()));

  const int dataPerFace = getBndFaceSize();
  const int dofPerFace  = getBndFluxSize();

  logDebug("solveRiemannProblemAtInterface(...)",
      "cell-description=" << cellDescription.toString());

  // @todo Doku im Header warum wir das hier brauchen,
  const int orientation = faceIndex % 2;
  const int direction   = (faceIndex-orientation)/2;
  if ( orientation==0 ) {
    const double* const QL = lQhbnd;
    const double* const QR = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
        (faceIndex * dataPerFace);
    double* FL = const_cast<double*>(lFhbnd); // TODO const-correct kernels
    double* FR = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
        (faceIndex * dofPerFace); // TODO const-correct kernels
    
    riemannSolver(
        FL, FR, QL, QR,
        cellDescription.getCorrectorTimeStepSize(),direction,false,faceIndex);
    
    #if defined(Debug) || defined(Asserts)
    for (int ii = 0; ii<dataPerFace; ii++) {
      assertion8(std::isfinite(QR[ii]), cellDescription.toString(),
          faceIndex, direction, ii, QR[ii], QL[ii], FR[ii], FL[ii]);
      assertion8(std::isfinite(QL[ii]), cellDescription.toString(),
          faceIndex, direction, ii, QR[ii], QL[ii], FR[ii], FL[ii]);
    }
    
    for (int ii = 0; ii<dofPerFace; ii++) {
      assertion8(std::isfinite(FL[ii]), cellDescription.toString(),
          faceIndex, ii, QR[ii], QL[ii], FR[ii], FL[ii],fromRank);
      assertion8(std::isfinite(FR[ii]), cellDescription.toString(),
          faceIndex, ii, QR[ii], QL[ii], FR[ii], FL[ii],fromRank);
    }
    #endif
  } else {
    const double* const QR = lQhbnd;
    const double* const QL = DataHeap::getInstance().getData(cellDescription.getExtrapolatedPredictor()).data() +
        (faceIndex * dataPerFace);
    double* FR = const_cast<double*>(lFhbnd); // TODO const-correct kernels
    double* FL = DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
        (faceIndex * dofPerFace); // TODO const-correct kernels
    
    riemannSolver(
        FL, FR, QL, QR,
        cellDescription.getCorrectorTimeStepSize(),direction,false,faceIndex);
    
    #if defined(Debug) || defined(Asserts)
    for (int ii = 0; ii<dataPerFace; ii++) {
      assertion8(std::isfinite(QR[ii]), cellDescription.toString(),
          faceIndex, direction, ii, QR[ii], QL[ii], FR[ii], FL[ii]);
      assertion8(std::isfinite(QL[ii]), cellDescription.toString(),
          faceIndex, direction, ii, QR[ii], QL[ii], FR[ii], FL[ii]);
    }

    for (int ii = 0; ii<dofPerFace; ii++) {
      assertion8(std::isfinite(FL[ii]), cellDescription.toString(),
          faceIndex, ii, QR[ii], QL[ii], FR[ii], FL[ii],fromRank);
      assertion8(std::isfinite(FR[ii]), cellDescription.toString(),
          faceIndex, ii, QR[ii], QL[ii], FR[ii], FL[ii],fromRank);
    }
    #endif
  }
  
  // directly perform the face integral
  int levelDelta= 0;
  tarch::la::Vector<DIMENSIONS-1,int> subfaceIndex(0);

  if ( cellDescription.getType()==CellDescription::Type::Descendant ) {
    levelDelta   = cellDescription.getLevel() - cellDescription.getParentCellLevel();

    const tarch::la::Vector<DIMENSIONS,int> subcellIndex =
        exahype::amr::computeSubcellIndex(
            cellDescription.getOffset(),cellDescription.getSize(),
              cellDescription.getParentOffset());

    subfaceIndex = exahype::amr::getSubfaceIndex(subcellIndex,direction);
  }
  double* update = DataHeap::getInstance().getData(cellDescription.getUpdate()).data();
  const double* const boundaryFlux =
      DataHeap::getInstance().getData(cellDescription.getFluctuation()).data() +
      (faceIndex * dofPerFace);
  faceIntegral(update,boundaryFlux,direction,orientation,subfaceIndex,levelDelta,cellDescription.getSize());
}

void exahype::solvers::ADERDGSolver::dropNeighbourData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  logDebug(
      "dropNeighbourData(...)", "drop "<<DataMessagesPerNeighbourCommunication<<" arrays from rank " <<
      fromRank << " for vertex x=" << x << ", level=" << level <<
      ", src=" << src << ", dest=" << dest
  );

  for(int receives=0; receives<DataMessagesPerNeighbourCommunication; ++receives)
    DataHeap::getInstance().receiveData(
        fromRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
}

///////////////////////////////////
// WORKER->MASTER
///////////////////////////////////
exahype::DataHeap::HeapEntries
exahype::solvers::ADERDGSolver::compileMessageForMaster(const int capacity) const {
  DataHeap::HeapEntries messageForMaster(0,std::max(3,capacity));
  messageForMaster.push_back(_minPredictorTimeStepSize);
  messageForMaster.push_back(_maxLevel);
  messageForMaster.push_back(_meshUpdateRequest ? 1.0 : -1.0);
  return messageForMaster;
}

/*
 * At the time of sending data to the master,
 * we have already performed a time step update locally
 * on the rank. We thus need to communicate the
 * current min predictor time step size to the master.
 * The next min predictor time step size is
 * already reset locally to the maximum double value.
 *
 * However on the master's side, we need to
 * merge the received time step size with
 * the next min predictor time step size since
 * the master has not yet performed a time step update
 * (i.e. called TimeStepSizeComputation::endIteration()).
 */
void exahype::solvers::ADERDGSolver::sendDataToMaster(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  DataHeap::HeapEntries messageForMaster = compileMessageForMaster();

  assertion1(messageForMaster.size()==3,messageForMaster.size());
  assertion1(std::isfinite(messageForMaster[0]),messageForMaster[0]);
  if (_timeStepping==TimeStepping::Global) {
    assertionNumericalEquals1(_minNextTimeStepSize,std::numeric_limits<double>::max(),
                                tarch::parallel::Node::getInstance().getRank());
  }

  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToMaster(...)","Sending time step data: " <<
             "data[0]=" << messageForMaster[0] <<
             ",data[1]=" << messageForMaster[1] <<
             ",data[2]=" << messageForMaster[2] <<
             " to rank " << masterRank <<
             ", message size="<<messageForMaster.size()
    );
  }

  DataHeap::getInstance().sendData(
      messageForMaster.data(), messageForMaster.size(),
      masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::ADERDGSolver::mergeWithWorkerData(const DataHeap::HeapEntries& message) {
  assertion1(message[0]>=0,message[0]);
  assertion1(std::isfinite(message[0]),message[0]);
  // The master solver has not yet updated its minNextPredictorTimeStepSize.
  // Thus it does not equal MAX_DOUBLE.

  int index=0;
  _minNextTimeStepSize    = std::min( _minNextTimeStepSize, message[index++] );
  _nextMaxLevel           = std::max( _nextMaxLevel,        static_cast<int>(message[index++]) );
  _nextMeshUpdateRequest |= (message[index++]) > 0 ? true : false;

  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","[post] Receiving time step data: " <<
        "data[0]=" << message[0] <<
        ",data[1]=" << message[1] <<
        ",data[2]=" << message[2] );
    logDebug("mergeWithWorkerData(...)","[post] Updated time step fields: " <<
        ",_minNextPredictorTimeStepSize=" << _minNextTimeStepSize <<
        ",_nextMeshUpdateRequest=" << _nextMeshUpdateRequest <<
        ",_nextMaxLevel=" << _nextMaxLevel);
  }
}

/**
 * At the time of the merging,
 * the workers and the master have already performed
 * at local update of the next predictor time step size
 * and of the predictor time stamp.
 * We thus need to minimise over both quantities.
 */
void exahype::solvers::ADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  DataHeap::HeapEntries messageFromWorker(3);

  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","Receiving time step data [pre] from rank " << workerRank);
  }

  DataHeap::getInstance().receiveData(
      messageFromWorker.data(),messageFromWorker.size(),workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);

  assertion1(messageFromWorker.size()==3,messageFromWorker.size());
  mergeWithWorkerData(messageFromWorker);
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
exahype::DataHeap::HeapEntries
exahype::solvers::ADERDGSolver::compileMessageForWorker(const int capacity) const {
  DataHeap::HeapEntries messageForWorker(0,std::max(7,capacity));
  messageForWorker.push_back(_minCorrectorTimeStamp);
  messageForWorker.push_back(_minCorrectorTimeStepSize);
  messageForWorker.push_back(_minPredictorTimeStamp);
  messageForWorker.push_back(_minPredictorTimeStepSize);

  messageForWorker.push_back(_maxLevel);

  messageForWorker.push_back(_meshUpdateRequest ? 1.0 : -1.0);

  messageForWorker.push_back(_stabilityConditionWasViolated ? 1.0 : -1.0);

  assertion1(messageForWorker.size()==7,messageForWorker.size());
  assertion1(std::isfinite(messageForWorker[0]),messageForWorker[0]);
  assertion1(std::isfinite(messageForWorker[1]),messageForWorker[1]);
  assertion1(std::isfinite(messageForWorker[2]),messageForWorker[2]);
  assertion1(std::isfinite(messageForWorker[3]),messageForWorker[3]);

  if (_timeStepping==TimeStepping::Global) {
    assertionEquals1(_minNextTimeStepSize,std::numeric_limits<double>::max(),
                     tarch::parallel::Node::getInstance().getRank());
  }

  return messageForWorker;
}

void exahype::solvers::ADERDGSolver::sendDataToWorker(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  DataHeap::HeapEntries messageForWorker = compileMessageForWorker();

  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug(
        "sendDataToWorker(...)","Broadcasting time step data: " <<
        " data[0]=" << messageForWorker[0] <<
        ",data[1]=" << messageForWorker[1] <<
        ",data[2]=" << messageForWorker[2] <<
        ",data[3]=" << messageForWorker[3] <<
        ",data[4]=" << messageForWorker[4] <<
        ",data[5]=" << messageForWorker[5] <<
        ",data[6]=" << messageForWorker[6]);
    logDebug("sendDataWorker(...)","_minNextPredictorTimeStepSize="<<_minNextTimeStepSize);
  }

  DataHeap::getInstance().sendData(
      messageForWorker.data(), messageForWorker.size(),
      workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::ADERDGSolver::mergeWithMasterData(const DataHeap::HeapEntries& message) {
  int index=0;
  _minCorrectorTimeStamp         = message[index++];
  _minCorrectorTimeStepSize      = message[index++];
  _minPredictorTimeStamp         = message[index++];
  _minPredictorTimeStepSize      = message[index++];

  _maxLevel                      = message[index++];

  _meshUpdateRequest             = (message[index++] > 0.0) ? true : false;
  _stabilityConditionWasViolated = (message[index++] > 0.0) ? true : false;

  logDebug("mergeWithMasterData(...)",
      "_stabilityConditionWasViolated="<< _stabilityConditionWasViolated);
}

void exahype::solvers::ADERDGSolver::mergeWithMasterData(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  DataHeap::HeapEntries messageFromMaster(7);
  DataHeap::getInstance().receiveData(
      messageFromMaster.data(),messageFromMaster.size(),masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);

  assertion1(messageFromMaster.size()==7,messageFromMaster.size());
  mergeWithMasterData(messageFromMaster);

  if (_timeStepping==TimeStepping::Global) {
    assertionNumericalEquals1(_minNextTimeStepSize,std::numeric_limits<double>::max(),
                                  _minNextTimeStepSize);
  }

  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug(
        "mergeWithMasterData(...)","Received time step data: " <<
        "data[0]="  << messageFromMaster[0] <<
        ",data[1]=" << messageFromMaster[1] <<
        ",data[2]=" << messageFromMaster[2] <<
        ",data[3]=" << messageFromMaster[3] <<
        ",data[4]=" << messageFromMaster[4] <<
        ",data[5]=" << messageFromMaster[5] <<
        ",data[6]=" << messageFromMaster[6]);
  }
}
#endif

std::string exahype::solvers::ADERDGSolver::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::solvers::ADERDGSolver::toString (std::ostream& out) const {
  out << "(";
  out << "_identifier:" << _identifier;
  out << ",";
  out << "_type:" << exahype::solvers::Solver::toString(_type);
  out << ",";
  out << "_numberOfVariables:" << _numberOfVariables;
  out << ",";
  out << "_numberOfParameters:" << _numberOfParameters;
  out << ",";
  out << "_nodesPerCoordinateAxis:" << _nodesPerCoordinateAxis;
  out << ",";
  out << "_maximumMeshSize:" << _maximumMeshSize;
  out << ",";
  out << "_timeStepping:" << exahype::solvers::Solver::toString(_timeStepping); // only solver attributes
  out << ",";
  out << "_unknownsPerFace:" << getUnknownsPerFace();
  out << ",";
  out << "_unknownsPerCellBoundary:" << getUnknownsPerCellBoundary();
  out << ",";
  out << "_unknownsPerCell:" << getUnknownsPerCell();
  out << ",";
  out << "_fluxUnknownsPerCell:" << getFluxUnknownsPerCell();
  out << ",";
  out << "_spaceTimeUnknownsPerCell:" << getSpaceTimeUnknownsPerCell();
  out << ",";
  out << "_spaceTimeFluxUnknownsPerCell:" << getSpaceTimeFluxUnknownsPerCell();
  out << ",";
  out << "_previousMinCorrectorTimeStepSize:" << _previousMinCorrectorTimeStepSize;
  out << ",";
  out << "_minCorrectorTimeStamp:" << _minCorrectorTimeStamp;
  out << ",";
  out << "_minCorrectorTimeStepSize:" << _minCorrectorTimeStepSize;
  out << ",";
  out << "_minPredictorTimeStepSize:" << _minPredictorTimeStepSize;
  out << ",";
  out << "_minNextPredictorTimeStepSize:" << _minNextTimeStepSize;
  out <<  ")";
}


exahype::solvers::ADERDGSolver::PredictionJob::PredictionJob(
  ADERDGSolver&     solver,
  const int         cellDescriptionsIndex,
  const int         element,
  const double      predictorTimeStamp,
  const double      predictorTimeStepSize,
  const bool        uncompressBefore,
  const bool        isSkeletonJob):
  _solver(solver),
  _cellDescriptionsIndex(cellDescriptionsIndex),
  _element(element),
  _predictorTimeStamp(predictorTimeStamp),
  _predictorTimeStepSize(predictorTimeStepSize),
  _uncompressBefore(uncompressBefore),
  _isSkeletonJob(isSkeletonJob) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter++;
  }
  lock.free();
}


bool exahype::solvers::ADERDGSolver::PredictionJob::operator()() {
  _solver.performPredictionAndVolumeIntegralBody(
      _cellDescriptionsIndex,_element,
      _predictorTimeStamp,_predictorTimeStepSize,
      _uncompressBefore,_isSkeletonJob); // ignore return value

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter--;
    assertion( jobCounter>=0 );
  }
  lock.free();
  return false;
}


exahype::solvers::ADERDGSolver::FusedTimeStepJob::FusedTimeStepJob(
  ADERDGSolver& solver,
  const int     cellDescriptionsIndex,
  const int     element,
  const bool    isSkeletonJob):
  _solver(solver),
  _cellDescriptionsIndex(cellDescriptionsIndex),
  _element(element),
  _isSkeletonJob(isSkeletonJob) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter++;
  }
  lock.free();
}

bool exahype::solvers::ADERDGSolver::FusedTimeStepJob::operator()() {
  _solver.fusedTimeStepBody(
      _cellDescriptionsIndex,_element, false, false, _isSkeletonJob, false );

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter--;
    assertion( jobCounter>=0 );
  }
  lock.free();
  return false;
}


exahype::solvers::ADERDGSolver::CompressionJob::CompressionJob(
  const ADERDGSolver& solver,
  CellDescription&    cellDescription,
  const bool          isSkeletonJob):
  _solver(solver),
  _cellDescription(cellDescription),
  _isSkeletonJob(isSkeletonJob) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter++;
  }
  lock.free();
}


bool exahype::solvers::ADERDGSolver::CompressionJob::operator()() {
  _solver.determineUnknownAverages(_cellDescription);
  _solver.computeHierarchicalTransform(_cellDescription,-1.0);
  _solver.putUnknownsIntoByteStream(_cellDescription);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter--;
    assertion( jobCounter>=0 );
  }
  lock.free();
  return false;
}


void exahype::solvers::ADERDGSolver::compress( CellDescription& cellDescription, const bool isSkeletonCell ) const {
  assertion1( cellDescription.getCompressionState() ==  CellDescription::Uncompressed, cellDescription.toString() );
  if (CompressionAccuracy>0.0) {
    if ( SpawnCompressionAsBackgroundJob ) {
      int& jobCounter = ( isSkeletonCell ) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
      cellDescription.setCompressionState(CellDescription::CurrentlyProcessed);
      CompressionJob compressionJob( *this, cellDescription, jobCounter );
      if ( isSkeletonCell ) {
        peano::datatraversal::TaskSet spawnedSet( compressionJob, peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible  );
      } else {
        peano::datatraversal::TaskSet spawnedSet( compressionJob, peano::datatraversal::TaskSet::TaskType::Background  );
      }
    }
    else {
      determineUnknownAverages(cellDescription);
      computeHierarchicalTransform(cellDescription,-1.0);
      putUnknownsIntoByteStream(cellDescription);
      cellDescription.setCompressionState(CellDescription::Compressed);
    }
  }
}


void exahype::solvers::ADERDGSolver::uncompress(CellDescription& cellDescription) const {
  #ifdef SharedMemoryParallelisation
  bool madeDecision = CompressionAccuracy<=0.0;
  bool uncompress   = false;

  while (!madeDecision) {
    peano::datatraversal::TaskSet::finishToProcessBackgroundJobs();

    tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
    madeDecision = cellDescription.getCompressionState() != CellDescription::CurrentlyProcessed;
    uncompress   = cellDescription.getCompressionState() == CellDescription::Compressed;
    if (uncompress) {
      cellDescription.setCompressionState( CellDescription::CurrentlyProcessed );
    }
    lock.free();
  }
  #else
  bool uncompress = CompressionAccuracy>0.0
      && cellDescription.getCompressionState() == CellDescription::Compressed;
  #endif

/*
  #ifdef Parallel
  assertion1(!cellDescription.getAdjacentToRemoteRank() || cellDescription.getCompressionState() == CellDescription::Compressed,
             cellDescription.toString());
  #endif
*/

  if (uncompress) {
    pullUnknownsFromByteStream(cellDescription);
    computeHierarchicalTransform(cellDescription,1.0);

    tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
      cellDescription.setCompressionState(CellDescription::Uncompressed);
    lock.free();
  }
}


void exahype::solvers::ADERDGSolver::determineUnknownAverages(
  CellDescription& cellDescription) const {
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getSolution()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolution()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getUpdate()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictor()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getFluctuation()), cellDescription.toString() );

  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getSolutionAverages()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionAverages()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getUpdateAverages()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorAverages()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationAverages()), cellDescription.toString() );

  const int dataPerNode  = getNumberOfParameters()+getNumberOfVariables();
  const int nodesPerCell = getDataPerCell()/ dataPerNode;
  const int nodesPerFace = getDataPerFace() / dataPerNode;

  auto& solutionAverages              = DataHeap::getInstance().getData( cellDescription.getSolutionAverages() );
  auto& previousSolutionAverage       = DataHeap::getInstance().getData( cellDescription.getPreviousSolutionAverages() );
  auto& updateAverages                = DataHeap::getInstance().getData( cellDescription.getUpdateAverages() );
  auto& extrapolatedPredictorAverages = DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictorAverages() );
  auto& fluctuationAverages           = DataHeap::getInstance().getData( cellDescription.getFluctuationAverages() );

  // patch data
  kernels::idx2 idx_cellData    (nodesPerCell,dataPerNode);
  kernels::idx2 idx_cellUnknowns(nodesPerCell,getNumberOfVariables());
  for (int i=0; i<nodesPerCell; i++) {
    for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
      solutionAverages[variableNumber]        += DataHeap::getInstance().getData( cellDescription.getSolution() )        [idx_cellData(i,variableNumber)];
      previousSolutionAverage[variableNumber] += DataHeap::getInstance().getData( cellDescription.getPreviousSolution() )[idx_cellData(i,variableNumber)];
    }
    for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
      updateAverages[variableNumber]          += DataHeap::getInstance().getData( cellDescription.getUpdate() )[idx_cellUnknowns(i,variableNumber)];
    }
  }
  for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
    solutionAverages[variableNumber]        = solutionAverages[variableNumber]        / (double) nodesPerCell;
    previousSolutionAverage[variableNumber] = previousSolutionAverage[variableNumber] / (double) nodesPerCell;
  }
  for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
    updateAverages[variableNumber]          = updateAverages[variableNumber]          / (double) nodesPerCell;
  }

  // face data
  kernels::idx2 idx_faceDataAvg    (DIMENSIONS_TIMES_TWO,dataPerNode);
  kernels::idx2 idx_faceUnknownsAvg(DIMENSIONS_TIMES_TWO,getNumberOfVariables());
  kernels::idx3 idx_faceData       (DIMENSIONS_TIMES_TWO,nodesPerFace,dataPerNode);
  kernels::idx3 idx_faceUnknowns   (DIMENSIONS_TIMES_TWO,nodesPerFace,getNumberOfVariables());
  for (int face=0; face<2*DIMENSIONS; face++) {
    for (int i=0; i<nodesPerFace; i++) {
      for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
        extrapolatedPredictorAverages[idx_faceDataAvg(face,variableNumber)] +=
            DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictor() )[idx_faceData(face,i,variableNumber)];
      }
      for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
        fluctuationAverages[idx_faceUnknownsAvg(face,variableNumber)] +=
            DataHeap::getInstance().getData( cellDescription.getFluctuation() )[idx_faceUnknowns(face,i,variableNumber)];
      }
    }
    for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
      extrapolatedPredictorAverages[idx_faceDataAvg(face,variableNumber)] =
          extrapolatedPredictorAverages[idx_faceDataAvg(face,variableNumber)] / (double) nodesPerFace;
    }
    for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
      fluctuationAverages[idx_faceUnknownsAvg(face,variableNumber)] =
          fluctuationAverages[idx_faceUnknownsAvg(face,variableNumber)]       / (double) nodesPerFace;
    }
  }
}


void exahype::solvers::ADERDGSolver::computeHierarchicalTransform(
    CellDescription& cellDescription, double sign) const {
  const int dataPerNode  = getNumberOfParameters()+getNumberOfVariables();
  const int nodesPerCell = getDataPerCell()/ dataPerNode;
  const int nodesPerFace = getDataPerFace() / dataPerNode;

  // patch data
  kernels::idx2 idx_cellData    (nodesPerCell,dataPerNode);
  kernels::idx2 idx_cellUnknowns(nodesPerCell,getNumberOfVariables());
  for (int i=0; i<nodesPerCell; i++) {
    for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) { // variables+parameters
      DataHeap::getInstance().getData( cellDescription.getSolution() )
                [idx_cellData(i,variableNumber)] += sign * DataHeap::getInstance().getData( cellDescription.getSolutionAverages() )[variableNumber];
      DataHeap::getInstance().getData( cellDescription.getPreviousSolution() )
                [idx_cellData(i,variableNumber)] += sign * DataHeap::getInstance().getData( cellDescription.getPreviousSolutionAverages() )[variableNumber];
    }
    for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) { // variables
      DataHeap::getInstance().getData( cellDescription.getUpdate() )
                [idx_cellUnknowns(i,variableNumber)] += sign * DataHeap::getInstance().getData( cellDescription.getUpdateAverages() )[variableNumber];
    }
  }

  // face data
  kernels::idx2 idx_faceDataAvg    (DIMENSIONS_TIMES_TWO,dataPerNode);
  kernels::idx2 idx_faceUnknownsAvg(DIMENSIONS_TIMES_TWO,getNumberOfVariables());
  kernels::idx3 idx_faceData       (DIMENSIONS_TIMES_TWO,nodesPerFace,dataPerNode);
  kernels::idx3 idx_faceUnknowns   (DIMENSIONS_TIMES_TWO,nodesPerFace,getNumberOfVariables());
  for (int face=0; face<DIMENSIONS_TIMES_TWO; face++) {
    for (int i=0; i<nodesPerFace; i++) {
      for (int variableNumber=0; variableNumber<dataPerNode; variableNumber++) {  // variables+parameters
        DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictor() )                           [idx_faceData(face,i,variableNumber)] +=
                    sign * DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictorAverages() )[idx_faceDataAvg(face,variableNumber)];
      }
      for (int variableNumber=0; variableNumber<getNumberOfVariables(); variableNumber++) {  // variables
        DataHeap::getInstance().getData( cellDescription.getFluctuation() )                                     [idx_faceUnknowns(face,i,variableNumber)] +=
                              sign * DataHeap::getInstance().getData( cellDescription.getFluctuationAverages() )[idx_faceUnknownsAvg(face,variableNumber)];
      }
    }
  }
}

void exahype::solvers::ADERDGSolver::putUnknownsIntoByteStream(
    CellDescription& cellDescription) const {
  assertion(CompressionAccuracy>0.0);

  assertion( cellDescription.getPreviousSolutionCompressed()==-1 );
  assertion( cellDescription.getSolutionCompressed()==-1 );
  assertion( cellDescription.getUpdateCompressed()==-1 );
  assertion( cellDescription.getExtrapolatedPredictorCompressed()==-1 );
  assertion( cellDescription.getFluctuationCompressed()==-1 );

  int compressionOfPreviousSolution;
  int compressionOfSolution;
  int compressionOfUpdate;
  int compressionOfExtrapolatedPredictor;
  int compressionOfFluctuation;

  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getUpdate() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictor() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getFluctuation() ));

  peano::datatraversal::TaskSet compressionFactorIdentification(
    [&]() -> bool { compressionOfPreviousSolution = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getPreviousSolution() ).data(),
      getDataPerCell(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&] () -> bool  { compressionOfSolution = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getSolution() ).data(),
      getDataPerCell(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&]() -> bool  { compressionOfUpdate = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getUpdate() ).data(),
      getUnknownsPerCell(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&]() -> bool  { compressionOfExtrapolatedPredictor = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictor() ).data(),
      getDataPerCellBoundary(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&]() -> bool  { compressionOfFluctuation = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getFluctuation() ).data(),
      getUnknownsPerCellBoundary(),
      CompressionAccuracy,true
      );
      return false;
      },
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
    true
  );

  assertion(1<=compressionOfPreviousSolution);
  assertion(1<=compressionOfSolution);
  assertion(1<=compressionOfUpdate);
  assertion(1<=compressionOfExtrapolatedPredictor);
  assertion(1<=compressionOfFluctuation);

  assertion(compressionOfPreviousSolution<=7);
  assertion(compressionOfSolution<=7);
  assertion(compressionOfUpdate<=7);
  assertion(compressionOfExtrapolatedPredictor<=7);
  assertion(compressionOfFluctuation<=7);

  peano::datatraversal::TaskSet runParallelTasks(
    [&]() -> bool {
      cellDescription.setBytesPerDoFInPreviousSolution(compressionOfPreviousSolution);
      if (compressionOfPreviousSolution<7) {
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          cellDescription.setPreviousSolutionCompressed( CompressedDataHeap::getInstance().createData(0,0) );
          assertion( cellDescription.getPreviousSolutionCompressed()>=0 );
        lock.free();

        const int numberOfEntries = getDataPerCell();
        tearApart(numberOfEntries, cellDescription.getPreviousSolution(), cellDescription.getPreviousSolutionCompressed(), compressionOfPreviousSolution);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getPreviousSolution() ).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getPreviousSolutionCompressed() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getPreviousSolution(), true );
          cellDescription.setPreviousSolution( -1 );
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getPreviousSolution() ).size() * 8.0;
          PipedCompressedBytes   += DataHeap::getInstance().getData( cellDescription.getPreviousSolution() ).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    [&]() -> bool {
      cellDescription.setBytesPerDoFInSolution(compressionOfSolution);
      if (compressionOfSolution<7) {
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          cellDescription.setSolutionCompressed(CompressedDataHeap::getInstance().createData(0,0));
          assertion1( cellDescription.getSolutionCompressed()>=0, cellDescription.getSolutionCompressed() );
        lock.free();

        const int numberOfEntries = getDataPerCell();

        tearApart(numberOfEntries, cellDescription.getSolution(), cellDescription.getSolutionCompressed(), compressionOfSolution);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getSolution() ).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getSolutionCompressed() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getSolution(), true );
          cellDescription.setSolution( -1 );
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getSolution() ).size() * 8.0;
          PipedCompressedBytes   += DataHeap::getInstance().getData( cellDescription.getSolution() ).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    [&]() -> bool {
      cellDescription.setBytesPerDoFInUpdate(compressionOfUpdate);
      if (compressionOfUpdate<7) {
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          cellDescription.setUpdateCompressed( CompressedDataHeap::getInstance().createData(0,0) );
          assertion( cellDescription.getUpdateCompressed()>=0 );
        lock.free();

        const int numberOfEntries = getUnknownsPerCell();
        tearApart(numberOfEntries, cellDescription.getUpdate(), cellDescription.getUpdateCompressed(), compressionOfUpdate);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getUpdate() ).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getUpdateCompressed() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getUpdate(), true );
          cellDescription.setUpdate( -1 );
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getUpdate() ).size() * 8.0;
          PipedCompressedBytes   += DataHeap::getInstance().getData( cellDescription.getUpdate() ).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    [&]() -> bool {
      cellDescription.setBytesPerDoFInExtrapolatedPredictor(compressionOfExtrapolatedPredictor);
      if (compressionOfExtrapolatedPredictor<7) {
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          cellDescription.setExtrapolatedPredictorCompressed( CompressedDataHeap::getInstance().createData(0,0) );
          assertion( cellDescription.getExtrapolatedPredictorCompressed()>=0 );
        lock.free();

        const int numberOfEntries = getDataPerCellBoundary();
        tearApart(numberOfEntries, cellDescription.getExtrapolatedPredictor(), cellDescription.getExtrapolatedPredictorCompressed(), compressionOfExtrapolatedPredictor);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictor() ).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictorCompressed() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getExtrapolatedPredictor(), true );
          cellDescription.setExtrapolatedPredictor( -1 );
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictor() ).size() * 8.0;
          PipedCompressedBytes   += DataHeap::getInstance().getData( cellDescription.getExtrapolatedPredictor() ).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    [&]() -> bool {
      cellDescription.setBytesPerDoFInFluctuation(compressionOfFluctuation);
      if (compressionOfFluctuation<7) {
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          cellDescription.setFluctuationCompressed( CompressedDataHeap::getInstance().createData(0,0) );
          assertion( cellDescription.getFluctuationCompressed()>=0 );
        lock.free();

        const int numberOfEntries = getUnknownsPerCellBoundary();
        tearApart(numberOfEntries, cellDescription.getFluctuation(), cellDescription.getFluctuationCompressed(), compressionOfFluctuation);

        #if defined(TrackGridStatistics)
        lock.lock();
          PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getFluctuation() ).size() * 8.0;
          PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getFluctuationCompressed() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
          DataHeap::getInstance().deleteData( cellDescription.getFluctuation(), true );
          cellDescription.setFluctuation( -1 );
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getFluctuation() ).size() * 8.0;
          PipedCompressedBytes   += DataHeap::getInstance().getData( cellDescription.getFluctuation() ).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
    true
  );
}


void exahype::solvers::ADERDGSolver::pullUnknownsFromByteStream(
    CellDescription& cellDescription) const {
  assertion(CompressionAccuracy>0.0);

  #if !defined(ValidateCompressedVsUncompressedData)
  const int dataPointsPerCell       = getDataPerCell();
  const int unknownsPerCellBoundary = getUnknownsPerCellBoundary();

  {
    tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
      cellDescription.setPreviousSolution( DataHeap::getInstance().createData( dataPointsPerCell, dataPointsPerCell ) );
      cellDescription.setSolution( DataHeap::getInstance().createData(         dataPointsPerCell, dataPointsPerCell ) );
      cellDescription.setUpdate( DataHeap::getInstance().createData(           getUpdateSize(),   getUpdateSize() ) );

      cellDescription.setExtrapolatedPredictor( DataHeap::getInstance().createData( unknownsPerCellBoundary, unknownsPerCellBoundary ) );
      cellDescription.setFluctuation( DataHeap::getInstance().createData(           unknownsPerCellBoundary, unknownsPerCellBoundary ) );
    lock.free();

    if (cellDescription.getPreviousSolution()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setPreviousSolution( DataHeap::getInstance().createData( dataPointsPerCell, dataPointsPerCell ) );
      lock.free();
    }
    if (cellDescription.getSolution()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setSolution( DataHeap::getInstance().createData( dataPointsPerCell, dataPointsPerCell ) );
      lock.free();
    }
    if (cellDescription.getUpdate()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setUpdate( DataHeap::getInstance().createData( getUpdateSize(), getUpdateSize() ) );
      lock.free();
    }
    if (cellDescription.getExtrapolatedPredictor()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setExtrapolatedPredictor( DataHeap::getInstance().createData(unknownsPerCellBoundary ) );
      lock.free();
    }
    if (cellDescription.getFluctuation()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
        cellDescription.setFluctuation( DataHeap::getInstance().createData( unknownsPerCellBoundary, unknownsPerCellBoundary ) );
      lock.free();
    }
  }
  #else
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getUpdate() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictor() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getFluctuation() ));
  #endif

  assertion1(
      CompressedDataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionCompressed() ),
      cellDescription.getPreviousSolutionCompressed()
    );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getSolutionCompressed() ),
    cellDescription.getSolutionCompressed()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getUpdateCompressed() ),
    cellDescription.getUpdateCompressed()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictorCompressed() ),
    cellDescription.getExtrapolatedPredictorCompressed()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getFluctuationCompressed() ),
    cellDescription.getFluctuationCompressed()
  );

  peano::datatraversal::TaskSet glueTasks(
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInPreviousSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolution() ), cellDescription.getPreviousSolution());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionCompressed() ));
        const int numberOfEntries = getDataPerCell();
        glueTogether(numberOfEntries, cellDescription.getPreviousSolution(), cellDescription.getPreviousSolutionCompressed(), cellDescription.getBytesPerDoFInPreviousSolution());
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getPreviousSolutionCompressed(), true );
          cellDescription.setPreviousSolutionCompressed( -1 );
        lock.free();
      }
      return false;
    },
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getSolution() ), cellDescription.getSolution() );
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getSolutionCompressed() ));
        const int numberOfEntries = getDataPerCell();
        glueTogether(numberOfEntries, cellDescription.getSolution(), cellDescription.getSolutionCompressed(), cellDescription.getBytesPerDoFInSolution());
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getSolutionCompressed(), true );
          cellDescription.setSolutionCompressed( -1 );
        lock.free();
      }
      return false;
    },
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInUpdate()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getUpdate() ), cellDescription.getUpdate());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getUpdateCompressed() ));
        const int numberOfEntries = getUnknownsPerCell();
        glueTogether(numberOfEntries, cellDescription.getUpdate(), cellDescription.getUpdateCompressed(), cellDescription.getBytesPerDoFInUpdate());
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getUpdateCompressed(), true );
          cellDescription.setUpdateCompressed( -1 );
        lock.free();
      }
      return false;
    },
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInExtrapolatedPredictor()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictor() ), cellDescription.getExtrapolatedPredictor());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedPredictorCompressed() ));
        const int numberOfEntries = getDataPerCellBoundary();
        glueTogether(numberOfEntries, cellDescription.getExtrapolatedPredictor(), cellDescription.getExtrapolatedPredictorCompressed(), cellDescription.getBytesPerDoFInExtrapolatedPredictor());
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getExtrapolatedPredictorCompressed(), true );
          cellDescription.setExtrapolatedPredictorCompressed( -1 );
        lock.free();
      }
      return false;
    },
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInFluctuation()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getFluctuation() ), cellDescription.getFluctuation());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getFluctuationCompressed() ));
        const int numberOfEntries = getUnknownsPerCellBoundary();
        glueTogether(numberOfEntries, cellDescription.getFluctuation(), cellDescription.getFluctuationCompressed(), cellDescription.getBytesPerDoFInFluctuation());
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
          CompressedDataHeap::getInstance().deleteData( cellDescription.getFluctuationCompressed(), true );
          cellDescription.setFluctuationCompressed( -1 );
        lock.free();
      }
      return false;
    },
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
    true
  );
}

