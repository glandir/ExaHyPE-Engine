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
 * \author Dominic E. Charrier, Tobias Weinzierl, Jean-Matthieu Gallard, Fabian Guera, Leonhard Rannabauer
 **/

// Definition of AMR-specific routines of the ADERDGSolver

#include "exahype/solvers/ADERDGSolver.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/amr/AdaptiveMeshRefinement.h"

// Allocate / deallocate

void exahype::solvers::ADERDGSolver::ensureNecessaryMemoryIsAllocated(
    CellDescription& cellDescription) const {
  // allocate solution
  if (
      isLeaf(cellDescription) &&
      !DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex())
  ) {
    assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()));

    tarch::multicore::Lock lock(exahype::HeapSemaphore);
    // Allocate volume DoF
    const int dataPerCell = getDataPerCell(); // Only the solution and previousSolution store material parameters
    cellDescription.setPreviousSolutionIndex( DataHeap::getInstance().createData( dataPerCell, dataPerCell ) );
    cellDescription.setSolutionIndex        ( DataHeap::getInstance().createData( dataPerCell, dataPerCell ) );
    checkDataHeapIndex(cellDescription,cellDescription.getPreviousSolutionIndex(),"getPreviousSolutionIndex()");
    checkDataHeapIndex(cellDescription,cellDescription.getSolutionIndex(),"getSolutionIndex()");
    cellDescription.setPreviousSolution( getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).data() ) ;
    cellDescription.setSolution        ( getDataHeapEntries(cellDescription.getSolutionIndex()).data() );
    std::fill_n(static_cast<double*>(cellDescription.getPreviousSolution()),getDataPerCell(),std::numeric_limits<double>::quiet_NaN());
    std::fill_n(static_cast<double*>(cellDescription.getSolution()),getDataPerCell(),std::numeric_limits<double>::quiet_NaN());

    cellDescription.setSolutionCompressedIndex(-1);
    cellDescription.setSolutionCompressed(nullptr);
    cellDescription.setPreviousSolutionCompressedIndex(-1);
    cellDescription.setPreviousSolutionCompressed(nullptr);

    const int dataPerNode = getNumberOfVariables()+getNumberOfParameters();
    cellDescription.setPreviousSolutionAveragesIndex( DataHeap::getInstance().createData( dataPerNode, dataPerNode ) );
    cellDescription.setSolutionAveragesIndex(         DataHeap::getInstance().createData( dataPerNode, dataPerNode ) );
    checkDataHeapIndex(cellDescription,cellDescription.getPreviousSolutionAveragesIndex(),"getPreviousSolutionAveragesIndex()");
    checkDataHeapIndex(cellDescription,cellDescription.getSolutionAveragesIndex(),"getSolutionAveragesIndex()");
    cellDescription.setPreviousSolutionAverages( getDataHeapEntries(cellDescription.getPreviousSolutionAveragesIndex()).data() ) ;
    cellDescription.setSolutionAverages        ( getDataHeapEntries(cellDescription.getSolutionAveragesIndex()).data() ) ;
    std::fill_n(static_cast<double*>(cellDescription.getPreviousSolutionAverages()),dataPerNode,std::numeric_limits<double>::quiet_NaN());
    std::fill_n(static_cast<double*>(cellDescription.getSolutionAverages()),dataPerNode,std::numeric_limits<double>::quiet_NaN());

    cellDescription.setCompressionState(CellDescription::Uncompressed);

    lock.free();
  }

  // allocate update and boundary arrays
  if (
      holdsFaceData(cellDescription) &&
      !DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex())
  ) {
    assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationIndex()));

    tarch::multicore::Lock lock(exahype::HeapSemaphore);

    // allocate update dof
    cellDescription.setUpdateIndex        ( DataHeap::getInstance().createData( getUpdateSize(), getUpdateSize() ) );
    cellDescription.setUpdateAveragesIndex( DataHeap::getInstance().createData( getNumberOfVariables(), getNumberOfVariables() ) );
    cellDescription.setUpdateCompressedIndex(-1);
    cellDescription.setUpdateCompressed(nullptr);
    checkDataHeapIndex(cellDescription,cellDescription.getUpdateIndex(),"getUpdate()");
    checkDataHeapIndex(cellDescription,cellDescription.getUpdateAveragesIndex(),"getUpdateAverages()");
    cellDescription.setUpdate        ( getDataHeapEntries(cellDescription.getUpdateIndex()).data() ) ;
    cellDescription.setUpdateAverages( getDataHeapEntries(cellDescription.getUpdateAveragesIndex()).data() ) ;
    // touch the memory
    std::fill_n(static_cast<double*>(cellDescription.getUpdate()),getUpdateSize(),std::numeric_limits<double>::quiet_NaN());
    std::fill_n(static_cast<double*>(cellDescription.getUpdateAverages()),getNumberOfVariables(),std::numeric_limits<double>::quiet_NaN());

    // extrapolated predictor
    const int dataPerBnd = getBndTotalSize();
    cellDescription.setExtrapolatedPredictorIndex( DataHeap::getInstance().createData(dataPerBnd, dataPerBnd) );
    cellDescription.setExtrapolatedPredictorCompressedIndex(-1);
    cellDescription.setExtrapolatedPredictorCompressed(nullptr);
    const int boundaryData = (getNumberOfParameters()+getNumberOfVariables()) * DIMENSIONS_TIMES_TWO; //TODO JMG / Dominic adapt for padding with optimized kernels //TODO Tobias: Does it make sense to pad these arrays.
    cellDescription.setExtrapolatedPredictorAveragesIndex( DataHeap::getInstance().createData( boundaryData,  boundaryData  ) );
    checkDataHeapIndex(cellDescription,cellDescription.getExtrapolatedPredictorIndex(),"getExtrapolatedPredictor()");
    checkDataHeapIndex(cellDescription,cellDescription.getExtrapolatedPredictorAveragesIndex(),"getExtrapolatedPredictorAverages()");
    cellDescription.setExtrapolatedPredictor        ( getDataHeapEntries(cellDescription.getExtrapolatedPredictorIndex()).data() ) ;
    cellDescription.setExtrapolatedPredictorAverages( getDataHeapEntries(cellDescription.getExtrapolatedPredictorAveragesIndex()).data() ) ;
    // touch the memory
    std::fill_n(static_cast<double*>(cellDescription.getExtrapolatedPredictor()),dataPerBnd,std::numeric_limits<double>::quiet_NaN());
    std::fill_n(static_cast<double*>(cellDescription.getExtrapolatedPredictorAverages()),boundaryData,std::numeric_limits<double>::quiet_NaN());

    if (isUseViscousFlux()) {
      // gradients of extrapolated predictor
      const int gradientSizePerBnd = _numberOfVariables * power(_nodesPerCoordinateAxis, DIMENSIONS - 1) * DIMENSIONS_TIMES_TWO * DIMENSIONS;
      cellDescription.setExtrapolatedPredictorGradientIndex( DataHeap::getInstance().createData(gradientSizePerBnd, gradientSizePerBnd) );
      cellDescription.setExtrapolatedPredictorGradient( getDataHeapEntries(cellDescription.getExtrapolatedPredictorGradientIndex()).data() ) ;
      // touch the memory
      std::fill_n(static_cast<double*>(cellDescription.getExtrapolatedPredictorGradient()),gradientSizePerBnd,std::numeric_limits<double>::quiet_NaN());
    } else {
      cellDescription.setExtrapolatedPredictorGradientIndex(-1);
    }

    //
    // fluctuations
    //
    const int dofPerBnd  = getBndFluxTotalSize();
    cellDescription.setFluctuationIndex( DataHeap::getInstance().createData(dofPerBnd,  dofPerBnd) );
    cellDescription.setFluctuationCompressedIndex(-1);
    cellDescription.setFluctuationCompressed(nullptr);
    const int boundaryUnknowns = getNumberOfVariables() * DIMENSIONS_TIMES_TWO;     //TODO JMG / Dominic adapt for padding with optimized kernels //TODO Tobias: Does it make sense to pad these arrays.
    cellDescription.setFluctuationAveragesIndex( DataHeap::getInstance().createData( boundaryUnknowns, boundaryUnknowns ) );
    checkDataHeapIndex(cellDescription,cellDescription.getFluctuationIndex(),"getFluctuation()");
    checkDataHeapIndex(cellDescription,cellDescription.getFluctuationAveragesIndex(),"getFluctuationAverages()");
    cellDescription.setFluctuation        ( getDataHeapEntries(cellDescription.getFluctuationIndex()).data() ) ;
    cellDescription.setFluctuationAverages( getDataHeapEntries(cellDescription.getFluctuationAveragesIndex()).data() ) ;
    // touch the memory
    std::fill_n(static_cast<double*>(cellDescription.getFluctuation()),dofPerBnd,std::numeric_limits<double>::quiet_NaN());
    std::fill_n(static_cast<double*>(cellDescription.getFluctuationAverages()),boundaryUnknowns,std::numeric_limits<double>::quiet_NaN());

    // Allocate volume DoF for limiter (we need for every of the 2*DIMENSIONS faces an array of min values
    // and array of max values of the neighbour at this face).
    const int numberOfObservables = getDMPObservables();
    if ( numberOfObservables>0 ) {
      cellDescription.setSolutionMinIndex(DataHeap::getInstance().createData(
          numberOfObservables * DIMENSIONS_TIMES_TWO, numberOfObservables * DIMENSIONS_TIMES_TWO ));
      cellDescription.setSolutionMaxIndex(DataHeap::getInstance().createData(
          numberOfObservables * DIMENSIONS_TIMES_TWO, numberOfObservables * DIMENSIONS_TIMES_TWO ));
      checkDataHeapIndex(cellDescription,cellDescription.getSolutionMinIndex(),"getSolutionMinIndex()");
      checkDataHeapIndex(cellDescription,cellDescription.getSolutionMaxIndex(),"getSolutionMaxIndex()");
      cellDescription.setSolutionMin( getDataHeapEntries(cellDescription.getSolutionMinIndex()).data() ) ;
      cellDescription.setSolutionMax( getDataHeapEntries(cellDescription.getSolutionMaxIndex()).data() ) ;
      // touch the memory
      std::fill_n(static_cast<double*>(cellDescription.getSolutionMin()),numberOfObservables * DIMENSIONS_TIMES_TWO,std::numeric_limits<double>::quiet_NaN());
      std::fill_n(static_cast<double*>(cellDescription.getSolutionMax()),numberOfObservables * DIMENSIONS_TIMES_TWO,std::numeric_limits<double>::quiet_NaN());
    }

    lock.free();
  }
}

void exahype::solvers::ADERDGSolver::ensureNoUnnecessaryMemoryIsAllocated(
    CellDescription& cellDescription) const {

  if (
      cellDescription.getType()!=CellDescription::Type::Leaf &&
      DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex())
  ) {
    tarch::multicore::Lock lock(exahype::HeapSemaphore);

    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()));
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()));

    if ( cellDescription.getSolutionIndex()>=0 ) {
      DataHeap::getInstance().deleteData(cellDescription.getSolutionIndex());
      assertion(cellDescription.getSolutionCompressedIndex()==-1);
    }
    else {
      assertion(CompressionAccuracy>0.0);
      assertion(cellDescription.getSolutionIndex()==-1);
      CompressedDataHeap::getInstance().deleteData(cellDescription.getSolutionCompressedIndex());
    }

    DataHeap::getInstance().deleteData(cellDescription.getSolutionAveragesIndex());
    DataHeap::getInstance().deleteData(cellDescription.getPreviousSolutionAveragesIndex());

    cellDescription.setPreviousSolutionIndex(-1);
    cellDescription.setPreviousSolution(nullptr);
    cellDescription.setSolutionIndex(-1);
    cellDescription.setSolution(nullptr);

    cellDescription.setPreviousSolutionAveragesIndex(-1);
    cellDescription.setPreviousSolutionAverages(nullptr);
    cellDescription.setSolutionAveragesIndex(-1);
    cellDescription.setSolutionAverages(nullptr);

    cellDescription.setPreviousSolutionCompressedIndex(-1);
    cellDescription.setPreviousSolutionCompressed(nullptr);
    cellDescription.setSolutionCompressedIndex(-1);
    cellDescription.setSolutionCompressed(nullptr);

    lock.free();
  }

  // deallocate update and boundary arrays
  if (
      !holdsFaceData(cellDescription) &&
      DataHeap::getInstance().isValidIndex(cellDescription.getUpdateIndex())
  ) {
    // update
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getUpdateIndex()));
    if ( cellDescription.getUpdateIndex()>=0 ) {
      assertion(cellDescription.getUpdateCompressedIndex()==-1);

      DataHeap::getInstance().deleteData(cellDescription.getUpdateIndex());
      cellDescription.setUpdateIndex(-1);
      cellDescription.setUpdate(nullptr);
    }
    else {
      assertion(CompressionAccuracy>0.0);
      assertion(cellDescription.getUpdateIndex()==-1);

      CompressedDataHeap::getInstance().deleteData(cellDescription.getUpdateCompressedIndex());
      cellDescription.setUpdateCompressedIndex(-1);
      cellDescription.setUpdateCompressed(nullptr);
    }
    DataHeap::getInstance().deleteData(cellDescription.getUpdateAveragesIndex());
    cellDescription.setUpdateAveragesIndex(-1);
    cellDescription.setUpdateAverages(nullptr);

    // extrapolated predictor
    tarch::multicore::Lock lock(exahype::HeapSemaphore);
    if ( cellDescription.getExtrapolatedPredictorIndex()>=0 ) {
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorIndex()));
      assertion(cellDescription.getExtrapolatedPredictorCompressedIndex()==-1);

      DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictorIndex());
      cellDescription.setExtrapolatedPredictorIndex(-1);
      cellDescription.setExtrapolatedPredictor(nullptr);
    }
    else {
      assertion(CompressionAccuracy>0.0);
      assertion(cellDescription.getExtrapolatedPredictorIndex()==-1);

      CompressedDataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictorCompressedIndex());
      cellDescription.setExtrapolatedPredictorCompressedIndex(-1);
      cellDescription.setExtrapolatedPredictorCompressed(nullptr);
    }
    DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictorAveragesIndex());
    cellDescription.setExtrapolatedPredictorAveragesIndex(-1);
    cellDescription.setExtrapolatedPredictorAverages(nullptr);

    // gradient of extrapolated predictor
    if (isUseViscousFlux() && cellDescription.getExtrapolatedPredictorGradientIndex() >= 0) {
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedPredictorGradientIndex()));

      DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedPredictorGradientIndex());
      cellDescription.setExtrapolatedPredictorGradientIndex(-1);
      cellDescription.setExtrapolatedPredictorGradient(nullptr);
    }

    // fluctuations
    if ( cellDescription.getFluctuationIndex()>=0 ) {
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getFluctuationIndex()));
      assertion(cellDescription.getFluctuationCompressedIndex()==-1);

      DataHeap::getInstance().deleteData(cellDescription.getFluctuationIndex());
      cellDescription.setFluctuationIndex(-1);
      cellDescription.setFluctuation(nullptr);
    }
    else {
      assertion(CompressionAccuracy>0.0);
      assertion(cellDescription.getFluctuationIndex()==-1);

      CompressedDataHeap::getInstance().deleteData(cellDescription.getFluctuationCompressedIndex());
      cellDescription.setFluctuationCompressedIndex(-1);
      cellDescription.setFluctuationCompressed(nullptr);
    }
    DataHeap::getInstance().deleteData(cellDescription.getFluctuationAveragesIndex());
    cellDescription.setFluctuationAveragesIndex(-1);
    cellDescription.setFluctuationAverages(nullptr);

    if ( getDMPObservables()>0 ) {
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionMinIndex()));
      assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionMaxIndex()));
      DataHeap::getInstance().deleteData(cellDescription.getSolutionMinIndex());
      DataHeap::getInstance().deleteData(cellDescription.getSolutionMaxIndex());

      cellDescription.setSolutionMinIndex(-1);
      cellDescription.setSolutionMin(nullptr);
      cellDescription.setSolutionMaxIndex(-1);
      cellDescription.setSolutionMax(nullptr);
    }

    lock.free();
  }
}

void exahype::solvers::ADERDGSolver::checkDataHeapIndex(const CellDescription& cellDescription, const int arrayIndex,const std::string arrayName) {
  assertion1(DataHeap::getInstance().isValidIndex(arrayIndex),cellDescription.toString());
  if ( arrayIndex < 0 ) {
    logError("checkDataHeapIndex(...)","The data heap array 'cellDescription."<<arrayName<<"' could not be allocated! Likely reason: Not enough memory available." <<
             " CellDescription="<<cellDescription.toString());
    std::abort();
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

// Operations

bool exahype::solvers::ADERDGSolver::progressMeshRefinementInEnterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const int  solverNumber,
    const bool inRefiningMode) {
  bool newComputeCell = false;

  // Fine grid cell based uniform mesh refinement.
  const int fineGridElement   = tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  const int coarseGridElement = tryGetElement(coarseGridCell.getCellDescriptionsIndex(),solverNumber);
  if (
      fineGridElement==exahype::solvers::Solver::NotFound &&
      fineGridVerticesEnumerator.getLevel()==_coarsestMeshLevel
  ) {
    logDebug("progressMeshRefinementInEnterCell(...)","Add new uniform grid cell at centre="<<fineGridVerticesEnumerator.getCellCenter() <<", level="<<fineGridVerticesEnumerator.getLevel()
        << " for solver=" << solverNumber);

    addNewCell(
        fineGridCell,fineGridVerticesEnumerator,
        multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex,
        solverNumber);
    newComputeCell = true;
  }
  else if ( fineGridElement!=exahype::solvers::Solver::NotFound ) {
    CellDescription& fineGridCellDescription = getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridElement);

    // wait for background jobs to complete
    waitUntilCompletedLastStep(fineGridCellDescription,false,false);

    #ifdef Asserts
    const tarch::la::Vector<DIMENSIONS,double> center = fineGridCellDescription.getOffset()+0.5*fineGridCellDescription.getSize();
    #endif
    assertion5(Vertex::equalUpToRelativeTolerance(fineGridVerticesEnumerator.getCellCenter(),center),
               fineGridVerticesEnumerator.getCellCenter(),center,fineGridVerticesEnumerator.getLevel(),fineGridCellDescription.getLevel(),tarch::parallel::Node::getInstance().getRank());
    assertionEquals3(fineGridVerticesEnumerator.getLevel(),fineGridCellDescription.getLevel(),fineGridVerticesEnumerator.getCellCenter(),fineGridCellDescription.getOffset()+0.5*fineGridCellDescription.getSize(),tarch::parallel::Node::getInstance().getRank());

    // Update the status flagging
    updateCommunicationStatus(fineGridCellDescription);
    ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
    ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);

    updateAugmentationStatus(fineGridCellDescription);

    updateRefinementStatus(fineGridCellDescription,fineGridCellDescription.getNeighbourMergePerformed());
    if ( coarseGridElement != exahype::solvers::Solver::NotFound ) {
      CellDescription& coarseGridCellDescription = getCellDescription(
          coarseGridCell.getCellDescriptionsIndex(),coarseGridElement);
      updateCoarseGridAncestorRefinementStatus(fineGridCellDescription,coarseGridCellDescription);
    }

    progressCollectiveRefinementOperationsInEnterCell(fineGridCellDescription);

    decideOnRefinement(fineGridCellDescription,inRefiningMode);
    decideOnVirtualRefinement(fineGridCellDescription);

    ensureFineGridCoarseGridConsistency(fineGridCellDescription,coarseGridCell.getCellDescriptionsIndex()); // must come after refinement status update
  }

  // Coarse grid cell based adaptive mesh refinement operations.
  // Add new cells to the grid and veto erasing or erasing virtual children
  // requests if there are cells on the fine level.
  if (coarseGridElement!=exahype::solvers::Solver::NotFound) {
    CellDescription& coarseGridCellDescription = getCellDescription(
        coarseGridCell.getCellDescriptionsIndex(),coarseGridElement);

    addNewVirtualCellIfVirtualRefiningRequested(
            fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
            coarseGridCellDescription,coarseGridCell.getCellDescriptionsIndex());
    newComputeCell |=
        addNewLeafIfRefinementRequested(
            fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
            coarseGridCellDescription,coarseGridCell.getCellDescriptionsIndex());
  }

  return newComputeCell;
}

int exahype::solvers::ADERDGSolver::evaluateRefinementCriterion(
    const CellDescription& cellDescription, const double* const solution, const double& timeStamp) {
  assertion1(isOfType(cellDescription,CellDescription::Type::Leaf),cellDescription.toString());

  Solver::RefinementControl refinementControl =
      refinementCriterion(
          solution,cellDescription.getOffset()+0.5*cellDescription.getSize(),
          cellDescription.getSize(),
          timeStamp, // is done after the update
          cellDescription.getLevel());

  switch ( refinementControl ) {
  case exahype::solvers::Solver::RefinementControl::Refine:
      return _refineOrKeepOnFineGrid;
  case exahype::solvers::Solver::RefinementControl::Keep:
    return ( cellDescription.getLevel()==getMaximumAdaptiveMeshLevel() ) ? _refineOrKeepOnFineGrid : Keep;
  case exahype::solvers::Solver::RefinementControl::Erase:
    return Erase;
  default:
    logError("adjustSolutionDuringMeshRefinementBody(...)",
        "unknown refinement control value=" << static_cast<int>(refinementControl) <<
        ". Please check the return values of your refinement criterion.");
    std::abort();
    return Pending-1;
  }
}

void exahype::solvers::ADERDGSolver::markForRefinement(CellDescription& cellDescription) {
  const double* const solution = static_cast<double*>(cellDescription.getSolution());
  const int refinementStatus = evaluateRefinementCriterion(
      cellDescription,solution,cellDescription.getTimeStamp());
  if ( refinementStatus==_refineOrKeepOnFineGrid ) {
    cellDescription.setRefinementFlag(true);
  }
  cellDescription.setRefinementStatus( std::max(cellDescription.getRefinementStatus(),refinementStatus) );
}

void exahype::solvers::ADERDGSolver::decideOnRefinement(
    CellDescription& fineGridCellDescription,const bool inRefiningMode) {
  // top-down refining
  if (
     inRefiningMode &&
     fineGridCellDescription.getType()==CellDescription::Type::Leaf &&
     && fineGridCellDescription.getLevel()<getMaximumAdaptiveMeshLevel()
     && fineGridCellDescription.getRefinementStatus() > Keep
  ) {
    fineGridCellDescription.setType(CellDescription::Type::LeafRefineInitiated);
  }
  // keep
  else if (
      inRefiningMode &&
      fineGridCellDescription.getType()==CellDescription::Type::Leaf
      &&
      (fineGridCellDescription.getRefinementStatus() == Keep
      ||
      (fineGridCellDescription.getLevel()==getMaximumAdaptiveMeshLevel()
      && fineGridCellDescription.getRefinementStatus() == _refineOrKeepOnFineGrid))
   ) {
    fineGridCellDescription.setType(CellDescription::Type::LeafKeep);
  }
  // coarsen
  else if (
      inRefiningMode &&
      fineGridCellDescription.getType()==CellDescription::Type::Leaf &&
      fineGridCellDescription.getRefinementStatus() == Keep ||
  ) {
    fineGridCellDescription.setType(CellDescription::Type::LeafKeep);
  }


  // bottom-up refining (halo refinement)
  if (
     inRefiningMode &&
     fineGridCellDescription.getType()==CellDescription::Type::Virtual &&
     fineGridCellDescription.getLevel()==getMaximumAdaptiveMeshLevel() &&
     (fineGridCellDescription.getRefinementStatus() > Keep ||
     fineGridCellDescription.getPreviousRefinementStatus() > Keep)
  ) {
    Solver::SubcellPosition subcellPosition = amr::computeSubcellPositionOfDescendant<CellDescription,Heap>(fineGridCellDescription);
    CellDescription& topMostParent =
      getCellDescription(subcellPosition.parentCellDescriptionsIndex,subcellPosition.parentElement);
    tarch::multicore::Lock lock(ADERDGSolver::CoarseGridSemaphore);
    if ( topMostParent.getType()==CellDescription::Type::Leaf ) {
      topMostParent.setRefinementStatus(_refineOrKeepOnFineGrid);
    }
    lock.free();
  }
}

void exahype::solvers::ADERDGSolver::addNewCell(
    exahype::Cell& fineGridCell,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    const int coarseGridCellDescriptionsIndex,
    const int solverNumber) {
  logDebug("addNewCell(...)","Add new grid cell with center "<<fineGridVerticesEnumerator.getCellCenter() <<
              " at level "<<fineGridVerticesEnumerator.getLevel());

  CellInfo cellInfo = fineGridCell.addNewCellDescription(
      solverNumber,
      CellDescription::Type::Leaf,
      fineGridVerticesEnumerator.getLevel(),
      coarseGridCellDescriptionsIndex,
      fineGridVerticesEnumerator.getCellSize(),
      fineGridVerticesEnumerator.getVertexPosition());

  const int fineGridElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  CellDescription& fineGridCellDescription = cellInfo._ADERDGCellDescriptions[fineGridElement]; //TODO(Dominic): Multi-solvers: Might need to lock this?
  if ( !exahype::solvers::Solver::SpawnAMRBackgroundJobs ) {
    ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
  }

  fineGridCellDescription.setPreviousRefinementStatus(Erase); // reasonable state after rollback
  fineGridCellDescription.setRefinementStatus(Pending);

  #ifdef Asserts
  fineGridCellDescription.setCreation(CellDescription::Creation::UniformRefinement);
  #endif
}

void exahype::solvers::ADERDGSolver::addNewVirtualCellIfVirtualRefiningRequested(
     exahype::Cell& fineGridCell,
     exahype::Vertex* const fineGridVertices,
     const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
     CellDescription& coarseGridCellDescription,
     const int coarseGridCellDescriptionsIndex) {
  const int fineGridElement = tryGetElement(
      fineGridCell.getCellDescriptionsIndex(),coarseGridCellDescription.getSolverNumber());
  // read from coarse grid
  tarch::multicore::Lock lock(CoarseGridSemaphore);
  const bool refineVirtualCell =
      (isLeaf(coarseGridCellDescription) || isVirtual(coarseGridCellDescription)) &&
      coarseGridCellDescription.getAugmentationStatus() >= MinimumAugmentationStatusForVirtualRefining;
  lock.free();

  // work on fine grid
  if ( refineVirtualCell &&
      fineGridElement==exahype::solvers::Solver::NotFound
  ) {
    fineGridCell.addNewCellDescription( // (EmptyDescendant),None
        coarseGridCellDescription.getSolverNumber(),
        CellDescription::Type::Virtual,
        fineGridVerticesEnumerator.getLevel(),
        coarseGridCellDescriptionsIndex,
        fineGridVerticesEnumerator.getCellSize(),
        fineGridVerticesEnumerator.getVertexPosition());

    #ifdef Asserts
    const int fineGridElement = tryGetElement(
        fineGridCell.getCellDescriptionsIndex(),coarseGridCellDescription.getSolverNumber());
    #endif
    assertion(fineGridElement!=exahype::solvers::Solver::NotFound);
  }
}

// TODO(Dominic): Cannot be called multiple times. Need extra state?
bool exahype::solvers::ADERDGSolver::addNewLeafIfRefinementRequested(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    CellDescription& coarseGridCellDescription,
    const int coarseGridCellDescriptionsIndex) {
  // read and modify coarse grid
  tarch::multicore::Lock lock(CoarseGridSemaphore);
  bool refiningOrRefiningRequested =
      coarseGridCellDescription.getType()==CellDescription::Type::LeafRefiningInitated ||
      coarseGridCellDescription.getType()==CellDescription::Type::LeafRefining;
  if ( refiningOrRefiningRequested ) {
    assertion1(coarseGridCellDescription.getType()==CellDescription::Type::Leaf,
               coarseGridCellDescription.toString());
    coarseGridCellDescription.setRefinementEvent(CellDescription::Type::LeafRefining);
  }
  lock.free();

  // work on fine grid
  if ( refiningOrRefiningRequested ) {
    const int fineGridElement = tryGetElement(
        fineGridCell.getCellDescriptionsIndex(),coarseGridCellDescription.getSolverNumber());

    if ( fineGridElement==exahype::solvers::Solver::NotFound ) {
      addNewCell(fineGridCell,fineGridVerticesEnumerator,
                 coarseGridCellDescriptionsIndex,
                 coarseGridCellDescription.getSolverNumber());
      CellDescription& fineGridCellDescription =
          getCellDescriptions(fineGridCell.getCellDescriptionsIndex()).back();
      fineGridCellDescription.setRefinementEvent(CellDescription::LeafProlongating);
      #ifdef Asserts
      fineGridCellDescription.setCreation(CellDescription::Creation::AdaptiveRefinement);
      #endif
    } else {
      CellDescription& fineGridCellDescription = getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridElement);
      #ifdef Parallel
      assertion4(fineGridCellDescription.getType()==CellDescription::Type::Virtual ||
                 fineGridCellDescription.getType()==CellDescription::Type::Leaf,
                 fineGridCellDescription.toString(),coarseGridCellDescription.toString(),
                 coarseGridCellDescriptionsIndex,
                 tarch::parallel::Node::getInstance().getRank());
      #else
      assertion2(fineGridCellDescription.getType()==CellDescription::Type::Virtual,
                 fineGridCellDescription.toString(),coarseGridCellDescription.toString());
      #endif
      assertion2(fineGridCellDescription.getParentIndex()==coarseGridCellDescriptionsIndex,
                 fineGridCellDescription.toString(),coarseGridCellDescriptionsIndex);

      fineGridCellDescription.setType(CellDescription::Type::LeafProlongatings);
      fineGridCellDescription.setCommunicationStatus(CellCommunicationStatus);
      fineGridCellDescription.setFacewiseCommunicationStatus(0); // implicit conversion
      ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
      #ifdef Asserts
      fineGridCellDescription.setCreation(CellDescription::Creation::AdaptiveRefinement);
      #endif
      fineGridCellDescription.setPreviousRefinementStatus(Erase); // reasonable state after rollback
      fineGridCellDescription.setRefinementStatus(Pending);
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
  double* solutionFine   = static_cast<double*>(fineGridCellDescription.getSolution());
  double* solutionCoarse = static_cast<double*>(coarseGridCellDescription.getSolution());
  volumeUnknownsProlongation(
      solutionFine,solutionCoarse,
      levelCoarse,levelFine,
      subcellIndex);

  // previous solution
  assertion(DataHeap::getInstance().isValidIndex(fineGridCellDescription.getPreviousSolutionIndex()));
  double* previousSolutionFine   = static_cast<double*>(fineGridCellDescription.getPreviousSolution());
  double* previousSolutionCoarse = static_cast<double*>(coarseGridCellDescription.getPreviousSolution());
  volumeUnknownsProlongation(
      previousSolutionFine,previousSolutionCoarse,
      levelCoarse,levelFine,
      subcellIndex);

  fineGridCellDescription.setTimeStamp(coarseGridCellDescription.getTimeStamp());
  fineGridCellDescription.setTimeStepSize(coarseGridCellDescription.getTimeStepSize());

  // TODO Dominic: This is a little inconsistent since I orignially tried to hide
  // the limiting from the pure ADER-DG scheme
  fineGridCellDescription.setPreviousRefinementStatus(Erase); // TODO(Dominic): NEW CELLS
  fineGridCellDescription.setRefinementStatus(Pending);

  // TODO Dominic:
  // During the inital mesh build where we only refine
  // according to the PAD, we don't want to have a too broad refined area.
  // We thus do not flag children with troubled
  if (
      !initialGrid &&
      coarseGridCellDescription.getRefinementStatus()>=_minRefinementStatusForTroubledCell
  ) {
    fineGridCellDescription.setRefinementStatus(_minRefinementStatusForTroubledCell);
    fineGridCellDescription.setIterationsToCureTroubledCell(coarseGridCellDescription.getIterationsToCureTroubledCell());
  }
  fineGridCellDescription.setFacewiseRefinementStatus(Pending);
}

bool exahype::solvers::ADERDGSolver::attainedStableState(
        exahype::Cell&                       fineGridCell,
        exahype::Vertex* const               fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
        const int                            solverNumber,
        const bool                           inRefiningMode) const {
  const int element = tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if ( element!=exahype::solvers::Solver::NotFound ) {
    CellDescription& cellDescription = getCellDescription(fineGridCell.getCellDescriptionsIndex(),element);
    // compute flagging gradients in inside cells
    bool flaggingHasConverged = true;
    if ( !peano::grid::aspects::VertexStateAnalysis::isOneVertexBoundary(fineGridVertices,fineGridVerticesEnumerator) ) { // no check on boundary
      if ( cellDescription.getType()==CellDescription::Type::Leaf || cellDescription.getType()==CellDescription::Type::Parent ) {
        for (int d=0; d<DIMENSIONS; d++) {
          flaggingHasConverged &=
              std::abs(cellDescription.getFacewiseAugmentationStatus(2*d+1)  - cellDescription.getFacewiseAugmentationStatus(2*d+0)) <= 2;
          flaggingHasConverged &=
              std::abs(cellDescription.getFacewiseCommunicationStatus(2*d+1) - cellDescription.getFacewiseCommunicationStatus(2*d+0)) <= 2;
        }
      }
      // refinement status is only spread on finest level
      if (
          cellDescription.getType()  == CellDescription::Type::Leaf &&
          cellDescription.getLevel() == getMaximumAdaptiveMeshLevel()
      ) {
        for (int d=0; d<DIMENSIONS; d++) {
          flaggingHasConverged &=
              std::abs(cellDescription.getFacewiseRefinementStatus(2*d+1) - cellDescription.getFacewiseRefinementStatus(2*d+0)) <= 2;
        }
      }
    }

    bool stable =
      flaggingHasConverged
      &&
      (cellDescription.getType()==CellDescription::Type::LeafKeep   ||
      cellDescription.getType()==CellDescription::Type::ParentKeep ||
      cellDescription.getType()==CellDescription::Type::Virtual)
      &&
      (!isLeaf(cellDescription) || // cell must not have pending refinement status and must not require refinement on coarser grids
      (cellDescription.getRefinementStatus()!=Pending &&
      (cellDescription.getLevel() == getMaximumAdaptiveMeshLevel() ||
      cellDescription.getRefinementStatus()<=0)))
      &&
      (isOfType(cellDescription,CellDescription::Type::Virtual) || // descendant must not have refinement status > 0 on finest level
      cellDescription.getLevel() != getMaximumAdaptiveMeshLevel() ||
      cellDescription.getRefinementStatus()<=0);

      if ( !stable ) {
        #ifdef MonitorMeshRefinement
        logInfo("attainedStableState(...)","cell has not attained stable state (yet):");
        logInfo("attainedStableState(...)","type="<<cellDescription.toString(cellDescription.getType()));
        logInfo("attainedStableState(...)","x="<<cellDescription.getOffset());
        logInfo("attainedStableState(...)","level="<<cellDescription.getLevel());
        logInfo("attainedStableState(...)","flaggingHasConverged="<<flaggingHasConverged);
        logInfo("attainedStableState(...)","refinementEvent="<<cellDescription.toString(cellDescription.getRefinementEvent()));
        logInfo("attainedStableState(...)","refinementStatus="<<cellDescription.getRefinementStatus());
        logInfo("attainedStableState(...)","getFacewiseAugmentationStatus="<<cellDescription.getFacewiseAugmentationStatus());
        logInfo("attainedStableState(...)","getFacewiseCommunicationStatus="<<cellDescription.getFacewiseCommunicationStatus());
        logInfo("attainedStableState(...)","getFacewiseRefinementStatus="<<cellDescription.getFacewiseRefinementStatus());
        #ifdef Asserts
        logInfo("attainedStableState(...)","cellDescription.getCreation()="<<cellDescription.toString(cellDescription.getCreation()));
        #endif
        logInfo("attainedStableState(...)","fineGridCell.getCellDescriptionsIndex()="<<fineGridCell.getCellDescriptionsIndex());
        logInfo("attainedStableState(...)","solver.getCoarsestMeshLevel="<<getCoarsestMeshLevel());
        logInfo("attainedStableState(...)","solver.getMaximumAdaptiveMeshLevel="<<getMaximumAdaptiveMeshLevel());
        logInfo("attainedStableState(...)","solver.getMaximumAdaptiveMeshDepth="<<getMaximumAdaptiveMeshDepth());
        #endif
        if (
            !inRefiningMode &&
            isOfType(cellDescription,CellDescription::Type::Leaf) &&
            cellDescription.getLevel() == getMaximumAdaptiveMeshLevel() &&
            cellDescription.getRefinementStatus() > 0
        ) {
          logError("attainedStableState(...)","Virtual subcell requests refining of coarse grid parent but mesh refinement is not in refining mode anymore. Inform Dominic.");
          logError("attainedStableState(...)","cell="<<cellDescription.toString());
          std::terminate();
        }
      }

    return stable;
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
    const int solverNumber,
    const bool inRefiningMode) {
  bool newComputeCell = false;

  const int fineGridElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if ( fineGridElement!=exahype::solvers::Solver::NotFound ) {
    CellDescription& fineGridCellDescription = getCellDescription(
        fineGridCell.getCellDescriptionsIndex(),fineGridElement);

    // start or finish collective operations
    progressCollectiveRefinementOperationsInLeaveCell(fineGridCellDescription,inRefiningMode);

    // skip remainder if the refinement criterion has not been evaluated yet for a Cell
    // Reading the refinement request might result into data race but this is accepted at this point
    // as we only read and not write
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

      restrictVolumeDataIfErasingRequested(
          fineGridCellDescription,coarseGridCellDescription);

      eraseCellDescriptionIfNecessary(
              fineGridCell.getCellDescriptionsIndex(),
              fineGridElement,
              coarseGridCellDescription);
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
    const int level,
    const bool checkThoroughly) const {
  if ( level < _coarsestMeshLevel ) {
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

        if (
            !checkThoroughly ||
            Vertex::equalUpToRelativeTolerance(cellDescription.getOffset(), cellOffset)
        ) {
          bool refineAdjacentVertices =
              cellDescription.getType()==CellDescription::Type::Parent ||
              cellDescription.getType()==CellDescription::Type::LeafRefineInitiated ||
              cellDescription.getType()==CellDescription::Type::LeafRefine ||
              cellDescription.getHasVirtualChildren() ||
              cellDescription.getAugmentationStatus() > 0 ||
              cellDescription.getRefinementStatus() > 0;

          bool eraseAdjacentVertices =
              (cellDescription.getType()==CellDescription::Type::LeafKeep ||
              cellDescription.getType()==CellDescription::Type::Virtual)
              &&
              !cellDescription.getHasVirtualChildren()
              &&
              cellDescription.getAugmentationStatus()==0 // TODO(Dominic): Probably can tune here. This is chosen to large
              &&
              cellDescription.getRefinementStatus()==0 &&
              cellDescription.getPreviousRefinementStatus()==0;

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
      static_cast<double*>(cellDescription.getSolution());
  std::fill_n(solution,getDataPerCell(),0.0);
  double* previousSolution =
      static_cast<double*>(cellDescription.getPreviousSolution());
  std::fill_n(previousSolution,getDataPerCell(),0.0);
}

void exahype::solvers::ADERDGSolver::changeLeafToParent(CellDescription& cellDescription) {
  assertion1(
      cellDescription.getType()==CellDescription::Type::LeafRefine ||
      cellDescription.getType()==CellDescription::Type::ParentCoarseningRequested
      ,cellDescription.toString());
  cellDescription.setType(CellDescription::Type::ParentKeep); // TODO(Dominic): Keep?
  cellDescription.setAugmentationStatus(MaximumAugmentationStatus);
  cellDescription.setHasVirtualChildren(false); // since we might replace descendants with cells
  cellDescription.setRefinementStatus(Keep);
  cellDescription.setCommunicationStatus(0);
  cellDescription.setFacewiseAugmentationStatus(0); // implicit conversion
  cellDescription.setFacewiseRefinementStatus(Pending);
  cellDescription.setFacewiseCommunicationStatus(0); // implicit conversion
  ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
}

void exahype::solvers::ADERDGSolver::progressCollectiveRefinementOperationsInEnterCell(
     CellDescription& fineGridCellDescription) {
  if (fineGridCellDescription.getType() == CellDescription::Type::LeafRefine ) {
    changeLeafToParent(fineGridCellDescription);
  }
}

bool exahype::solvers::ADERDGSolver::markPreviousAncestorForRefinement(CellDescription& cellDescription) {
  assertion1(cellDescription.getType()==CellDescription::Type::ParentCoarsening,cellDescription.toString());
    double* solution = static_cast<double*>(cellDescription.getSolution());
    adjustSolution(solution,
          cellDescription.getOffset()+0.5*cellDescription.getSize(),
          cellDescription.getSize(),
          cellDescription.getTimeStamp(),
          cellDescription.getTimeStepSize());

    double* previousSolution = static_cast<double*>(cellDescription.getPreviousSolution());
    adjustSolution(previousSolution,
          cellDescription.getOffset()+0.5*cellDescription.getSize(),
          cellDescription.getSize(),
          cellDescription.getPreviousTimeStamp(),
          cellDescription.getPreviousTimeStepSize());

    // TODO(Dominic): Probably not a good idea to overwrite these values
    cellDescription.setRefinementStatus( evaluateRefinementCriterion(
            cellDescription,solution,cellDescription.getTimeStamp())
    );
    cellDescription.setPreviousRefinementStatus( evaluateRefinementCriterion(
            cellDescription,previousSolution,cellDescription.getPreviousTimeStamp())
    );

    return cellDescription.getRefinementStatus()        !=_refineOrKeepOnFineGrid &&
           cellDescription.getPreviousRefinementStatus()!=_refineOrKeepOnFineGrid;
}

void exahype::solvers::ADERDGSolver::progressCollectiveRefinementOperationsInLeaveCell(
     CellDescription& fineGridCellDescription,
     const bool inRefiningMode) {
  switch ( fineGridCellDescription.getType() ) {
    case CellDescription::Type::ParentCoarseningRequested:
      // evaluate refinement criterion now that fine grid cells have restricted their data
      if ( markPreviousAncestorForRefinement(fineGridCellDescription) ) {
        fineGridCellDescription.setType(CellDescription::Type::ParentCoarsening);
      } else { // veto erasing request
        changeLeafToParent(fineGridCellDescription);
      }
      break;
    case CellDescription::Type::ParentCoarsening:
      //logInfo("progressCollectiveRefinementOperationsInLeaveCell(...)","ErasingChildren done");
      fineGridCellDescription.setType(CellDescription::Type::Leaf);
      break;
    default:
      break;
  }

  if (  // The children of this cell description have all flagged themselves and their parent (this cell) with Erase
      !inRefiningMode &&
      fineGridCellDescription.getType()==CellDescription::Type::ParentCoarseningRequested &&
      fineGridCellDescription.getVetoErasingChildren()==false &&
      fineGridCellDescription.getRefinementStatus()==Erase &&
      fineGridCellDescription.getPreviousRefinementStatus()==Erase
  ) {
    //logInfo("progressCollectiveRefinementOperationsInLeaveCell(...)","ErasingChildren requested: "<<fineGridCellDescription.getRefinementStatus()<< ", "<<fineGridCellDescription.getPreviousRefinementStatus());
    fineGridCellDescription.setType(CellDescription::Type::Leaf);
    fineGridCellDescription.setAugmentationStatus(0);
    fineGridCellDescription.setFacewiseAugmentationStatus(0); // implicit conversion
    fineGridCellDescription.setCommunicationStatus(CellCommunicationStatus);
    fineGridCellDescription.setFacewiseCommunicationStatus(CellCommunicationStatus); // implicit conversion
    ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
    prepareVolumeDataRestriction(fineGridCellDescription);
    fineGridCellDescription.setType(CellDescription::Type::ErasingChildrenRequested);
  }
}

void exahype::solvers::ADERDGSolver::eraseCellDescriptionIfNecessary(
    const int cellDescriptionsIndex,
    const int fineGridElement,
    CellDescription& coarseGridCellDescription) {
  tarch::multicore::Lock lock(CoarseGridSemaphore);
  const bool erasingChildren =
      coarseGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ErasingChildren;
  const bool changeChildrenToDescendants =
      coarseGridCellDescription.getRefinementEvent()==CellDescription::ChangeChildrenToVirtualChildren;
  const bool deaugmentingChildren =
      coarseGridCellDescription.getRefinementEvent()==CellDescription::RefinementEvent::ErasingVirtualChildren;
  lock.free();

  if ( changeChildrenToDescendants ) {
    CellDescription& fineGridCellDescription = getCellDescription(
        cellDescriptionsIndex,fineGridElement);

    // erase cell description // or change to descendant
    fineGridCellDescription.setType(CellDescription::Type::Virtual);
    fineGridCellDescription.setCommunicationStatus(0);
    fineGridCellDescription.setFacewiseCommunicationStatus(0); // implicit conversion
    ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);
  }
  else if ( erasingChildren || deaugmentingChildren ) {
    CellDescription& fineGridCellDescription = getCellDescription(
          cellDescriptionsIndex,fineGridElement);

    fineGridCellDescription.setType(CellDescription::Erased);
    fineGridCellDescription.setCommunicationStatus(0);
    fineGridCellDescription.setFacewiseCommunicationStatus(0); // implicit conversion
    ensureNoUnnecessaryMemoryIsAllocated(fineGridCellDescription);

    getCellDescriptions(cellDescriptionsIndex).erase(
        getCellDescriptions(cellDescriptionsIndex).begin()+fineGridElement);
  }
}

void exahype::solvers::ADERDGSolver::restrictVolumeDataIfErasingRequested(
    const CellDescription& fineGridCellDescription,
    const CellDescription& coarseGridCellDescription) {
//  assertion1(coarseGridCellDescription.getLimiterStatus()==CellDescription::LimiterStatus::Ok,
//      coarseGridCellDescription.toString()); // TODO(Dominic): Does not always apply see veto
  tarch::multicore::Lock lock(CoarseGridSemaphore);
  const bool restrictVolumeData =
      coarseGridCellDescription.getType()==CellDescription::Type::ParentCoarseningRequested;
  lock.free();

  if ( restrictVolumeData ) {
    //logInfo("restrictVolumeData(..)","restricting solution");

    tarch::la::Vector<DIMENSIONS,int> subcellIndex =
        exahype::amr::computeSubcellIndex(
            fineGridCellDescription.getOffset(),
            fineGridCellDescription.getSize(),coarseGridCellDescription.getOffset());

    // restrict values.
    tarch::multicore::Lock lock(RestrictionSemaphore);
    assertion1(fineGridCellDescription.getRefinementStatus()==-1,fineGridCellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(fineGridCellDescription.getSolutionIndex()),fineGridCellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(coarseGridCellDescription.getSolutionIndex()),coarseGridCellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(fineGridCellDescription.getPreviousSolutionIndex()),fineGridCellDescription.toString());
    assertion1(DataHeap::getInstance().isValidIndex(coarseGridCellDescription.getPreviousSolutionIndex()),coarseGridCellDescription.toString());

    const int levelFine   = fineGridCellDescription.getLevel();
    const int levelCoarse = coarseGridCellDescription.getLevel();
    assertion(levelCoarse < levelFine);

    if ( !DataHeap::getInstance().isValidIndex(fineGridCellDescription.getSolutionIndex()) ) {
      logError("restrictVolumeData(..)","solution not valid for cell="<<fineGridCellDescription.toString());
      std::abort();
    }

    // restrict current solution
    double* solutionFine   = static_cast<double*>(fineGridCellDescription.getSolution());
    double* solutionCoarse = static_cast<double*>(coarseGridCellDescription.getSolution());
    volumeUnknownsRestriction(
        solutionCoarse,solutionFine,
        levelCoarse,levelFine,
        subcellIndex);

    // restrict next solution
    double* previousSolutionFine   = static_cast<double*>(fineGridCellDescription.getPreviousSolution());
    double* previousSolutionCoarse = static_cast<double*>(coarseGridCellDescription.getPreviousSolution());
    volumeUnknownsRestriction(
        previousSolutionCoarse,previousSolutionFine,
        levelCoarse,levelFine,
        subcellIndex);

    // TODO(Dominic): Do later, move out

    // Reset the min and max
    const int numberOfObservables = getDMPObservables();
    if ( numberOfObservables>0 ) {
      double* solutionMin = static_cast<double*>(coarseGridCellDescription.getSolutionMin());
      std::fill_n(solutionMin,DIMENSIONS_TIMES_TWO*numberOfObservables,
          std::numeric_limits<double>::infinity());
      double* solutionMax = static_cast<double*>(coarseGridCellDescription.getSolutionMax());
      std::fill_n(solutionMax,DIMENSIONS_TIMES_TWO*numberOfObservables,
          -std::numeric_limits<double>::infinity()); // Be aware of "-"
    }

    // TODO(Dominic): What to do with the time step data for anarchic time stepping?
    // Tobias proposed some waiting procedure. Until they all have reached
    // a certain time level.
    //  coarseGridCellDescription.setTimeStamp(fineGridCellDescription.getTimeStamp());
    //  coarseGridCellDescription.setPredictorTimeStamp(fineGridCellDescription.getPredictorTimeStamp());
    //  coarseGridCellDescription.setTimeStepSize(fineGridCellDescription.getTimeStepSize());
    //  coarseGridCellDescription.setPredictorTimeStepSize(fineGridCellDescription.getPredictorTimeStepSize());
    // TODO(Dominic): Reconsider for anarchic time stepping.
    // coarseGridCellDescription.setTimeStamp(fineGridCellDescription.getTimeStamp());
    // coarseGridCellDescription.setPredictorTimeStamp(fineGridCellDescription.getPredictorTimeStamp());
    lock.free();
  }
}

void exahype::solvers::ADERDGSolver::ensureFineGridCoarseGridConsistency(
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

    #if defined(Asserts) || defined(Debug)
    assertion2(coarseGridElement==exahype::solvers::Solver::NotFound || fineGridCellDescription.getParentIndex()==coarseGridCellDescriptionsIndex,
        fineGridCellDescription.toString(), getCellDescription(coarseGridCellDescriptionsIndex,coarseGridElement).toString());
    #endif

  } else {
    fineGridCellDescription.setParentIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
  }
}

void exahype::solvers::ADERDGSolver::finaliseStateUpdates(
      const int solverNumber,
      CellInfo& cellInfo) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element!=exahype::solvers::Solver::NotFound ) {
    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];
    if (cellDescription.getType()==CellDescription::Type::Leaf) {
      validateCellDescriptionData(cellDescription,cellDescription.getTimeStamp()>0,false,true,"finaliseStateUpdates");
    }
    cellDescription.setRefinementFlag(false);

    if ( getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested ) {
      cellDescription.setPreviousRefinementStatus(cellDescription.getRefinementStatus());
    }
  }
}

//////////////////////////
/// Forks and Joins
//////////////////////////


#ifdef Parallel

const int exahype::solvers::ADERDGSolver::DataMessagesPerForkOrJoinCommunication   = 2;

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
      // wait for background jobs to complete
      if ( !cellDescription.getHasCompletedLastStep() ) {
        peano::datatraversal::TaskSet::startToProcessBackgroundJobs();
      }
      int numberOfBackgroundJobsToProcess = 1;
      while ( !cellDescription.getHasCompletedLastStep() ) {
        tarch::multicore::jobs::processBackgroundJobs(numberOfBackgroundJobsToProcess);
        numberOfBackgroundJobsToProcess++;
      }
      oneSolverRequiresVerticalCommunication &=
          cellDescription.getType()==CellDescription::Type::Virtual && cellDescription.getHasVirtualChildren();
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

  // Default field data indices
  cellDescription.setSolutionIndex(-1);
  cellDescription.setSolution(nullptr);
  cellDescription.setPreviousSolutionIndex(-1);
  cellDescription.setPreviousSolution(nullptr);
  cellDescription.setUpdateIndex(-1);
  cellDescription.setUpdate(nullptr);
  cellDescription.setExtrapolatedPredictorIndex(-1);
  cellDescription.setExtrapolatedPredictor(nullptr);
  cellDescription.setFluctuationIndex(-1);
  cellDescription.setFluctuation(nullptr);

  // Limiter meta data (oscillations identificator)
  cellDescription.setSolutionMinIndex(-1);
  cellDescription.setSolutionMin(nullptr);
  cellDescription.setSolutionMaxIndex(-1);
  cellDescription.setSolutionMax(nullptr);

  // compression
  cellDescription.setCompressionState(CellDescription::CompressionState::Uncompressed);

  cellDescription.setExtrapolatedPredictorCompressedIndex(-1);
  cellDescription.setExtrapolatedPredictorCompressed(nullptr);
  cellDescription.setFluctuationCompressedIndex(-1);
  cellDescription.setFluctuationCompressed(nullptr);
  cellDescription.setSolutionCompressedIndex(-1);
  cellDescription.setSolutionCompressed(nullptr);
  cellDescription.setPreviousSolutionCompressedIndex(-1);
  cellDescription.setPreviousSolutionCompressed(nullptr);
  cellDescription.setUpdateCompressedIndex(-1);
  cellDescription.setUpdateCompressed(nullptr);

  cellDescription.setSolutionAveragesIndex(-1);
  cellDescription.setSolutionAverages(nullptr);
  cellDescription.setSolutionAveragesIndex(-1);
  cellDescription.setSolutionAverages(nullptr);
  cellDescription.setUpdateAveragesIndex(-1);
  cellDescription.setUpdateAverages(nullptr);
  cellDescription.setExtrapolatedPredictorAveragesIndex(-1);
  cellDescription.setExtrapolatedPredictorAverages(nullptr);
  cellDescription.setFluctuationAveragesIndex(-1);
  cellDescription.setFluctuationAverages(nullptr);
  cellDescription.setBytesPerDoFInPreviousSolution(-1);
  cellDescription.setBytesPerDoFInSolution(-1);
  cellDescription.setBytesPerDoFInUpdate(-1);
  cellDescription.setBytesPerDoFInExtrapolatedPredictor(-1);
  cellDescription.setBytesPerDoFInFluctuation(-1);

  // reset the facewise flags
  cellDescription.setFacewiseAugmentationStatus(0);
  cellDescription.setFacewiseCommunicationStatus(0);

  cellDescription.setFacewiseRefinementStatus(Pending);

  // limiter flagging
  cellDescription.setIterationsToCureTroubledCell(-1);

  #ifdef Asserts
  cellDescription.setCreation(CellDescription::Creation::ReceivedDueToForkOrJoin);
  #endif

  // background jobs
  cellDescription.setHasCompletedLastStep(true);
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
    metadata.push_back(cellDescription.getRefinementStatus());
    metadata.push_back( (cellDescription.getHasVirtualChildren()) ? 1 : 0 );
  } else {
    for (int i = 0; i < exahype::MasterWorkerCommunicationMetadataPerSolver; ++i) {
      metadata.push_back(exahype::InvalidMetadataEntry); // implicit conversion
    }
  }
}

void exahype::solvers::ADERDGSolver::deduceChildCellErasingEvents(CellDescription& cellDescription) const {
  const int coarseGridElement = tryGetElement(cellDescription.getParentIndex(),cellDescription.getSolverNumber());
  if ( coarseGridElement!=exahype::solvers::Solver::NotFound ) {
    CellDescription& coarseGridCellDescription =
        getCellDescription(cellDescription.getParentIndex(),coarseGridElement);

    tarch::multicore::Lock lock(CoarseGridSemaphore);
    switch (coarseGridCellDescription.getType()) {
    case CellDescription::Type::ParentCoarseningRequested: {
      assertion1(coarseGridCellDescription.getType()==CellDescription::Type::Parent,
          coarseGridCellDescription.toString());
      cellDescription.setType(CellDescription::Type::LeafParentRequestedCoarsening); // TODO(Dominic): Introduce these states
    } break;
    case CellDescription::Type::ParentCoarsening: {
      assertion1(coarseGridCellDescription.getType()==CellDescription::Type::Parent,
          coarseGridCellDescription.toString());
      cellDescription.setType(CellDescription::Type::LeafParentRequestedCoarsening); // TODO(Dominic): Introduce these states
    } break;
    //
    lock.free();
    }
  }
}

void exahype::solvers::ADERDGSolver::progressMeshRefinementInPrepareSendToWorker(
    const int workerRank,
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const int solverNumber) {
  // coarse grid based operations
  const int coarseGridCellDescriptionsIndex = coarseGridCell.getCellDescriptionsIndex();
  const int coarseGridCellElement = tryGetElement(coarseGridCellDescriptionsIndex,solverNumber);
  if (coarseGridCellElement!=exahype::solvers::Solver::NotFound) {
    CellDescription& coarseGridCellDescription = getCellDescription(
        coarseGridCell.getCellDescriptionsIndex(),coarseGridCellElement);

    addNewVirtualCellIfVirtualRefiningRequested(
        fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
        coarseGridCellDescription,coarseGridCell.getCellDescriptionsIndex());

    bool addedNewCell =
        addNewLeafIfRefinementRequested(
            fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
            coarseGridCellDescription,coarseGridCell.getCellDescriptionsIndex());

    const int cellDescriptionsIndex = fineGridCell.getCellDescriptionsIndex();
    const int element = tryGetElement(cellDescriptionsIndex,solverNumber);

    if ( element!=exahype::solvers::Solver::NotFound ) {
      CellDescription& fineGridCellDescription = getCellDescription(cellDescriptionsIndex,element);
      if (
        fineGridCellDescription.getType()==CellDescription::Type::Virtual &&
        fineGridCellDescription.getHasVirtualChildren()
      ) {
        Solver::SubcellPosition subcellPosition =  amr::computeSubcellPositionOfVirtualCell<CellDescription,Heap>(fineGridCellDescription);
        CellDescription& topMostParentCellDescription =
            getCellDescription(subcellPosition.parentCellDescriptionsIndex,subcellPosition.parentElement);
        if ( topMostParentCellDescription.getType()==CellDescription::Type::Leaf ) {
           logDebug( "progressMeshRefinementInPrepareSendToWorker(...)"," try to refine parent " << topMostParentCellDescription.toString());
           topMostParentCellDescription.setRefinementStatus(_refineOrKeepOnFineGrid);
        }
      }
    }

    if ( addedNewCell ) {
      CellDescription& fineGridCellDescription = getCellDescription(cellDescriptionsIndex,element);
      prolongateVolumeData(fineGridCellDescription,getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested);
      assertion1( fineGridCellDescription.getType()==CellDescription::Type::LeafProlongating,fineGridCellDescription.toString());
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
  if ( fineGridCellDescription.getType()==CellDescription::Type::LeafProlongating ) {
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

  if ( receivedCellDescription.getType()==CellDescription::Type::LeafProlongating ) {
    logDebug( "receiveDataFromMasterIfProlongating(...)","receiving prolongated solution from rank "<<masterRank<< " at x="<<x.toString()<< ",level="<<level);

    mergeWithWorkerOrMasterDataDueToForkOrJoin(
      masterRank,receivedCellDescriptionsIndex,receivedElement,
      peano::heap::MessageType::MasterWorkerCommunication,x,level);
  }
}

bool exahype::solvers::ADERDGSolver::progressMeshRefinementInMergeWithWorker(
    const int localCellDescriptionsIndex,
    const int receivedCellDescriptionsIndex, const int receivedElement) {
  auto& receivedCellDescriptions = getCellDescriptions(receivedCellDescriptionsIndex);
  assertion1( isValidCellDescriptionIndex(localCellDescriptionsIndex), localCellDescriptionsIndex );
  auto& localCellDescriptions = getCellDescriptions(localCellDescriptionsIndex);

  int localElement = NotFound;
  for (unsigned int element = 0; element < localCellDescriptions.size(); ++element) {
    if ( localCellDescriptions[element].getSolverNumber()==receivedCellDescriptions[receivedElement].getSolverNumber() ) {
      localElement = element;
    }
    element++;
  }
  if ( localElement==NotFound ) { // We have already received data and allocated all memory in this case
    CellDescription& receivedCellDescription = receivedCellDescriptions[receivedElement];
    localCellDescriptions.push_back(receivedCellDescription);
    assertion1(receivedCellDescription.getType()==CellDescription::Type::Virtual,receivedCellDescription.toString());
    receivedCellDescriptions.erase(receivedCellDescriptions.begin()+receivedElement);
    return false;
  } else {
    CellDescription& receivedCellDescription = receivedCellDescriptions[receivedElement];
    if ( receivedCellDescriptions[receivedElement].getType()==CellDescription::Type::LeafProlongating ) {
      CellDescription& localCellDescription = localCellDescriptions[localElement];

      logDebug( "progressMeshRefinementInMergeWithWorker(...)","merging prolongated solution");

      assertion( localCellDescription.getType()==CellDescription::Type::Leaf ||
          localCellDescription.getType()==CellDescription::Type::Virtual);
      assertion(receivedCellDescription.getType()==CellDescription::Type::Leaf);
      assertion(DataHeap::getInstance().isValidIndex(receivedCellDescription.getSolutionIndex()));
      assertion(DataHeap::getInstance().isValidIndex(receivedCellDescription.getPreviousSolutionIndex()));

      // we know we have received data in this case
      localCellDescription.setType(CellDescription::Type::LeafProlongating);
      localCellDescription.setPreviousRefinementStatus(Pending);
      localCellDescription.setCommunicationStatus(CellCommunicationStatus);
      localCellDescription.setFacewiseCommunicationStatus(0); // implicit conversion

      ensureNecessaryMemoryIsAllocated(localCellDescription); // copy indices
      localCellDescription.setSolution(receivedCellDescription.getSolution());
      localCellDescription.setPreviousSolution(receivedCellDescription.getPreviousSolution());
      receivedCellDescription.setSolutionIndex(-1);
      receivedCellDescription.setPreviousSolutionIndex(-1);
      receivedCellDescription.setSolution(nullptr);
      receivedCellDescription.setPreviousSolution(nullptr);

      // adjust solution
      // TODO(Dominic): do not understand this part
      localCellDescription.setType(CellDescription::Type::Leaf);
      return true;
    } else {
      return false;
    }
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
      cellDescription.getType()==CellDescription::Type::LeafParentCoarseningRequested
  ) {
    sendDataToWorkerOrMasterDueToForkOrJoin(masterRank,cellDescriptionsIndex,element,
        peano::heap::MessageType::MasterWorkerCommunication,x,level); // assumes blocking/copy
  }

  // erase or change type of cell descriptions
  if ( cellDescription.getType()==CellDescription::Type::LeafParentCoarsening ) {
    cellDescription.setType(CellDescription::Type::Erased);
    ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
    getCellDescriptions(cellDescriptionsIndex).erase(
        getCellDescriptions(cellDescriptionsIndex).begin()+element);
  }
}

bool exahype::solvers::ADERDGSolver::progressMeshRefinementInMergeWithMaster(
    const int worker,
    const int localCellDescriptionsIndex,
    const int localElement,
    const int coarseGridCellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level,
    const bool                                   stillInRefiningMode) {
  CellDescription& cellDescription = getCellDescription(localCellDescriptionsIndex,localElement);
  ensureFineGridCoarseGridConsistency(cellDescription,coarseGridCellDescriptionsIndex);
  #ifdef Asserts
  cellDescription.setCreation(CellDescription::Creation::ReceivedFromWorker);
  #endif

  const int coarseGridElement = tryGetElement(
      cellDescription.getParentIndex(),cellDescription.getSolverNumber());
  if ( // receive restricted data
      cellDescription.getType()==CellDescription::Type::LeafParentRequestedCoarsening
  ) {
    assertion1(coarseGridElement!=exahype::solvers::Solver::NotFound,cellDescription.toString());
    mergeWithWorkerOrMasterDataDueToForkOrJoin(worker,localCellDescriptionsIndex,localElement,
        peano::heap::MessageType::MasterWorkerCommunication,x,level); // assumes blocking/copy

    // use the received data
    CellDescription& coarseGridCellDescription =
        getCellDescription(cellDescription.getParentIndex(),coarseGridElement); // TODO(Dominic): Have helper function for that

    restrictVolumeDataIfErasingRequested(cellDescription,coarseGridCellDescription);
  }

  // work with the data
  if ( cellDescription.getType()==CellDescription::Type::LeafParentCoarsens ) {
    CellDescription& coarseGridCellDescription =
        getCellDescription(cellDescription.getParentIndex(),coarseGridElement); // TODO(Dominic): Have helper function for that
    assertion2( coarseGridCellDescription.getType()==CellDescription::Type::ParentCoarsening,
        cellDescription.toString(),coarseGridCellDescription.toString());

    eraseCellDescriptionIfNecessary(localCellDescriptionsIndex,localElement,coarseGridCellDescription);
  }

  // potentially veto
  if ( coarseGridElement != exahype::solvers::Solver::NotFound ) {
    CellDescription& coarseGridCellDescription = getCellDescription(
        cellDescription.getParentIndex(),coarseGridElement);
    updateCoarseGridAncestorRefinementStatus(cellDescription,coarseGridCellDescription);

    if (
      coarseGridCellDescription.getType()==CellDescription::Type::ParentCoarseningRequested &&
      cellDescription.getHasVirtualChildren()
    ) {            // no lock required; serial context
       coarseGridCellDescription.setType(CellDescription::Type::ParentKeep);
    }
  }

  progressCollectiveRefinementOperationsInLeaveCell(cellDescription,stillInRefiningMode);
  // ignore return value as responsibiliy is still on fine grid.

  // check if any cell description requires vertical communication
  bool solverRequiresVerticalCommunication = false;

  // block erasing request of coarse grid cell description if deployed cell
  // does not want to be erased
  decideOnRefinement(cellDescription,stillInRefiningMode);

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

  #ifdef Asserts
  const tarch::la::Vector<DIMENSIONS,double> center = cellDescription.getOffset()+0.5*cellDescription.getSize();
  #endif
  assertion5(Vertex::equalUpToRelativeTolerance(x,center),x,center,level,cellDescription.getLevel(),tarch::parallel::Node::getInstance().getRank());
  assertion2(cellDescription.getLevel()==level,cellDescription.getLevel(),level);

  if ( cellDescription.getType()==CellDescription::Type::Leaf ) {
    logDebug("sendDataToWorkerOrMasterDueToForkOrJoin(...)",""
            "solution of solver " << cellDescription.getSolverNumber() << " sent to rank "<<toRank<<
                 ", cell: "<< x << ", level: " << level);

    assertion2(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()),
        cellDescriptionsIndex,cellDescription.toString());
    assertion2(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()),
            cellDescriptionsIndex,cellDescription.toString());
    DataHeap::getInstance().sendData(
        static_cast<double*>(cellDescription.getSolution()),
        getDataPerCell(), toRank, x, level,messageType);
    DataHeap::getInstance().sendData(
        static_cast<double*>(cellDescription.getPreviousSolution()),
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
  #ifdef Asserts
  const tarch::la::Vector<DIMENSIONS,double> center = cellDescription.getOffset()+0.5*cellDescription.getSize();
  #endif
  assertion5(Vertex::equalUpToRelativeTolerance(x,center),x,center,level,cellDescription.getLevel(),tarch::parallel::Node::getInstance().getRank());
  assertion2(cellDescription.getLevel()==level,cellDescription.getLevel(),level);

  // allocate memory
  if ( messageType==peano::heap::MessageType::ForkOrJoinCommunication ) {
    ensureNecessaryMemoryIsAllocated(cellDescription);
    ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
  } else if ( cellDescription.getType()==CellDescription::Type::Leaf ) {
    ensureNecessaryMemoryIsAllocated(cellDescription);
    ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
  }

  // receive data
  if ( cellDescription.getType()==CellDescription::Type::Leaf ) {
    logDebug("mergeWithRemoteDataDueToForkOrJoin(...)","[solution] receive from rank "<<fromRank<<
             ", cell: "<< x << ", level: " << level);

    getDataHeapEntries(cellDescription.getSolutionIndex()).clear();
    getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).clear();
    DataHeap::getInstance().receiveData(cellDescription.getSolutionIndex(),
        fromRank,x,level,messageType);
    DataHeap::getInstance().receiveData(cellDescription.getPreviousSolutionIndex(),
        fromRank,x,level,messageType);
  }
}

#endif
