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
 * \author Dominic E. Charrier, Tobias Weinzierl
 **/

#include "LimitingADERDGSolver.h"

#include "kernels/limiter/generic/Limiter.h"

#include "exahype/VertexOperations.h"

#include "exahype/amr/AdaptiveMeshRefinement.h"


namespace exahype {
namespace solvers {

} /* namespace solvers */
} /* namespace exahype */

tarch::logging::Log exahype::solvers::LimitingADERDGSolver::_log("exahype::solvers::LimitingADERDGSolver");

int exahype::solvers::LimitingADERDGSolver::getMaxMinimumHelperStatusForTroubledCell() {
  int result = 0;
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
    ) {
      result = std::max(
          result,
          static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver()->
          getMinimumLimiterStatusForTroubledCell());
    }
  }
  return result;
}

bool exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLimiterStatusSpreading() {
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
        &&
        (static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainChange()
        !=exahype::solvers::LimiterDomainChange::Regular
        ||
        solver->getMeshUpdateRequest())
    ) {
      return true;
    }
  }
  return false;
}

bool exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalOrGlobalRecomputation() {
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG &&
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainChange()
        !=exahype::solvers::LimiterDomainChange::Regular
    ) {
      return true;
    }
  }
  return false;
}

bool exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalRecomputation(){
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG &&
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainChange()
        ==exahype::solvers::LimiterDomainChange::Irregular
    ) {
      return true;
    }
  }
  return false;
}

bool exahype::solvers::LimitingADERDGSolver::oneSolverRequestedGlobalRecomputation(){
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG &&
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainChange()
        ==exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate
    ) {
      return true;
    }
  }
  return false;
}

bool exahype::solvers::LimitingADERDGSolver::isValidCellDescriptionIndex(
    const int cellDescriptionsIndex) const  {
  return _solver->isValidCellDescriptionIndex(cellDescriptionsIndex);
}

exahype::solvers::Solver::SubcellPosition exahype::solvers::LimitingADERDGSolver::computeSubcellPositionOfCellOrAncestor(
      const int cellDescriptionsIndex,
      const int element) const {
  return _solver->computeSubcellPositionOfCellOrAncestor(cellDescriptionsIndex,element);
}

exahype::solvers::LimitingADERDGSolver::LimitingADERDGSolver(
    const std::string& identifier,
    std::unique_ptr<exahype::solvers::ADERDGSolver> solver,
    std::unique_ptr<exahype::solvers::FiniteVolumesSolver> limiter,
    const double DMPRelaxationParameter,
    const double DMPDifferenceScaling,
    const int iterationsToCureTroubledCell)
    :
    exahype::solvers::Solver(identifier, Solver::Type::LimitingADERDG, solver->getNumberOfVariables(),
        solver->getNumberOfParameters(), solver->getNodesPerCoordinateAxis(), solver->getMaximumMeshSize(),
          solver->getMaximumAdaptiveMeshDepth(),
          solver->getTimeStepping()),
          _solver(std::move(solver)),
          _limiter(std::move(limiter)),
          _limiterDomainChange(LimiterDomainChange::Regular),
          _nextLimiterDomainChange(LimiterDomainChange::Regular),
          _DMPMaximumRelaxationParameter(DMPRelaxationParameter),
          _DMPDifferenceScaling(DMPDifferenceScaling),
          _iterationsToCureTroubledCell(iterationsToCureTroubledCell)
{
  assertion(_solver->getNumberOfParameters() == 0);
  assertion(_solver->getTimeStepping()==_limiter->getTimeStepping());
}

void exahype::solvers::LimitingADERDGSolver::updateNextMeshUpdateRequest(const bool& meshUpdateRequest)  {
  _solver->updateNextMeshUpdateRequest(meshUpdateRequest);
}
bool exahype::solvers::LimitingADERDGSolver::getNextMeshUpdateRequest() const  {
  return _solver->getNextMeshUpdateRequest();
}
bool exahype::solvers::LimitingADERDGSolver::getMeshUpdateRequest() const  {
  return _solver->getMeshUpdateRequest();
}
void exahype::solvers::LimitingADERDGSolver::setNextMeshUpdateRequest()  {
  _solver->setNextMeshUpdateRequest();
}

void exahype::solvers::LimitingADERDGSolver::updateNextAttainedStableState(const bool& attainedStableState)  {
  _solver->updateNextAttainedStableState(attainedStableState);
}
bool exahype::solvers::LimitingADERDGSolver::getNextAttainedStableState() const  {
  return _solver->getNextAttainedStableState();
}
bool exahype::solvers::LimitingADERDGSolver::getAttainedStableState() const  {
  return _solver->getAttainedStableState();
}
void exahype::solvers::LimitingADERDGSolver::setNextAttainedStableState()  {
  _solver->setNextAttainedStableState();
}

double exahype::solvers::LimitingADERDGSolver::getMinTimeStamp() const {
  return _solver->getMinTimeStamp();
}

double exahype::solvers::LimitingADERDGSolver::getMinTimeStepSize() const {
  return _solver->getMinTimeStepSize();
}

double exahype::solvers::LimitingADERDGSolver::getMinNextTimeStepSize() const {
  return _solver->getMinNextTimeStepSize();
}

void exahype::solvers::LimitingADERDGSolver::updateMinNextTimeStepSize(double value) {
  _solver->updateMinNextTimeStepSize(value);
}

void exahype::solvers::LimitingADERDGSolver::initSolver(
    const double timeStamp,
    const tarch::la::Vector<DIMENSIONS,double>& domainOffset,
    const tarch::la::Vector<DIMENSIONS,double>& domainSize,
    const tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize
) {
  _domainOffset=domainOffset;
  _domainSize=domainSize;
  _coarsestMeshLevel =
      exahype::solvers::Solver::computeMeshLevel(_maximumMeshSize,boundingBoxSize[0]);

  _limiterDomainChange=LimiterDomainChange::Regular;

  _solver->initSolver(timeStamp, domainOffset, domainSize, boundingBoxSize);
}

bool exahype::solvers::LimitingADERDGSolver::isSending(
    const exahype::records::State::AlgorithmSection& section) const {
  bool isSending = false;

  switch (section) {
    case exahype::records::State::AlgorithmSection::TimeStepping:
      isSending = true;
      break;
    case exahype::records::State::AlgorithmSection::LimiterStatusSpreading:
      isSending = false; // value doesn't actually matter
      break;
    case exahype::records::State::AlgorithmSection::MeshRefinement:
      isSending |= getMeshUpdateRequest();
      break;
    case exahype::records::State::AlgorithmSection::MeshRefinementOrGlobalRecomputation:
      isSending |= getMeshUpdateRequest();
      isSending |= getLimiterDomainChange()==exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate;
      break;
    case exahype::records::State::AlgorithmSection::MeshRefinementOrGlobalRecomputationAllSend:
      isSending = true;
      break;
    case exahype::records::State::AlgorithmSection::MeshRefinementOrLocalOrGlobalRecomputation:
      isSending |= getMeshUpdateRequest();
      isSending |= getLimiterDomainChange()==exahype::solvers::LimiterDomainChange::Irregular;
      isSending |= getLimiterDomainChange()==exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate;
      break;
    case exahype::records::State::AlgorithmSection::LocalRecomputationAllSend:
      isSending = true;
      break;
    case exahype::records::State::AlgorithmSection::PredictionRerunAllSend:
      isSending = true;
  }

  return isSending;
}

bool exahype::solvers::LimitingADERDGSolver::isComputing(
    const exahype::records::State::AlgorithmSection& section) const {
  bool isComputing = false;

  switch (section) {
    case exahype::records::State::AlgorithmSection::TimeStepping:
      isComputing = true;
      break;
    case exahype::records::State::AlgorithmSection::LimiterStatusSpreading:
      isComputing |= getMeshUpdateRequest();
      isComputing |= getLimiterDomainChange()==exahype::solvers::LimiterDomainChange::Irregular;
      isComputing |= getLimiterDomainChange()==exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate;
      break;
    case exahype::records::State::AlgorithmSection::MeshRefinement:
      isComputing = getMeshUpdateRequest();
      break;
    case exahype::records::State::AlgorithmSection::MeshRefinementOrGlobalRecomputationAllSend:
      isComputing |= getMeshUpdateRequest();
      isComputing |= getLimiterDomainChange()==exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate;
      break;
    case exahype::records::State::AlgorithmSection::MeshRefinementOrGlobalRecomputation:
      isComputing |= getMeshUpdateRequest();
      isComputing |= getLimiterDomainChange()==exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate;
      break;
    case exahype::records::State::AlgorithmSection::MeshRefinementOrLocalOrGlobalRecomputation:
      isComputing |= getMeshUpdateRequest();
      isComputing |= getLimiterDomainChange()==exahype::solvers::LimiterDomainChange::Irregular;
      isComputing |= getLimiterDomainChange()==exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate;
      break;
    case exahype::records::State::AlgorithmSection::LocalRecomputationAllSend:
      isComputing = getLimiterDomainChange()==exahype::solvers::LimiterDomainChange::Irregular;
      break;
    case exahype::records::State::AlgorithmSection::PredictionRerunAllSend:
      isComputing = _solver->getStabilityConditionWasViolated();
  }

  return isComputing;
}

void exahype::solvers::LimitingADERDGSolver::synchroniseTimeStepping(
    const int cellDescriptionsIndex,
    const int solverElement) {
  _solver->synchroniseTimeStepping(cellDescriptionsIndex,solverElement);
  ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::startNewTimeStep() {
  _solver->startNewTimeStep();

  _minCellSize     = _nextMinCellSize;
  _maxCellSize     = _nextMaxCellSize;
  _nextMinCellSize = std::numeric_limits<double>::max();
  _nextMaxCellSize = -std::numeric_limits<double>::max(); // "-", min

  logDebug("startNewTimeStep()","_limiterDomainHasChanged="<<static_cast<int>(_limiterDomainChange)<<
           ",nextLimiterDomainChange="<<static_cast<int>(_nextLimiterDomainChange));
}

void exahype::solvers::LimitingADERDGSolver::updateTimeStepSizesFused()  {
  _solver->updateTimeStepSizesFused();
}

void exahype::solvers::LimitingADERDGSolver::updateTimeStepSizes()  {
  _solver->updateTimeStepSizes();
}

void exahype::solvers::LimitingADERDGSolver::zeroTimeStepSizes() {
  _solver->zeroTimeStepSizes();
}

exahype::solvers::LimiterDomainChange
exahype::solvers::LimitingADERDGSolver::getNextLimiterDomainChange() const {
  return _nextLimiterDomainChange;
}

void exahype::solvers::LimitingADERDGSolver::setNextLimiterDomainChange() {
  _limiterDomainChange     = _nextLimiterDomainChange;
  _nextLimiterDomainChange = LimiterDomainChange::Regular;
}
void exahype::solvers::LimitingADERDGSolver::updateNextLimiterDomainChange(
    exahype::solvers::LimiterDomainChange limiterDomainChange) {
  _nextLimiterDomainChange =
      std::max( _nextLimiterDomainChange, limiterDomainChange );
}
exahype::solvers::LimiterDomainChange
exahype::solvers::LimitingADERDGSolver::getLimiterDomainChange() const {
  return _limiterDomainChange;
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep() {
  _solver->rollbackToPreviousTimeStep();
}

void exahype::solvers::LimitingADERDGSolver::reconstructStandardTimeSteppingDataAfterRollback() {
  _solver->reconstructStandardTimeSteppingDataAfterRollback();
}

void exahype::solvers::LimitingADERDGSolver::updateNextMinCellSize(double minCellSize) {
  _solver->updateNextMinCellSize(minCellSize);
  _nextMinCellSize = _solver->getNextMinCellSize();
}

void exahype::solvers::LimitingADERDGSolver::updateNextMaxCellSize(double maxCellSize) {
  _solver->updateNextMaxCellSize(maxCellSize);
  _nextMaxCellSize = _solver->getNextMaxCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getNextMinCellSize() const {
  return _solver->getNextMinCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getNextMaxCellSize() const {
  return _solver->getNextMaxCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getMinCellSize() const {
  return _solver->getMinCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getMaxCellSize() const {
  return _solver->getMaxCellSize();
}

///////////////////////////////////
// MODIFY CELL DESCRIPTION
///////////////////////////////////
bool exahype::solvers::LimitingADERDGSolver::markForRefinementBasedOnLimiterStatus(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const bool initialGrid,
    const int solverNumber) {
  const int solverElement =
      _solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);

  if (solverElement!=exahype::solvers::Solver::NotFound) {
    SolverPatch& solverPatch =
        ADERDGSolver::getCellDescription(fineGridCell.getCellDescriptionsIndex(),solverElement);
    bool refineFineGridCell =
        evaluateLimiterStatusRefinementCriterion(
            fineGridCell.getCellDescriptionsIndex(),solverElement);
//    refineFineGridCell |=
//        evaluateDiscreteMaximumPrincipleRefinementCriterion(
//            fineGridCell.getCellDescriptionsIndex(),solverElement); // TODO(Dominic): Reenable

    if (refineFineGridCell) {
      solverPatch.setRefinementEvent(SolverPatch::RefinementEvent::RefiningRequested);
      return true;
    }

    if (solverPatch.getLimiterStatus()<
        computeMinimumLimiterStatusForRefinement(solverPatch.getLevel())) {
      return _solver->markForRefinement(
              fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
              coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
              fineGridPositionOfCell,
              initialGrid,
              solverNumber);
    }
  }
  return false;
}

bool exahype::solvers::LimitingADERDGSolver::updateLimiterStatusDuringLimiterStatusSpreading(
    const int cellDescriptionsIndex, const int solverElement) const {
  SolverPatch& solverPatch =
      ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  if (solverPatch.getLimiterStatus()<_solver->getMinimumLimiterStatusForTroubledCell()) {
    solverPatch.setLimiterStatus(ADERDGSolver::determineLimiterStatus(solverPatch));
  } else {
    solverPatch.setLimiterStatus(_solver->getMinimumLimiterStatusForTroubledCell());
    solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1);
  }
  solverPatch.setFacewiseLimiterStatus(0);
  deallocateLimiterPatchOnHelperCell(cellDescriptionsIndex,solverElement);
  if (solverPatch.getPreviousLimiterStatus()==0) { // We might throw away old valuable information if we do not perform this check TODO(Dominic): Add to docu
    ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(cellDescriptionsIndex,solverElement);
  }
  return ensureRequiredLimiterPatchIsAllocated(cellDescriptionsIndex,solverElement);
}

bool exahype::solvers::LimitingADERDGSolver::markForRefinement(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const bool initialGrid,
    const int solverNumber)  {
  // update limiter status if neighbours are in a consistent state
  const int fineGridSolverElement =
      _solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (fineGridSolverElement!=exahype::solvers::Solver::NotFound) {
    const tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D,int>&
    indicesAdjacentToFineGridVertices =
        VertexOperations::readCellDescriptionsIndex(
            fineGridVerticesEnumerator,fineGridVertices);
    if (multiscalelinkedcell::adjacencyInformationIsConsistent(
        indicesAdjacentToFineGridVertices)) {
      bool refineFineGridCell =
          markForRefinementBasedOnLimiterStatus(
              fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
              coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
              fineGridPositionOfCell,
              initialGrid,
              solverNumber);

      vetoErasingChildrenRequestBasedOnLimiterStatus(
          fineGridCell.getCellDescriptionsIndex(),fineGridSolverElement,
          coarseGridCell.getCellDescriptionsIndex());

      return refineFineGridCell;
    }
  }

  return false;
}

exahype::solvers::Solver::UpdateStateInEnterCellResult exahype::solvers::LimitingADERDGSolver::updateStateInEnterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const bool initialGrid,
    const int solverNumber)  {
  UpdateStateInEnterCellResult result =
      _solver->updateStateInEnterCell(
          fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
          coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
          fineGridPositionOfCell,initialGrid,solverNumber);

  return result;
}

void exahype::solvers::LimitingADERDGSolver::vetoErasingChildrenRequestBasedOnLimiterStatus(
    const int fineGridCellDescriptionsIndex,
    const int fineGridSolverElement,
    const int coarseGridCellDescriptionsIndex) const {
  SolverPatch& fineGridSolverPatch = ADERDGSolver::getCellDescription(
          fineGridCellDescriptionsIndex,fineGridSolverElement);

  if (
      fineGridSolverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::ErasingChildrenRequested
      ||
      fineGridSolverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::ChangeChildrenToDescendantsRequested
  ) {
    fineGridSolverPatch.setRefinementEvent(SolverPatch::RefinementEvent::None); // TODO(Dominic): Reenable erasing again later on
//    if (fineGridSolverPatch.getLimiterStatus()>=
//        computeMinimumLimiterStatusForRefinement(fineGridSolverPatch.getLevel())
//        ||
//        fineGridSolverPatch.getPreviousLimiterStatus()>=
//        computeMinimumLimiterStatusForRefinement(fineGridSolverPatch.getLevel()) // TODO(Dominic): Add to docu: This is necessary for not erasing cells in global recomputation and it further adds some laziness in erasing.
//    ) {
//      fineGridSolverPatch.setRefinementEvent(SolverPatch::RefinementEvent::None);
//    }
  }
}

bool exahype::solvers::LimitingADERDGSolver::updateStateInLeaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber)  {
  const int solverElement =
      _solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  const int parentSolverElement =
      _solver->tryGetElement(coarseGridCell.getCellDescriptionsIndex(),solverNumber);
  if (solverElement!=exahype::solvers::Solver::NotFound
      &&
      parentSolverElement!=exahype::solvers::Solver::NotFound
  ) {
    _solver->restrictLimiterStatus(
        fineGridCell.getCellDescriptionsIndex(),solverElement,
        coarseGridCell.getCellDescriptionsIndex(),parentSolverElement);
  }

  return _solver->updateStateInLeaveCell(
      fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
      coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
      fineGridPositionOfCell,solverNumber);
}

bool exahype::solvers::LimitingADERDGSolver::attainedStableState(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    const int solverNumber) const {
  return _solver->attainedStableState(
      fineGridCell,fineGridVertices,fineGridVerticesEnumerator,solverNumber);
}

void exahype::solvers::LimitingADERDGSolver::finaliseStateUpdates(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) {
  _solver->finaliseStateUpdates(
      fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
      coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
      fineGridPositionOfCell,solverNumber);
}

///////////////////////////////////
// CELL-LOCAL
//////////////////////////////////
int exahype::solvers::LimitingADERDGSolver::computeMinimumLimiterStatusForRefinement(
    int level) const {
  const int levelDelta = getMaximumAdaptiveMeshLevel() - level;

  assertion(levelDelta>=0);
  const int status = _solver->getMinimumLimiterStatusForActiveFVPatch()-1 + (levelDelta-1);
  return (_solver->getMinimumLimiterStatusForTroubledCell()-status>=2) ?
      status : _solver->getMinimumLimiterStatusForTroubledCell()-2;
}

bool exahype::solvers::LimitingADERDGSolver::evaluateLimiterStatusRefinementCriterion(
    const int cellDescriptionsIndex,const int solverElement) const {
  if (solverElement!=exahype::solvers::Solver::NotFound) {
    SolverPatch& solverPatch =
        ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
    if (
        solverPatch.getType()==SolverPatch::Type::Cell
        &&
        solverPatch.getLevel() < getMaximumAdaptiveMeshLevel()
        &&
        (solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None ||
         solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::DeaugmentingChildrenRequested ||
         solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::AugmentingRequested)
    ) {
      return solverPatch.getLimiterStatus() >=
                computeMinimumLimiterStatusForRefinement(solverPatch.getLevel());
    }
  }
  return false;
}

bool exahype::solvers::LimitingADERDGSolver::evaluateDiscreteMaximumPrincipleRefinementCriterion(
    const int cellDescriptionsIndex,
    const int solverElement) const {
  if (solverElement!=exahype::solvers::Solver::NotFound) {
    SolverPatch& solverPatch =
        ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
    const int numberOfObservables = _solver->getDMPObservables();
    if (solverPatch.getType()==SolverPatch::Type::Cell
        &&
        numberOfObservables > 0
    ) {
      double* observablesMin = DataHeap::getInstance().getData(
          solverPatch.getSolutionMin()).data();
      double* observablesMax = DataHeap::getInstance().getData(
          solverPatch.getSolutionMax()).data();

      bool discreteMaximumPrincipleSatisfied = true;
      if (
          solverPatch.getLevel() < getMaximumAdaptiveMeshLevel()
          &&
          (solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None ||
         solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::DeaugmentingChildrenRequested ||
         solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::AugmentingRequested)
      ) {
        const double* localMinPerObservable    = observablesMin+0;
        const double* localMaxPerObservable    = observablesMin+0;
        const double* boundaryMinPerObservable = observablesMin+numberOfObservables;
        const double* boundaryMaxPerObservable = observablesMax+numberOfObservables;

        const double differenceScaling   = _DMPDifferenceScaling;
        const double relaxationParameter = _DMPMaximumRelaxationParameter * 0.1; // TODO(Dominic): Should be a parameter in the spec file (optional)

        for(int v = 0; v < numberOfObservables; v++) {
          double boundaryMin = boundaryMinPerObservable[v];
          double boundaryMax = boundaryMaxPerObservable[v];
          double scaledDifference = (boundaryMax - boundaryMin) * differenceScaling;

          assertion5(tarch::la::greaterEquals(scaledDifference,0.0),scaledDifference,boundaryMin,boundaryMax,localMinPerObservable[v],localMaxPerObservable[v]);
          scaledDifference = std::max( scaledDifference, relaxationParameter );

          if((localMinPerObservable[v] < (boundaryMin - scaledDifference)) ||
              (localMaxPerObservable[v] > (boundaryMax + scaledDifference))) {
            discreteMaximumPrincipleSatisfied=false;
          }
        }

        return false; // TODO(Dominic): Reenable
//        return !discreteMaximumPrincipleSatisfied;
      }

      // copy local min and max onto boundary
      //  TODO(Dominic):
      // 2. Copy the result on the other faces as well
      for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
        std::copy_n(
            observablesMin,numberOfObservables, // past-the-end element
            observablesMin+i*numberOfObservables);
        std::copy_n(
            observablesMax,numberOfObservables, // past-the-end element
            observablesMax+i*numberOfObservables);
      }
    }
  }
  return false;
}

bool
exahype::solvers::LimitingADERDGSolver::evaluateRefinementCriterionAfterSolutionUpdate(
      const int cellDescriptionsIndex,const int element) {
  // First evaluate the limiter status based refinement criterion
  bool refinementRequested =
      evaluateLimiterStatusRefinementCriterion(cellDescriptionsIndex,element);

  refinementRequested |=
      evaluateDiscreteMaximumPrincipleRefinementCriterion(cellDescriptionsIndex,element);

  // If no refinement was requested, evaluate the user's refinement criterion
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);
  if (!refinementRequested &&
      solverPatch.getLimiterStatus() <
      computeMinimumLimiterStatusForRefinement(solverPatch.getLevel())
  ) {
    refinementRequested =
        _solver->evaluateRefinementCriterionAfterSolutionUpdate(
            cellDescriptionsIndex,element);
  }

  return refinementRequested;
}


double exahype::solvers::LimitingADERDGSolver::startNewTimeStep(
    const int cellDescriptionsIndex,
    const int solverElement)  {
  double admissibleTimeStepSize =
      _solver->startNewTimeStep(cellDescriptionsIndex,solverElement);
  ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);

  return admissibleTimeStepSize;
}


void exahype::solvers::LimitingADERDGSolver::reconstructStandardTimeSteppingData(
    const int cellDescriptionsIndex,
    const int element) const {
  _solver->reconstructStandardTimeSteppingData(cellDescriptionsIndex,element);

  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  const int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
    limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
    limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());
  }
}

double exahype::solvers::LimitingADERDGSolver::updateTimeStepSizesFused(
      const int cellDescriptionsIndex,
      const int solverElement) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    const double admissibleTimeStepSize =
        _solver->updateTimeStepSizesFused(cellDescriptionsIndex,solverElement);

    ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);

    return admissibleTimeStepSize;
  }

  return std::numeric_limits<double>::max();
}

double exahype::solvers::LimitingADERDGSolver::updateTimeStepSizes(
      const int cellDescriptionsIndex,
      const int solverElement) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    const double admissibleTimeStepSize =
        _solver->updateTimeStepSizes(cellDescriptionsIndex,solverElement);

    ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);

    return admissibleTimeStepSize;
  }

  return std::numeric_limits<double>::max();
}

void exahype::solvers::LimitingADERDGSolver::zeroTimeStepSizes(
    const int cellDescriptionsIndex, const int solverElement) const {
  _solver->zeroTimeStepSizes(cellDescriptionsIndex,solverElement);
  ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep(
    const int cellDescriptionsIndex,
    const int solverElement) {
  _solver->rollbackToPreviousTimeStep(cellDescriptionsIndex,solverElement);
  ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::reconstructStandardTimeSteppingDataAfterRollback(
      const int cellDescriptionsIndex,
      const int solverElement) const {
  _solver->reconstructStandardTimeSteppingDataAfterRollback(cellDescriptionsIndex,solverElement);
  ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::adjustSolution(
    const int cellDescriptionsIndex,
    const int solverElement) {
  _solver->adjustSolution(
      cellDescriptionsIndex,solverElement);

  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    SolverPatch& solverPatch   = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
    LimiterPatch& limiterPatch = FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex,limiterElement);
    copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);

    _limiter->adjustSolution(
        cellDescriptionsIndex,limiterElement);
  }
}

exahype::solvers::LimitingADERDGSolver::LimiterPatch& exahype::solvers::LimitingADERDGSolver::getLimiterPatchForSolverPatch(
    const int cellDescriptionsIndex, const SolverPatch& solverPatch) const {
  assertion(solverPatch.getLimiterStatus()>=0);
  const int limiterElement =
      _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  return getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex,limiterElement);
}

void exahype::solvers::LimitingADERDGSolver::ensureLimiterPatchTimeStepDataIsConsistent(
    const int cellDescriptionsIndex, const int solverElement) const {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    assertion(solverPatch.getPreviousLimiterStatus()>0 || solverPatch.getLimiterStatus()>0);
    copyTimeStepDataFromSolverPatch(solverPatch,cellDescriptionsIndex,limiterElement);
  }
}

void exahype::solvers::LimitingADERDGSolver::copyTimeStepDataFromSolverPatch(
    const SolverPatch& solverPatch, const int cellDescriptionsIndex, const int limiterElement) {
  LimiterPatch& limiterPatch =
      FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex,limiterElement);
  copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
}

void exahype::solvers::LimitingADERDGSolver::copyTimeStepDataFromSolverPatch(
    const SolverPatch& solverPatch, LimiterPatch& limiterPatch) {
  limiterPatch.setPreviousTimeStepSize(solverPatch.getPreviousCorrectorTimeStepSize());
  limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
  limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());
}

exahype::solvers::LimitingADERDGSolver::LimiterPatch& exahype::solvers::LimitingADERDGSolver::getLimiterPatchForSolverPatch(
    const SolverPatch& solverPatch, const int cellDescriptionsIndex, const int limiterElement) const {
  assertion(limiterElement!=exahype::solvers::Solver::NotFound);
  LimiterPatch& limiterPatch =
      FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex,limiterElement);
  // Ensure time stamps and step sizes are consistent
  copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
  return limiterPatch;
}

exahype::solvers::Solver::UpdateResult exahype::solvers::LimitingADERDGSolver::fusedTimeStep(
    const int cellDescriptionsIndex,
    const int element,
    double** tempSpaceTimeUnknowns,
    double** tempSpaceTimeFluxUnknowns,
    double*  tempUnknowns,
    double*  tempFluxUnknowns,
    double** tempPointForceSources
) {
  SolverPatch& cellDescription =
        _solver->getCellDescription(cellDescriptionsIndex,element);
  // solver->synchroniseTimeStepping(cellDescription); // assumes this was done in neighbour merge

  updateSolution(cellDescriptionsIndex,element);

  UpdateResult result;
  result._limiterDomainChange =
      updateLimiterStatusAndMinAndMaxAfterSolutionUpdate(cellDescriptionsIndex,element);
  result._refinementRequested=
      evaluateRefinementCriterionAfterSolutionUpdate(cellDescriptionsIndex,element);

  if (cellDescription.getLimiterStatus()<_solver->getMinimumLimiterStatusForTroubledCell()) {
    _solver->performPredictionAndVolumeIntegral(
        cellDescription,
        tempSpaceTimeUnknowns,tempSpaceTimeFluxUnknowns,
        tempUnknowns,tempFluxUnknowns,tempPointForceSources);
  }

  result._timeStepSize=startNewTimeStep(cellDescriptionsIndex,element);
  return result;
}


/**
 * This method assumes the ADERDG solver's cell-local limiter status has
 * already been determined.
 */
void exahype::solvers::LimitingADERDGSolver::updateSolution(
    const int cellDescriptionsIndex,
    const int element)  {
  SolverPatch& solverPatch =
      ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  // 0. Erase old cells; now it's safe (TODO(Dominic): Add to docu)
  deallocateLimiterPatchOnHelperCell(cellDescriptionsIndex,element);

  // 1. Write back the limiter status to the previous limiter status field
  solverPatch.setPreviousLimiterStatus(solverPatch.getLimiterStatus()); // TODO(Dominic): This might cause problems
                                                                        // for the global recomputation

  // 2. Update the solution in the cells
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
      const int limiterElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(solverPatch.getLimiterStatus()>=0);

      if (solverPatch.getLimiterStatus()==0) {
        _solver->updateSolution(solverPatch);
      }
      else if (solverPatch.getLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch()) {
        _solver->updateSolution(solverPatch);

        LimiterPatch& limiterPatch =
            getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);

        double* solverSolution = DataHeap::getInstance().getData(
            solverPatch.getSolution()).data();
        _limiter->swapSolutionAndPreviousSolution(limiterPatch);
        double* limiterSolution = DataHeap::getInstance().getData(
            limiterPatch.getSolution()).data();
        kernels::limiter::generic::c::projectOnFVLimiterSpace( // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
            solverSolution,_solver->getNumberOfVariables(),
            _solver->getNodesPerCoordinateAxis(),
            _limiter->getGhostLayerWidth(),
            limiterSolution);
      }
      else { // solverPatch.getLimiterStatus()>=ADERDGSolver::MinimumLimiterStatusForActiveFVPatch
        assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());

        LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
        _limiter->updateSolution(cellDescriptionsIndex,limiterElement);

        double* solverSolution = DataHeap::getInstance().getData(
            solverPatch.getSolution()).data();
        double* limiterSolution = DataHeap::getInstance().getData(
            limiterPatch.getSolution()).data();

        kernels::limiter::generic::c::projectOnDGSpace(
            limiterSolution,_solver->getNumberOfVariables(),
            _solver->getNodesPerCoordinateAxis(), // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
            _limiter->getGhostLayerWidth(),
            solverSolution);
      }

      // 3. Only after the solution update, we are allowed to remove limiter patches.
      ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(cellDescriptionsIndex,element);
    } else {
      _solver->updateSolution(solverPatch);
    }
  }
}

exahype::solvers::LimiterDomainChange
exahype::solvers::LimitingADERDGSolver::updateLimiterStatusAndMinAndMaxAfterSolutionUpdate(
    const int cellDescriptionsIndex,
    const int solverElement) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  LimiterDomainChange limiterDomainChange = LimiterDomainChange::Regular;
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    bool solutionIsValid =
        evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(solverPatch)
        && evaluatePhysicalAdmissibilityCriterion(solverPatch); // after min and max was found
    if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
        solverPatch.getLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch()) {
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
      determineLimiterMinAndMax(solverPatch,limiterPatch);
    } // else: Keep the previously computed min and max values

    limiterDomainChange =
        determineLimiterStatusAfterSolutionUpdate(solverPatch,!solutionIsValid);

    allocateLimiterPatchAfterSolutionUpdate(cellDescriptionsIndex,solverElement);
  } else {
    solverPatch.setLimiterStatus(ADERDGSolver::determineLimiterStatus(solverPatch));
    deallocateLimiterPatchOnHelperCell(cellDescriptionsIndex,solverElement);
  }
  solverPatch.setFacewiseLimiterStatus(0);

  return limiterDomainChange;
}

void exahype::solvers::LimitingADERDGSolver::allocateLimiterPatchAfterSolutionUpdate(
    const int cellDescriptionsIndex,const int solverElement) const {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  assertion(solverPatch.getType()==SolverPatch::Type::Cell);

  bool limiterPatchAllocated =
      ensureRequiredLimiterPatchIsAllocated(cellDescriptionsIndex,solverPatch.getSolverNumber());
  if (limiterPatchAllocated &&
      solverPatch.getLimiterStatus()>=0 &&
      solverPatch.getLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch()) {
    assertion(solverPatch.getLevel()==getMaximumAdaptiveMeshLevel());
    LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
    projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
  }
}

exahype::solvers::LimiterDomainChange
exahype::solvers::LimitingADERDGSolver::determineLimiterStatusAfterSolutionUpdate(
    SolverPatch& solverPatch,const bool isTroubled) const {
  assertion1(solverPatch.getType()==SolverPatch::Type::Cell,solverPatch.toString());

  LimiterDomainChange limiterDomainChange = LimiterDomainChange::Regular;
  if (isTroubled) {
    solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1);
    if (solverPatch.getLimiterStatus()>_solver->getMinimumLimiterStatusForActiveFVPatch()) {
      // TODO(Dominic): Add to docu: width=2: 5->5 and 4->5 are okay. width=1: 3->3 is okay
      limiterDomainChange = LimiterDomainChange::Regular;
      solverPatch.setLimiterStatus(_solver->getMinimumLimiterStatusForTroubledCell());
    } else {
      limiterDomainChange = LimiterDomainChange::Irregular;
      solverPatch.setLimiterStatus(_solver->getMinimumLimiterStatusForTroubledCell());
    }
    if (solverPatch.getLevel()<getMaximumAdaptiveMeshLevel()) {
      limiterDomainChange = LimiterDomainChange::IrregularRequiringMeshUpdate;
    }
  } else {
    if (solverPatch.getPreviousLimiterStatus()>=_solver->getMinimumLimiterStatusForTroubledCell()) {
      solverPatch.setLimiterStatus(_solver->getMinimumLimiterStatusForTroubledCell());
      solverPatch.setIterationsToCureTroubledCell(
          solverPatch.getIterationsToCureTroubledCell()-1);
      if (solverPatch.getIterationsToCureTroubledCell()==0) {
        solverPatch.setLimiterStatus(_solver->getMinimumLimiterStatusForTroubledCell()-1);
        solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1); // TODO(Dominic): Probably not necessary
      }
    } else { // merge limiter status normally
      solverPatch.setLimiterStatus(ADERDGSolver::determineLimiterStatus(solverPatch));
    }
  }

  assertion2(
      !isTroubled ||
      solverPatch.getLimiterStatus()>=_solver->getMinimumLimiterStatusForTroubledCell(),
      isTroubled,
      solverPatch.getLimiterStatus());
  return limiterDomainChange;
}

bool exahype::solvers::LimitingADERDGSolver::evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(SolverPatch& solverPatch) {
  double* solution = DataHeap::getInstance().getData(
      solverPatch.getSolution()).data();

  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    double* observablesMin = DataHeap::getInstance().getData(
        solverPatch.getSolutionMin()).data();
    double* observablesMax = DataHeap::getInstance().getData(
        solverPatch.getSolutionMax()).data();

    // 1. Check if the DMP is satisfied and search for the min and max
    // Write the new min and max to the storage reserved for face 0
    bool dmpIsSatisfied = kernels::limiter::generic::c::discreteMaximumPrincipleAndMinAndMaxSearch(
          solution,_solver.get(),
          _DMPMaximumRelaxationParameter, _DMPDifferenceScaling,
          observablesMin,observablesMax);

    // TODO(Dominic):
//    // 2. Copy the result on the other faces as well
//    for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
//      std::copy_n(
//          observablesMin,numberOfObservables, // past-the-end element
//          observablesMin+i*numberOfObservables);
//      std::copy_n(
//          observablesMax,numberOfObservables, // past-the-end element
//          observablesMax+i*numberOfObservables);
//    }

    return dmpIsSatisfied;
  } else {
    return true;
  }
}

void
exahype::solvers::LimitingADERDGSolver::updateLimiterStatusAndMinAndMaxAfterAdjustSolution(
    const int cellDescriptionsIndex,
    const int solverElement) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  if (
      solverPatch.getType()==SolverPatch::Type::Cell
  ) {
    determineSolverMinAndMax(solverPatch);
    if (!evaluatePhysicalAdmissibilityCriterion(solverPatch)) {
      solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1);
      solverPatch.setLimiterStatus(_solver->getMinimumLimiterStatusForTroubledCell());
    }
  }
}

bool exahype::solvers::LimitingADERDGSolver::evaluatePhysicalAdmissibilityCriterion(SolverPatch& solverPatch) {
  double* observablesMin = nullptr;
  double* observablesMax = nullptr;

  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables > 0) {
    observablesMin = DataHeap::getInstance().getData(
        solverPatch.getSolutionMin()).data();
    observablesMax = DataHeap::getInstance().getData(
        solverPatch.getSolutionMax()).data();
  }

  const double* const solution = DataHeap::getInstance().getData(
        solverPatch.getSolution()).data();

  return _solver->isPhysicallyAdmissible(
      solution,
      observablesMin,observablesMax,numberOfObservables,
      solverPatch.getOffset()+0.5*solverPatch.getSize(),solverPatch.getSize(),
      solverPatch.getCorrectorTimeStamp(),solverPatch.getCorrectorTimeStepSize());
}

void exahype::solvers::LimitingADERDGSolver::determineMinAndMax(
    const int cellDescriptionsIndex,
    const int solverElement) {
  SolverPatch& solverPatch =
      ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
      assertion(solverPatch.getLimiterStatus()>=0);
      if (solverPatch.getLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch()) {
        determineSolverMinAndMax(solverPatch);
      } else { // solverPatch.getLimiterStatus()>=ADERDGSolver::MinimumLimiterStatusForActiveFVPatch
        LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
        determineLimiterMinAndMax(solverPatch,limiterPatch);
      }
    } else {
      determineSolverMinAndMax(solverPatch);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::determineSolverMinAndMax(SolverPatch& solverPatch) {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolution()),
            solverPatch.toString());
    assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMin()),
            solverPatch.toString());

    const double* const solution = DataHeap::getInstance().getData(
        solverPatch.getSolution()).data();

    double* observablesMin = DataHeap::getInstance().getData(
        solverPatch.getSolutionMin()).data();
    double* observablesMax = DataHeap::getInstance().getData(
        solverPatch.getSolutionMax()).data();

    // Write the result to the face with index "0"
    kernels::limiter::generic::c::findCellLocalMinAndMax(
        solution,_solver.get(),
        observablesMin,observablesMax);

    // Copy the result on the other faces as well
    for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
      std::copy_n(
          observablesMin,numberOfObservables, // past-the-end element
          observablesMin+i*numberOfObservables);
      std::copy_n(
          observablesMax,numberOfObservables, // past-the-end element
          observablesMax+i*numberOfObservables);
    }

    for (int i=0; i<DIMENSIONS_TIMES_TWO*numberOfObservables; ++i) {
      assertion2(*(observablesMin+i)<std::numeric_limits<double>::max(),i,solverPatch.toString());
      assertion2(*(observablesMax+i)>-std::numeric_limits<double>::max(),i,solverPatch.toString());
    } // Dead code elimination will get rid of this loop
  }
}

void exahype::solvers::LimitingADERDGSolver::determineLimiterMinAndMax(SolverPatch& solverPatch,LimiterPatch& limiterPatch) {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    double* limiterSolution = DataHeap::getInstance().getData(
        limiterPatch.getSolution()).data();

    double* observablesMin = DataHeap::getInstance().getData(
        solverPatch.getSolutionMin()).data();
    double* observablesMax = DataHeap::getInstance().getData(
        solverPatch.getSolutionMax()).data();

    // Write the result to the face with index "0"
    kernels::limiter::generic::c::findCellLocalLimiterMinAndMax(
        limiterSolution,_solver.get(),
        _limiter->getGhostLayerWidth(),observablesMin,observablesMax);

    // Copy the result on the other faces as well
    const int numberOfObservables = _solver->getDMPObservables();
    for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
      std::copy_n(
          observablesMin,numberOfObservables, // past-the-end element
          observablesMin+i*numberOfObservables);
      std::copy_n(
          observablesMax,numberOfObservables, // past-the-end element
          observablesMax+i*numberOfObservables);
    }

    for (int i=0; i<DIMENSIONS_TIMES_TWO*numberOfObservables; ++i) {
      assertion(*(observablesMin+i)<std::numeric_limits<double>::max());
      assertion(*(observablesMax+i)>-std::numeric_limits<double>::max());
    } // Dead code elimination will get rid of this loop
  }
}

void exahype::solvers::LimitingADERDGSolver::deallocateLimiterPatch(
        const int cellDescriptionsIndex,
        const int solverElement) const {
  assertion(solverElement!=exahype::solvers::Solver::NotFound);
  const int limiterElement =
        tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  assertion(limiterElement!=exahype::solvers::Solver::NotFound);
  LimiterPatch& limiterPatch = FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex,limiterElement);
  limiterPatch.setType(LimiterPatch::Type::Erased);
  _limiter->ensureNoUnnecessaryMemoryIsAllocated(limiterPatch);

  tarch::multicore::Lock lock(exahype::HeapSemaphore);
  LimiterHeap::getInstance().getData(cellDescriptionsIndex).erase(
      LimiterHeap::getInstance().getData(cellDescriptionsIndex).begin()+limiterElement);
  lock.free();
}

void exahype::solvers::LimitingADERDGSolver::deallocateLimiterPatchOnHelperCell(
    const int cellDescriptionsIndex,
    const int solverElement) const {
  assertion(solverElement!=exahype::solvers::Solver::NotFound);
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);

  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    switch(solverPatch.getType()) {
    case SolverPatch::Type::Ancestor:
    case SolverPatch::Type::Descendant:
    case SolverPatch::Type::Erased: {
      deallocateLimiterPatch(cellDescriptionsIndex,solverElement);
    } break;
    case SolverPatch::Type::Cell: {
      // do nothing
    } break;
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(
    const int cellDescriptionsIndex,
    const int solverElement) const {
  assertion(solverElement!=exahype::solvers::Solver::NotFound);
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  if (limiterElement!=exahype::solvers::Solver::NotFound
      && solverPatch.getLimiterStatus()==0) {
    deallocateLimiterPatch(cellDescriptionsIndex,solverElement);
  }
}

int exahype::solvers::LimitingADERDGSolver::allocateLimiterPatch(
        const int cellDescriptionsIndex,
        const int solverElement) const {
  assertion(solverElement!=exahype::solvers::Solver::NotFound);
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  #if defined(Asserts)
  const int previouslimiterElement =
          tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  #endif
  assertion(previouslimiterElement==exahype::solvers::Solver::NotFound);

  tarch::multicore::Lock lock(exahype::HeapSemaphore);
  exahype::solvers::FiniteVolumesSolver::addNewCellDescription(
      cellDescriptionsIndex,
      solverPatch.getSolverNumber(),
      LimiterPatch::Type::Cell,
      LimiterPatch::RefinementEvent::None,
      solverPatch.getLevel(),
      solverPatch.getParentIndex(),
      solverPatch.getSize(),
      solverPatch.getOffset());
  lock.free();

  assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getPreviousSolution()),solverPatch.toString());
  assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolution()),solverPatch.toString());

  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  assertion(limiterElement!=exahype::solvers::Solver::NotFound);
  LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(
      solverPatch,cellDescriptionsIndex,limiterElement);
  _limiter->ensureNecessaryMemoryIsAllocated(limiterPatch);

  // Copy more geometry information from the solver patch
  limiterPatch.setIsInside(solverPatch.getIsInside());

  return limiterElement;
}

bool exahype::solvers::LimitingADERDGSolver::ensureRequiredLimiterPatchIsAllocated(
        const int cellDescriptionsIndex,
        const int solverElement) const {
  assertion(solverElement!=exahype::solvers::Solver::NotFound);
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      limiterElement==exahype::solvers::Solver::NotFound    &&
      solverPatch.getType()==SolverPatch::Type::Cell        &&
      solverPatch.getLimiterStatus()>0
  ) {
    assertion1(solverPatch.getPreviousLimiterStatus()==0,solverPatch.toString());
    allocateLimiterPatch(cellDescriptionsIndex,solverElement);
    return true;
  }
  return false;
}


void exahype::solvers::LimitingADERDGSolver::projectDGSolutionOnFVSpace(
    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
  const double* solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();
  double*       limiterSolution = DataHeap::getInstance().getData(limiterPatch.getSolution()).data();

  // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
  kernels::limiter::generic::c::projectOnFVLimiterSpace(
      solverSolution,_solver->getNumberOfVariables(),
      _solver->getNodesPerCoordinateAxis(),
      _limiter->getGhostLayerWidth(),
      limiterSolution);
}

// TODO(Dominic): Check that we have rolled back in time as well
void exahype::solvers::LimitingADERDGSolver::rollbackSolverSolutionsGlobally(
    const int cellDescriptionsIndex, const int solverElement) const {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  // 1. Rollback solution to previous time step
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
      assertion(solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None);

      if (solverPatch.getPreviousLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch()) {
        LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
        _limiter->swapSolutionAndPreviousSolution(limiterPatch); // roll back limiter

        if (solverPatch.getPreviousLimiterStatus()<_solver->getMinimumLimiterStatusForTroubledCell()) {
          LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
          projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
        }
      }
      else if (solverPatch.getPreviousLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch()){
        _solver->swapSolutionAndPreviousSolution(solverPatch); // roll back solver

        if (solverPatch.getPreviousLimiterStatus() > 0) {
          LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
          projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
        }
      }
    }
    else { // solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel()
      _solver->swapSolutionAndPreviousSolution(solverPatch);
    }
  }
  // 2. Update the limiter status (do not overwrite the previous limiter status)
  //solverPatch.setLimiterStatus(ADERDGSolver::determineLimiterStatus(solverPatch));
  //solverPatch.setFacewiseLimiterStatus(solverPatch.getLimiterStatus());
}

void exahype::solvers::LimitingADERDGSolver::reinitialiseSolversGlobally(
    const int cellDescriptionsIndex,
        const int solverElement) const {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  // 1. Overwrite the limiter status with the previous one
  solverPatch.setLimiterStatus(solverPatch.getPreviousLimiterStatus());
  assertion1(tarch::la::max(solverPatch.getFacewiseLimiterStatus())==0,solverPatch.toString());

  // 2. Reset the iterationsToCure on all troubled cells to maximum value if cell is troubled
  if (solverPatch.getLimiterStatus()>=_solver->getMinimumLimiterStatusForTroubledCell()) {
    solverPatch.setIterationsToCureTroubledCell(1+_iterationsToCureTroubledCell);
  }

  // 3. Only after the reinitialisation, it is safe to deallocate the limiter patch
  deallocateLimiterPatchOnHelperCell(cellDescriptionsIndex,solverElement);
  ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::reinitialiseSolversLocally(
    const int cellDescriptionsIndex, const int solverElement) const {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  // TODO(Dominic): Add to docu: Assumes that no merge is performed in adapter FinaliseMeshRefinementAndReinitialisation
  // 0. Update the limiter status (do not overwrite the previous limiter status)
  // solverPatch.setLimiterStatus(ADERDGSolver::determineLimiterStatus(solverPatch));
  // solverPatch.setFacewiseLimiterStatus(0);

  // 1. Ensure limiter patch is allocated
  ensureRequiredLimiterPatchIsAllocated(cellDescriptionsIndex,solverElement);

  // 2. Now roll back to the last valid solution
  if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getLimiterStatus()>0) { // this is one of the important differences to the global recomputation where we rollback also cells with limiter status == 0
    assertion(solverPatch.getType()==SolverPatch::Type::Cell);
    assertion(solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None);

    if (solverPatch.getPreviousLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch()) {
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
      _limiter->swapSolutionAndPreviousSolution(limiterPatch);
      projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
    }
    else {
      _solver->swapSolutionAndPreviousSolution(solverPatch);
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
      projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
    }

    // 2.1 Reset the iterationsToCure on all troubled cells to maximum value if cell is troubled
    if (solverPatch.getLimiterStatus()>=_solver->getMinimumLimiterStatusForTroubledCell()) {
      solverPatch.setIterationsToCureTroubledCell(1+_iterationsToCureTroubledCell);
    }
  }

  // 3. Only after the reinitialisation, it is safe to deallocate the limiter patch
  deallocateLimiterPatchOnHelperCell(cellDescriptionsIndex,solverElement);
  ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::projectFVSolutionOnDGSpace(
    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
  const double* limiterSolution = DataHeap::getInstance().getData(limiterPatch.getSolution()).data();
  double*       solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

  // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
  kernels::limiter::generic::c::projectOnDGSpace(
      limiterSolution,_solver->getNumberOfVariables(),
      _solver->getNodesPerCoordinateAxis(),
      _limiter->getGhostLayerWidth(),
      solverSolution);
}

void exahype::solvers::LimitingADERDGSolver::recomputeSolutionLocally(
        const int cellDescriptionsIndex, const int solverElement) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getType()==SolverPatch::Type::Cell &&
      solverPatch.getLimiterStatus()>0
  ) {
    if (solverPatch.getLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch()) {
      const int limiterElement =
          tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(
          solverPatch,cellDescriptionsIndex,limiterElement);

      // 1. Evolve solution to desired  time step again
      _limiter->updateSolution(cellDescriptionsIndex,limiterElement);
      // 2. Project FV solution on ADER-DG space
      projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
    }
    else { // solverPatch.getLimiterStatus()<ADERDGSolver::MinimumLimiterStatusForActiveFVPatch
      if (solverPatch.getPreviousLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch()) {
        _solver->swapSolutionAndPreviousSolution(solverPatch);
        LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
        projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
        // TODO(Dominic): Add to docu: Here, were just went back one time step to supply the NT neighbours with old limiter unknowns.
      }
      else {
        LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
        _limiter->swapSolutionAndPreviousSolution(limiterPatch);
        projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
      }
    }
  } else { // limiterStatus=0 or not on finest level
      #if defined(Asserts)
      const int limiterElement =
          tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
      #endif
      assertion(limiterElement==exahype::solvers::Solver::NotFound);
  }
}

void exahype::solvers::LimitingADERDGSolver::recomputePredictorLocally(
    const int cellDescriptionsIndex,
    const int element,
    exahype::solvers::PredictionTemporaryVariables& predictionTemporaryVariables) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getType()==SolverPatch::Type::Cell &&
      solverPatch.getLimiterStatus() < _solver->getMinimumLimiterStatusForTroubledCell()  &&
      solverPatch.getPredictorTimeStepSize() > 0) {
    assertion(solverPatch.getLimiterStatus()>=0);
    if (
        solverPatch.getLimiterStatus() >= _solver->getMinimumLimiterStatusForActiveFVPatch()
        ||
        (solverPatch.getLimiterStatus() < _solver->getMinimumLimiterStatusForActiveFVPatch() &&
         solverPatch.getPreviousLimiterStatus()>=_solver->getMinimumLimiterStatusForTroubledCell())
    ) {
      _solver->performPredictionAndVolumeIntegral(
          solverPatch,
          predictionTemporaryVariables._tempSpaceTimeUnknowns    [solverPatch.getSolverNumber()],
          predictionTemporaryVariables._tempSpaceTimeFluxUnknowns[solverPatch.getSolverNumber()],
          predictionTemporaryVariables._tempUnknowns             [solverPatch.getSolverNumber()],
          predictionTemporaryVariables._tempFluxUnknowns         [solverPatch.getSolverNumber()],
          predictionTemporaryVariables._tempPointForceSources    [solverPatch.getSolverNumber()]);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::preProcess(
    const int cellDescriptionsIndex,
    const int element) const {
  _solver->preProcess(cellDescriptionsIndex,element);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=Solver::NotFound) {
    _limiter->preProcess(cellDescriptionsIndex,limiterElement);
  }
}

void exahype::solvers::LimitingADERDGSolver::postProcess(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->postProcess(cellDescriptionsIndex,element);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=Solver::NotFound) {
    _limiter->postProcess(cellDescriptionsIndex,limiterElement);
  }
}

void exahype::solvers::LimitingADERDGSolver::prolongateDataAndPrepareDataRestriction(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->prolongateDataAndPrepareDataRestriction(cellDescriptionsIndex,element);

  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);
  if (
      solverPatch.getType()==SolverPatch::Type::Ancestor &&
      solverPatch.getPreviousLimiterStatus()>=_solver->getMinimumLimiterStatusForTroubledCell()) {
    solverPatch.setLimiterStatus(_solver->getMinimumLimiterStatusForTroubledCell()-1);
    solverPatch.setFacewiseLimiterStatus(0);
  }
}

void exahype::solvers::LimitingADERDGSolver::restrictToNextParent(
        const int fineGridCellDescriptionsIndex,
        const int fineGridElement,
        const int coarseGridCellDescriptionsIndex,
        const int coarseGridElement) const {
  _solver->restrictToNextParent(
      fineGridCellDescriptionsIndex,
      fineGridElement,coarseGridCellDescriptionsIndex,coarseGridElement);
}

void exahype::solvers::LimitingADERDGSolver::restrictToTopMostParent(
    const int cellDescriptionsIndex,
    const int solverElement,
    const int parentCellDescriptionsIndex,
    const int parentSolverElement,
    const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  _solver->restrictToTopMostParent(cellDescriptionsIndex,solverElement,parentCellDescriptionsIndex,parentSolverElement,subcellIndex);
}


///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::mergeNeighboursMetadata(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) const {
  _solver->mergeNeighboursMetadata(cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);
}

void exahype::solvers::LimitingADERDGSolver::mergeNeighbours(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    double**                                  tempFaceUnknowns) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  // 1. Solve the riemann problems
  mergeNeighboursBasedOnLimiterStatus(
      cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,
      pos1,pos2,
      false, /* isRecomputation */
      tempFaceUnknowns);

  // 2. Merge the min and max of both cell description's solver's
  // solution value.
  mergeSolutionMinMaxOnFace(
      cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);
}

// TODO(Dominic): Remove limiterstatus1 and limiterStatus2 argument.
// They depend on the isRecomputation value
void exahype::solvers::LimitingADERDGSolver::mergeNeighboursBasedOnLimiterStatus(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    const bool                                isRecomputation,
    double**                                  tempFaceUnknowns) const {
  SolverPatch& solverPatch1 = ADERDGSolver::getCellDescription(cellDescriptionsIndex1,element1);
  SolverPatch& solverPatch2 = ADERDGSolver::getCellDescription(cellDescriptionsIndex2,element2);
  const int limiterElement1 = tryGetLimiterElement(cellDescriptionsIndex1,solverPatch1.getSolverNumber());
  const int limiterElement2 = tryGetLimiterElement(cellDescriptionsIndex2,solverPatch2.getSolverNumber());

  // We only limit on the finest mesh level
  if (solverPatch1.getLevel()==getMaximumAdaptiveMeshLevel()) {
    assertion2(solverPatch1.getLevel()==getMaximumAdaptiveMeshLevel(),solverPatch1.toString(),solverPatch2.toString());
    // 1. Merge solver solution or limiter solution values in
    // non-overlapping parts of solver and limiter domain:
    if (solverPatch1.getLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch() &&
        solverPatch2.getLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch()) {
      if (!isRecomputation) {
        _solver->mergeNeighbours(
            cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
            tempFaceUnknowns); // Which one is left and right is checked internally again.
      }
    }
    else if (solverPatch1.getLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch() &&
             solverPatch2.getLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch()) {
      _limiter->mergeNeighbours(
          cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2,
          tempFaceUnknowns); // Which one is left and right is checked internally again.
    }
    // 2. Merge limiter solution values in overlapping part
    // of solver and limiter domain:
    else if (
        (solverPatch1.getLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch() &&
         solverPatch2.getLimiterStatus() <_solver->getMinimumLimiterStatusForActiveFVPatch())
        ||
        (solverPatch1.getLimiterStatus() <_solver->getMinimumLimiterStatusForActiveFVPatch() &&
         solverPatch2.getLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch())
    ) {
      assertion2(limiterElement1!=exahype::solvers::Solver::NotFound,solverPatch1.toString(),solverPatch2.toString());
      assertion2(limiterElement2!=exahype::solvers::Solver::NotFound,solverPatch2.toString(),solverPatch1.toString());
      _limiter->mergeNeighbours(
          cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2,
          tempFaceUnknowns); // Which one is left and right is checked internally again.
      if (!isRecomputation) {
        _solver->mergeNeighbours(
            cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
            tempFaceUnknowns); // Which one is left and right is checked internally again.
      }
    }
    else {
      logError("mergeNeighboursBasedOnLimiterStatus(...)","Neighbours cannot communicate. " <<
          std::endl << "cell1=" << solverPatch1.toString() <<
          std::endl << ".cell2=" << solverPatch1.toString());
      abort();
    }

  // On the other levels, we work with the ADER-DG solver only
  } else { // solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel()
    if (!isRecomputation) {
      _solver->mergeNeighbours(
          cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
          tempFaceUnknowns); // Which one is left and right is checked internally again.
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) const {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  SolverPatch& solverPatch1 = ADERDGSolver::getCellDescription(cellDescriptionsIndex1,element1);
  SolverPatch& solverPatch2 = ADERDGSolver::getCellDescription(cellDescriptionsIndex2,element2);

  // 1.2. Merge min/max of both solver patches
  const int direction    = tarch::la::equalsReturnIndex(pos1,pos2);
  const int orientation1 = (1 + pos2(direction) - pos1(direction))/2;
  const int orientation2 = 1-orientation1;

  const int faceIndex1 = 2*direction+orientation1;
  const int faceIndex2 = 2*direction+orientation2;

  mergeSolutionMinMaxOnFace(
      solverPatch1,solverPatch2,faceIndex1,faceIndex2);
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch& solverPatch1,
  SolverPatch& solverPatch2,
  const int faceIndex1,
  const int faceIndex2
) const {
  if (
      _solver->getDMPObservables() > 0
      &&
      (solverPatch1.getType()==SolverPatch::Cell ||
      solverPatch2.getType()==SolverPatch::Cell)
  ) {
    assertion( solverPatch1.getSolverNumber() == solverPatch2.getSolverNumber() );
    const int numberOfObservables = _solver->getDMPObservables();
    double* min1 = DataHeap::getInstance().getData( solverPatch1.getSolutionMin()  ).data() + faceIndex1 * numberOfObservables;
    double* min2 = DataHeap::getInstance().getData( solverPatch2.getSolutionMin()  ).data() + faceIndex2 * numberOfObservables;
    double* max1 = DataHeap::getInstance().getData( solverPatch1.getSolutionMax()  ).data() + faceIndex1 * numberOfObservables;
    double* max2 = DataHeap::getInstance().getData( solverPatch2.getSolutionMax()  ).data() + faceIndex2 * numberOfObservables;

    for (int i=0; i<numberOfObservables; i++) {
      const double min = std::min(
          *(min1+i),
          *(min2+i)
      );
      const double max = std::max(
          *(max1+i),
          *(max2+i)
      );

      *(min1+i) = min;
      *(min2+i) = min;

      *(max1+i) = max;
      *(max2+i) = max;
    }
  } // else do nothing
}

void exahype::solvers::LimitingADERDGSolver::mergeWithBoundaryData(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
      double**                                  tempFaceUnknowns) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  mergeWithBoundaryDataBasedOnLimiterStatus(
      cellDescriptionsIndex,element,
      solverPatch.getLimiterStatus(),
      posCell,posBoundary,
      false, /* isRecomputation */
      tempFaceUnknowns);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithBoundaryDataBasedOnLimiterStatus(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const int                                 limiterStatus, //TODO(Dominic): Still necessary?
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
      const bool                                isRecomputation,
      double**                                  tempFaceUnknowns) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());

      if (solverPatch.getLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch()) {
        assertion(limiterStatus==0 || limiterElement!=exahype::solvers::Solver::NotFound);
        if (!isRecomputation) {
          _solver->mergeWithBoundaryData(cellDescriptionsIndex,element,posCell,posBoundary,
              tempFaceUnknowns);
        }
      }
      else {
        assertion(limiterElement!=exahype::solvers::Solver::NotFound);
        _limiter->mergeWithBoundaryData(cellDescriptionsIndex,limiterElement,posCell,posBoundary,
            tempFaceUnknowns);
      }
    }
    else { // solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel()
      if (!isRecomputation) {
        _solver->mergeWithBoundaryData(cellDescriptionsIndex,element,posCell,posBoundary,
            tempFaceUnknowns);
      }
    }
  }
}

#ifdef Parallel
///////////////////////////////////
// NEIGHBOUR - Mesh refinement
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::appendNeighbourCommunicationMetadata(
    exahype::MetadataHeap::HeapEntries& metadata,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  _solver->appendNeighbourCommunicationMetadata(
      metadata,src,dest,
      cellDescriptionsIndex,solverNumber);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourMetadata(
      const exahype::MetadataHeap::HeapEntries& metadata,
      const tarch::la::Vector<DIMENSIONS, int>& src,
      const tarch::la::Vector<DIMENSIONS, int>& dest,
      const int cellDescriptionsIndex,
      const int element) const {
  _solver->mergeWithNeighbourMetadata(metadata,src,dest,cellDescriptionsIndex,element);
}

///////////////////////////////////
// NEIGHBOUR - Time marching
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  sendMinAndMaxToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);

  sendDataToNeighbourBasedOnLimiterStatus(
      toRank,cellDescriptionsIndex,element,src,dest,
      false,/* isRecomputation */
      x,level);

  // send order:   minAndMax,solver,limiter
  // receive order limiter,solver,minAndMax
}

void exahype::solvers::LimitingADERDGSolver::sendMinAndMaxToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);
    if (ADERDGSolver::holdsFaceData(solverPatch)) {
      const int direction   = tarch::la::equalsReturnIndex(src, dest);
      const int orientation = (1 + dest(direction) - src(direction))/2;
      const int faceIndex   = 2*direction+orientation;

      assertion(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMin()));
      assertion(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMax()));
      const double* observablesMin = DataHeap::getInstance().getData(
          solverPatch.getSolutionMin()).data() +
          (faceIndex * numberOfObservables);
      const double* observablesMax = DataHeap::getInstance().getData(
          solverPatch.getSolutionMax()).data() +
          (faceIndex * numberOfObservables);

      DataHeap::getInstance().sendData(
          observablesMin, numberOfObservables, toRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
      DataHeap::getInstance().sendData(
          observablesMax, numberOfObservables, toRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
    } else {
      for (int sends = 0; sends < 2; ++sends) {
        DataHeap::getInstance().sendData(
            EmptyDataHeapMessage, toRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
      }
    }
  }
}


void exahype::solvers::LimitingADERDGSolver::sendDataToNeighbourBasedOnLimiterStatus(
        const int                                    toRank,
        const int                                    cellDescriptionsIndex,
        const int                                    element,
        const tarch::la::Vector<DIMENSIONS, int>&    src,
        const tarch::la::Vector<DIMENSIONS, int>&    dest,
        const bool                                   isRecomputation,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int                                    level) const {
  if (level==getMaximumAdaptiveMeshLevel()) {
    SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);
    logDebug("sendDataToNeighbourBasedOnLimiterStatus(...)", "send data for solver " << _identifier << " from rank " <<
                 toRank << " at vertex x=" << x << ", level=" << level <<
                 ", source=" << src << ", destination=" << dest <<", limiterStatus="<<solverPatch.getLimiterStatus());

    if (solverPatch.getLimiterStatus()==0) {
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
      _limiter->sendEmptyDataToNeighbour(toRank,x,level); // !!! Receive order must be inverted in neighbour comm.
    }
    else if (solverPatch.getLimiterStatus()>0 &&
             solverPatch.getLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch()) {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
      _limiter->sendDataToNeighbour(toRank,cellDescriptionsIndex,limiterElement,src,dest,x,level);
    }
    else if (solverPatch.getLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch() &&
             solverPatch.getLimiterStatus()<_solver->getMinimumLimiterStatusForTroubledCell()) {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
      _limiter->sendDataToNeighbour(toRank,cellDescriptionsIndex,limiterElement,src,dest,x,level);
    }
    else { // solverPatch.getLimiterStatus()>=ADERDGSolver::MinimumLimiterStatusForTroubledCell
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());
      _solver->sendEmptyDataToNeighbour(toRank,x,level);
      _limiter->sendDataToNeighbour(toRank,cellDescriptionsIndex,limiterElement,src,dest,x,level);
    }
  } else {
    _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  // send an empty minAndMax message
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables) {
    for(int sends=0; sends<2; ++sends)
      DataHeap::getInstance().sendData(
          exahype::EmptyDataHeapMessage, toRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
  }

  _solver->sendEmptyDataToNeighbour(toRank,x,level);
  if (level==getMaximumAdaptiveMeshLevel()) {
    _limiter->sendEmptyDataToNeighbour(toRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourData(
    const int                                    fromRank,
    const exahype::MetadataHeap::HeapEntries&    neighbourMetadata,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    double**                                     tempFaceUnknowns,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  logDebug("mergeWithNeighbourDataBasedOnLimiterStatus(...)", "receive data for solver " << _identifier << " from rank " <<
          fromRank << " at vertex x=" << x << ", level=" << level <<
          ", source=" << src << ", destination=" << dest);

  mergeWithNeighbourDataBasedOnLimiterStatus(
      fromRank,neighbourMetadata,cellDescriptionsIndex,element,src,dest,
      false,/*isRecomputation*/
      tempFaceUnknowns,x,level);

  mergeWithNeighbourMinAndMax(fromRank,cellDescriptionsIndex,element,src,dest,x,level);

  // send order:   minAndMax,solver,limiter
  // receive order limiter,solver,minAndMax
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourDataBasedOnLimiterStatus(
    const int                                    fromRank,
    const exahype::MetadataHeap::HeapEntries&    neighbourMetadata,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const bool                                   isRecomputation,
    double**                                     tempFaceUnknowns,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  if (level==getMaximumAdaptiveMeshLevel()) {
    SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);
    logDebug("mergeWithNeighbourDataBasedOnLimiterStatus(...)", "receive data for solver " << _identifier << " from rank " <<
        fromRank << " at vertex x=" << x << ", level=" << level <<
        ", source=" << src << ", destination=" << dest << ",limiterStatus=" << solverPatch.getLimiterStatus());
    assertion1(solverPatch.getLimiterStatus()>=0,solverPatch.toString());
    if (solverPatch.getLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch()) {
      _limiter->dropNeighbourData(fromRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
      if (!isRecomputation) {
        _solver->mergeWithNeighbourData(
            fromRank,neighbourMetadata,cellDescriptionsIndex,element,
            src,dest,tempFaceUnknowns,x,level);
      } else {
        _solver->dropNeighbourData(fromRank,src,dest,x,level);
      }
    }
    else { // solverPatch.getLimiterStatus()>=ADERDGSolver::MinimumLimiterStatusForActiveFVPatch) {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _limiter->mergeWithNeighbourData(
          fromRank,neighbourMetadata,cellDescriptionsIndex,limiterElement,
          src,dest,tempFaceUnknowns,x,level);
      _solver->dropNeighbourData(fromRank,src,dest,x,level);
    }
  } else {
    logDebug("mergeWithNeighbourDataBasedOnLimiterStatus(...)", "receive data for solver " << _identifier << " from rank " <<
        fromRank << " at vertex x=" << x << ", level=" << level <<
        ", source=" << src << ", destination=" << dest);

    _solver->mergeWithNeighbourData(
        fromRank,neighbourMetadata,cellDescriptionsIndex,element,
        src,dest,tempFaceUnknowns,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourMinAndMax(
    const int                                    fromRank,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

    if (ADERDGSolver::holdsFaceData(solverPatch)) {
      const int direction   = tarch::la::equalsReturnIndex(src, dest);
      const int orientation = (1 + src(direction) - dest(direction))/2;
      const int faceIndex   = 2*direction+orientation;

      const int receivedMaxIndex = DataHeap::getInstance().createData(numberOfObservables, numberOfObservables);
      const int receivedMinIndex = DataHeap::getInstance().createData(numberOfObservables, numberOfObservables);
      DataHeap::HeapEntries& receivedMax = DataHeap::getInstance().getData(receivedMaxIndex);
      DataHeap::HeapEntries& receivedMin = DataHeap::getInstance().getData(receivedMinIndex);
      assertionEquals(DataHeap::getInstance().getData(receivedMaxIndex).size(),static_cast<size_t>(numberOfObservables));
      assertionEquals(DataHeap::getInstance().getData(receivedMinIndex).size(),static_cast<size_t>(numberOfObservables));

      // Inverted send-receive order: TODO(Dominic): Add to docu
      // Send order:    min,max
      // Receive order; max,min
      DataHeap::getInstance().receiveData(
          receivedMax.data(), numberOfObservables, fromRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
      DataHeap::getInstance().receiveData(
          receivedMin.data(), numberOfObservables, fromRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);

      mergeSolutionMinMaxOnFace(
          solverPatch,faceIndex,receivedMin.data(),receivedMax.data());

      DataHeap::getInstance().deleteData(receivedMinIndex,true);
      DataHeap::getInstance().deleteData(receivedMaxIndex,true);
    } else {
      for(int receives=0; receives<2; ++receives)
        DataHeap::getInstance().receiveData(
            fromRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch&  solverPatch,
  const int           faceIndex,
  const double* const min, const double* const max) const {
  assertion1(ADERDGSolver::holdsFaceData(solverPatch),solverPatch.toString());

  double* solutionMin = DataHeap::getInstance().getData( solverPatch.getSolutionMin()  ).data();
  double* solutionMax = DataHeap::getInstance().getData( solverPatch.getSolutionMax()  ).data();

  const int numberOfObservables = _solver->getDMPObservables();
  for (int i=0; i<numberOfObservables; i++) {
    solutionMin[i+faceIndex*numberOfObservables]  = std::min( solutionMin[i+faceIndex*numberOfObservables], min[i] );
    solutionMax[i+faceIndex*numberOfObservables]  = std::max( solutionMax[i+faceIndex*numberOfObservables], max[i] );
  }
}

void exahype::solvers::LimitingADERDGSolver::dropNeighbourData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  // send order:   minAndMax,solver,limiter
  // receive order limiter,solver,minAndMax
  if (level==getMaximumAdaptiveMeshLevel()) {
    _limiter->dropNeighbourData(fromRank,src,dest,x,level);
  }
  _solver->dropNeighbourData(fromRank,src,dest,x,level);

  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    for(int receives=0; receives<2; ++receives)
      DataHeap::getInstance().receiveData(
          fromRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
  }
}

///////////////////////////////////////////////////////////
// NEIGHBOUR - Solution recomputation
///////////////////////////////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendEmptySolverAndLimiterDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendEmptyDataToNeighbour(toRank,x,level);
  if (level==getMaximumAdaptiveMeshLevel()) {
    _limiter->sendEmptyDataToNeighbour(toRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::dropNeighbourSolverAndLimiterData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  logDebug("dropNeighbourSolverAndLimiterData(...)", "drop data for solver " << _identifier << " from rank " <<
            fromRank << " at vertex x=" << x << ", level=" << level <<
            ", src=" << src << ", dest=" << dest);

  _limiter->dropNeighbourData(fromRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
  if (level==getMaximumAdaptiveMeshLevel()) {
    _solver->dropNeighbourData(fromRank,src,dest,x,level);
  }
}

/////////////////////////////////////
// MASTER<=>WORKER
/////////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::prepareMasterCellDescriptionAtMasterWorkerBoundary(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->prepareMasterCellDescriptionAtMasterWorkerBoundary(cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::prepareWorkerCellDescriptionAtMasterWorkerBoundary(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->prepareWorkerCellDescriptionAtMasterWorkerBoundary(cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::appendMasterWorkerCommunicationMetadata(
    exahype::MetadataHeap::HeapEntries& metadata,
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  _solver->appendMasterWorkerCommunicationMetadata(
      metadata,cellDescriptionsIndex,solverNumber);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithMasterMetadata(
    const exahype::MetadataHeap::HeapEntries& metadata,
    const int                                 cellDescriptionsIndex,
    const int                                 element) {
  _solver->mergeWithMasterMetadata(
      metadata,cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerMetadata(
    const exahype::MetadataHeap::HeapEntries& receivedMetadata,
    const int                                 cellDescriptionsIndex,
    const int                                 element) {
  _solver->mergeWithWorkerMetadata(
      receivedMetadata,cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::sendDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendDataToWorkerOrMasterDueToForkOrJoin(
      toRank,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->sendDataToWorkerOrMasterDueToForkOrJoin(
          toRank,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
          toRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
        toRank,x,level);
  _limiter->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
        toRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->mergeWithWorkerOrMasterDataDueToForkOrJoin(
      fromRank,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->mergeWithWorkerOrMasterDataDueToForkOrJoin(
        fromRank,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->dropWorkerOrMasterDataDueToForkOrJoin(
        fromRank,x,level);
  } // !!! Receive order must be the same in master<->worker comm.
}

void exahype::solvers::LimitingADERDGSolver::dropWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->dropWorkerOrMasterDataDueToForkOrJoin(
          fromRank,x,level);
  _limiter->dropWorkerOrMasterDataDueToForkOrJoin(
          fromRank,x,level);
}

///////////////////////////////////
// WORKER->MASTER
///////////////////////////////////

void exahype::solvers::LimitingADERDGSolver::sendDataToMaster(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  DataHeap::HeapEntries messageForMaster =
      _solver->compileMessageForMaster(5);

  // Send additional data to master
  messageForMaster.push_back(
      exahype::solvers::convertToDouble(_limiterDomainChange));

  assertion1(messageForMaster.size()==5,messageForMaster.size());
  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToMaster(...)","Sending time step data: " <<
        "data[0]=" << messageForMaster[0] <<
        ",data[1]=" << messageForMaster[1] <<
        ",data[2]=" << messageForMaster[2] <<
        ",data[3]=" << messageForMaster[3] <<
        ",data[4]=" << messageForMaster[4]);
  }
  DataHeap::getInstance().sendData(
      messageForMaster.data(), messageForMaster.size(),
      masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  DataHeap::HeapEntries messageFromWorker(5); // !!! Creates and fills the vector
  DataHeap::getInstance().receiveData(
      messageFromWorker.data(),messageFromWorker.size(),workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);

  // pass message to the ADER-DG solver
  _solver->mergeWithWorkerData(messageFromWorker);

  // merge own flags
  const int firstEntry=4;
  LimiterDomainChange workerLimiterDomainChange =
      exahype::solvers::convertToLimiterDomainChange(messageFromWorker[firstEntry]);
  updateNextLimiterDomainChange(workerLimiterDomainChange); // !!! It is important that we merge with the "next" field here

  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","Received data from worker:" <<
             " messageFromWorker[0]=" << messageFromWorker[0] <<
             " messageFromWorker[1]=" << messageFromWorker[1] <<
             " messageFromWorker[2]=" << messageFromWorker[2] <<
             " messageFromWorker[3]=" << messageFromWorker[3] <<
             " messageFromWorker[4]=" << messageFromWorker[4]);
    logDebug("mergeWithWorkerData(...)","nextLimiterDomainChange=" << static_cast<int>(_nextLimiterDomainChange));
  }
}

bool exahype::solvers::LimitingADERDGSolver::hasToSendDataToMaster(
      const int cellDescriptionsIndex,
      const int element) const {
  #if defined(Asserts) || defined(Debug)
  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    assertion(_limiter->hasToSendDataToMaster(cellDescriptionsIndex,limiterElement));
  }
  #endif

  return _solver->hasToSendDataToMaster(cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::sendDataToMaster(
    const int                                     masterRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendDataToMaster(
      masterRank,cellDescriptionsIndex,element,x,level);
  // limiter is only active on the finest mesh level
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToMaster(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendEmptyDataToMaster(masterRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const exahype::MetadataHeap::HeapEntries&    workerMetadata,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->mergeWithWorkerData(
      workerRank,workerMetadata,cellDescriptionsIndex,element,x,level);
}

void exahype::solvers::LimitingADERDGSolver::dropWorkerData(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->dropWorkerData(workerRank,x,level);
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToWorker(
    const                                        int workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  DataHeap::HeapEntries messageForWorker = _solver->compileMessageForWorker(9);

  // append additional data
  messageForWorker.push_back(
      exahype::solvers::convertToDouble(_limiterDomainChange));

  assertion1(messageForWorker.size()==9,messageForWorker.size());
  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToWorker(...)","sending data to worker: " <<
             "data[0]=" << messageForWorker[0] <<
             ",data[1]=" << messageForWorker[1] <<
             ",data[2]=" << messageForWorker[2] <<
             ",data[3]=" << messageForWorker[3] <<
             ",data[4]=" << messageForWorker[4] <<
             ",data[5]=" << messageForWorker[5] <<
             ",data[6]=" << messageForWorker[6] <<
             ",data[7]=" << messageForWorker[7] <<
             ",data[8]=" << messageForWorker[8]);
  }
  DataHeap::getInstance().sendData(
      messageForWorker.data(), messageForWorker.size(),
      workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithMasterData(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  DataHeap::HeapEntries messageFromMaster(9); // !!! Creates and fills the vector
  DataHeap::getInstance().receiveData(
      messageFromMaster.data(),messageFromMaster.size(),masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);

  // pass message to the ADER-DG solver
  _solver->mergeWithMasterData(messageFromMaster);

  // merge own data
  const int firstEntry=8;
  LimiterDomainChange workerLimiterDomainChange =
      exahype::solvers::convertToLimiterDomainChange(messageFromMaster[firstEntry]);
  updateNextLimiterDomainChange(workerLimiterDomainChange); // !!! It is important that we merge with the "next" field here

  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","received data from master:" <<
             " messageFromMaster[0]=" << messageFromMaster[0] <<
             " messageFromMaster[1]=" << messageFromMaster[1] <<
             " messageFromMaster[2]=" << messageFromMaster[2] <<
             " messageFromMaster[3]=" << messageFromMaster[3] <<
             " messageFromMaster[4]=" << messageFromMaster[4] <<
             " messageFromMaster[5]=" << messageFromMaster[5] <<
             " messageFromMaster[6]=" << messageFromMaster[6] <<
             " messageFromMaster[7]=" << messageFromMaster[7] <<
             " messageFromMaster[8]=" << messageFromMaster[8]);
    logDebug("mergeWithWorkerData(...)","nextLimiterDomainChange=" << static_cast<int>(_nextLimiterDomainChange));
  }
}

void exahype::solvers::LimitingADERDGSolver::sendDataToWorker(
    const int                                     workerRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendDataToWorker(
      workerRank,cellDescriptionsIndex,element,x,level);
  // limiter is only active on the finest mesh level
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToWorker(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendEmptyDataToWorker(workerRank,x,level);
  // limiter is only active on the finest mesh level
}

void exahype::solvers::LimitingADERDGSolver::mergeWithMasterData(
    const int                                     masterRank,
    const exahype::MetadataHeap::HeapEntries&     masterMetadata,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->mergeWithMasterData(
      masterRank,masterMetadata,cellDescriptionsIndex,element,x,level);
}

void exahype::solvers::LimitingADERDGSolver::dropMasterData(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->dropMasterData(masterRank,x,level);
}
#endif

std::string exahype::solvers::LimitingADERDGSolver::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::solvers::LimitingADERDGSolver::toString (std::ostream& out) const {
  out << getIdentifier() << "{_ADERDG: ";
  out << _solver->toString() << "}\n";
  out << getIdentifier() << "{_FV: ";
  out << _limiter->toString() << "}";
}
