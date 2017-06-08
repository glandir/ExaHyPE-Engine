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

#include "tarch/parallel/NodePool.h"

#include <limits>

bool exahype::State::FuseADERDGPhases           = false;
double exahype::State::WeightForPredictionRerun = 0.9;

bool exahype::State::VirtuallyExpandBoundingBox = false;

bool exahype::State::EnableMasterWorkerCommunication = true;
bool exahype::State::EnableNeighbourCommunication    = true;

exahype::State::State() : Base() {
  // @todo Guidebook

  _stateData.setMaxRefinementLevelAllowed(0);
}

exahype::State::State(const Base::PersistentState& argument) : Base(argument) {
  // do nothing
}

void exahype::State::merge(const exahype::State& anotherState) {
  // do nothing
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

  if (nodePoolHasGivenOutRankSizeLastQuery && _stateData.getMaxRefinementLevelAllowed()>=3) {
    _stateData.setMaxRefinementLevelAllowed( _stateData.getMaxRefinementLevelAllowed()-3);
  }
  // We already thought that we are done but then suddenly the load balancing
  // decided to balance once more. So reset counter.
  else if (nodePoolHasGivenOutRankSizeLastQuery && _stateData.getMaxRefinementLevelAllowed()<0) {
    _stateData.setMaxRefinementLevelAllowed(0);
  }
  // No more nodes left. Start to enforce refinement
  else if (!idleNodesLeft && _stateData.getMaxRefinementLevelAllowed()>=0) {
    _stateData.setMaxRefinementLevelAllowed(-1);
  }
  // Refinement is enforced. So we decrease counter. Once we underrun -2, grid
  // construction can terminate as all enforced refined went through.
  else if (_stateData.getMaxRefinementLevelAllowed()<=-1) {
    _stateData.setMaxRefinementLevelAllowed( _stateData.getMaxRefinementLevelAllowed()-1 );
  }
  // Seems that max permitted level has exceeded max grid level. We may assume
  // that there are more MPI ranks than available trees.
  else if (isGridStationary() && _stateData.getMaxRefinementLevelAllowed()>finestGridLevelPossible) {
    _stateData.setMaxRefinementLevelAllowed( -1 );
  }
  // Nothing has changed in this grid iteration in the grid and we haven't
  // given out new workers. So increase the permitted maximum grid level by
  // one and give another try whether the grid adds more vertices.
  else if ( !(nodePoolHasGivenOutRankSizeLastQuery) && isGridStationary() && _stateData.getMaxRefinementLevelAllowed()>=0 ) {
    _stateData.setMaxRefinementLevelAllowed( _stateData.getMaxRefinementLevelAllowed()+1);
  }
}


exahype::State::RefinementAnswer exahype::State::mayRefine(
  bool isCreationalEvent, int level) const
{
  #ifdef Parallel
  if (
    _stateData.getMaxRefinementLevelAllowed()<=-2
    &&
    isCreationalEvent
  ) {
    // The simple refine works, but we want enforceRefine to work
    //return RefinementAnswer::Refine;
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
  if (isCreationalEvent) {
    return RefinementAnswer::Refine;
  }
  else {
    return RefinementAnswer::DontRefineYet;
  }
  #endif
}


bool exahype::State::continueToConstructGrid() const {
  #ifdef Parallel
  return _stateData.getMaxRefinementLevelAllowed()>=-3;
  #else
  return !isGridBalanced();
  #endif
}

void exahype::State::setAlgorithmSection(const records::State::AlgorithmSection& algorithmSection) {
  _stateData.setAlgorithmSection(algorithmSection);
}

exahype::records::State::AlgorithmSection exahype::State::getAlgorithmSection() const {
  return _stateData.getAlgorithmSection();
}

 exahype::records::State::MergeMode exahype::State::getMergeMode() const {
   return _stateData.getMergeMode();
 }

 exahype::records::State::SendMode exahype::State::getSendMode() const {
   return _stateData.getSendMode();
 }

 void exahype::State::switchToInitialConditionAndTimeStepSizeComputationContext() {
   _stateData.setMergeMode(records::State::MergeMode::BroadcastAndMergeTimeStepData);
   _stateData.setSendMode (records::State::SendMode::ReduceAndMergeTimeStepData);
 }

 void exahype::State::switchToPredictionAndFusedTimeSteppingInitialisationContext() {
   _stateData.setMergeMode(records::State::MergeMode::BroadcastAndMergeTimeStepData);
   _stateData.setSendMode (records::State::SendMode::SendFaceData);
 }

 void exahype::State::switchToADERDGTimeStepContext() {
   if (EnableMasterWorkerCommunication) {
     if (EnableNeighbourCommunication) {
       _stateData.setMergeMode(records::State::MergeMode::BroadcastAndMergeTimeStepDataAndMergeFaceData);
       _stateData.setSendMode (records::State::SendMode::ReduceAndMergeTimeStepDataAndSendFaceData);
     } else { //!EnableNeighbourCommunication
       _stateData.setMergeMode(records::State::MergeMode::BroadcastAndMergeTimeStepData);
       _stateData.setSendMode (records::State::SendMode::ReduceAndMergeTimeStepData);
     }
   } else { // !EnableMasterWorkerCommunication
     if (EnableNeighbourCommunication) {
       _stateData.setMergeMode(records::State::MergeMode::MergeFaceData);
       _stateData.setSendMode (records::State::SendMode::SendFaceData);
     } else { //!EnableNeighbourCommunication
       _stateData.setMergeMode(records::State::MergeMode::MergeNothing);
       _stateData.setSendMode (records::State::SendMode::SendNothing);
     }
   }
 }

 void exahype::State::switchToPredictionRerunContext() {
   _stateData.setMergeMode(records::State::MergeMode::BroadcastAndMergeTimeStepDataAndDropFaceData);
   _stateData.setSendMode (records::State::SendMode::SendFaceData);
 }

 void exahype::State::switchToNeighbourDataMergingContext() {
   _stateData.setMergeMode(records::State::MergeMode::MergeFaceData);
   _stateData.setSendMode (records::State::SendMode::SendNothing);
 }

 void exahype::State::switchToPredictionContext() {
   _stateData.setMergeMode(records::State::MergeMode::BroadcastAndMergeTimeStepData);
   _stateData.setSendMode (records::State::SendMode::SendFaceData);
 }

 void exahype::State::switchToSolutionUpdateContext() {
   _stateData.setMergeMode(records::State::MergeMode::MergeNothing);
   _stateData.setSendMode (records::State::SendMode::SendNothing);
 }

 void exahype::State::switchToTimeStepSizeComputationContext() {
   _stateData.setMergeMode(records::State::MergeMode::MergeNothing);
   _stateData.setSendMode (records::State::SendMode::ReduceAndMergeTimeStepData);
 }

 void exahype::State::switchToUpdateMeshContext() {
   _stateData.setMergeMode(records::State::MergeMode::BroadcastAndMergeTimeStepData);
   _stateData.setSendMode (records::State::SendMode::SendNothing);
 }

 void exahype::State::switchToPostAMRContext() {
   _stateData.setMergeMode(records::State::MergeMode::MergeNothing);
   _stateData.setSendMode (records::State::SendMode::ReduceAndMergeTimeStepData);
 }

 void exahype::State::switchToLimiterStatusSpreadingContext() {
   _stateData.setMergeMode(records::State::MergeMode::BroadcastAndMergeTimeStepData);
   _stateData.setSendMode (records::State::SendMode::SendNothing);
 }

 void exahype::State::switchToLimiterStatusSpreadingFusedTimeSteppingContext() {
   _stateData.setMergeMode(records::State::MergeMode::BroadcastAndMergeTimeStepDataAndDropFaceData);
   _stateData.setSendMode (records::State::SendMode::SendNothing);
 }

 void exahype::State::switchToReinitialisationContext() {
   // We are merging a limiter status but we do not use the merging and sending mappings. So, we can use any value here.
   _stateData.setMergeMode(records::State::MergeMode::MergeNothing);
   _stateData.setSendMode (records::State::SendMode::SendNothing);
 }

 void exahype::State::switchToRecomputeSolutionAndTimeStepSizeComputationContext() {
   _stateData.setMergeMode(records::State::MergeMode::MergeNothing);
   _stateData.setSendMode (records::State::SendMode::ReduceAndMergeTimeStepData);
 }

 void exahype::State::switchToLocalRecomputationAndTimeStepSizeComputationFusedTimeSteppingContext() {
   _stateData.setMergeMode(records::State::MergeMode::MergeNothing);
   _stateData.setSendMode (records::State::SendMode::ReduceAndMergeTimeStepDataAndSendFaceData);
 }

 void exahype::State::switchToNeighbourDataDroppingContext() {
   _stateData.setMergeMode(records::State::MergeMode::DropFaceData);
   _stateData.setSendMode (records::State::SendMode::SendNothing);
 }

 bool exahype::State::fuseADERDGPhases() {
   return FuseADERDGPhases;
 }

 double exahype::State::getTimeStepSizeWeightForPredictionRerun() {
   return WeightForPredictionRerun;
 }
