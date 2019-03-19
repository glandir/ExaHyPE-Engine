// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "TimingStatistics_FV.h"

#include "TimingStatisticsWriter.cpph"
#include "AbstractGRMHDSolver_FV.h"

GRMHD::TimingStatistics_FV::TimingStatistics_FV() : exahype::plotters::FiniteVolumes2UserDefined::FiniteVolumes2UserDefined(){
  const int basisSize = GRMHD::AbstractGRMHDSolver_FV::PatchSize;
  const int basisSizeD =  basisSize * basisSize * (DIMENSIONS == 3 ? basisSize : 1);
  const int numberOfVariables = GRMHD::AbstractGRMHDSolver_FV::NumberOfVariables;
  const int unknownsPerPatch = basisSizeD * numberOfVariables;
  
  writer = new TimingStatisticsWriter("output/timing-statistics.asc", unknownsPerPatch);
}


void GRMHD::TimingStatistics_FV::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
    double timeStamp) {
  writer->numberOfPatches++;
}


void GRMHD::TimingStatistics_FV::startPlotting( double time) {
  writer->newTimeStep(time);
  writer->numberOfPatches = 0;
}


void GRMHD::TimingStatistics_FV::finishPlotting() {
  // @TODO Please insert your code here.
}


