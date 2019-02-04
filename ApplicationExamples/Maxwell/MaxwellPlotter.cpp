// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "MaxwellPlotter.h"

Maxwell::MaxwellPlotter::MaxwellPlotter(Maxwell::MaxwellSolver& solver) {
  // @TODO Please insert your code here.
}

Maxwell::MaxwellPlotter::~MaxwellPlotter() {
}

void Maxwell::MaxwellPlotter::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void Maxwell::MaxwellPlotter::finishPlotting() {
  // @TODO Please insert your code here.
}

void Maxwell::MaxwellPlotter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  const int writtenUnknowns = 8;
  for (int i=0; i<writtenUnknowns; i++){ 
    outputQuantities[i] = Q[i];
  }
}