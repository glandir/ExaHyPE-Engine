// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ProbeWriter14.h"

SWE::ProbeWriter14::ProbeWriter14(SWE::MySWESolver& solver) {
  // @TODO Please insert your code here.
}

SWE::ProbeWriter14::~ProbeWriter14() {
}

void SWE::ProbeWriter14::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void SWE::ProbeWriter14::finishPlotting() {
  // @TODO Please insert your code here.
}

void SWE::ProbeWriter14::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
        outputQuantities[0] = Q[3];
    outputQuantities[1] = Q[3] + Q[4];
}
