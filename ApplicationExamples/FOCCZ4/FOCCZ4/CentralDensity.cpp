// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "CentralDensity.h"
#include "PDE.h"

FOCCZ4::CentralDensity::CentralDensity(FOCCZ4::FOCCZ4Solver& solver) {
  // @TODO Please insert your code here.
}

FOCCZ4::CentralDensity::~CentralDensity() {
}

void FOCCZ4::CentralDensity::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void FOCCZ4::CentralDensity::finishPlotting() {
  // @TODO Please insert your code here.
}

void FOCCZ4::CentralDensity::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  double V[96];
  pdecons2prim_(V,Q);
  outputQuantities[0] = V[59];
}
