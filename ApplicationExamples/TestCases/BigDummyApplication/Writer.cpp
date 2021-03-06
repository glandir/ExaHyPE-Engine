// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "Writer.h"
#include "tarch/logging/Log.h"


Dummy::Writer::Writer(Dummy::DummySolver& solver) {
  // @TODO Please insert your code here.
}

Dummy::Writer::~Writer() {
}

void Dummy::Writer::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void Dummy::Writer::finishPlotting() {
  // @TODO Please insert your code here.
}

void Dummy::Writer::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  const int writtenUnknowns = 60;
  for (int i=0; i<writtenUnknowns; i++){ 
    outputQuantities[i] = Q[i];
    // Check for consistency
    const double eps = 1e-8;
    if(Q[i] > 1 + eps || Q[i] < 1 - eps) {
        tarch::logging::Log _log("Dummy::Writer");
        logError("mapQuantities()", "Consistency check failed, Q["<<i<<"]="<<Q[i]<<" which is too far away from the ID, despite the PDE Q_t=0 is solved");
    }
  }
}
