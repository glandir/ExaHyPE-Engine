// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "PrimitiveWriter.h"
#include "C2P-DIM.h"

DIM::PrimitiveWriter::PrimitiveWriter(exahype::solvers::LimitingADERDGSolver&  solver) {
  // @TODO Please insert your code here.
}

DIM::PrimitiveWriter::~PrimitiveWriter() {
}

void DIM::PrimitiveWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void DIM::PrimitiveWriter::finishPlotting() {
  // @TODO Please insert your code here.
}

void DIM::PrimitiveWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
/*  const int writtenUnknowns = 14;
  for (int i=0; i<writtenUnknowns; i++){ 
    outputQuantities[i] = Q[i];
  }
  */
  pdecons2prim_(outputQuantities, Q);
}