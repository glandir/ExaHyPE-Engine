// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "PrimitiveWriter.h"

GRMHDb::PrimitiveWriter::PrimitiveWriter(GRMHDb::GRMHDbSolver& solver) {
  // @TODO Please insert your code here.
}

GRMHDb::PrimitiveWriter::~PrimitiveWriter() {
}

void GRMHDb::PrimitiveWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void GRMHDb::PrimitiveWriter::finishPlotting() {
  // @TODO Please insert your code here.
}

void GRMHDb::PrimitiveWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  const int writtenUnknowns = 19;
  for (int i=0; i<writtenUnknowns; i++){ 
    outputQuantities[i] = Q[i];
  }
}