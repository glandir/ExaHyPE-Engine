// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "PrimitiveWriter.h"
#include "C2P-GPR.h"

GPR::PrimitiveWriter::PrimitiveWriter(exahype::solvers::LimitingADERDGSolver&  solver) {
  // @TODO Please insert your code here.
}

GPR::PrimitiveWriter::~PrimitiveWriter() {
}

void GPR::PrimitiveWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void GPR::PrimitiveWriter::finishPlotting() {
  // @TODO Please insert your code here.
}

void GPR::PrimitiveWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  pdecons2prim_(outputQuantities, Q);
}