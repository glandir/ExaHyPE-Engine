// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "PrimitiveWriterSubcellVTK.h"

#include "PDE.h"

GRMHDb::PrimitiveWriterSubcellVTK::PrimitiveWriterSubcellVTK(GRMHDb::GRMHDbSolver& solver) {
  // @TODO Please insert your code here.
}

GRMHDb::PrimitiveWriterSubcellVTK::~PrimitiveWriterSubcellVTK() {
}

void GRMHDb::PrimitiveWriterSubcellVTK::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void GRMHDb::PrimitiveWriterSubcellVTK::finishPlotting() {
  // @TODO Please insert your code here.
}

void GRMHDb::PrimitiveWriterSubcellVTK::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
	int iErr = 0;
	pdecons2prim_(outputQuantities, Q, &iErr);
  //const int writtenUnknowns = 19;
  //for (int i=0; i<writtenUnknowns; i++){ 
  //  outputQuantities[i] = Q[i];
  //}
}