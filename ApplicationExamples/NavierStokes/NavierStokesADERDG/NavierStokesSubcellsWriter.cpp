// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "NavierStokesSubcellsWriter.h"


NavierStokesADERDG::NavierStokesSubcellsWriter::NavierStokesSubcellsWriter(exahype::solvers::LimitingADERDGSolver&  solver) {
  // @TODO Please insert your code here.
}

NavierStokesADERDG::NavierStokesSubcellsWriter::~NavierStokesSubcellsWriter() {
}

void NavierStokesADERDG::NavierStokesSubcellsWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void NavierStokesADERDG::NavierStokesSubcellsWriter::finishPlotting() {
  // @TODO Please insert your code here.
}

void NavierStokesADERDG::NavierStokesSubcellsWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  const int writtenUnknowns = 5;
  for (int i=0; i<writtenUnknowns; i++){ 
    outputQuantities[i] = Q[i];
  }
}