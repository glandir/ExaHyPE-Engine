// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ProbeWriter.h"


SWE::ProbeWriter::ProbeWriter(MySWESolver_FV&  solver) {
  // @TODO Please insert your code here.
}

SWE::ProbeWriter::~ProbeWriter() {
}

void SWE::ProbeWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void SWE::ProbeWriter::finishPlotting() {
  // @TODO Please insert your code here.
}

void SWE::ProbeWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  const int writtenUnknowns = 6;
  for (int i=0; i<writtenUnknowns - 2; i++){
    outputQuantities[i] = Q[i];
  }
  outputQuantities[4] = Q[3] + Q[0];

  //calculate error for OscillationLake
  //  outputQuantities[5] = std::abs(Q[0] - std::max(0.0, 0.05 * (2 * x[0] * std::cos(omega * timeStamp) + 2 * x[1] * std::sin(omega * timeStamp)) + 0.075 - Q[3]));
  if(Q[3] > 0.0){
    outputQuantities[5] = 0.0;
  }else{
    outputQuantities[5] = outputQuantities[4];
  }
}
