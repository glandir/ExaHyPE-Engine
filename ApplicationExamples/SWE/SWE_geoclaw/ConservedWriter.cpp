// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ConservedWriter.h"

extern double grav;

SWE::ConservedWriter::ConservedWriter(SWE::MySWESolver& solver) {
  // @TODO Please insert your code here.
}

SWE::ConservedWriter::~ConservedWriter() {}

void SWE::ConservedWriter::startPlotting(double time) {
  // @TODO Please insert your code here.
}

void SWE::ConservedWriter::finishPlotting() {
  // @TODO Please insert your code here.
}

void SWE::ConservedWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>& pos, double* const Q,
    double* const outputQuantities, double timeStamp) {
  const int writtenUnknowns = 6;
  for (int i = 0; i < writtenUnknowns; i++) {
    outputQuantities[i] = Q[i];
  }
  outputQuantities[4] = Q[3] + Q[0];

  // calculate error for OscillationLake
  double omega = sqrt(0.2 * grav);
  outputQuantities[5] = std::abs(
      Q[0] - std::max(0.0, 0.05 * (2 * x[0] * std::cos(omega * timeStamp) +
                                   2 * x[1] * std::sin(omega * timeStamp)) +
                               0.075 - Q[3]));
}