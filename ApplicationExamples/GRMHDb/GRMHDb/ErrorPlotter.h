// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
// ========================
//   www.exahype.eu
// ========================
#ifndef POSTPROCESSING_ErrorPlotter_CLASS_HEADER_
#define POSTPROCESSING_ErrorPlotter_CLASS_HEADER_

#include "exahype/plotters/Plotter.h"

namespace GRMHDb {
  class GRMHDbSolver_ADERDG;
  class ErrorPlotter;
}

class GRMHDb::ErrorPlotter : public exahype::plotters::Plotter::UserOnTheFlyPostProcessing {
public:
  ErrorPlotter(GRMHDb::GRMHDbSolver_ADERDG& solver);
  virtual ~ErrorPlotter();

  void startPlotting(double time) override;
  void finishPlotting() override;
  void mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp) override;
};

#endif /* POSTPROCESSING_ErrorPlotter_CLASS_HEADER_ */