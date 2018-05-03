// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
// ========================
//   www.exahype.eu
// ========================
#ifndef POSTPROCESSING_ConservedWriter_CLASS_HEADER_
#define POSTPROCESSING_ConservedWriter_CLASS_HEADER_

#include "exahype/plotters/Plotter.h"
#include "exahype/solvers/LimitingADERDGSolver.h"


namespace SWE {
  
  class ConservedWriter;
}

class SWE::ConservedWriter : public exahype::plotters::Plotter::UserOnTheFlyPostProcessing {
public:
  ConservedWriter(exahype::solvers::LimitingADERDGSolver&  solver);
  virtual ~ConservedWriter();

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

#endif /* POSTPROCESSING_ConservedWriter_CLASS_HEADER_ */