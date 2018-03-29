// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "exahype/plotters/Plotter.h"
namespace Euler{
  class ErrorPlotter;
}

/**
 * Forward declaration
 */
namespace exahype {
  namespace solvers {
    class LimitingADERDGSolver;
  }
}



class Euler::ErrorPlotter: public exahype::plotters::Plotter::UserOnTheFlyPostProcessing{
  public:
  ErrorPlotter(exahype::solvers::LimitingADERDGSolver&  solver);
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
