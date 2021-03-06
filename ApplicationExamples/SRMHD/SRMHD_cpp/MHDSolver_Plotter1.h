// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "exahype/plotters/Plotter.h"
namespace MHDSolver{
  class MHDSolver_Plotter1;

  /**
   * Forward declaration
   */
  class MHDSolver;
}


#include "TimeSeriesReductions.h"
#include "GeneratedConstants.h"
#define nVars MY_NUMBER_OF_VARIABLES /* MHD */

class MHDSolver::MHDSolver_Plotter1: public exahype::plotters::Plotter::UserOnTheFlyPostProcessing{
  public:
  TimeSeriesReductions* conserved[nVars];
  TimeSeriesReductions* primitives[nVars];
  TimeSeriesReductions* errors[nVars];
  TimeSeriesReductions* statistics;

  MHDSolver_Plotter1(MHDSolver&  solver);
  virtual ~MHDSolver_Plotter1();
  void startPlotting(double time) override;
  void finishPlotting() override;
  void mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp) override;
};
