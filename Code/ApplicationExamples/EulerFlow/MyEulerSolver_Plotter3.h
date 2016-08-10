// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "exahype/plotters/Plotter.h"
/**
 * Forward declaration
 */
class MyEulerSolver;


class MyEulerSolver_Plotter3: public exahype::plotters::Plotter::UserOnTheFlyPostProcessing{
  private:
    double time;
  public:
  MyEulerSolver_Plotter3(MyEulerSolver& solver);
  virtual ~MyEulerSolver_Plotter3();
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
