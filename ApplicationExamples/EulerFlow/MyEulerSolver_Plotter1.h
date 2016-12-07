// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "exahype/plotters/Plotter.h"

namespace Euler{
  class MyEulerSolver_Plotter1;

  /**
   * Forward declaration
   */
  class MyEulerSolver;
}

/* I hope these modifications are not overwritten... */
#include "TimeSeriesReductions.h"
#include "MyEulerSolver.h"
static const int nVar = Euler::MyEulerSolver::nVar; // shortcut

namespace Euler{
  class MyEulerSolver_Plotter1;

  /**
   * Forward declaration
   */
  class MyEulerSolver;
}


class Euler::MyEulerSolver_Plotter1: public exahype::plotters::Plotter::UserOnTheFlyPostProcessing{
  TimeSeriesReductions* conserved[nVar];
  TimeSeriesReductions* primitives[nVar];
  TimeSeriesReductions* errors[nVar];
  TimeSeriesReductions* statistics;
  double time;
  
  // Auto generated:
  public:
  MyEulerSolver_Plotter1(MyEulerSolver& solver);
  virtual ~MyEulerSolver_Plotter1();
  virtual void startPlotting(double time);
  virtual void finishPlotting();
  virtual void mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
  double timeStamp);
};