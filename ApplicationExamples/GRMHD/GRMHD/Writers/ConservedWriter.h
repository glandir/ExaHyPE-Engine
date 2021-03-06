#ifndef __CONSERVED_WRITER_FV_CLASS_HEADER__
#define __CONSERVED_WRITER_FV_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "exahype/plotters/Plotter.h"
#include "exahype/solvers/LimitingADERDGSolver.h"

namespace GRMHD{
  class ConservedWriter;

  /**
   * Forward declaration
   */
  class GRMHDSolver_ADERDG;
  class GRMHDSolver_FV;
}




class GRMHD::ConservedWriter: public exahype::plotters::Plotter::UserOnTheFlyPostProcessing{
  public:
  ConservedWriter(GRMHDSolver_ADERDG&  solver);
  ConservedWriter(GRMHDSolver_FV&  solver);
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

#endif /* __CONSERVED_WRITER_FV_CLASS_HEADER__ */
