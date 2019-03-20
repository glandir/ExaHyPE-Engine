#ifndef __H5_WRITER_FV_CLASS_HEADER__
#define __H5_WRITER_FV_CLASS_HEADER__
// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "exahype/plotters/Plotter.h"
namespace GRMHD{
  class H5Writer;

  /**
   * Forward declaration
   */
  class GRMHDSolver_FV;
  class GRMHDSolver_ADERDG;
}




class GRMHD::H5Writer: public exahype::plotters::Plotter::UserOnTheFlyPostProcessing{
  public:
  H5Writer(GRMHDSolver_FV&  solver);
  H5Writer(GRMHDSolver_ADERDG&  solver);
  virtual ~H5Writer();
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

#endif