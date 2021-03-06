// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
// ========================
//   www.exahype.eu
// ========================
#ifndef TecplotWriter_CLASS_HEADER_
#define TecplotWriter_CLASS_HEADER_

#include "exahype/plotters/LimitingADERDG2UserDefined.h"
#include "exahype/plotters/ADERDG2UserDefined.h"
#include "exahype/plotters/ascii/MultipleReductionsWriter.h"
#include "exahype/solvers/LimitingADERDGSolver.h"
#include "GRMHDbSolver_ADERDG_Variables.h"
#include "GRMHDbSolver_ADERDG.h"
#include "GRMHDbSolver_FV.h"

namespace GRMHDb {
  class TecplotWriter;

  //class GRMHDSolver_ADERDG;
  //class GRMHDSolver_FV;
}

class GRMHDb::TecplotWriter : public exahype::plotters::LimitingADERDG2UserDefined {
private:
	int plotForADERSolver;
	int plotForFVSolver;
	int mpirank;
	//int counterloc;
	//int counterloc2;
	//static constexpr int numberOfVariables = AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
	//static constexpr int numberOfParameters = AbstractGRMHDbSolver_ADERDG::NumberOfParameters;
	//static constexpr int nVar = numberOfVariables + numberOfParameters;
 public:
  //static constexpr int nVar = GRMHDb::AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
 //static constexpr int order = GRMHDb::AbstractGRMHDbSolver_ADERDG::Order;
 //static constexpr int basisSize = order + 1;
  /**
   * Constructor.
   * 
   * \note ExaHyPE does not increment file counters for
   * you if you use user defined plotting. You have
   * to declare and manage such member variables yourself. 
   */
  TecplotWriter();
  //TecplotWriter(GRMHDb::GRMHDbSolver_FV& solver);
  //TecplotWriter(GRMHDb::GRMHDbSolver_ADERDG& solver);
  //TecplotWriter(exahype::solvers::LimitingADERDGSolver&  solver);

  /**
   * This method is invoked every time a cell 
   * is touched by the plotting plotter.
   *
   * \note Use the protected variables _order, _variables to
   * determine the size of u. 
   * The array u has the size _variables * (_order+1)^DIMENSIONS.
   * 
   * \param[in] offsetOfPatch the offset of the cell/patch.
   * \param[in] sizeOfPatch the offset of the cell/patch.
   * \param[in] u the degrees of freedom "living" inside of the patch.
   */
  /* void plotPatch(
	  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
	  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
	  double timeStamp) override;*/
  void plotADERDGPatch(
	  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* const u,
	  double timeStamp) override;

  /**
   * This method is invoked every time a Finite Volume cell 
   * is touched by the plotting plotter.
   *
   * \note Use the protected variables _order, _variables to
   * determine the size of u. 
   * The array u has the size _variables * (_order+1)^DIMENSIONS.
   * 
   * \param[in] offsetOfPatch the offset of the cell/patch.
   * \param[in] sizeOfPatch the offset of the cell/patch.
   * \param[in] u the degrees of freedom "living" inside of the patch.
   */
  void plotFiniteVolumesPatch(
	  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* const u,
	  double timeStamp) override;

  /** 
   * This method is called at the beginning of the plotting.
   * You can use it to reset member variables, e.g., those
   * used for calculations, or to increment file counters.
   *
   * \param[in] time a characteristic solver time stamp.
   *            Usually the global minimum.
   */
  void startPlotting( double time) override;
  
  /** 
   * This method is called at the end of the plotting.
   * You can use it to reset member variables, finalise calculations (compute square roots etc.),
   * or to increment file counters
   */
  void finishPlotting() override;
};

#endif /* TecplotWriter_CLASS_HEADER_ */