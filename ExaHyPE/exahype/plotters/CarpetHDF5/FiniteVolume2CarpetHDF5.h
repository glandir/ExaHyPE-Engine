/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
#ifndef _EXAHYPE_PLOTTERS_FV_2_CARPETHDF5_H_
#define _EXAHYPE_PLOTTERS_FV_2_CARPETHDF5_H_

#include "exahype/plotters/Plotter.h"

namespace exahype {
  namespace plotters {
    class FiniteVolume2CarpetHDF5;
    class CarpetHDF5Writer;
  }
}

namespace kernels {
  class index; // instead #include "kernels/KernelUtils.h"
}

/**
 * Writing CarpetHDF5 files from FiniteVolume solvers.
 * 
 **/
class exahype::plotters::FiniteVolume2CarpetHDF5 : public exahype::plotters::Plotter::Device {
 public:
  /**
   * Pimpl idiom: In order to avoid any HDF5 dependency all HDF5 logic is hidden inside this
   * class (instance).
   **/
  CarpetHDF5Writer* writer;

  static std::string getIdentifier();

  FiniteVolume2CarpetHDF5(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);

  virtual ~FiniteVolume2CarpetHDF5();

  virtual void init(const std::string& filename, int orderPlusOne, int solverUnknowns, int writtenUnknowns, const std::string& select);

  virtual void plotPatch(
        const int cellDescriptionsIndex,
        const int element);
  
  virtual void plotPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
      double timeStamp);
  
  virtual void startPlotting( double time );
  virtual void finishPlotting();

  // TODO: These FV interpolating routines should be in some generic library
  void interpolateCartesianPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    double *u,
    double *mappedCell,
    double timeStamp
  );
};

#endif/* _EXAHYPE_PLOTTERS_ADERDG_2_CARPETHDF5_H_ */
