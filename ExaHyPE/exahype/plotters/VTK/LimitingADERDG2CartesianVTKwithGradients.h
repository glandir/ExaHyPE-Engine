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

#ifndef _EXAHYPE_PLOTTERS_LIMITING_ADERDG_2_CARTESIAN_VTK_WG_H_
#define _EXAHYPE_PLOTTERS_LIMITING_ADERDG_2_CARTESIAN_VTK_WG_H_

#include "tarch/plotter/griddata/blockstructured/PatchWriterUnstructured.h"
#include "tarch/plotter/griddata/VTUTimeSeriesWriter.h"

#include "exahype/plotters/Plotter.h"
#include "exahype/plotters/slicing/Slicer.h"

namespace exahype {
  namespace plotters {
    class LimitingADERDG2CartesianVTKwithGradients;

    class LimitingADERDG2CartesianVerticesVTKAsciiwithGradients;
    class LimitingADERDG2CartesianVerticesVTKBinarywithGradients;
    class LimitingADERDG2CartesianCellsVTKAsciiwithGradients;
    class LimitingADERDG2CartesianCellsVTKBinarywithGradients;

    class LimitingADERDG2CartesianVerticesVTUAsciiwithGradients;
    class LimitingADERDG2CartesianVerticesVTUBinarywithGradients;
    class LimitingADERDG2CartesianCellsVTUAsciiwithGradients;
    class LimitingADERDG2CartesianCellsVTUBinarywithGradients;

    class Slicer; // external forward decl, #include exahype/plotters/slicing/Slicer.h
  }
}

/**
 * Common VTK class. Usually not used directly but through one of the subclasses.
 */
class exahype::plotters::LimitingADERDG2CartesianVTKwithGradients: public exahype::plotters::Plotter::Device {
protected:
   enum class PlotterType {
     BinaryVTK,
     ASCIIVTK,
     BinaryVTU,
     ASCIIVTU
   };
private:
  int           _fileCounter     = -1;
  const PlotterType _plotterType = PlotterType::BinaryVTU;
  const bool    _plotCells       = false;
  std::string   _filename        = "";
  int           _order           = -1;
  int           _solverUnknowns  = -1;
  int           _writtenUnknowns = -1;
  exahype::parser::ParserView   _plotterParameters ;

  /**
   * To resolve features of a high order polynomial solution
   * in a VTK compatible way, we apply the volume interpolation operator
   * to the DG solution polynomial.
   *
   * @return number of times the prolongation operators should be applied.
   */
  int                          _resolution = 0;


  /**
   * Temporary solution and gradient arrays used for interpolating the solution
   * onto finer grids.
   */
  std::vector<double>          _tempSolution;
  std::vector<double>          _tempGradient;
  /**
   * The ghost layer width the finite volumes patch is using.
   */
  const int     _ghostLayerWidth    = -1;

  exahype::plotters::Slicer *_slicer = nullptr;

  static tarch::logging::Log _log;

  /**
   * To memorise the time argument from startPlotter(). We need it when we close the plotter for the time series.
   */
  double _time = 0;

  /**
   * Is obviously only used if we use vtu instead of the vtk legacy format.
   */
  tarch::plotter::griddata::VTUTimeSeriesWriter _timeSeriesWriter;

  tarch::plotter::griddata::blockstructured::PatchWriter::SinglePatchWriter* _gridWriter  = nullptr;
  tarch::plotter::griddata::blockstructured::PatchWriterUnstructured*        _patchWriter = nullptr;

  tarch::plotter::griddata::Writer::VertexDataWriter*                        _vertexDataWriter = nullptr;
  tarch::plotter::griddata::Writer::CellDataWriter*                          _cellDataWriter   = nullptr;
  tarch::plotter::griddata::Writer::VertexDataWriter*                        _timeStampVertexDataWriter = nullptr;
  tarch::plotter::griddata::Writer::CellDataWriter*                          _timeStampCellDataWriter   = nullptr;
  tarch::plotter::griddata::Writer::CellDataWriter*                          _cellRefinementStatusWriter   = nullptr;
  tarch::plotter::griddata::Writer::VertexDataWriter*                        _vertexRefinementStatusWriter         = nullptr;
  tarch::plotter::griddata::Writer::CellDataWriter*                          _cellPreviousRefinementStatusWriter   = nullptr;
  tarch::plotter::griddata::Writer::VertexDataWriter*                        _vertexPreviousRefinementStatusWriter = nullptr;

  void writeTimeStampDataToADERDGPatch( double timeStamp, int vertexIndex );

  void plotVertexData(
      int firstVertexIndex,
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
      double* u,
      double timeStamp,
      const int RefinementStatus,
      const int previousRefinementStatusAsInt
  );

  void plotCellData(
      int firstCellIndex,
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
      double* u,
      double timeStamp,
      const int RefinementStatusAsInt,
      const int previousRefinementStatusAsInt
  );
public:
  LimitingADERDG2CartesianVTKwithGradients(
      exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
      const int ghostLayerWidth,
      PlotterType plotterType,
      bool plotCells);

  virtual ~LimitingADERDG2CartesianVTKwithGradients();

  virtual void init(const std::string& filename, int orderPlusOne, int solverUnknowns, int writtenUnknowns, exahype::parser::ParserView plotterParameters);

  void plotPatch(const int solverNumber,solvers::Solver::CellInfo& cellInfo) override;

  /**
   * Plot an ADER-DG solution.
   */
  void plotADERDGPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
      double timeStamp,
      const int RefinementStatusAsInt,
      const int previousRefinementStatusAsInt);

  /**
   * Plot a finite volumes solution.
   */
  void plotFiniteVolumesPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
      double timeStamp,
      const int numberOfCellsPerAxis,
      const int RefinementStatusAsInt,
      const int previousRefinementStatusAsInt);

  void startPlotting( double time ) override;
  void finishPlotting() override;
};

// VTK subclasses
class exahype::plotters::LimitingADERDG2CartesianVerticesVTKAsciiwithGradients: public exahype::plotters::LimitingADERDG2CartesianVTKwithGradients {
  public:
    static std::string getIdentifier();
    LimitingADERDG2CartesianVerticesVTKAsciiwithGradients(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
                                             const int ghostLayerWidth);
};
class exahype::plotters::LimitingADERDG2CartesianVerticesVTKBinarywithGradients: public exahype::plotters::LimitingADERDG2CartesianVTKwithGradients {
  public:
    static std::string getIdentifier();
    LimitingADERDG2CartesianVerticesVTKBinarywithGradients(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
                                              const int ghostLayerWidth);
};
class exahype::plotters::LimitingADERDG2CartesianCellsVTKAsciiwithGradients: public exahype::plotters::LimitingADERDG2CartesianVTKwithGradients {
  public:
    static std::string getIdentifier();
    LimitingADERDG2CartesianCellsVTKAsciiwithGradients(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
                                          const int ghostLayerWidth);
};
class exahype::plotters::LimitingADERDG2CartesianCellsVTKBinarywithGradients: public exahype::plotters::LimitingADERDG2CartesianVTKwithGradients {
  public:
    static std::string getIdentifier();
    LimitingADERDG2CartesianCellsVTKBinarywithGradients(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
                                           const int ghostLayerWidth);
};

// VTU subclasses
class exahype::plotters::LimitingADERDG2CartesianVerticesVTUAsciiwithGradients: public exahype::plotters::LimitingADERDG2CartesianVTKwithGradients {
  public:
    static std::string getIdentifier();
    LimitingADERDG2CartesianVerticesVTUAsciiwithGradients(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,const int ghostLayerWidth);
};


class exahype::plotters::LimitingADERDG2CartesianVerticesVTUBinarywithGradients: public exahype::plotters::LimitingADERDG2CartesianVTKwithGradients {
  public:
    static std::string getIdentifier();
    LimitingADERDG2CartesianVerticesVTUBinarywithGradients(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,const int ghostLayerWidth);
};


class exahype::plotters::LimitingADERDG2CartesianCellsVTUAsciiwithGradients: public exahype::plotters::LimitingADERDG2CartesianVTKwithGradients {
  public:
    static std::string getIdentifier();
    LimitingADERDG2CartesianCellsVTUAsciiwithGradients(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,const int ghostLayerWidth);
};


class exahype::plotters::LimitingADERDG2CartesianCellsVTUBinarywithGradients: public exahype::plotters::LimitingADERDG2CartesianVTKwithGradients {
  public:
    static std::string getIdentifier();
    LimitingADERDG2CartesianCellsVTUBinarywithGradients(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,const int ghostLayerWidth);
};

#endif // _EXAHYPE_PLOTTERS_LIMITING_ADERDG_2_CARTESIAN_VTK_H_
