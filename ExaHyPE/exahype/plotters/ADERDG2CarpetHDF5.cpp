#include "exahype/plotters/ADERDG2CarpetHDF5.h"

#include <cstdlib>
#include <stdio.h>
#include <sstream>


std::string exahype::plotters::ADERDG2CarpetHDF5::getIdentifier() {
	return std::string("experimental::CarpetHDF5");
}


#ifndef HDF5
/*************************************************************************************************
 * ADERDG2CarpetHDF5 Dummy implementation in case HDF5 support is skipped.
 *************************************************************************************************/

exahype::plotters::ADERDG2CarpetHDF5::ADERDG2CarpetHDF5(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
    Device(postProcessing) {
	printf("ERROR: Compile with HDF5, otherwise you cannot use the HDF5 plotter.\n");
	abort();
}

// all other methods are stubs
exahype::plotters::ADERDG2CarpetHDF5::~ADERDG2CarpetHDF5() {}
void exahype::plotters::ADERDG2CarpetHDF5::init(const std::string& filename, int orderPlusOne, int solverUnknowns, int writtenUnknowns, const std::string& select) {}
void exahype::plotters::ADERDG2CarpetHDF5::plotPatch(const int cellDescriptionsIndex, const int element) {}
void exahype::plotters::ADERDG2CarpetHDF5::startPlotting(double time) {}
void exahype::plotters::ADERDG2CarpetHDF5::finishPlotting() {}

#else
/*************************************************************************************************
 * The HDF5 Writer code actually starts here.
 * I embed the C++ class here mainly for laziness.
 *************************************************************************************************/

#include "kernels/KernelUtils.h" // indexing
#include "kernels/DGMatrices.h"
#include "peano/utils/Loop.h"
#include "exahype/solvers/ADERDGSolver.h"
#include "kernels/DGBasisFunctions.h"

// HDF5 stuff
#include "H5Cpp.h"

class exahype::plotters::ADERDG2CarpetHDF5Writer {
public:
  // The ADERDG2CarpetHDF5 instance
  exahype::plotters::ADERDG2CarpetHDF5* device;
  FILE* fh;
  H5::H5File* hf;
  
  std::string txtfilename;
  std::string h5filename;
  
  ADERDG2CarpetHDF5Writer(exahype::plotters::ADERDG2CarpetHDF5* _device) : device(_device) {
  }
  
  void init() {
	// txt for cross checking
	txtfilename = device->filename + ".txt";
	h5filename = device->filename + ".h5";
	
	fh = fopen(txtfilename.c_str(), "w");
	if(!fh) {
		printf("Could not open file '%s'\n", txtfilename.c_str());
		abort();
	}
	fprintf(fh, "# This is the ADERDG2CarpetHDF5Writer, just checking out...");
	
	using namespace H5;
	
	hf = new H5File( H5std_string(h5filename.c_str()), H5F_ACC_TRUNC  );
	hf->setComment("Created by ExaHyPE");
	
	// write basic group
	Group* parameters = new Group(hf->createGroup( "/Parameters and Global Attributes" ));
	
	int ranks = 1;
	
	DataSpace dspace(H5S_SCALAR);
	Attribute nioprocs = parameters->createAttribute("nioprocs", PredType::NATIVE_INT, dspace);
	nioprocs.write(PredType::NATIVE_INT, &ranks);
  }

  /**
   * Assumes u to be a (order+1,order+1) matrix of single values.
   * This is 2D.
   * This is only one writeOut quantity.
   * This is simple.
   **/
  void plotPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
      double* mappedCell,
      double timeStamp, int iteration, int component) {
        using namespace H5;
	
	const int order = device->order;
	const int rank = 2; // DIMENSIONS
	hsize_t  dims[rank] = {order+1, order+1};
	DataSpace patch_dataspace(rank, dims);
	kernels::idx2 mappedIdx(order+1, order+1);
	
	char component_name[100];
	sprintf(component_name, "EXA::rho it=%d tl=0 m=0 rl=0 c=%d", iteration, component);
	
	DataSet table = hf->createDataSet(component_name, PredType::NATIVE_DOUBLE, patch_dataspace);
	table.write(mappedCell, PredType::NATIVE_DOUBLE);
	
	// write the meat information about the new table
	
	hsize_t dtuple_len[rank] = {1,1};
	DataSpace dtuple(rank, dtuple_len);
	
	Attribute origin = table.createAttribute("origin", PredType::NATIVE_DOUBLE, dtuple);
	origin.write(PredType::NATIVE_DOUBLE, offsetOfPatch.data());
	
	double iorigin_data[rank] = {0.,0.};
	Attribute iorigin = table.createAttribute("iorigin", PredType::NATIVE_DOUBLE, dtuple);
	iorigin.write(PredType::NATIVE_DOUBLE, &iorigin_data);
	
	
	// alternatively, also write to text file for comparison
	dfor(i,order+1) {
		fprintf(fh, "%d %d %f\n", i(0), i(1), mappedCell[mappedIdx(i(0),i(1))]);
	}
  }

}; // class ADERDG2CarpetHDF5Impl

/*************************************************************************************************
 * ADERDG2CarpetHDF5 non-dummy implementation
 *************************************************************************************************/

exahype::plotters::ADERDG2CarpetHDF5::ADERDG2CarpetHDF5(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
    Device(postProcessing) {
	writer = new exahype::plotters::ADERDG2CarpetHDF5Writer(this);
}

// all other methods are stubs
exahype::plotters::ADERDG2CarpetHDF5::~ADERDG2CarpetHDF5() {
	delete writer;
}

void exahype::plotters::ADERDG2CarpetHDF5::init(const std::string& _filename, int _orderPlusOne, int _solverUnknowns, int _writtenUnknowns, const std::string& _select) {
	filename          = _filename;
	order             = _orderPlusOne-1;
	solverUnknowns    = _solverUnknowns;
	select            = _select;
	writtenUnknowns   = _writtenUnknowns;
	iteration         = 0;
	
	// todo at this place: Accept 3D slicing or so, cf. ADERDG2CartesianVTK
	
	writer->init(); // open files and so
}

void exahype::plotters::ADERDG2CarpetHDF5::plotPatch(
        const int cellDescriptionsIndex,
        const int element) {
  auto& aderdgCellDescription = exahype::solvers::ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  if (aderdgCellDescription.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Cell) {
    double* solverSolution = DataHeap::getInstance().getData(aderdgCellDescription.getSolution()).data();

    plotPatch(
        aderdgCellDescription.getOffset(),
        aderdgCellDescription.getSize(), solverSolution,
        aderdgCellDescription.getCorrectorTimeStamp());
  }
}

void exahype::plotters::ADERDG2CarpetHDF5::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    double* u,
    double timeStamp) {

    kernels::idx2 mappedIdx(order+1, order+1);
    double* mappedCell  = new double[mappedIdx.size];

    interpolateCartesianPatch(offsetOfPatch, sizeOfPatch, u, mappedCell, timeStamp);

    component++;
    writer->plotPatch(offsetOfPatch, sizeOfPatch, mappedCell, timeStamp, iteration, component);
}


void exahype::plotters::ADERDG2CarpetHDF5::startPlotting( double _time ) {
	time = _time;
	component = 0;
	_postProcessing->startPlotting(time);
	//writer->startPlotting(time);
}

void exahype::plotters::ADERDG2CarpetHDF5::finishPlotting() {
	_postProcessing->finishPlotting();
	iteration++;
}


void exahype::plotters::ADERDG2CarpetHDF5::interpolateCartesianPatch(
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
  double *u,
  double *mappedCell,
  double timeStamp
) {
  // the layout of an ADERDG cell, assuming 2D here for the time being.
  assert(DIMENSIONS==2);
  // make sure we only map to ONE written unknown, as this is how CarpetHDF5 works in the moment.
  assert(writtenUnknowns == 1);
  kernels::idx2 mappedIdx(order+1, order+1);

  double* interpoland = new double[solverUnknowns];
  //double* mappedCell  = new double[mappedIdx.size];
  //double* value       = writtenUnknowns==0 ? nullptr : new double[writtenUnknowns];
  
  dfor(i,order+1) {
    for (int unknown=0; unknown < solverUnknowns; unknown++) {
      interpoland[unknown] = 0.0;
      dfor(ii,order+1) { // Gauss-Legendre node indices
        int iGauss = peano::utils::dLinearisedWithoutLookup(ii,order + 1);
        interpoland[unknown] += kernels::equidistantGridProjector1d[order][ii(1)][i(1)] *
                 kernels::equidistantGridProjector1d[order][ii(0)][i(0)] *
                 #if DIMENSIONS==3
                 kernels::equidistantGridProjector1d[order][ii(2)][i(2)] *
                 #endif
                 u[iGauss * solverUnknowns + unknown];
        assertion3(interpoland[unknown] == interpoland[unknown], offsetOfPatch, sizeOfPatch, iGauss);
      }
    }

    double* value = mappedCell + mappedIdx(i(0), i(1));
    assertion(sizeOfPatch(0)==sizeOfPatch(1));
    _postProcessing->mapQuantities(
      offsetOfPatch,
      sizeOfPatch,
      offsetOfPatch + i.convertScalar<double>()* (sizeOfPatch(0)/(order)),
      i,
      interpoland,
      value,
      timeStamp
    );
  }

  delete[] interpoland;
  //if (value!=nullptr)        delete[] value;
}



#endif /* HDF5 */