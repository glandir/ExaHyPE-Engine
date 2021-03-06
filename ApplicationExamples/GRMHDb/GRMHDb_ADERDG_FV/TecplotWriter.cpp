// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "TecplotWriter.h"

#include "TECPLOTinterface.h"
#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include <stdio.h>
#include "GRMHDbSolver_FV.h"
#include "GRMHDbSolver_ADERDG.h"

#include "kernels/GaussLegendreBasis.h"
#include "kernels/KernelUtils.h"

#include "peano/utils/Loop.h"

#include "tarch/la/VectorOperations.h"

#include <algorithm>

#include <iomanip>

//GRMHDb::TecplotWriter::TecplotWriter() : exahype::plotters::ADERDG2UserDefined::ADERDG2UserDefined(){
GRMHDb::TecplotWriter::TecplotWriter() : exahype::plotters::LimitingADERDG2UserDefined::LimitingADERDG2UserDefined(){
  // @TODO Please insert your code here.
}




void GRMHDb::TecplotWriter::plotADERDGPatch(
	const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* const u,
	double timeStamp) {

	plotForADERSolver = 1;
	plotForFVSolver = 0;
	elementcalltecplotaderdgplotter_(u, &offsetOfPatch[0], &sizeOfPatch[0], &plotForADERSolver);
}


void GRMHDb::TecplotWriter::plotFiniteVolumesPatch(
	const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* const u,
	double timeStamp) {
	// @TODO Please insert your code here.
	plotForADERSolver = 0;
	plotForFVSolver = 1;
	elementcalltecplotfvplotter_(u, &offsetOfPatch[0], &sizeOfPatch[0], &plotForADERSolver);
}

void GRMHDb::TecplotWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
  mpirank = tarch::parallel::Node::getInstance().getRank();
	// printf("\n******************************************************************");
	// printf("\n********* initializetecplotplotter *************");
	// printf("\n******************************************************************");
    initializetecplotplotter_(&time);
}


void GRMHDb::TecplotWriter::finishPlotting() {
  // @TODO Please insert your code here.

		finishtecplotplotter_(&mpirank);
  //printf("2: finishtecplotplotter: DONE! %d ", mpirank);
  //fflush(stdout);
}