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



GRMHDb::TecplotWriter::TecplotWriter() : exahype::plotters::ADERDG2UserDefined::ADERDG2UserDefined(){
	// @TODO Please insert your code here.
}

/*GRMHDb::TecplotWriter::TecplotWriter(GRMHDb::GRMHDbSolver_ADERDG& solver) : TecplotWriter(){
	plotForADERSolver = 1;
}*/



   void GRMHDb::TecplotWriter::plotPatch(
   const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
   const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
   double timeStamp) {
// @TODO Please insert your code here.
//counterloc++;
//#ifdef Parallel
//  if (!tarch::parallel::NodePool::getInstance().isIdleNode(mpirank)) { 
//
//#endif
        
        //elementcalltecplotplotter_(u, &offsetOfPatch[0], &sizeOfPatch[0],&plotForADERSolver);

        elementcalltecplotaderdgplotter_(u, &offsetOfPatch[0], &sizeOfPatch[0], &plotForADERSolver);
        
//#ifdef Parallel
//	}
//#endif
//elementcalltecplotplotter_(u,&offsetOfPatch[0],&sizeOfPatch[0],&counterloc2);
}



/*
void GRMHDb::TecplotWriter::plotADERDGPatch(
		const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
		const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
		double timeStamp) {
	//constexpr int numberOfVariables = AbstractGRMHDbSolver_ADERDG::NumberOfVariables;
	//constexpr int numberOfParameters = AbstractGRMHDbSolver_ADERDG::NumberOfParameters;
	//constexpr int numberOfData = numberOfVariables + numberOfParameters;
	//constexpr int basisSize = AbstractGRMHDbSolver_ADERDG::Order + 1;
	//constexpr int order = basisSize - 1;

	// @TODO Please insert your code here.
	//counterloc++;
	//#ifdef Parallel
	//  if (!tarch::parallel::NodePool::getInstance().isIdleNode(mpirank)) {
	//
	plotForADERSolver = 1;
	plotForFVSolver = 0;
	//#endif
	//elementcalltecplotplotter_(u, &offsetOfPatch[0], &sizeOfPatch[0], &plotForADERSolver);
	elementcalltecplotaderdgplotter_(u, &offsetOfPatch[0], &sizeOfPatch[0], &plotForADERSolver);
	//#ifdef Parallel
	//      }
	//#endif
	//elementcalltecplotplotter_(u,&offsetOfPatch[0],&sizeOfPatch[0],&counterloc2);
}


void GRMHDb::TecplotWriter::plotFiniteVolumesPatch(
		const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
		const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
		double timeStamp) {
	// @TODO Please insert your code here.
	//constexpr int numberOfVariables = AbstractGRMHDbSolver_FV::NumberOfVariables;
	//constexpr int numberOfParameters = AbstractGRMHDbSolver_FV::NumberOfParameters;
	//constexpr int numberOfData = numberOfVariables + numberOfParameters;
	//constexpr int basisSize = AbstractGRMHDbSolver_FV::PatchSize;
	//constexpr int ghostLayerWidth = AbstractGRMHDbSolver_FV::GhostLayerWidth;
	//constexpr int order = basisSize - 1;

	plotForADERSolver = 0;
	plotForFVSolver = 1;
	//counterloc++;
	//#ifdef Parallel
	//  if (!tarch::parallel::NodePool::getInstance().isIdleNode(mpirank)) {
	//
	//#endif
	//elementcalltecplotplotter_(u, &offsetOfPatch[0], &sizeOfPatch[0], &plotForADERSolver);
	//printf("\n******************************************************************");
	//printf("\n********* initializetecplotplotter *************");
	//printf("\n**                      plotForADERSolver :  %d ", plotForADERSolver);
	elementcalltecplotfvplotter_(u, &offsetOfPatch[0], &sizeOfPatch[0], &plotForADERSolver);
	//#ifdef Parallel
	//      }
	//#endif
	//elementcalltecplotplotter_(u,&offsetOfPatch[0],&sizeOfPatch[0],&counterloc2);
}*/


void GRMHDb::TecplotWriter::startPlotting( double time) {
	// @TODO Please insert your code here.
	mpirank = tarch::parallel::Node::getInstance().getRank();
	//counterloc = 0;
	//counterloc2 = 0;
	//#ifdef Parallel
	//  if (!tarch::parallel::NodePool::getInstance().isIdleNode(mpirank)) {
	//	  printf("\n******************************************************************");
	//	  printf("\n********* initializetecplotplotter *************");
	//	  printf("\n******************************************************************");
	//#endif
	initializetecplotplotter_(&time);
	//#ifdef Parallel
	//  }
	//#endif
}


void GRMHDb::TecplotWriter::finishPlotting() {
	// @TODO Please insert your code here.
	//printf("1: finishtecplotplotter: DONE! %d ", mpirank);
	//fflush(stdout);
	//printf("  counterloc=%d", counterloc);
	//printf("  counterloc2=%d", counterloc2);

	//#ifdef Parallel
	//  if (!tarch::parallel::NodePool::getInstance().isIdleNode(mpirank)) {
	//#endif
	finishtecplotplotter_(&mpirank);
	//#ifdef Parallel
	//}
	//#endif
	//printf("2: finishtecplotplotter: DONE! %d ", mpirank);
	//fflush(stdout);
}
