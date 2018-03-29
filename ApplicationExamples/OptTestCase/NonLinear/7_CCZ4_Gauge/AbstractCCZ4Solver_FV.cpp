// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "AbstractCCZ4Solver_FV.h"

#include "kernels/finitevolumes/commons/c/commons.h"
#include "kernels/finitevolumes/musclhancock/c/musclhancock.h"
#include "kernels/finitevolumes/riemannsolvers/c/riemannsolvers.h"

#include "CCZ4Solver_FV.h" // Have to include a proper declaration. Cannot use forward declared classes in static_cast.

#include <stdio.h>
#include <cstdlib> // abort()
#include "kernels/KernelUtils.h" // idx

#include "exahype/disableOptimization.h" // we experience compiler bugs sometimes.

CCZ4::CCZ4Solver_FV::CCZ4Solver_FV(double maximumMeshSize,int maximumAdaptiveMeshDepth,exahype::solvers::Solver::TimeStepping timeStepping):
  AbstractCCZ4Solver_FV::AbstractCCZ4Solver_FV(maximumMeshSize,maximumAdaptiveMeshDepth,timeStepping) {
}

CCZ4::AbstractCCZ4Solver_FV::AbstractCCZ4Solver_FV(double maximumMeshSize,int maximumAdaptiveMeshDepth,exahype::solvers::Solver::TimeStepping timeStepping):
  exahype::solvers::FiniteVolumesSolver("CCZ4Solver_FV",NumberOfVariables,NumberOfParameters,PatchSize,
                                        GhostLayerWidth,maximumMeshSize,maximumAdaptiveMeshDepth,timeStepping) {
}

void CCZ4::AbstractCCZ4Solver_FV::constantsToString(std::ostream& os) {
	// This string is used in the --version output to identify compile time constants
	os << "CCZ4::AbstractCCZ4Solver_FV("
	   << "nVar=" << NumberOfVariables << ", "
	   << "nParam=" << NumberOfParameters << ", "
	   << "PatchSize=" << PatchSize << ", "
	   << "GhostLayerWidth=" << GhostLayerWidth
	   << ")";
}

void CCZ4::AbstractCCZ4Solver_FV::abortWithMsg(const char* const msg) {
	// verbosily fail even without assertions turned on
	puts(msg);
	abort();
}

void CCZ4::AbstractCCZ4Solver_FV::solutionUpdate(double* luhNew,const double* luh,const tarch::la::Vector<DIMENSIONS,double>& dx,const double dt,double& maxAdmissibleDt) {
  maxAdmissibleDt = kernels::finitevolumes::musclhancock::c::solutionUpdate<true, true, false, CCZ4Solver_FV>(*static_cast<CCZ4Solver_FV*>(this),luhNew,luh,dx,dt);
}


double CCZ4::AbstractCCZ4Solver_FV::stableTimeStepSize(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx) {
  double maxAdmissibleDt = kernels::finitevolumes::commons::c::stableTimeStepSize<CCZ4Solver_FV>(*static_cast<CCZ4Solver_FV*>(this),luh,dx);
  return maxAdmissibleDt;
}

void CCZ4::AbstractCCZ4Solver_FV::adjustSolution(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) {
  kernels::finitevolumes::commons::c::solutionAdjustment<CCZ4Solver_FV>(*static_cast<CCZ4Solver_FV*>(this),luh,center,dx,t,dt);
}

void CCZ4::AbstractCCZ4Solver_FV::boundaryConditions(double* luh,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const tarch::la::Vector<DIMENSIONS, int>& posCell,const tarch::la::Vector<DIMENSIONS, int>& posBoundary) {
  constexpr int cellsPerFace = PatchSize*PatchSize*GhostLayerWidth;
  constexpr int sizeLuhbnd = (NumberOfVariables+NumberOfParameters)*cellsPerFace;
  
  const int indexLuhbndInside        = exahype::DataHeap::getInstance().createData(
       sizeLuhbnd,sizeLuhbnd,exahype::DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired);
  const int indexLuhbndOutside       = exahype::DataHeap::getInstance().createData(
       sizeLuhbnd,sizeLuhbnd,exahype::DataHeap::Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired);
  double* luhbndInside               = exahype::DataHeap::getInstance().getData(indexLuhbndInside).data();
  double* luhbndOutside              = exahype::DataHeap::getInstance().getData(indexLuhbndOutside).data();
  
  const int direction   = tarch::la::equalsReturnIndex(posCell, posBoundary);
  const int orientation = (1 + posBoundary(direction) - posCell(direction))/2;
  const int faceIndex   = 2*direction+orientation;
  
  boundaryLayerExtraction(luhbndInside,luh,posBoundary-posCell);
  kernels::finitevolumes::commons::c::boundaryConditions<CCZ4Solver_FV>(*static_cast<CCZ4Solver_FV*>(this),luhbndOutside,luhbndInside,cellCentre,cellSize,t,dt,faceIndex,direction);
  ghostLayerFillingAtBoundary(luh,luhbndOutside,posBoundary-posCell);
  
  exahype::DataHeap::getInstance().deleteData(indexLuhbndInside,true);
  exahype::DataHeap::getInstance().deleteData(indexLuhbndOutside,true);
}


void CCZ4::AbstractCCZ4Solver_FV::ghostLayerFilling(double* luh,const double* luhNeighbour,const tarch::la::Vector<DIMENSIONS,int>& neighbourPosition) {
  kernels::finitevolumes::commons::c::ghostLayerFilling<CCZ4Solver_FV>(*static_cast<CCZ4Solver_FV*>(this),luh,luhNeighbour,neighbourPosition);
}

void CCZ4::AbstractCCZ4Solver_FV::ghostLayerFillingAtBoundary(double* luh,const double* luhbnd,const tarch::la::Vector<DIMENSIONS,int>& boundaryPosition) {
  kernels::finitevolumes::commons::c::ghostLayerFillingAtBoundary<CCZ4Solver_FV>(*static_cast<CCZ4Solver_FV*>(this),luh,luhbnd,boundaryPosition);
}

void CCZ4::AbstractCCZ4Solver_FV::boundaryLayerExtraction(double* luhbnd,const double* luh,const tarch::la::Vector<DIMENSIONS,int>& boundaryPosition) {
  kernels::finitevolumes::commons::c::boundaryLayerExtraction<CCZ4Solver_FV>(*static_cast<CCZ4Solver_FV*>(this),luhbnd,luh,boundaryPosition);
}


double CCZ4::AbstractCCZ4Solver_FV::riemannSolver(double* fL, double *fR, const double* qL, const double* qR, int direction) {
  // Default FV Riemann Solver
  return kernels::finitevolumes::riemannsolvers::c::rusanov<true, false, CCZ4Solver_FV>(*static_cast<CCZ4Solver_FV*>(this), fL,fR,qL,qR,direction);
}


/**
 * Fallback implementations of joined functions. Users can either safely ignore this
 * or overwrite with their own implementations.
 **/
 // @todo Can we remove this one?
#include "kernels/fusedMethods.cpph"
