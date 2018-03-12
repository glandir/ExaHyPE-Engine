// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "AbstractGRMHDSolver_FV.h"

#include "kernels/finitevolumes/commons/c/commons.h"
#include "kernels/finitevolumes/godunov/c/godunov.h"
#include "kernels/finitevolumes/riemannsolvers/c/riemannsolvers.h"

#include "GRMHDSolver_FV.h" // Have to include a proper declaration. Cannot use forward declared classes in static_cast.

#include <stdio.h>
#include <cstdlib> // abort()
#include "kernels/KernelUtils.h" // idx

#include "exahype/disableOptimization.h" // we experience compiler bugs sometimes.

GRMHD::GRMHDSolver_FV::GRMHDSolver_FV(double maximumMeshSize,int maximumAdaptiveMeshDepth,exahype::solvers::Solver::TimeStepping timeStepping):
  AbstractGRMHDSolver_FV::AbstractGRMHDSolver_FV(maximumMeshSize,maximumAdaptiveMeshDepth,timeStepping) {
}

GRMHD::AbstractGRMHDSolver_FV::AbstractGRMHDSolver_FV(double maximumMeshSize,int maximumAdaptiveMeshDepth,exahype::solvers::Solver::TimeStepping timeStepping):
  exahype::solvers::FiniteVolumesSolver("GRMHDSolver_FV",NumberOfVariables,NumberOfParameters,PatchSize,
                                        GhostLayerWidth,maximumMeshSize,maximumAdaptiveMeshDepth,timeStepping) {
}

void GRMHD::AbstractGRMHDSolver_FV::constantsToString(std::ostream& os) {
	// This string is used in the --version output to identify compile time constants
	os << "GRMHD::AbstractGRMHDSolver_FV("
	   << "nVar=" << NumberOfVariables << ", "
	   << "nParam=" << NumberOfParameters << ", "
	   << "PatchSize=" << PatchSize << ", "
	   << "GhostLayerWidth=" << GhostLayerWidth
	   << ")";
}

void GRMHD::AbstractGRMHDSolver_FV::abortWithMsg(const char* const msg) {
	// verbosily fail even without assertions turned on
	puts(msg);
	abort();
}

void GRMHD::AbstractGRMHDSolver_FV::solutionUpdate(double* luhNew,const double* luh,const tarch::la::Vector<DIMENSIONS,double>& dx,const double dt,double& maxAdmissibleDt) {
  maxAdmissibleDt = kernels::finitevolumes::godunov::c::solutionUpdate<true, true, true, GRMHDSolver_FV>(*static_cast<GRMHDSolver_FV*>(this),luhNew,luh,dx,dt);
}


double GRMHD::AbstractGRMHDSolver_FV::stableTimeStepSize(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx) {
  double maxAdmissibleDt = kernels::finitevolumes::commons::c::stableTimeStepSize<GRMHDSolver_FV>(*static_cast<GRMHDSolver_FV*>(this),luh,dx);
  return maxAdmissibleDt;
}

void GRMHD::AbstractGRMHDSolver_FV::adjustSolution(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) {
  kernels::finitevolumes::commons::c::solutionAdjustment<GRMHDSolver_FV>(*static_cast<GRMHDSolver_FV*>(this),luh,center,dx,t,dt);
}

void GRMHD::AbstractGRMHDSolver_FV::boundaryConditions(double* luh,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const tarch::la::Vector<DIMENSIONS, int>& posCell,const tarch::la::Vector<DIMENSIONS, int>& posBoundary) {
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
  kernels::finitevolumes::commons::c::boundaryConditions<GRMHDSolver_FV>(*static_cast<GRMHDSolver_FV*>(this),luhbndOutside,luhbndInside,cellCentre,cellSize,t,dt,faceIndex,direction);
  ghostLayerFillingAtBoundary(luh,luhbndOutside,posBoundary-posCell);
  
  exahype::DataHeap::getInstance().deleteData(indexLuhbndInside,true);
  exahype::DataHeap::getInstance().deleteData(indexLuhbndOutside,true);
}


void GRMHD::AbstractGRMHDSolver_FV::ghostLayerFilling(double* luh,const double* luhNeighbour,const tarch::la::Vector<DIMENSIONS,int>& neighbourPosition) {
  kernels::finitevolumes::commons::c::ghostLayerFilling<GRMHDSolver_FV>(*static_cast<GRMHDSolver_FV*>(this),luh,luhNeighbour,neighbourPosition);
}

void GRMHD::AbstractGRMHDSolver_FV::ghostLayerFillingAtBoundary(double* luh,const double* luhbnd,const tarch::la::Vector<DIMENSIONS,int>& boundaryPosition) {
  kernels::finitevolumes::commons::c::ghostLayerFillingAtBoundary<GRMHDSolver_FV>(*static_cast<GRMHDSolver_FV*>(this),luh,luhbnd,boundaryPosition);
}

void GRMHD::AbstractGRMHDSolver_FV::boundaryLayerExtraction(double* luhbnd,const double* luh,const tarch::la::Vector<DIMENSIONS,int>& boundaryPosition) {
  kernels::finitevolumes::commons::c::boundaryLayerExtraction<GRMHDSolver_FV>(*static_cast<GRMHDSolver_FV*>(this),luhbnd,luh,boundaryPosition);
}


double GRMHD::AbstractGRMHDSolver_FV::riemannSolver(double* fL, double *fR, const double* qL, const double* qR, int direction) {
  // Default FV Riemann Solver
  return kernels::finitevolumes::riemannsolvers::c::rusanov<true, true, GRMHDSolver_FV>(*static_cast<GRMHDSolver_FV*>(this), fL,fR,qL,qR,direction);
}


/**
 * Fallback implementations of joined functions. Users can either safely ignore this
 * or overwrite with their own implementations.
 **/
 // @todo Can we remove this one?
#include "kernels/fusedMethods.cpph"
