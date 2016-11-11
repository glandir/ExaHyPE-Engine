#include "MHDSolver_Plotter5.h"

#include "MHDSolver.h"
#include "InitialDataAdapter.h"

MHDSolver::MHDSolver_Plotter5::MHDSolver_Plotter5(MHDSolver&  solver) {
  // @todo Please insert your code here
}


MHDSolver::MHDSolver_Plotter5::~MHDSolver_Plotter5() {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter5::startPlotting(double time) {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter5::finishPlotting() {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter5::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  const int nDim = MHDSolver::MHDSolver::nDim;
  const int nVar = MHDSolver::MHDSolver::nVar;
  
  // convert tarch::la::Vector to double*  
  double xpos[nDim];
  double exact[nVar];
  for(int i=0; i<nDim; i++) xpos[i] = x[i];

  alfenwave_(xpos, exact, &timeStamp);
  
  // compute the difference
  for (int i=0; i<nVar; i++){ 
    outputQuantities[i] = ( Q[i] - exact[i] ) / exact[i];
  }
}


