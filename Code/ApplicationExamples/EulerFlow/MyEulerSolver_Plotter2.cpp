#include "MyEulerSolver_Plotter2.h"
#include "GeneratedConstants.h"
#include "Primitives.h"


Euler::MyEulerSolver_Plotter2::MyEulerSolver_Plotter2(MyEulerSolver& solver) {
  // @todo Please insert your code here
}


Euler::MyEulerSolver_Plotter2::~MyEulerSolver_Plotter2() {
  // @todo Please insert your code here
}


void Euler::MyEulerSolver_Plotter2::startPlotting(double time) {
  // @todo Please insert your code here
}


void Euler::MyEulerSolver_Plotter2::finishPlotting() {
  // @todo Please insert your code here
}


void Euler::MyEulerSolver_Plotter2::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {

/**
 * This is the Primitive Variable plotter
 **/
  //double V[Euler::MyEulerSolver::nVar];
  cons2prim(outputQuantities, Q);

  //for (int i=0; i<0; i++){ 
  //  outputQuantities[i] = V[i];
  //}
}


