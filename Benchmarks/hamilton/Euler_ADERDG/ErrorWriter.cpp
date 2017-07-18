// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ErrorWriter.h"

#include "EulerSolver_ADERDG.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h"

#include "peano/utils/Loop.h"

#include "tarch/la/VectorOperations.h"

#include <algorithm>

#include <iomanip>

Euler::ErrorWriter::ErrorWriter() : exahype::plotters::ADERDG2UserDefined::ADERDG2UserDefined(){
  // @TODO Please insert your code here.
}


void Euler::ErrorWriter::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
    double timeStamp) {
  constexpr int numberOfVariables = AbstractEulerSolver_ADERDG::NumberOfVariables;
  constexpr int basisSize         = AbstractEulerSolver_ADERDG::Order+1;
  constexpr int order             = basisSize-1;

  double x[DIMENSIONS];

  kernels::idx4 idx(basisSize,basisSize,basisSize,numberOfVariables);
  dfor(i,basisSize) {
     double w_dV = 1.0;
     for (int d=0; d<DIMENSIONS; d++) {
       x[d]  = offsetOfPatch[d] + sizeOfPatch[d] * kernels::gaussLegendreNodes[order][i(d)];
       w_dV *= sizeOfPatch[d] * kernels::gaussLegendreWeights[order][i(d)];
     }

     double uAna[numberOfVariables];
     EulerSolver_ADERDG::referenceSolution(x,timeStamp,uAna);

     const double* uNum = u + idx ( (DIMENSIONS==3) ? i(2) : 0, i(1), i(0), 0);

     for (int v=0; v<numberOfVariables; v++) {
        const double uDiff = std::abs(uNum[v]-uAna[v]);
        errorL2[v]   += uDiff*uDiff * w_dV;
        errorL1[v]   += uDiff * w_dV;
        errorLInf[v]  = std::max( errorLInf[v], uDiff );

        normL1Ana[v]  += std::abs(uAna[v]) * w_dV;
        normL2Ana[v]  += uAna[v] * uAna[v] * w_dV;
        normLInfAna[v] = std::max( normLInfAna[v], std::abs(uAna[v]) );
     }
  }
}

void Euler::ErrorWriter::startPlotting( double time) {
  _timeStamp = time;

  std::fill_n(errorL1,  AbstractEulerSolver_ADERDG::NumberOfVariables, 0.0);
  std::fill_n(errorL2,  AbstractEulerSolver_ADERDG::NumberOfVariables, 0.0);
  std::fill_n(errorLInf,AbstractEulerSolver_ADERDG::NumberOfVariables, 0.0);
  
  std::fill_n(normL1Ana,  AbstractEulerSolver_ADERDG::NumberOfVariables, 0.0);
  std::fill_n(normL2Ana,  AbstractEulerSolver_ADERDG::NumberOfVariables, 0.0);
  std::fill_n(normLInfAna,AbstractEulerSolver_ADERDG::NumberOfVariables, 0.0);
}

void Euler::ErrorWriter::finishPlotting() {
  constexpr int numberOfVariables = AbstractEulerSolver_ADERDG::NumberOfVariables;

  for (int v=0; v<numberOfVariables; v++) {
    errorL2[v]   = sqrt(errorL2[v]);
    normL2Ana[v] = sqrt(normL2Ana[v]);
  }

  std::cout << "**Errors for ADER-DG solver with order="<<AbstractEulerSolver_ADERDG::Order<<"**" << std::endl;
  std::cout << "t_eval : "<<_timeStamp << std::endl;
  std::cout << "variable     : ";
  for (int v=0; v<numberOfVariables; v++) {
    std::cout << v << ", ";
  }
  std::cout << std::endl;

  std::cout << "absErrorL1   : ";
  for (int v=0; v<numberOfVariables; v++) {
    std::cout << std::setprecision(2) << errorL1[v] << ", ";
  }
  std::cout << std::endl;

  std::cout << "absErrorL2   : ";
  for (int v=0; v<numberOfVariables; v++) {
    std::cout << std::setprecision(2) << errorL2[v] << ", ";
  }
  std::cout << std::endl;

  std::cout << "absErrorLInf : ";
  for (int v=0; v<numberOfVariables; v++) {
    std::cout << std::setprecision(2) << errorLInf[v] << ", ";
  }
  std::cout << std::endl;

  std::cout << "relErrorL1   : ";
  for (int v=0; v<numberOfVariables; v++) {
    std::cout << std::setprecision(2) << errorL1[v]/normL1Ana[v] << ", ";
  }
  std::cout << std::endl;

  std::cout << "relErrorL2   : ";
  for (int v=0; v<numberOfVariables; v++) {
    std::cout << std::setprecision(2) << errorL2[v]/normL2Ana[v] << ", ";
  }
  std::cout << std::endl;

  std::cout << "relErrorLInf : ";
  for (int v=0; v<numberOfVariables; v++) {
    std::cout << std::setprecision(2) << errorLInf[v]/normLInfAna[v] << ", ";
  }
  std::cout << std::endl;
}
