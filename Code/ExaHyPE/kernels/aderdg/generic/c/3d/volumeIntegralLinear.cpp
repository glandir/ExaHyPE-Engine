/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released unter the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
#include "../../Kernels.h"

#include <cstring>

namespace kernels {
namespace aderdg {
namespace generic {
namespace c {

#if DIMENSIONS == 3

void volumeIntegralLinear(double* lduh, const double* const lFhi,
                          const tarch::la::Vector<DIMENSIONS, double>& dx,
                          const int numberOfVariables, const int basisSize) {
  // for linear non-conservative PDE, the volume integral is trivial, since it
  // only involves the element mass matrix, which later will cancel

  const int basisSize2 = basisSize * basisSize;
  const int basisSize3 = basisSize2 * basisSize;
  const int order = basisSize - 1;

  std::memset(lduh, 0, basisSize3 * numberOfVariables * sizeof(double));

  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      for (int k = 0; k < basisSize; k++) {
        double weight = kernels::gaussLegendreWeights[order][i] *
                        kernels::gaussLegendreWeights[order][j] *
                        kernels::gaussLegendreWeights[order][k];

        // Fortran: lduh(:,k,j,i) = -SUM(lFhi(:,k,j,i,1:nDim), dim = 5) * weight
        for (int l = 0; l < numberOfVariables; l++) {
          for (int m = 0; m < 3; m++) {
            lduh[i * basisSize2 * numberOfVariables +
                 j * basisSize * numberOfVariables + k * numberOfVariables +
                 l] -= weight * lFhi[m * basisSize3 * numberOfVariables +
                                     i * basisSize2 * numberOfVariables +
                                     j * basisSize * numberOfVariables +
                                     k * numberOfVariables + l];
          }
        }
      }
    }
  }
}

#endif  // DIMENSIONS == 3

}  // namespace c
}  // namespace generic
}  // namespace aderdg
}  // namespace kernels
