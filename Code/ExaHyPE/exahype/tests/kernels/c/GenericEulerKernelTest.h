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
 
// @todo ExaHyPE Lizenz
#ifndef _EXAHYPE_TESTS_GENERIC_EULER_KERNEL_TEST_H_
#define _EXAHYPE_TESTS_GENERIC_EULER_KERNEL_TEST_H_

#include "peano/utils/Globals.h"
#include "tarch/tests/TestCase.h"

namespace exahype {
namespace tests {
namespace c {

class GenericEulerKernelTest : public tarch::tests::TestCase {
 public:
  GenericEulerKernelTest();
  virtual ~GenericEulerKernelTest();

  void run() override;

 private:
  void testPDEFluxes();
  void testSpaceTimePredictorLinear();
  void testSpaceTimePredictorNonlinear();
  void testVolumeIntegralLinear();
  void testVolumeIntegralNonlinear();
  void testRiemannSolverLinear();
  void testRiemannSolverNonlinear();
  void testSurfaceIntegralLinear();
  void testSurfaceIntegralNonlinear();
  void testSolutionUpdate();
  void testVolumeUnknownsProjection();
  void testFaceUnknownsProjection();
  void testEquidistantGridProjection();

  static void testFlux(const double* const Q, double** F);

  static void testEigenvalues(const double* const Q,
                              const int normalNonZeroIndex, double* lambda);

  static void testNCP(const double* const Q, const double* const gradQ,
                      double* BgradQ);

  static void testMatrixB(const double* const Q, const int normalNonZero,
                          double* Bn);

  const double eps = 1.0e-10;  // for quick adaption of the test cases (say,
                               // switch to single precision)
};

}  // namespace c
}  // namespace tests
}  // namespace exahype

#endif
