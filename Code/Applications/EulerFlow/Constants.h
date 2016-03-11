/*
 * ! This file is autogenerated and should not be modified !
 * (This file should be generated according to a config file)
 *
 * Constants.h
 */

#ifndef EXAHYPE_CONSTANTS_H_
#define EXAHYPE_CONSTANTS_H_

#include "peano/utils/Globals.h"

/**
 * Mesh parameters.
 */
///@{
#define EXAHYPE_INITIAL_GLOBAL_REFINEMENT_LEVEL \
  3  //!< The refinement level of the initial grid.
#define EXAHYPE_INITIAL_ADAPTIVE_REFINEMENT_LEVEL \
  0  //!< The initial level of the adaptive refinement.
#define EXAHYPE_MAXIMUM_ADAPTIVE_REFINEMENT_LEVEL \
  0  //!< The maximum level of adaptive refinement.
#define EXAHYPE_MAXIMUM_REFINEMENT_LEVEL \
  3  //!< The sum of ::EXAHYPE_INITIAL_GLOBAL_REFINEMENT_LEVEL and
     //::EXAHYPE_MAXIMUM_ADAPTIVE_REFINEMENT_LEVEL.

#define EXAHYPE_PATCH_SIZE_X \
  1  //!< Number of patches in each coordinate direction (+ two ghost cells).
#define EXAHYPE_PATCH_SIZE_Y \
  1  //!< Number of patches in each coordinate direction (+ two ghost cells).
//#define EXAHYPE_PATCH_SIZE_Z                                1    //!< Number
//of patches in each coordinate direction (+ two ghost cells/0 ghost cells in
//2D).
#define EXAHYPE_PATCH_SIZE_TOTAL \
  (EXAHYPE_PATCH_SIZE_X + 2) * (EXAHYPE_PATCH_SIZE_Y + 2)
#define EXAHYPE_PATCH_SIZE_TOTAL_TIMES_DIMENSIONS_TIMES_TWO \
  EXAHYPE_PATCH_SIZE_TOTAL* DIMENSIONS_TIMES_TWO
///@}
/*
#define EXAHYPE_INITIAL_GLOBAL_REFINEMENT_LEVEL             5  //!< The
refinement level of the initial grid.
#define EXAHYPE_INITIAL_ADAPTIVE_REFINEMENT_LEVEL           0  //!< The initial
level of the adaptive refinement.
#define EXAHYPE_MAXIMUM_ADAPTIVE_REFINEMENT_LEVEL           0  //!< The maximum
level of adaptive refinement.
#define EXAHYPE_MAXIMUM_REFINEMENT_LEVEL                    3  //!< The sum of
::EXAHYPE_INITIAL_GLOBAL_REFINEMENT_LEVEL and
::EXAHYPE_MAXIMUM_ADAPTIVE_REFINEMENT_LEVEL.

#define EXAHYPE_PATCH_SIZE_X                                3 //!< Number of
patches in each coordinate direction (+ two ghost cells).
#define EXAHYPE_PATCH_SIZE_Y                                3 //!< Number of
patches in each coordinate direction (+ two ghost cells).
//#define EXAHYPE_PATCH_SIZE_Z                                1    //!< Number
of patches in each coordinate direction (+ two ghost cells/0 ghost cells in 2D).
#define EXAHYPE_PATCH_SIZE_TOTAL
(EXAHYPE_PATCH_SIZE_X+2) * (EXAHYPE_PATCH_SIZE_Y+2)
#define EXAHYPE_PATCH_SIZE_TOTAL_TIMES_DIMENSIONS_TIMES_TWO
EXAHYPE_PATCH_SIZE_TOTAL * DIMENSIONS_TIMES_TWO
#endif
*/

/**
 * Microarchitecture-specific alignment.
 *   Architecture   Alignment   Flag        Padding
 *   wsm            16          SSE3        2 (2 doubles, 4 floats)
 *   snb            32          AVX         4 (4 doubles, 8 floats)
 *   hsw            32          AVX2        4
 *   knl            64          AVX512F     8
 */
#if defined(__AVX512F__)
#define ALIGNMENT 64  // 512 bits, 64 byte aligned
#define SIMD_SIZE 8   // 8 doubles at a time can be processed
#elif defined(__AVX__) || defined(__AVX2__)
#define ALIGNMENT 32  // 256 bits, 32 byte aligned
#define SIMD_SIZE 4   // 4 doubles at a time
#elif defined(__SSE3__)
#define ALIGNMENT 16  // 16 byte aligned
#define SIMD_SIZE 2   // 2 doubles/4 floats at a time
#else                 // fallback solution
#define ALIGNMENT 16  // some alignment doesn't harm, 8 is also feasible
#define SIMD_SIZE 1   // no padding
#endif

/**
 *  Computes the eventual size of an array after padding has been applied. The
 * size
 *  of the vectors are known at compile time; hence we compute the size
 * incorporating
 *  padding also at compile time.
 *  \param[in]  vectorSize  the size of the vector in its bare form, i.e.
 * without padding
 *  \returns    the original vector size plus padding
 */
constexpr int sizeWithPadding(const int vectorSize) {
  return SIMD_SIZE * ((vectorSize + (SIMD_SIZE - 1)) / SIMD_SIZE);
}

/**
 * Simulation runtime parameters.
 */
///@{
#define EXAHYPE_SIMULATION_TIME 4  //!< The total simulation time.
                                   ///@}

/**
 * PDE parameters.
 */
///@{
#define EXAHYPE_NPROBLEMS \
  1  //!< The number of different problems we want to solve at once.

#define EXAHYPE_NVARS 5
#define EXAHYPE_ORDER 3

namespace exahype {
constexpr static const double nvars[1] = {5};
constexpr static const double order[1] = {3};

constexpr static const double patchSizeX[1] = {12};
constexpr static const double patchSizeY[1] = {12};
constexpr static const double patchSizeZ[1] = {1};
}  // namespace exahype

#define EXAHYPE_NBASIS_POWER_DIMENSIONS 16

#define EXAHYPE_CFL_FACTOR 0.9
///@}

/**
 * Visualization parameters.
 */
///@{
#define EXAHYPE_PLOTTING_STRIDE \
  5  //!< The number of snapshots of the solution we want to create during the
     //simulation.
///@}

#endif /* EXAHYPE_CONSTANTS_H_ */
