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
#define EXAHYPE_INITIAL_GLOBAL_REFINEMENT_LEVEL             2  //!< The refinement level of the initial grid.
#define EXAHYPE_INITIAL_ADAPTIVE_REFINEMENT_LEVEL           0  //!< The initial level of the adaptive refinement.
#define EXAHYPE_MAXIMUM_ADAPTIVE_REFINEMENT_LEVEL           0  //!< The maximum level of adaptive refinement.
#define EXAHYPE_MAXIMUM_REFINEMENT_LEVEL                    2  //!< The sum of ::EXAHYPE_INITIAL_GLOBAL_REFINEMENT_LEVEL and ::EXAHYPE_MAXIMUM_ADAPTIVE_REFINEMENT_LEVEL.

#define EXAHYPE_PATCH_SIZE_X                                12 //!< Number of patches in each coordinate direction (+ two ghost cells).
#define EXAHYPE_PATCH_SIZE_Y                                12 //!< Number of patches in each coordinate direction (+ two ghost cells).
//#define EXAHYPE_PATCH_SIZE_Z                                1    //!< Number of patches in each coordinate direction (+ two ghost cells/0 ghost cells in 2D).
#define EXAHYPE_PATCH_SIZE_TOTAL                            (EXAHYPE_PATCH_SIZE_X+2) * (EXAHYPE_PATCH_SIZE_Y+2)
#define EXAHYPE_PATCH_SIZE_TOTAL_TIMES_DIMENSIONS_TIMES_TWO EXAHYPE_PATCH_SIZE_TOTAL * DIMENSIONS_TIMES_TWO
///@}

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
#define EXAHYPE_NUMBER_OF_PDES 1  //!< The number of snapshots of the solution we want to write to disk during the simulation.

namespace exahype {
static const int numberOfVariables[EXAHYPE_NUMBER_OF_PDES] = { //!< The number of variables for each PDE.
    5
};
static const int order[EXAHYPE_NUMBER_OF_PDES] = {             //!< The order of approximation for each PDE.
    3
};
}

#define EXAHYPE_CFL_FACTOR 0.9
///@}


/**
 * Visualization parameters.
 */
///@{
#define EXAHYPE_NUMBER_OF_PLOTS 10 //!< The number of snapshots of the solution we want to create during the simulation.
///@}

#endif /* EXAHYPE_CONSTANTS_H_ */
