#ifndef __DIMSolver_ADERDG_CLASS_HEADER__
#define __DIMSolver_ADERDG_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include <ostream>

#include "AbstractDIMSolver_ADERDG.h"

/**
 * We use Peano's logging
 */
#include "tarch/logging/Log.h"

namespace DIM{
  class DIMSolver_ADERDG;
}

class DIM::DIMSolver_ADERDG : public DIM::AbstractDIMSolver_ADERDG {
  private:
    /**
     * Log device
     */
    static tarch::logging::Log _log;
  public:
    DIMSolver_ADERDG(const double maximumMeshSize,const int maximumMeshDepth,const int haloCells,const int regularisedFineGridLevels,const exahype::solvers::Solver::TimeStepping timeStepping,const int limiterHelperLayers,const int DMPObservables);

    /**
     * Initialise the solver.
     *
     * \param[in] cmdlineargs the command line arguments.
     */
    void init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) final override;

    /**
     * Adjust the conserved variables and parameters (together: Q) at a given time t at the (quadrature) point x.
     *
     * \note Please overwrite function adjustSolution(...) if you want to
     * adjust the solution degrees of freedom in a cellwise manner.
     *
     * \param[in]    x         the physical coordinate on the face.
     * \param[in]    t         the start of the time interval.
     * \param[in]    dt        the width of the time interval.
     * \param[inout] Q         the conserved variables (and parameters) associated with a quadrature point
     *                         as C array (already allocated).
     */
    void adjustPointSolution(const double* const x,const double t,const double dt,double* Q) final override;

    /**
     * Compute the eigenvalues of the flux tensor per coordinate direction \p d.
     *
     * \param[in] Q  the conserved variables associated with a quadrature node
     *               as C array (already allocated).
     * \param[in] d  the column of the flux vector (d=0,1,...,DIMENSIONS).
     * \param[inout] lambda the eigenvalues as C array (already allocated).
     */
    void eigenvalues(const double* const Q,const int d,double* lambda) final override;
    
    /**
     * Impose boundary conditions at a point on a boundary face
     * within the time interval [t,t+dt].
     *
     * \param[in]    x         the physical coordinate on the face.
     * \param[in]    t         the start of the time interval.
     * \param[in]    dt        the width of the time interval.
     * \param[in]    faceIndex indexing of the face (0 -- {x[0]=xmin}, 1 -- {x[1]=xmax}, 2 -- {x[1]=ymin}, 3 -- {x[2]=ymax}, and so on,
     *                         where xmin,xmax,ymin,ymax are the bounds of the cell containing point x.
     * \param[in]    d         the coordinate direction the face normal is pointing to.
     * \param[in]    QIn       the conserved variables at point x from inside of the domain
     *                         and time-averaged (over [t,t+dt]) as C array (already allocated).
     * \param[in]    FIn       the normal fluxes at point x from inside of the domain
     *                         and time-averaged (over [t,t+dt]) as C array (already allocated).
     * \param[inout] QOut      the conserved variables at point x from outside of the domain
     *                         and time-averaged (over [t,t+dt]) as C array (already allocated).
     * \param[inout] FOut      the normal fluxes at point x from outside of the domain
     *                         and time-averaged (over [t,t+dt]) as C array (already allocated).
     */
    void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double * const fluxIn,const double* const stateIn,double *fluxOut,double* stateOut) final override;
    
    /**
     * Evaluate the refinement criterion within a cell.
     *
     * \note Instead of a variables array at a single quadrature point we give
     * you all NumberOfVariables*(Order+1)^DIMENSIONS solution degrees of freedom.
     *
     * \note Use this function and ::adjustSolution to set initial conditions.
     *
     * \param[in]    centre    The centre of the cell.
     * \param[in]    dx        The extent of the cell.
     * \param[in]    t         the start of the time interval.
     * \param[in]    dt        the width of the time interval.
     * \return One of exahype::solvers::Solver::RefinementControl::{Erase,Keep,Refine}.
     */
    exahype::solvers::Solver::RefinementControl refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& centre,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) override;
    
    //PDE

    /**
     * Compute the flux tensor.
     *
     * \param[in]    Q the conserved variables (and parameters) associated with a quadrature point
     *                 as C array (already allocated).
     * \param[inout] F the fluxes at that point as C array (already allocated).
     */
    void flux(const double* const Q,double** F) final override;


    /**
     * Compute the nonconservative term $B(Q) \nabla Q$.
     * 
     * This function shall return a vector BgradQ which holds the result
     * of the full term. To do so, it gets the vector Q and the matrix
     * gradQ which holds the derivative of Q in each spatial direction.
     * Currently, the gradQ is a continous storage and users can use the
     * kernels::idx2 class in order to compute the positions inside gradQ.
     *
     * @TODO: Check if the following is still right:
     * 
     * !!! Warning: BgradQ is a vector of size NumberOfVariables if you
     * use the ADER-DG kernels for nonlinear PDEs. If you use
     * the kernels for linear PDEs, it is a tensor with dimensions
     * Dim x NumberOfVariables.
     * 
     * \param[in]   Q   the vector of unknowns at the given position
     * \param[in]   gradQ   the gradients of the vector of unknowns,
     *                  stored in a linearized array.
     * \param[inout]  The vector BgradQ (extends nVar), already allocated. 
     *
     **/
    void nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) final override;

    void mapDiscreteMaximumPrincipleObservables(double* observables, const double* const Q) const override;

    bool isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const override;

};

#endif // __DIMSolver_ADERDG_CLASS_HEADER__
