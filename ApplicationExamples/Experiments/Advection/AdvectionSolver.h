#ifndef __AdvectionSolver_CLASS_HEADER__
#define __AdvectionSolver_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include <ostream>

#include "AbstractAdvectionSolver.h"


/**
 * We use Peano's logging
 */
#include "tarch/logging/Log.h"





namespace Trivial{
  class AdvectionSolver;
}

class Trivial::AdvectionSolver: public Trivial::AbstractAdvectionSolver {
  private:
    /**
     * Log device
     */
    static tarch::logging::Log _log;
  public:
    AdvectionSolver(double maximumMeshSize,int maximumAdaptiveMeshDepth,int DMPObservables,int limiterHelperLayers,exahype::solvers::Solver::TimeStepping timeStepping,std::vector<std::string>& cmdlineargs);

    /**
     * Initialise the solver.
     *
     * \param[in] cmdlineargs the command line arguments.
     */
    void init(const std::vector<std::string>& cmdlineargs,const exahype::Parser::ParserView& constants) final override;
    
   
    
    /**
     * Adjust the conserved variables and parameters (together: Q) at a given time t at the (quadrature) point x.
     *
     *
     * \param[in]    x         the physical coordinate on the face.
     * \param[in]    w         (deprecated) the quadrature weight corresponding to the quadrature point w.
     * \param[in]    t         the start of the time interval.
     * \param[in]    dt        the width of the time interval.
     * \param[inout] Q         the conserved variables (and parameters) associated with a quadrature point
     *                         as C array (already allocated).
     */
    void adjustPointSolution(const double* const x,const double t,const double dt,double* Q) override;
    
    /**
     * Compute the flux tensor.
     *
     * \param[in]    Q the conserved variables (and parameters) associated with a quadrature point
     *                 as C array (already allocated).
     * \param[inout] F the fluxes at that point as C array (already allocated).
     */
    virtual void flux(const double* const Q,double** F);
    
    /**
     * Compute the eigenvalues of the flux tensor per coordinate direction \p d.
     *
     * \param[in] Q  the conserved variables associated with a quadrature node
     *               as C array (already allocated).
     * \param[in] d  the column of the flux vector (d=0,1,...,DIMENSIONS).
     * \param[inout] lambda the eigenvalues as C array (already allocated).
     */
    void eigenvalues(const double* const Q,const int d,double* lambda);
    
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
    void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double * const fluxIn,const double* const stateIn,double *fluxOut,double* stateOut);
    
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
};

#endif // __AdvectionSolver_CLASS_HEADER__
