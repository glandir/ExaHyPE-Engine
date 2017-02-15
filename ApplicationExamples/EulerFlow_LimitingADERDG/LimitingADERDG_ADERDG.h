#ifndef __LimitingADERDG_ADERDG_CLASS_HEADER__
#define __LimitingADERDG_ADERDG_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include <ostream>

#include "AbstractLimitingADERDG_ADERDG.h"




namespace Euler{
  class LimitingADERDG_ADERDG;
}

class Euler::LimitingADERDG_ADERDG: public Euler::AbstractLimitingADERDG_ADERDG {
  public:
    LimitingADERDG_ADERDG(double maximumMeshSize,exahype::solvers::Solver::TimeStepping timeStepping,std::vector<std::string>& cmdlineargs);

    /**
     * Initialise the solver.
     *
     * \param[in] cmdlineargs the command line arguments.
     */
    void init(std::vector<std::string>& cmdlineargs);
    
    /**
     * Check if we need to adjust the conserved variables and parameters (together: Q) in a cell
     * within the time interval [t,t+dt].
     *
     * \note Use this function and ::adjustedSolutionValues to set initial conditions.
     *
     * \param[in]    centre    The centre of the cell.
     * \param[in]    dx        The extent of the cell.
     * \param[in]    t         the start of the time interval.
     * \param[in]    dt        the width of the time interval.
     * \return true if the solution has to be adjusted.
     */
    bool hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& centre,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) override;
    
    /**
     * Adjust the conserved variables and parameters (together: Q) at a given time t at the (quadrature) point x.
     *
     * \note Use this function and ::hasToAdjustSolution to set initial conditions.
     *
     * \param[in]    x         the physical coordinate on the face.
     * \param[in]    w         (deprecated) the quadrature weight corresponding to the quadrature point w.
     * \param[in]    t         the start of the time interval.
     * \param[in]    dt        the width of the time interval.
     * \param[inout] Q         the conserved variables (and parameters) associated with a quadrature point
     *                         as C array (already allocated).
     */
    void adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q);
    
    /**
     * Compute the flux tensor.
     *
     * \param[in]    Q the conserved variables (and parameters) associated with a quadrature point
     *                 as C array (already allocated).
     * \param[inout] F the fluxes at that point as C array (already allocated).
     */
    void flux(const double* const Q,double** F);
    
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
     * Compute the source.
     *
     * \param[in]    Q the conserved variables (and parameters) associated with a quadrature point
     *                 as C array (already allocated).
     * \param[inout] S the source point as C array (already allocated).
     */
    void source(const double* const Q,double* S);
    
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
    void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int dir,const double * const fluxIn,const double* const stateIn,double *fluxOut,double* stateOut);
    
    void ncp(const double* const Q,const double* const gradQ,double* BgradQ);
    void matrixb(const double* const Q,const int d,double* Bn);
    bool isDummyKRequired() const override;
    void dummyK_Value(const double* const x,const double t,const double dt, double* forceVector, double* x0); //TODO KD
    
    /**
     * Evaluate the refinement criterion within a cell.
     *
     * \note Instead of a variables array at a single quadrature point we give
     * you all NumberOfVariables*(Order+1)^DIMENSIONS solution degrees of freedom.
     *
     * \note Use this function and ::adjustedSolutionValues to set initial conditions.
     *
     * \param[in]    centre    The centre of the cell.
     * \param[in]    dx        The extent of the cell.
     * \param[in]    t         the start of the time interval.
     * \param[in]    dt        the width of the time interval.
     * \return One of exahype::solvers::Solver::RefinementControl::{Erase,Keep,Refine}.
     */
    exahype::solvers::Solver::RefinementControl refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& centre,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) override;
    
    /**
     * Evaluate the physical admissibility detection (PAD) criterion within a cell.
     *
     *
     * \param[in] QMin the minimum value per variable.
     * \param[in] QMax the maximum value per variable.
     * \return true if the bounds are physically admissible.
     */
    bool physicalAdmissibilityDetection(const double* const QMin,const double* const QMax) override;
};

#endif // __LimitingADERDG_ADERDG_CLASS_HEADER__