#ifndef __MyElastodynamicsSolver_CLASS_HEADER__
#define __MyElastodynamicsSolver_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include <ostream>

#include "AbstractMyElastodynamicsSolver.h"


/**
 * We use Peano's logging
 */
#include "tarch/logging/Log.h"





namespace Elastodynamics{
  class MyElastodynamicsSolver;
}

class Elastodynamics::MyElastodynamicsSolver: public Elastodynamics::AbstractMyElastodynamicsSolver {
  private:
    /**
     * Log device
     */
    static tarch::logging::Log _log;
  public:
    MyElastodynamicsSolver(const double maximumMeshSize,const int maximumMeshDepth,const int haloCells,const int regularisedFineGridLevels,const exahype::solvers::Solver::TimeStepping timeStepping,const int DMPObservables);

    /**
     * Initialise the solver.
     *
     * \param[in] cmdlineargs the command line arguments.
     */
    void init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) final override;
    
   
    
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
    void adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) override;
    
    /**
     * Compute the flux tensor.
     *
     * \param[in]    Q the conserved variables (and parameters) associated with a quadrature point
     *                 as C array (already allocated).
     * \param[inout] F the fluxes at that point as C array (already allocated).
     */
    virtual void flux(const double* const Q,double** const F);
    
    /**
     * Compute the eigenvalues of the flux tensor per coordinate direction \p d.
     *
     * \param[in] Q  the conserved variables associated with a quadrature node
     *               as C array (already allocated).
     * \param[in] d  the column of the flux vector (d=0,1,...,DIMENSIONS).
     * \param[inout] lambda the eigenvalues as C array (already allocated).
     */
    void eigenvalues(const double* const Q,const int d,double* const lambda);
    
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
    void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut);
    
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
    exahype::solvers::Solver::RefinementControl refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& centre,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) override;

    virtual void nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ);

    virtual void coefficientMatrix(const double* const Q,const int d,double* const Bn);
    virtual void pointSource(const double* const x,const double t,const double dt, double* const forceVector, double* const x0);

    void riemannSolver(double* const FL,double* const FR,const double* const QL,const double* const QR,double* const tempFaceUnknownsArray,double** const tempStateSizedVectors,double** const tempStateSizedSquareMatrices,const double dt,const int normalNonZeroIndex,bool isBounaryFace) override;

    void riemannSolver_Nodal(double v_p,double v_m, double sigma_p, double sigma_m, double z_p , double z_m, double& v_hat_p , double& v_hat_m, double& sigma_hat_p, double& sigma_hat_m);
    void localBasis(double* const n, double* const m, double* const l, int d);
    void Gram_Schmidt(double* const y, double* const z);
    void riemannSolver_BC0(double v, double sigma, double z, double r, double& v_hat, double& sigma_hat);
    void riemannSolver_BCn(double v, double sigma, double z, double r, double& v_hat, double& sigma_hat);
};

#endif // __MyElastodynamicsSolver_CLASS_HEADER__
