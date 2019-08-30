#ifndef __MyElasticWaveSolver_CLASS_HEADER__
#define __MyElasticWaveSolver_CLASS_HEADER__

// This file was initially generated by the ExaHyPE toolkit.
// You can modify it in order to extend your solver with features.
// Whenever this file is present, a re-run of the ExaHyPE toolkit will
// not overwrite it. Delete it to get it regenerated.
//
// ========================
//   www.exahype.eu
// ========================

#include <ostream>

#include "AbstractMyElasticWaveSolver.h"
#include "exahype/parser/ParserView.h"

/**
 * We use Peano's logging
 */
#include "tarch/logging/Log.h"
#include "tarch/la/Vector.h"

namespace Elastic{
  class MyElasticWaveSolver;
}

class Elastic::MyElasticWaveSolver : public Elastic::AbstractMyElasticWaveSolver {
  private:
    /**
     * Log device
     */
    static tarch::logging::Log _log;
  public:
    MyElasticWaveSolver(const double maximumMeshSize,const int maximumMeshDepth,const int haloCells,const int haloBufferCells,const int limiterBufferCells,const int regularisedFineGridLevels,const exahype::solvers::Solver::TimeStepping timeStepping,const int DMPObservables);

    /**
     * Initialise the solver.
     *
     * \param[in] cmdlineargs the command line arguments.
     * \param[in] constants   access to the constants specified for the solver.
     */
    void init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) final override;

    /**
     * Patchwise adjust
     * @TODO LR : Document
     */
    void adjustSolution(double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) final override;

    /**
     * Compute the eigenvalues of the flux tensor per coordinate direction \p d.
     *
     * \param[in] Q  the conserved variables associated with a quadrature node
     *               as C array (already allocated).
     * \param[in] d  the column of the flux vector (d=0,1,...,DIMENSIONS).
     * \param[inout] lambda the eigenvalues as C array (already allocated).
     */
    void eigenvalues(const double* const Q,const int d,double* const lambda) final override;
    
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
  void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) override;
    
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
    
    //PDE

    /**
     * Compute the flux tensor.
     *
     * \param[in]    Q the conserved variables (and parameters) associated with a quadrature point
     *                 as C array (already allocated).
     * \param[inout] F the fluxes at that point as C array (already allocated).
     */
    void flux(const double* const Q,double** const F) final override;

/* algebraicSource() function not included, as requested by the specification file */

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
    void nonConservativeProduct(const double* const Q,const double* const * const gradQ,double** const BgradQ) final override;

    /**
     * Compute a pointSource contribution.
     * 
     * @TODO: Document me, please.
    **/
     void pointSource(const double* const Q,const double* const x,const double t,const double dt, double* const forceVector,int n) override;

    /**
     * @TODO LR : document
     */
    void multiplyMaterialParameterMatrix(const double* const Q, double** const rhs) final override;

  void riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,const double t, const double dt, const tarch::la::Vector<DIMENSIONS, double>& lengthScale, const int direction, bool isBoundaryFace, int faceIndex) final override;
    void riemannSolver_Nodal(double v_p,double v_m, double sigma_p, double sigma_m, double z_p , double z_m, double& v_hat_p , double& v_hat_m, double& sigma_hat_p, double& sigma_hat_m);
    void riemannSolver_boundary(int faceIndex,double r, double vn , double vm , double vl, double Tn , double Tm ,double Tl , double zp, double zs,  double& vn_hat , double& vm_hat ,double& vl_hat , double& Tn_hat , double& Tm_hat ,double& Tl_hat);
    
    void riemannSolver_BC0(double v, double sigma, double z,  double r, double& v_hat, double& sigma_hat);
    void riemannSolver_BCn(double v, double sigma, double z,  double r, double& v_hat, double& sigma_hat);

    void localBasis(double* const n, double* const m, double* const l, int d);
    void Gram_Schmidt(double* const y, double* const z);
    
    void generate_fluctuations_right(double z,  double T,double T_hat,double v, double v_hat, double& F);
    void generate_fluctuations_left(double z,  double T,double T_hat,double v, double v_hat, double& F);
    void rotate_into_physical_basis(double* const n,double* const m,double* const l, double Fn,double Fm,double Fl, double& Fx, double& Fy, double& Fz);
    void rotate_into_orthogonal_basis(double* const n,double* const m,double* const l, double Tx,double Ty,double Tz, double& Tn, double& Tm, double& Tl);
    void extract_tractions_and_particle_velocity(double* const n, const double* const Q, double& Tx,double& Ty,double& Tz,double& vx,double& vy,double& vz );


    // Vect PDE
    void flux_vect(const double* const * const restrict Q, double* const * const * const restrict F, const int size);
    void nonConservativeProduct_vect(const double* const * const Q,const double* const * const * const gradQ, double* const * const * const BgradQ, const int size);
    void multiplyMaterialParameterMatrix_vect(const double* const * const Q, double* const * const * const rhs, const int size);
};

#endif // __MyElasticWaveSolver_CLASS_HEADER__
