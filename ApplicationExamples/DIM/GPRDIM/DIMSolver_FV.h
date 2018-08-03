#ifndef __DIMSolver_FV_CLASS_HEADER__
#define __DIMSolver_FV_CLASS_HEADER__

// This file was initially generated by the ExaHyPE toolkit.
// You can modify it in order to extend your solver with features.
// Whenever this file is present, a re-run of the ExaHyPE toolkit will
// not overwrite it. Delete it to get it regenerated.
//
// ========================
//   www.exahype.eu
// ========================

#include <ostream>

#include "AbstractDIMSolver_FV.h"
#include "exahype/parser/ParserView.h"

/**
 * We use Peano's logging
 */
#include "tarch/logging/Log.h"


namespace GPRDIM{
  class DIMSolver_FV;
}

class GPRDIM::DIMSolver_FV : public GPRDIM::AbstractDIMSolver_FV {
  private:
    /**
     * Log device
     */
    static tarch::logging::Log _log;
  public:
    DIMSolver_FV(const double maximumMeshSize,const exahype::solvers::Solver::TimeStepping timeStepping);
    
    /**
     * Initialise the solver.
     *
     * \param[in] cmdlineargs the command line arguments.
     * \param[in] constants   access to the constants specified for the solver.
     */
    void init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) final override ;

    /**
     * @see FiniteVolumesSolver
     */    
    void adjustSolution(const double* const x,const double t,const double dt, double* Q) override; 
    
    /**
     * Compute the eigenvalues of the flux tensor per coordinate direction \p d.
     *
     * \param[in] Q  the conserved variables associated with a quadrature node
     *               as C array (already allocated).
     * \param[in] d  the column of the flux vector (d=0,1,...,DIMENSIONS).
     * \param[inout] lambda the eigenvalues as C array (already allocated).
     */
    void eigenvalues(const double* const Q,const int d,double* lambda) override;
        
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
     * \param[inout] QOut      the conserved variables at point x from outside of the domain
     *                         and time-averaged (over [t,t+dt]) as C array (already allocated).
     */
    void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const stateIn,double* stateOut) override;
    
    /**
     * Compute the flux tensor.
     *
     * \param[in]    Q the conserved variables (and parameters) associated with a quadrature point
     *                 as C array (already allocated).
     * \param[inout] F the fluxes at that point as C array (already allocated).
     */
    void flux(const double* const Q,double** F) override;


    /**
     * Compute the algebraic source (right hand side).
     *
     * \param[in]   Q      the conserved variables associated with a quadrature node as C array (already allocated).
     * \param[out]  S      the source term vector (already allocated, same shape as Q).
     **/
    void algebraicSource(const double* const Q,double* S) override;

    /**
     * Compute the nonconservative product at a given position.
     *
     * \param[in]   Q       the conserved variables associated with a quadrature node
     *                      as C array (already allocated).
     * \param[in]   gradQ   The linearized gradient matrix of derivatives of Q (see guidebook for details).
     * \param[out]  BgradQ  The nonconservative product (already allocated, same shape as Q).
     *
     **/
    void nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) override;
    
	double riemannSolver(double* fL, double *fR, const double* qL, const double* qR, int normalNonZero) override;
    /* pointSource() function not included, as requested in the specification file */
};


#endif // __DIMSolver_FV_CLASS_HEADER__
