#ifndef __LimitingADERDG_FV_CLASS_HEADER__
#define __LimitingADERDG_FV_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include <ostream>

#include "AbstractLimitingADERDG_FV.h"

/**
 * We use Peano's logging
 */
#include "tarch/logging/Log.h"


namespace Euler{
  class LimitingADERDG_FV;
}

class Euler::LimitingADERDG_FV : public Euler::AbstractLimitingADERDG_FV {
  private:
    /**
     * Log device
     */
    static tarch::logging::Log _log;
  public:
    LimitingADERDG_FV(double maximumMeshSize,int maximumAdaptiveMeshDepth,exahype::solvers::Solver::TimeStepping timeStepping,std::vector<std::string>& cmdlineargs);
    
    /**
     * Initialise the solver.
     *
     * \param[in] cmdlineargs the command line arguments.
     */
    void init(const std::vector<std::string>& cmdlineargs,const exahype::Parser::ParserView& constants) final override;

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
    
};


#endif // __LimitingADERDG_FV_CLASS_HEADER__
