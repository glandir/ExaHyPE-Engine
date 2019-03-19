#ifndef __PotentialEulerSolver_CLASS_HEADER__
#define __PotentialEulerSolver_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "exahype/parser/ParserView.h"


#include <ostream>

#include "AbstractPotentialEulerSolver.h"


/**
 * We use Peano's logging
 */
#include "tarch/logging/Log.h"


namespace PotentialHydro{
  class PotentialEulerSolver;
}

class PotentialHydro::PotentialEulerSolver : public PotentialHydro::AbstractPotentialEulerSolver {
  private:
    /**
     * Log device
     */
    static tarch::logging::Log _log;
  public:
    PotentialEulerSolver(const double maximumMeshSize,const exahype::solvers::Solver::TimeStepping timeStepping);
    
    /**
     * Initialise the solver.
     *
     * \param[in] cmdlineargs the command line arguments.
     */
    void init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) final override;

    /**
     * @see FiniteVolumesSolver
     */    
    bool useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const override;
    
    /**
     * @see FiniteVolumesSolver
     */    
    void adjustSolution(const double* const x,const double w,const double t,const double dt, double* const Q) override; 
    
    /**
     * Compute the flux tensor.
     *
     * \param[in]    Q the conserved variables (and parameters) associated with a quadrature point
     *                 as C array (already allocated).
     * \param[inout] F the fluxes at that point as C array (already allocated).
     */
    void flux(const double* const Q,double** const F);
    
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
     * \param[inout] QOut      the conserved variables at point x from outside of the domain
     *                         and time-averaged (over [t,t+dt]) as C array (already allocated).
     */
    void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut);
    
    /** Has currently no effect for the Finite Volumes Solver. */
    exahype::solvers::Solver::RefinementControl refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) override;
    
    
    /*** Implement a source potential as nonconservative product ***/
    //bool useAlgebraicSource() const override { return true; } // <-- ADERDG
    //     virtual bool useSource()                 const {return false;} // <-- FV
    //void algebraicSource(const double* const Q, double* const S) override;

    //bool useNonConservativeProduct() const override { return true; }
    void nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) override;
};


#endif // __PotentialEulerSolver_CLASS_HEADER__
