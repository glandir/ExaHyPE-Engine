#ifndef __EulerSolver_FV_CLASS_HEADER__
#define __EulerSolver_FV_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "exahype/parser/ParserView.h"

#include "AbstractEulerSolver_FV.h"

/**
 * We use Peano's logging
 */
#include "tarch/logging/Log.h"

#include <ostream>

namespace Euler{
class EulerSolver_FV;
}

class Euler::EulerSolver_FV : public Euler::AbstractEulerSolver_FV {
private:
  enum class Reference { EntropyWave=0, SodShockTube=1, SphericalExplosion=2, RarefactionWave=3, ShuOsher=4 };
  static Reference ReferenceChoice;

  /**
   * (Smooth solution)
   *
   * Entropy wave is a moving Gaussian matter distribution where it is simple
   * to give an analytic result.
   *
   * See also chapter 7.13.2 in "I do like CFD, VOL.1" by Katate Masatsuka.
   */
  static void entropyWave(const double* const x,double t, double* const Q);

  /**
   * (Discontinuous solution)
   *
   * Analytical solution to Sod's shock tube problem:
   * # reference 1: https://en.wikipedia.org/wiki/Sod_shock_tube
   * # reference 2: https://gitlab.com/fantaz/Riemann_exact/tree/master
   * # reference 3: Toro, E., Riemann Solvers and Numerical Methods for Fluid Dynamics (4.3.3 Numerical tests)
   *
   *   |             |         |       |
   *   |             |         |       |
   *   | rarefaction | contact | shock |
   *___|_____________|_________|_______|_________
   *   x1           x2   x0   x3      x4
   */
  static void sodShockTube(const double* const x,double t, double* const Q);
  
  /**
   * (Discontinuous solution)
   *
   * Copied from [2] and corrected ρ_R:
   *
   * "This test [1] is a 1D hydrodynamic shock tube. 
   * The downstream ﬂow has a sinusoidal density ﬂuctuation ρ = 1 − ε sin(λπ x) with a wave length of λ = 5 and an amplitude of ε = 0.2. 
   * A Mach 3 shock front is initially located at x = −4 on domain [−5; 5]. The left and the right states are given by 
   * ρ_L = 3.857143, u_L = 2.629369, p_L = 10.33333 and ρ_R = 1 + 0.2 sin(5 x), u_R = 0 and p_R = 1. 
   * The ﬁnal time is set to t ﬁnal = 0.18. 
   * This problem involves small scales after the shock has interacted with the sine wave that can be 
   * captured either with a ﬁne enough mesh or with high order accurate method."
   *  
   * [1] Chi-Wang Shu, Stanley Osher, 
   * Efficient implementation of essentially non-oscillatory shock-capturing schemes, II, Journal of Computational Physics, Volume 83, Issue 1, 1989, Pages 32-78, ISSN 0021-9991, 
   * https://doi.org/10.1016/0021-9991(89)90222-2.
   *
   * [2] M. Dumbser, O. Zanotti, R. Loubère, and S. Diot, “A posteriori subcell limiting of the discontinuous Galerkin finite element method for hyperbolic conservation laws,” 
   * Journal of Computational Physics, vol. 278, pp. 47–75, Dec. 2014.
   */
  static void shuOsher(const double* const x, double t, double* const Q);

  /**
   * Spherical explosion.
   */
  static void sphericalExplosion(const double* const x,double t, double* const Q);

  /**
   * (Smooth solution)
   *
   * The rarefaction wave is not a consistent solution of Euler's equations, but instead
   * a perturbation which immediately changes rho, vx, vy and vz. It corresponds
   * to consistent initial data with roughly
   *   rho = exp( - (r - v0*t)**2)
   *   vx  = v0*t * exp(-y**2)  // actually even more complicated, it is more
   *   vy  = v0*t * exp(-x**2)  // a velocity on a kind of ring of radius v0*t
   *   E   = p/(gamma-1) + rho/2 * v**2 = p/(gamma-1) + alpha*exp(- (r-v0*t)**2)
   *   p   = 1
   * However, it is much more complicated to write the closed form solution
   * instead of the perturbation approach. However, the closed form solution allows to
   * specify the solution at any time while while the initial perturbation form
   * does *not* allow to specify the solution.
   */
  static void rarefactionWave(const double* const x,double t, double* const Q);

  /**
   * Log device
   */
  static tarch::logging::Log _log;

public:
  EulerSolver_FV(const double maximumMeshSize,const exahype::solvers::Solver::TimeStepping timeStepping);

  /**
   * Initialise the solver.
   *
   * \param[in] cmdlineargs the command line arguments.
   */
  void init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) final override;

  /**
   * Calls ::sodShockTube if constant 'reference' is set to 'sod'.
   * calls ::entropyWave if constant 'reference is set to 'entropywave'.
   * ErrorWriter, ErrorPlotter write errors of numerical solution for both choices.
   *
   * Calls ::sphericalExplosion if constant 'reference' is set to 'explosion'.
   * calls ::rarefactionWave if constant 'reference is set to 'rarefactionwave'.
   * Computes errors for both choices.
   * ErrorWriter, ErrorPlotter write norms of numerical solution for both choices.
   */
  static void referenceSolution(const double* const x, const double t, double* const Q);

  /**
   * @see FiniteVolumesSolver
   */
  void adjustSolution(const double* const x,const double t,const double dt, double* const Q) override;

  /**
   * Use generalised Osher Solomon flux.
   */
  double riemannSolver(double* const fL, double* const fR, const double* const qL, const double* const qR, int direction)  final override;

  /**
   * Compute the eigenvectors of the jacobian matrix (coefficient matrix of derivative in direction @direction).
   *
   * This function also returns eigenvalues as numerical packages often compute eigenvalues and 
   * eigenvectors together. 
   *
   * see: https://www3.nd.edu/~dbalsara/Numerical-PDE-Course/Appendix_LesHouches/LesHouches_Lecture_5_Approx_RS.pdf
   *
   * @param[in]    Q       state variables and parameters; range: [0,nVar+nData].
   * @param[in]    in      the normal direction; index of the x-axis in reference space.
   * @param[in]    is      index of the y-axis in reference space.
   * @param[in]    it      index of the z-axis in reference space.
   * @param[inout] R       the right eigenvectors (each column stores an eigenvector); initialised to zero; range: [0,nVar]^2.
   * @param[inout] eigvals the right eigenvectors (each column stores an eigenvector); range: [0,nVar].
   * @param[inout] iR      the left eigenvectors  (inverse of RM); initialised to zero; range: [0,nVar]^2.
   */
   void eigenvectors(const double* const Q,const int in,const int is,const int it,double (&RM)[NumberOfVariables][NumberOfVariables],double (&eigvals)[NumberOfVariables], double (&iRM)[NumberOfVariables][NumberOfVariables]);

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
  void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const stateIn,double* const stateOut);
};


#endif // __EulerSolver_FV_CLASS_HEADER__
