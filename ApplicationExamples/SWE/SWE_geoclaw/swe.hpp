#pragma once

/// \file swe.hpp
///
/// This file was created in an effort to extract the riemann solvers into pure
/// functions for ease of understanding and testing.
///
/// Given the quantities Q in two neighboring cells (here left and right),
/// the riemann solver calculates the resulting fluxes from the left and right
/// cells across the shared wall.
/// TODO: Must these not always be equal and opposite?.
///
/// Q is defined as the four numbers:
///
/// h:  The height of water column (distance from the seabed to the surface).
///     TODO: Is this correct? It seems to be, based on descriptions like "water
///     depth" used in ExaHyPE documentation and the calculations the riemann
///     solver does, but different formulations exist.
///
/// hu: The momentum in the x direction.
///
/// hv: The momentum in the y direction.
///
/// b:  The Absolute height of the sea bed.

#define SUPPRESS_SOLVER_DEBUG_OUTPUT
#include "extern/AugRie.hpp"
#include <cmath>
#include <ostream>

extern "C" {

enum SolverType {
  GEOCLAW_FWAVE = 1,
  GEOCLAW_SSQ_FWAVE,
  GEOCLAW_AUG_RIEMANN,
};

void c_bind_geoclaw_solver_dp(int* i_solver, int* i_maxIter,
                              int* i_numberOfFWaves, double* i_hL, double* i_hR,
                              double* i_huL, double* i_huR, double* i_hvL,
                              double* i_hvR, double* i_bL, double* i_bR,
                              double* i_dryTol, double* i_g,
                              double* o_netUpdatesLeft,
                              double* o_netUpdatesRight,
                              double* o_maxWaveSpeed);

}  // extern "C"

namespace swe {

/// Helper to print a 4d vector behind a `double*`.
struct PrintVec4 {
  const double* v;

  friend std::ostream& operator<<(std::ostream& os, const PrintVec4& p) {
    return os << "{" << p.v[0] << ",\t" << p.v[1] << ",\t" << p.v[2] << ",\t"
              << p.v[3] << "}";
  }
};

/// eigenvalues: h, hu -> lambda[4]
/// \param Q The values of Q in a given cell; only `h, hu` matter.
/// \param dIndex Direction of `u`.
/// \param[out] lambda Eigenvector (`float[4]`).
inline auto eigenvalues(const double* const Q, const int dIndex,
                        double* const lambda, double grav, double epsilon)
    -> void {
  // Dimensions             = 2
  // Number of variables    = 4 + #parameters
  /// Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0

  const double c = std::sqrt(grav * Q[0]);
  const double ih = 1. / Q[0];
  double u_n = Q[dIndex + 1] * ih;

  if (Q[0] < epsilon) {
    lambda[0] = 0.0;
    lambda[1] = 0.0;
    lambda[2] = 0.0;
    lambda[3] = 0.0;
    //    std::cout << 0.0 << std::endl;
  } else {
    lambda[0] = u_n + c;
    lambda[1] = u_n - c;
    lambda[2] = u_n;
    lambda[3] = 0.0;
    //    std::cout << eigs.h() + std::abs(c) << std::endl;
  }
}

/// \param[in] Q The values of Q in a given cell.
/// \param[out] F The resulting fluxes in that cell for both directions
/// (`{.flux_x = {dh, dhu, dhv, db}, .flux_y = {dh, dhu, dhv, db}}`).
inline auto flux(const double* Q, double** F, double epsilon, double grav)
    -> void {
  const double ih = 1. / Q[0];

  double* f = F[0];
  double* g = F[1];

  if (Q[0] < epsilon) {
    f[0] = 0.0;
    f[1] = 0.0;
    f[2] = 0.0;
    f[3] = 0.0;

    g[0] = 0.0;
    g[1] = 0.0;
    g[2] = 0.0;
    g[3] = 0.0;
  } else {
    f[0] = Q[1];
    f[1] = Q[1] * Q[1] * ih + 0.5 * grav * Q[0] * Q[0];
    f[2] = Q[1] * Q[2] * ih;
    f[3] = 0.0;

    g[0] = Q[2];
    g[1] = Q[1] * Q[2] * ih;
    g[2] = Q[2] * Q[2] * ih + 0.5 * grav * Q[0] * Q[0];
    g[3] = 0.0;
  }
}

/// \param[out] fL Fluxes from the left cell at the cell wall.
/// \param[out] fR Fluxes from the right cell at the cell wall.
/// \param qL Quantities in the left cell (`{h, hu, hv, b}`).
/// \param qR Quantities in the right cell (`{h, hu, hv, b}`).
/// \param direction Direction to consider: 0 => use u, 1 => use v.
///
/// \return The maximum speed (exact definition?); used to maintain CFL
/// condition.
inline auto originalRiemannSolver(double* fL, double* fR, const double* qL,
                                  const double* qR, int direction, double grav,
                                  double epsilon) -> double {
  constexpr auto NumberOfVariables = 4;
  constexpr auto Dimensions = 2;
  /// Eigenvalues for left cell.
  double LL[NumberOfVariables] = {0.0};
  /// Eigenvalues for right cell.
  double LR[NumberOfVariables] = {0.0};

  eigenvalues(qL, direction, LL, grav, epsilon);
  eigenvalues(qR, direction, LR, grav, epsilon);

  /// Maximum speed, somehow determined by eigenvalues.
  /// Used to maintain CFL condition.
  double smax = 0.0;
  for (int i = 0; i < NumberOfVariables; i++) {
    const double abs_sL_i = std::abs(LL[i]);
    smax = std::max(abs_sL_i, smax);
  }
  for (int i = 0; i < NumberOfVariables; i++) {
    const double abs_sR_i = std::abs(LR[i]);
    smax = std::max(abs_sR_i, smax);
  }

  double FL2[Dimensions][NumberOfVariables] = {{0.0}};
  double FR2[Dimensions][NumberOfVariables] = {{0.0}};
  /// "Internal" fluxes in left cell (`double[2][4]`).
  double* FL[Dimensions] = {FL2[0], FL2[1]};
  /// "Internal" fluxes in right cell (`double[2][4]`).
  double* FR[Dimensions] = {FR2[0], FR2[1]};
  flux(qL, FL, epsilon, grav);
  flux(qR, FR, epsilon, grav);

  /// Flux across cell wall (`double[4]`).
  double flux[NumberOfVariables] = {
      0.5 * (FL[direction][0] + FR[direction][0]) -
          0.5 * smax * (qR[0] + qR[3] - qL[0] - qL[3]),
      0.5 * (FL[direction][1] + FR[direction][1]) -
          0.5 * smax * (qR[1] - qL[1]),
      0.5 * (FL[direction][2] + FR[direction][2]) -
          0.5 * smax * (qR[2] - qL[2]),
      0.5 * (FL[direction][3] + FR[direction][3]),
  };

  /// Average water column height (at cell wall).
  const double hRoe = 0.5 * (qL[0] + qR[0]);

  /// Maximum bottom height of the two cells.
  const double bm = std::max(qL[3], qR[3]);

  /// ?
  const double Deta =
      std::max(qR[0] + qR[3] - bm, 0.0) - std::max(qL[0] + qL[3] - bm, 0.0);

  /// Flux due to height differential and gravity; only nonzero for hu.
  double djump[NumberOfVariables] = {0.0};

  djump[direction + 1] = 0.5 * grav * hRoe * Deta;

  flux[0] = 0.5 * (FL[direction][0] + FR[direction][0]) - 0.5 * smax * Deta;
  for (int i = 0; i < NumberOfVariables; i++) {
    fL[i] = flux[i] + djump[i];
    fR[i] = flux[i] - djump[i];
  }

  return smax;
}

/// \param[out] fL Fluxes from the left cell at the cell wall.
/// \param[out] fR Fluxes from the right cell at the cell wall.
/// \param qL Quantities in the left cell (`{h, hu, hv, b}`).
/// \param qR Quantities in the right cell (`{h, hu, hv, b}`).
/// \param direction Direction to consider: 0 => use u, 1 => use v.
///
/// \return The maximum speed (exact definition?); used to maintain CFL
/// condition.
inline auto riemannSolver(double* fL, double* fR, const double* qL,
                          const double* qR, int direction, double grav,
                          double epsilon) -> double {
  const double& i_hLeft = qL[0];
  const double& i_hRight = qR[0];
  const double& i_huLeft = qL[direction + 1];
  const double& i_huRight = qR[direction + 1];
  const double& i_bLeft = qL[3];
  const double& i_bRight = qR[3];

  double& o_hUpdateLeft = fL[0];
  double& o_hUpdateRight = fR[0];
  double& o_huUpdateLeft = fL[direction + 1];
  double& o_huUpdateRight = fR[direction + 1];
  double o_maxWaveSpeed;

  const auto dryTolerance = epsilon;
  const auto newtonTolerance = 1e-6;
  const auto newtonIterations = 10;
  const auto zeroTolerance = std::numeric_limits<double>::epsilon();

  solver::AugRie<double> solver(dryTolerance, grav, newtonTolerance,
                                newtonIterations, zeroTolerance);
  // clang-format off
  solver.computeNetUpdates(
            i_hLeft,
            i_hRight,
            i_huLeft,
            i_huRight,
            i_bLeft,
            i_bRight,

            o_hUpdateLeft,
            o_hUpdateRight,
            o_huUpdateLeft,
            o_huUpdateRight,
            o_maxWaveSpeed);
  // clang-format on

  // godunov.ccph expects fluxes directed towards the wall instead of in x
  // direction, unfortunately.
  // fR[0] = -fR[0];
  // fR[1 + direction] = -fR[1 + direction];

  // Fluxes for bathymetry and momentum in other direction are always zero.
  // fL[2 - direction] = 0;
  // fL[3] = 0;
  // fR[2 - direction] = 0;
  // fR[3] = 0;

  return o_maxWaveSpeed;
}

inline auto samoaRiemannSolver(double* fL, double* fR, const double* qL,
                               const double* qR, int direction, double grav,
                               double epsilon) -> double {
  int solverType = GEOCLAW_AUG_RIEMANN;
  int maxIter = 10;
  int numberOfFWaves = 3;
  double dryTolerance = epsilon;
  // double o_netUpdatesLeft[3];
  // double o_netUpdatesRight[3];
  double* o_netUpdatesLeft = fL;
  double* o_netUpdatesRight = fR;
  double o_maxWaveSpeed;

  int x_direction = 1 + direction;
  int y_direction = 2 - direction;

  double qR_[4] = {qR[0], qR[1], qR[2], qR[3]};
  double qL_[4] = {qL[0], qL[1], qL[2], qL[3]};

  c_bind_geoclaw_solver_dp(&solverType, &maxIter, &numberOfFWaves, qL_ + 0,
                           qR_ + 0, qL_ + x_direction, qR_ + x_direction,
                           qL_ + y_direction, qR_ + y_direction, qL_ + 3,
                           qR_ + 3, &dryTolerance, &grav, o_netUpdatesLeft,
                           o_netUpdatesRight, &o_maxWaveSpeed);
  // godunov.ccph expects fluxes directed towards the wall instead of in x
  // direction, unfortunately.
  // fR[0] = -fR[0];
  // fR[x_direction] = -fR[x_direction];

  // Fluxes for bathymetry and momentum in other direction are always zero.
  // fL[y_direction] = 0;
  // fL[3] = 0;
  // fR[y_direction] = 0;
  // fR[3] = 0;
  // TODO
  return o_maxWaveSpeed;
}

}  // namespace swe
