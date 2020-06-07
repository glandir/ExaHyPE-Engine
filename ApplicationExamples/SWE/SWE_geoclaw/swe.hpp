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

  const double c = std::sqrt(grav * Q[0]);
  double u_n =
      Q[dIndex + 1] * Q[0] * std::sqrt(2) /
      std::sqrt(std::pow(Q[0], 4) + std::pow(std::max(Q[0], epsilon), 4));

  lambda[0] = u_n + c;
  lambda[1] = u_n - c;
  lambda[2] = u_n;
  lambda[3] = 0.0;
}

/// \param[in] Q The values of Q in a given cell.
/// \param[out] F The resulting fluxes in that cell for both directions
/// (`{.flux_x = {dh, dhu, dhv, db}, .flux_y = {dh, dhu, dhv, db}}`).
inline auto flux(const double* Q, double** F, double epsilon) -> void {
  // Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0

  double* const f = F[0];
  double* const g = F[1];

  double u_n =
      Q[1] * Q[0] * std::sqrt(2) /
      std::sqrt(std::pow(Q[0], 4) + std::pow(std::max(Q[0], epsilon), 4));
  double v_n =
      Q[2] * Q[0] * std::sqrt(2) /
      std::sqrt(std::pow(Q[0], 4) + std::pow(std::max(Q[0], epsilon), 4));

  f[0] = Q[0] * u_n;
  f[1] = Q[0] * u_n * u_n;  // 0.5 * grav * Q[0] * Q[0];
  f[2] = Q[0] * u_n * v_n;
  f[3] = 0.0;

  g[0] = Q[0] * v_n;
  g[1] = Q[0] * u_n * v_n;
  g[2] = Q[0] * v_n * v_n;  // 0.5 * grav * Q[0] * Q[0];
  g[3] = 0.0;
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
  flux(qL, FL, epsilon);
  flux(qR, FR, epsilon);

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
  fR[0] = -fR[0];
  fR[1 + direction] = -fR[1 + direction];

  // Fluxes for bathymetry and momentum in other direction are always zero.
  fL[2 - direction] = 0;
  fL[3] = 0;
  fR[2 - direction] = 0;
  fR[3] = 0;

  return o_maxWaveSpeed;
}

}  // namespace swe
