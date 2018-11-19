// This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
// 
// ========================
//   www.exahype.eu
// ========================

#include <cmath>
#include <map>
#include "NavierStokesSolverDG.h"
#include "NavierStokesSolverDG_Variables.h"

#include "totalVariation.h"
#include "stableDiffusiveTimeStepSize.h"
#if DIMENSIONS == 2
#include "diffusiveRiemannSolver2d.h"
#elif DIMENSIONS == 3
#include "diffusiveRiemannSolver3d.h"
#endif

#include "kernels/aderdg/generic/Kernels.h"
#include "kernels/KernelUtils.h"

#include "Scenarios/ScenarioFactory.h"
#include "Scenarios/Atmosphere.h"


tarch::logging::Log NavierStokes::NavierStokesSolverDG::_log( "NavierStokes::NavierStokesSolverDG" );

void NavierStokes::NavierStokesSolverDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  assert(constants.isValueValidString("scenario"));

  double referenceViscosity;
  if (constants.isValueValidString("viscosity") &&
      constants.getValueAsString("viscosity") == "default") {
   throw -1;
  } else {
    assert(constants.isValueValidDouble("viscosity"));
    referenceViscosity = constants.getValueAsDouble("viscosity");
  }

  scenarioName = constants.getValueAsString("scenario");
  scenario = ScenarioFactory::createScenario(scenarioName);

  std::cout << referenceViscosity << " " << scenario->getGasConstant() << std::endl;
  ns = PDE(referenceViscosity, *scenario);
}

void NavierStokes::NavierStokesSolverDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  if (tarch::la::equals(t, 0.0)) {
    Variables vars(Q);
    scenario->initialValues(x, ns, vars);
    for (int i = 0; i < vars.variables(); ++i) {
      assertion2(std::isfinite(Q[i]), i, Q[i]);

    }
  }
}

void NavierStokes::NavierStokesSolverDG::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
   scenario->source(x, t, ns, Q, S);
}

void NavierStokes::NavierStokesSolverDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
					const double * const fluxIn,const double* const stateIn, const double* const gradStateIn,
  double *fluxOut,double* stateOut) {
  constexpr auto basisSize = Order + 1;
  constexpr auto gradSize = NumberOfVariables * DIMENSIONS;

  auto gradStateOut = std::array<double, gradSize>{{0.0}};
  kernels::idx2 idxGradQ(DIMENSIONS,NumberOfVariables);

  std::fill_n(fluxOut, NumberOfVariables, 0.0);
  std::fill_n(stateOut, NumberOfVariables, 0.0);

  double _F[DIMENSIONS][NumberOfVariables]={0.0};
#if DIMENSIONS == 2
  double* F[2] = {_F[0], _F[1]};
#elif DIMENSIONS == 3
  double* F[3] = {_F[0], _F[1], _F[2]};
#endif

  for (int i = 0; i < NumberOfVariables; i++) {
    assertion2(std::isfinite(stateIn[i]), stateIn[i], i);
  }
  assertion1(stateIn[0] > 0, stateIn[0]);

  if (scenario->getBoundaryType(faceIndex) == NavierStokes::BoundaryType::analytical) {
    // Integrate over time.
    auto curStateOut = std::array<double, NumberOfVariables>{0.0};
    Variables curVarsOut(curStateOut.data());
    for (int i = 0; i < basisSize; ++i) {
      // TODO(Lukas): Check if we need to reset this data here.
      std::fill(curStateOut.begin(), curStateOut.end(), 0.0);
      std::fill(gradStateOut.begin(), gradStateOut.end(), 0.0);

      const double weight = kernels::gaussLegendreWeights[Order][i];
      const double xi = kernels::gaussLegendreNodes[Order][i];
      const double ti = t + xi * dt;

      scenario->analyticalSolution(x, ti, ns, curVarsOut, gradStateOut.data());

      //flux(curStateOut.data(), gradStateOut.data(), F);
      ns.evaluateFlux(curStateOut.data(), gradStateOut.data(), F, true);

      for (int j = 0; j < NumberOfVariables; ++j) {
        stateOut[j] += weight * curStateOut[j];
        fluxOut[j] += weight * F[normalNonZero][j];
      }

    }
    return;
  }

  assertion(scenario->getBoundaryType(faceIndex) == BoundaryType::wall ||
                 scenario->getBoundaryType(faceIndex) == BoundaryType::hydrostaticWall);

  // Set no slip wall boundary conditions.
  ReadOnlyVariables varsIn(stateIn);
  Variables varsOut(stateOut);

  // Rho/E extrapolated, velocity mirrored.
  // Leads to zero velocity after Riemann solver.
  std::copy_n(stateIn, NumberOfVariables, stateOut);
  varsOut.j(0) = -varsIn.j(0);
  varsOut.j(1) = -varsIn.j(1);
#if DIMENSIONS == 3
  varsOut.j(2) = -varsIn.j(2);
#endif

  // Extrapolate gradient.
  std::copy_n(gradStateIn, gradSize, gradStateOut.data());

  // We deal with heat conduction by computing the flux at the boundary without heat conduction,
  // To do this, we reconstruct the incoming flux using the extrapolated/time-averaged state/gradient.
  // Note that this incurs an error.
  // The incoming flux is reconstructed in boundaryConditions.

  // Then compute the outgoing flux.
  if (scenario->getBoundaryType(faceIndex) == BoundaryType::hydrostaticWall) {
    // We need to reconstruct the temperature gradient here.
    // TODO(Lukas) Put those constants into the scenario.
    const auto g = 9.81;
    const auto backgroundPotTemperature = 300;
#if DIMENSIONS == 2
    const auto posZ = x[1];
#else
    const auto posZ = x[2];
#endif

    const double equilibriumTemperatureGradient = computeHydrostaticTemperatureGradient(ns, g, posZ,
            backgroundPotTemperature);

    // Use no viscous effects and use equilibrium temperature gradient.
    ns.evaluateFlux(stateOut, gradStateOut.data(), F, false, true, equilibriumTemperatureGradient);
  } else {
    ns.evaluateFlux(stateOut, gradStateOut.data(), F);
  }

  std::copy_n(F[normalNonZero], NumberOfVariables, fluxOut);

}

exahype::solvers::Solver::RefinementControl NavierStokes::NavierStokesSolverDG::refinementCriterion(
    const double* luh,
    const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx,
    const double t,
    const int level) {
  const bool isAmrScenario =
          scenarioName == "two-bubbles" ||
          scenarioName == "density-current";
  if (!isAmrScenario || DIMENSIONS != 2) {
    return exahype::solvers::Solver::RefinementControl::Keep;
  }

  if (t == 0) {
    // Global observables are reduced after first timestep!
    return exahype::solvers::Solver::RefinementControl::Keep;
  }

  const auto maxGlobal = _globalObservables[0];
  const auto minGlobal = _globalObservables[1];
  const auto maxGlobalDiff = maxGlobal - minGlobal;

  const auto curObservables = mapGlobalObservables(luh, dx);
  const auto maxDiff = curObservables[0] - minGlobal;

  if (maxDiff > maxGlobalDiff/10) {
    return exahype::solvers::Solver::RefinementControl::Refine;
  }
  if (maxDiff < maxGlobalDiff/15.) {
    return exahype::solvers::Solver::RefinementControl::Erase;
  }

  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void NavierStokes::NavierStokesSolverDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  ns.evaluateEigenvalues(Q, d, lambda);
}

void NavierStokes::NavierStokesSolverDG::viscousEigenvalues(const double* const Q,const int d,double* lambda) {
  ns.evaluateDiffusiveEigenvalues(Q, d, lambda);
}

void NavierStokes::NavierStokesSolverDG::viscousFlux(const double *const Q, const double *const gradQ, double **F) {
  ns.evaluateFlux(Q, gradQ, F, true);
}

double NavierStokes::NavierStokesSolverDG::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {
  return (0.7/0.9) * stableDiffusiveTimeStepSize<NavierStokesSolverDG>(*static_cast<NavierStokesSolverDG*>(this),luh,dx);
  //return kernels::aderdg::generic::c::stableTimeStepSize<NavierStokesSolverDG, true>(*static_cast<NavierStokesSolverDG*>(this),luh,dx);
}

void NavierStokes::NavierStokesSolverDG::riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,const double dt,const tarch::la::Vector<DIMENSIONS, double>& lengthScale, const int direction, bool isBoundaryFace, int faceIndex) {
  assertion2(direction>=0,dt,direction);
  assertion2(direction<DIMENSIONS,dt,direction);
  //kernels::aderdg::generic::c::riemannSolverNonlinear<false, NavierStokesSolverDG>(*static_cast<NavierStokesSolverDG*>(this),FL,FR,QL,QR,dt,direction);
  riemannSolverNonlinear<false,NavierStokesSolverDG>(*static_cast<NavierStokesSolverDG*>(this),FL,FR,QL,QR,lengthScale, dt,direction);

}

void NavierStokes::NavierStokesSolverDG::boundaryConditions( double* const fluxIn, const double* const stateIn, const double* const gradStateIn, const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& cellCentre, const tarch::la::Vector<DIMENSIONS,double>&  cellSize, const double t,const double dt, const int direction, const int orientation) {
  constexpr int basisSize     = (Order+1);
  constexpr int sizeStateOut = (NumberOfVariables+NumberOfParameters)*basisSize;
  constexpr int sizeFluxOut  = NumberOfVariables*basisSize;

  constexpr int totalSize = sizeStateOut + sizeFluxOut;
  double* block = new double[totalSize];

  double* memory = block;

  double* stateOut = memory; memory+=sizeStateOut;
  double* fluxOut  = memory; memory+=sizeFluxOut;

  const int faceIndex = 2*direction+orientation;

  kernels::aderdg::generic::c::boundaryConditions<true, NavierStokesSolverDG>(*static_cast<NavierStokesSolverDG*>(this),fluxOut,stateOut,fluxIn,stateIn,gradStateIn, cellCentre,cellSize,t,dt,faceIndex,direction);

  if (orientation==0) {
    double* FL = fluxOut; const double* const QL = stateOut;
    double* FR =  fluxIn;  const double* const QR = stateIn;

    riemannSolverNonlinear<false,NavierStokesSolverDG>(*static_cast<NavierStokesSolverDG*>(this),FL,FR,QL,QR,cellSize, dt,direction);
  }
  else {
    double* FL =  fluxIn;  const double* const QL = stateIn;
    double* FR = fluxOut; const double* const QR = stateOut;

    riemannSolverNonlinear<false,NavierStokesSolverDG>(*static_cast<NavierStokesSolverDG*>(this),FL,FR,QL,QR,cellSize, dt,direction);
  }

  if (scenario->getBoundaryType(faceIndex) == NavierStokes::BoundaryType::wall) {
    static_assert(DIMENSIONS == 2, "BC only implemented for 2D!"); // TODO(Lukas) Implement for 3D
    kernels::idx2 idx_F(Order + 1, NumberOfVariables);
    for (int i = 0; i < (Order + 1); ++i) {
      // Set energy flux to zero!
      fluxIn[idx_F(i, NavierStokesSolverDG_Variables::shortcuts::E)] = 0.0; // TODO(Lukas) Is fluxIn later reused?
    }
  }

  delete[] block;
}

std::vector<double> NavierStokes::NavierStokesSolverDG::mapGlobalObservables(const double *const Q,
        const tarch::la::Vector<DIMENSIONS,double>& dx) const {
 auto observables = resetGlobalObservables();
 // TODO(Lukas) Implement global observables for 3D!
 const auto idxQ = kernels::idx3(Order+1,Order+1,NumberOfVariables + NumberOfParameters);

 auto computePotT = [this](const double *const Q) {
   const auto vars = ReadOnlyVariables{Q};
   const auto pressure = ns.evaluatePressure(vars.E(), vars.rho(), vars.j());
   const auto temperature = ns.evaluateTemperature(vars.rho(), pressure);
   return ns.evaluatePotentialTemperature(temperature, pressure);
 };

 const auto tv = totalVariation(Q, Order, NumberOfVariables, NumberOfParameters, dx, false, computePotT);

 /*
 dfor(i,Order+1) {
   const auto vars = ReadOnlyVariables{Q + idxQ(i(1), i(0), 0)};
   const auto pressure = ns.evaluatePressure(vars.E(), vars.rho(), vars.j());
   const auto temperature = ns.evaluateTemperature(vars.rho(), pressure);
   const auto potT = ns.evaluatePotentialTemperature(temperature, pressure);
   reduceGlobalObservables(observables, std::vector<double>{potT, potT});
 }
  */
 observables = {tv, tv};

 return observables;
}

std::vector<double> NavierStokes::NavierStokesSolverDG::resetGlobalObservables() const {
  return {std::numeric_limits<double>::min(), std::numeric_limits<double>::max()};
}

void NavierStokes::NavierStokesSolverDG::reduceGlobalObservables(
        std::vector<double> &reducedGlobalObservables,
        const std::vector<double> &curGlobalObservables) const {
  assertion2(reducedGlobalObservables.size() == curGlobalObservables.size(),
          reducedGlobalObservables.size(),
          curGlobalObservables.size());
  reducedGlobalObservables[0] = std::max(reducedGlobalObservables[0],
                                         std::abs(curGlobalObservables[0]));

  reducedGlobalObservables[1] = std::min(reducedGlobalObservables[1],
                                         std::abs(curGlobalObservables[1]));
}
