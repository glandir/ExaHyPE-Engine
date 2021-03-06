#include "SodShockTube.h"
void NavierStokes::SodShockTube::initialValues(const double* const x,
                                               const PDE& ns, Variables& vars) {
  double pressure = 0.0;
  if (x[1] < 0.45) {
    vars.rho() = 0.125;
    pressure = 0.1;
  } else {
    vars.rho() = 1.0;
    pressure = 1.0;
  }
#if DIMENSIONS == 2
  vars.j(0.0, 0.0);
#elif DIMENSIONS == 3
  vars.j(0.0, 0.0, 0.0);
#endif
  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());
}
