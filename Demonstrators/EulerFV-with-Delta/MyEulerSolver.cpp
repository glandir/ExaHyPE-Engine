#include "MyEulerSolver.h"

#include "MyEulerSolver_Variables.h"

#include "Logo.h"

#include "delta/ContactPoint.h"
#include "delta/primitives/Cylinder.h"
#include "delta/primitives/Sphere.h"
#include "delta/contactdetection/sphere.h"
#include "delta/io/vtk.h"


tarch::logging::Log EulerFV::MyEulerSolver::_log( "EulerFV::MyEulerSolver" );


void EulerFV::MyEulerSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  _embeddedGeometry = new delta::primitives::Cylinder(
    0.5, 0.5, 0.0,  // center
	0.2,       // rdius
	-0.1, 0.1, // min/maxZ
	0.1
  );

/*
  _embeddedGeometry = new delta::primitives::Sphere(
    0.5, 0.5, 0.0,  // center
	0.4,       // rdius
	0.1
  );
*/

  delta::io::writeVTK(
    _embeddedGeometry->getNumberOfTriangles(),
	_embeddedGeometry->getXCoordinates(),
	_embeddedGeometry->getYCoordinates(),
	_embeddedGeometry->getZCoordinates(),
    "embedded-geometry.vtk"
  );
}


void EulerFV::MyEulerSolver::adjustSolution(const double* const x,const double t,const double dt, double* Q) {
  Variables vars(Q);
  if ( tarch::la::equals( t,0.0 ) ) {
    tarch::la::Vector<DIMENSIONS,double> myX( x[0] - 0.06, 1.0-x[1] - 0.25 ); // translate
    myX *= static_cast<double>(Image.width);
    tarch::la::Vector<DIMENSIONS,int>    myIntX( 1.2*myX(0) , 1.2*myX(1) );  // scale

    double Energy = 0.1;

    if (
      myIntX(0) > 0 && myIntX(0) < static_cast<int>(Image.width)
      &&
      myIntX(1) > 0 && myIntX(1) < static_cast<int>(Image.height)
    ) {
      Energy += 1.0-Image.pixel_data[myIntX(1)*Image.width+myIntX(0)];
    }

    vars.rho() = 1.0;
    vars.E()   = Energy;
    vars.j(0,0,0);
  }

  // Boundary stuff
  const double maxDistance = 1e-2; // should be roughly sqrt(h)*d
  std::vector< delta::ContactPoint > contact =
    delta::contactdetection::sphereToSphere(
      _embeddedGeometry->getCentreX(),
      _embeddedGeometry->getCentreY(),
      _embeddedGeometry->getCentreZ(),
      _embeddedGeometry->getBoundingSphereRadius(),

      x[0], // voxel centre
      x[1], // voxel centre
      x[2], // voxel centre

	  1e-2,  // voxel size
	  maxDistance
    );
  if ( contact.empty() ) {
    vars.inside() = 2.0*maxDistance; // If maxDistance is the epsilon environment
                                     // we can twice this.
  }
  else {
    vars.inside() = contact[0].distance;
    assertion6(
      vars.inside()<=2.0 * maxDistance,
	  vars.inside(), maxDistance,
	  x[0], x[1], x[2],
	  contact[0].toString()
	);
  }
  logDebug( "adjustSolution(...)", "voxel at " << x[0] << "," << x[1] << "," << x[2] << " has value " << vars.inside() );
}


void EulerFV::MyEulerSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  double u_n = Q[normalNonZeroIndex + 1] * irho;
  double c  = std::sqrt(GAMMA * p * irho);

  eigs.rho()=u_n - c;
  eigs.E()  =u_n + c;
  eigs.j(u_n,u_n,u_n);
}


void EulerFV::MyEulerSolver::flux(const double* const Q, double** F) {
  ReadOnlyVariables vars(Q);
  Fluxes fluxes(F);

  tarch::la::Matrix<3,3,double> I;
  I = 1, 0, 0,
      0, 1, 0,
      0, 0, 1;

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  fluxes.rho ( vars.j()                                 );
  fluxes.j   ( irho * outerDot(vars.j(),vars.j()) + p*I );
  fluxes.E   ( irho * (vars.E() + p) * vars.j()         );
}


void EulerFV::MyEulerSolver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int normalNonZero,
    const double* const stateInside,
    double* stateOutside) {
  ReadOnlyVariables varsInside(stateInside);
  Variables         varsOutside(stateOutside);

  varsOutside = varsInside;
}
