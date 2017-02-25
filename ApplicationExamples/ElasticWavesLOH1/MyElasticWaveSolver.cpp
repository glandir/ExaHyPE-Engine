#include "MyElasticWaveSolver.h"

#include "MyElasticWaveSolver_Variables.h"

void ElasticWave::MyElasticWaveSolver::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

bool ElasticWave::MyElasticWaveSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) {
  // @todo Please implement/augment if required
    if (tarch::la::equals(t, 0.0, 1e-15)) {  // @todo precision
    return true;
  }
   return false;
}

void ElasticWave::MyElasticWaveSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 3
  // Number of variables    = 12 (#unknowns + #parameters)
  
  // @todo Please implement/augment if required
  // State variables:
  Variables vars(Q);
  vars.v(0,0,0);
  vars.s(0,0,0,0,0,0);

  vars.rho() = 2.6;     // gm/cm^3
  vars.cs() = 2.0;      // km/s
  vars.cp() = 4.0;      // km/s
  
  if (x[0] > 1.0) {
  vars.rho() = 2.7;     // gm/cm^3
  vars.cs() = 3.464;    // km/s
  vars.cp() = 6.0;      // km/s
  }
  // Q[ 0] = 0.0;
  // Q[ 1] = 0.0;
  // Q[ 2] = 0.0;
  // Q[ 3] = 0.0;
  // Q[ 4] = 0.0;
  // Q[ 5] = 0.0;
  // Q[ 6] = 0.0;
  // Q[ 7] = 0.0;
  // Q[ 8] = 0.0;  // Material parameters:
  // Q[ 9] = 0.0;
  // Q[10] = 0.0;
  // Q[11] = 0.0;
}

void ElasticWave::MyElasticWaveSolver::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions             = 3
  // Number of variables    = 12 (#unknowns + #parameters)
  
  // @todo Please implement/augment if required
   ReadOnlyVariables vars(Q);
   double cs  = vars.cs();   // km/s
  double cp  = vars.cp();   // km/s

  lambda[0] = -cp;
  lambda[1] = -cs;
  lambda[2] = -cs;
  lambda[3] = 0.0;
  lambda[4] = 0.0;
  lambda[5] = 0.0;
  lambda[6] = +cs;
  lambda[7] = +cs;
  lambda[8] = +cp;
}


void ElasticWave::MyElasticWaveSolver::flux(const double* const Q,double** F) {
  // Dimensions             = 3
  // Number of variables    = 12 (#unknowns + #parameters)
  
  // @todo Please implement/augment if required
  F[0][ 0] = 0.0;
  F[0][ 1] = 0.0;
  F[0][ 2] = 0.0;
  F[0][ 3] = 0.0;
  F[0][ 4] = 0.0;
  F[0][ 5] = 0.0;
  F[0][ 6] = 0.0;
  F[0][ 7] = 0.0;
  F[0][ 8] = 0.0;

  F[1][ 0] = 0.0;
  F[1][ 1] = 0.0;
  F[1][ 2] = 0.0;
  F[1][ 3] = 0.0;
  F[1][ 4] = 0.0;
  F[1][ 5] = 0.0;
  F[1][ 6] = 0.0;
  F[1][ 7] = 0.0;
  F[1][ 8] = 0.0;

  F[2][ 0] = 0.0;
  F[2][ 1] = 0.0;
  F[2][ 2] = 0.0;
  F[2][ 3] = 0.0;
  F[2][ 4] = 0.0;
  F[2][ 5] = 0.0;
  F[2][ 6] = 0.0;
  F[2][ 7] = 0.0;
  F[2][ 8] = 0.0;
}


void ElasticWave::MyElasticWaveSolver::source(const double* const Q,double* S) {
  // Dimensions             = 3
  // Number of variables    = 12 (#unknowns + #parameters)
  
  // @todo Please implement/augment if required
  S[ 0] = 0.0;
  S[ 1] = 0.0;
  S[ 2] = 0.0;
  S[ 3] = 0.0;
  S[ 4] = 0.0;
  S[ 5] = 0.0;
  S[ 6] = 0.0;
  S[ 7] = 0.0;
  S[ 8] = 0.0;
}


void ElasticWave::MyElasticWaveSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int d,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions             = 3
  // Number of variables    = 12 (#unknowns + #parameters)

  // @todo Please implement/augment if required
  stateOut[0] = 0.0;
  stateOut[1] = 0.0;
  stateOut[2] = 0.0;

  stateOut[3] = 0.0;
  stateOut[4] = 0.0;
  stateOut[5] = 0.0;

  stateOut[6] = 0.0;
  stateOut[7] = 0.0;
  stateOut[8] = 0.0;

  fluxOut[0] = 0.0;
  fluxOut[1] = 0.0;
  fluxOut[2] = 0.0;

  fluxOut[3] = 0.0;
  fluxOut[4] = 0.0;
  fluxOut[5] = 0.0;

  fluxOut[6] = 0.0;
  fluxOut[7] = 0.0;
  fluxOut[8] = 0.0;

  if (faceIndex == 0){

  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];

  stateOut[3] = -stateIn[3];
  stateOut[4] = stateIn[4];
  stateOut[5] = stateIn[5];

  stateOut[6] = -stateIn[6];
  stateOut[7] = -stateIn[7];
  stateOut[8] = stateIn[8];

  fluxOut[0] = fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = fluxIn[2];

  fluxOut[3] = fluxIn[3];
  fluxOut[4] = fluxIn[4];
  fluxOut[5] = fluxIn[5];

  fluxOut[6] = fluxIn[6];
  fluxOut[7] = fluxIn[7];
  fluxOut[8] = fluxIn[8];
  }

}


exahype::solvers::Solver::RefinementControl ElasticWave::MyElasticWaveSolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


bool ElasticWave::MyElasticWaveSolver::physicalAdmissibilityDetection(const double* const QMin,const double* const QMax) {
  // @todo Please implement/augment if required
  return true;
}


void ElasticWave::MyElasticWaveSolver::ncp(const double* const Q,const double* const gradQ,double* BgradQ) {
  // Dimensions             = 3
  // Number of variables    = 12 (#unknowns + #parameters)

  // @todo Please implement/augment if required
  ReadOnlyVariables vars(Q);

  const double* Qx = &gradQ[0]; // (:,1)
  const double* Qy = &gradQ[9]; // (:,2)
  const double* Qz = &gradQ[18]; // (:,3)

  //  std::cout << "ncp " << vars.rho() << "," << vars.cs() << "," << vars.cp() << std::endl;
  //  double cs = 3.464;  // km/s
  //  double cp = 6.0;    // km/s
  //  double rho = 2.67;  // gm/cm^3

  double cs  = vars.cs();   // km/s
  double cp  = vars.cp();   // km/s
  double rho = vars.rho();  // gm/cm^3

  double lam = rho * (cp * cp - 2. * cs * cs);  // lambda (MPa)
  double mu  = rho * cs * cs;                   // mu (MPa)

  BgradQ[0] = -1.0/rho*Qx[3];
  BgradQ[1] = -1.0/rho*Qx[6];
  BgradQ[2] = -1.0/rho*Qx[7];
  BgradQ[3] = -(2.0*mu + lam)*Qx[0];
  BgradQ[4] = -lam*Qx[0];
  BgradQ[5] = -lam*Qx[0];
  BgradQ[6] = -mu*Qx[1];
  BgradQ[7] = -mu*Qx[2];
  BgradQ[8] = 0.0;

  BgradQ[9]  = -1.0/rho*Qy[6];
  BgradQ[10] = -1.0/rho*Qy[4];
  BgradQ[11] = -1.0/rho*Qy[8];
  BgradQ[12] = -lam*Qy[1];
  BgradQ[13] = -(lam+2.0*mu)*Qy[1];
  BgradQ[14] = -lam*Qy[1];
  BgradQ[15] = -mu*Qy[0];
  BgradQ[16] = -0.0;
  BgradQ[17] = -mu*Qy[2];

  BgradQ[18] = -1.0/rho*Qz[7];
  BgradQ[19] = -1.0/rho*Qz[8];
  BgradQ[20] = -1.0/rho*Qz[5];
  BgradQ[21] = -lam*Qz[2];
  BgradQ[22] = -lam*Qz[2];
  BgradQ[23] = -(lam+2.0*mu)*Qz[2];
  BgradQ[24] = 0.0;
  BgradQ[25] = -mu*Qz[0];
  BgradQ[26] = -mu*Qz[1];
}


void ElasticWave::MyElasticWaveSolver::matrixb(const double* const Q,const int d,double* Bn) {
  // Dimensions             = 3
  // Number of variables    = 12 (#unknowns + #parameters)
  
  // @todo Please implement/augment if required
  Bn[ 0] = 0.0;
  Bn[ 1] = 0.0;
  Bn[ 2] = 0.0;
  Bn[ 3] = 0.0;
  Bn[ 4] = 0.0;
  Bn[ 5] = 0.0;
  Bn[ 6] = 0.0;
  Bn[ 7] = 0.0;
  Bn[ 8] = 0.0;
  Bn[ 9] = 0.0;
  Bn[10] = 0.0;
  Bn[11] = 0.0;
  Bn[12] = 0.0;
  Bn[13] = 0.0;
  Bn[14] = 0.0;
  Bn[15] = 0.0;
  Bn[16] = 0.0;
  Bn[17] = 0.0;
  Bn[18] = 0.0;
  Bn[19] = 0.0;
  Bn[20] = 0.0;
  Bn[21] = 0.0;
  Bn[22] = 0.0;
  Bn[23] = 0.0;
  Bn[24] = 0.0;
  Bn[25] = 0.0;
  Bn[26] = 0.0;
  Bn[27] = 0.0;
  Bn[28] = 0.0;
  Bn[29] = 0.0;
  Bn[30] = 0.0;
  Bn[31] = 0.0;
  Bn[32] = 0.0;
  Bn[33] = 0.0;
  Bn[34] = 0.0;
  Bn[35] = 0.0;
  Bn[36] = 0.0;
  Bn[37] = 0.0;
  Bn[38] = 0.0;
  Bn[39] = 0.0;
  Bn[40] = 0.0;
  Bn[41] = 0.0;
  Bn[42] = 0.0;
  Bn[43] = 0.0;
  Bn[44] = 0.0;
  Bn[45] = 0.0;
  Bn[46] = 0.0;
  Bn[47] = 0.0;
  Bn[48] = 0.0;
  Bn[49] = 0.0;
  Bn[50] = 0.0;
  Bn[51] = 0.0;
  Bn[52] = 0.0;
  Bn[53] = 0.0;
  Bn[54] = 0.0;
  Bn[55] = 0.0;
  Bn[56] = 0.0;
  Bn[57] = 0.0;
  Bn[58] = 0.0;
  Bn[59] = 0.0;
  Bn[60] = 0.0;
  Bn[61] = 0.0;
  Bn[62] = 0.0;
  Bn[63] = 0.0;
  Bn[64] = 0.0;
  Bn[65] = 0.0;
  Bn[66] = 0.0;
  Bn[67] = 0.0;
  Bn[68] = 0.0;
  Bn[69] = 0.0;
  Bn[70] = 0.0;
  Bn[71] = 0.0;
  Bn[72] = 0.0;
  Bn[73] = 0.0;
  Bn[74] = 0.0;
  Bn[75] = 0.0;
  Bn[76] = 0.0;
  Bn[77] = 0.0;
  Bn[78] = 0.0;
  Bn[79] = 0.0;
  Bn[80] = 0.0;

  ReadOnlyVariables vars(Q);
  double rho = vars.rho();
  double cs  = vars.cs();
  double cp  = vars.cp();
  
  double lam = rho * (cp * cp - 2 * cs * cs);  // lambda (MPa)
  double mu= rho * cs * cs;

  double nv[3] = {0.0};
  nv[d] = 1.0;
  
  double B1[9][9];
  double B2[9][9];
  double B3[9][9];

  for(int i=0; i<9; i++) {
   for(int j=0; j<9; j++) {
     B1[i][j] = 0.0;
     B2[i][j] = 0.0;
     B3[i][j] = 0.0;
   }
 }

  B1[3][0] = -1.0/rho; 
  B1[6][1] = -1.0/rho;
  B1[7][2] = -1.0/rho;

  B1[0][3] = -(2.0*mu + lam); 
  B1[0][4] = -lam;
  B1[0][5] = -lam;

  B1[1][6] = -mu; 
  B1[2][7] = -mu;

  B2[6][0] = -1.0/rho; 
  B2[4][1] = -1.0/rho;
  B2[8][2] = -1.0/rho;

  B2[1][3] = -lam;
  B2[1][4] = -(2.0*mu + lam);
  B2[1][5] = -lam;

  B2[0][6] = -mu; 
  B2[2][8] = -mu;

  B3[7][0] = -1.0/rho; 
  B3[8][1] = -1.0/rho;
  B3[5][2] = -1.0/rho;

  B3[2][3] = -lam;
  B3[2][4] = -lam;
  B3[2][5] = -(2.0*mu + lam);

  B3[0][7] = -mu; 
  B3[1][8] = -mu;


  for(int i=0; i<9; i++) {
   for(int j=0; j<9; j++) {
      Bn[i*9 + j] = nv[0]*B1[i][j] + nv[1]*B2[i][j] + nv[2]*B3[i][j];
   }
 }
}


//TODO KD
// tell the user what it is
bool ElasticWave::MyElasticWaveSolver::hasToApplyPointSource() const { 
  return true;
}

//TODO KD
void ElasticWave::MyElasticWaveSolver::pointSource(const double* const x,const double t,const double dt, double* forceVector, double* x0) {
  //TODO KD // @todo Please implement/augment if required and set bool function

   double pi = 3.14159265359;
  double sigma = 0.1149;
  //double t0 = 0.7;
  double t0 = 0.1;
  double f;
  double M0 = 1000.0;
  

  //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-pow((t-t0)/(std::sqrt(2)*sigma),2.0)));

  //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));


  f = M0*t/(t0*t0)*std::exp(-t/t0);

  x0[0] = 2.;
  x0[1] = 15.;
  x0[2] = 15.;

  forceVector[0] = 0.0;
  forceVector[1] = 0.0;
  forceVector[2] = 0.0;
  forceVector[3] = 0.0*f;
  forceVector[4] = 0.0*f;
  forceVector[5] = 0.0*f;
  forceVector[6] = 0.0;
  forceVector[7] = 0.0;
  forceVector[8] = 1.0*f;
}
