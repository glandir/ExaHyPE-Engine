#include "MyLinearSolver.h"
#include "kernels/KernelUtils.h"

#include "kernels/aderdg/generic/Kernels.h"

#include "MyLinearSolver_Variables.h"

#include <algorithm>

#include "CurvilinearTransformation.h"


tarch::logging::Log Linear::MyLinearSolver::_log( "Linear::MyLinearSolver" );


void Linear::MyLinearSolver::init(std::vector<std::string>& cmdlineargs) {
  static tarch::logging::Log _log("MyLinearSolver::init");
  // @todo Please implement/augment if required
}

exahype::solvers::ADERDGSolver::AdjustSolutionValue Linear::MyLinearSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  // @todo Please implement/augment if required
  static tarch::logging::Log _log("MyLinearSolver::AdjustSolutionValue");
  return tarch::la::equals(t,0.0) ? exahype::solvers::ADERDGSolver::AdjustSolutionValue::PointWisely : exahype::solvers::ADERDGSolver::AdjustSolutionValue::No;
}

void Linear::MyLinearSolver::adjustPointSolution(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 2
  // Number of variables    = 3 + #parameters
  
  // @todo Please implement/augment if required
  // State variables:
  static tarch::logging::Log _log("MyLinearSolver::adjustPointSolution");

  Variables vars(Q);

  vars.p() = std::exp(-((x[0]-0.5)*(x[0]-0.5)+(x[1]-0.5)*(x[1]-0.5))/0.01);
  //  vars.p()=0;
  vars.v(0,0);

  int num_pos_x = 4;
  int num_pos_y = 4;


  double left_bnd_x[4] = {0,0,0,0};
  double left_bnd_y[4] ={0,1.0/3.0,2.0/3.0,1};
  double right_bnd_x[4]  ={1,1,1,1};
  double right_bnd_y[4] ={0,1.0/3.0,2.0/3.0,1};

  double top_bnd_y[4] ={1,1+1.0/3.0*0.2,1+2.0/3.0*0.2,1.2};
  double top_bnd_x[4] ={0,1.0/3.0,2.0/3.0,1};

  double bottom_bnd_y[4] ={0,0,0,0};
  double bottom_bnd_x[4] ={0,1.0/3.0,2.0/3.0,1};

  double curvilinear_x[16];
  double curvilinear_y[16];


  transFiniteInterpolation( num_pos_x,  num_pos_y,  left_bnd_x,  left_bnd_y,  right_bnd_x,  right_bnd_y,  bottom_bnd_x,  bottom_bnd_y,  top_bnd_x,  top_bnd_y,  curvilinear_x ,  curvilinear_y );

  kernels::idx2 id_xy(num_pos_x,num_pos_y);

  for(int j = 0 ; j < num_pos_y ; j ++){
    //    printf("%f \n",top_bnd_y[j]);
    for(int i = 0 ; i< num_pos_x ; i ++){

      // printf("%d \n",i);
      // printf("%d \n",j);

      // printf("%f\n",curvilinear_x[id_xy(j,i)]);
      //       printf("%f\n",curvilinear_y[id_xy(i,j)]);
    } 
  } 

  for(int j = 0 ; j < num_pos_y ; j ++){
    for(int i = 0 ; i< num_pos_x ; i ++){

      // printf("%d \n",i);
      // printf("%d \n",j);

      // printf("%f\n",curvilinear_x[id_xy(j,i)]);
    } 
  } 


 double unif_mesh[4] ={0,1.0/3.0,2.0/3.0,1};

 double coeff=0;
 interpolate(0.11,0.12,unif_mesh,unif_mesh,curvilinear_x,num_pos_y,coeff);

 double coeff2=0;
 interpolate(0.11,0.12,unif_mesh,unif_mesh,curvilinear_y,num_pos_y,coeff2);
 double gl_vals[16];
 getValuesAtQuadNodes(unif_mesh,unif_mesh,curvilinear_y,num_pos_y,gl_vals);

 

 for(int j = 0 ; j < num_pos_y ; j ++){
    for(int i = 0 ; i< num_pos_x ; i ++){

      // printf("%d %d:  %f\n",i,j,gl_vals[id_xy(i,j)]);

    } 
  } 
 
 double der_x;
 computeDerivatives_x(2,3,gl_vals,num_pos_x,der_x);

 double der_y;
 computeDerivatives_y(2,3,gl_vals,num_pos_x,der_y);


 // printf("%f\n",der_x);
 // printf("%f\n",der_y);

   

 std::exit(-1);
}

void Linear::MyLinearSolver::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 3 + #parameters
  
  // @todo Please implement/augment if required

  static tarch::logging::Log _log("MyLinearSolver::eigenvalues");
  lambda[0] = 1.0;
  lambda[1] = 1.0;
  lambda[2] = 1.0;
}


void Linear::MyLinearSolver::flux(const double* const Q,double** F) {
  // Dimensions             = 2
  // Number of variables    = 3 + #parameters
  
  // @todo Please implement/augment if required
  static tarch::logging::Log _log("MyLinearSolver::flux");
  F[0][0] = 0.0;
  F[0][1] = 0.0;
  F[0][2] = 0.0;

  F[1][0] = 0.0;
  F[1][1] = 0.0;
  F[1][2] = 0.0;

}


void Linear::MyLinearSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
					      const double * const fluxIn,const double* const stateIn,
					      double *fluxOut,double* stateOut) {
  // Dimensions             = 2
  // Number of variables    = 3 + #parameters

  // @todo Please implement/augment if required
static tarch::logging::Log _log("MyLinearSolver::boundaryValues");

stateOut[0] = stateIn[0];
stateOut[1] = stateIn[1];
stateOut[2] = stateIn[2];
fluxOut[0] =  fluxIn[0];
fluxOut[1] = fluxIn[1];
fluxOut[2] = fluxIn[2];
}


exahype::solvers::Solver::RefinementControl Linear::MyLinearSolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  static tarch::logging::Log _log("MyLinearSolver::refinementCriterion");
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void Linear::MyLinearSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ){
  kernels::idx2 idx(DIMENSIONS, NumberOfVariables);

  static tarch::logging::Log _log("MyLinearSolver::nonConservativeProduct");

  BgradQ[0] = -gradQ[1];
  BgradQ[1] = -gradQ[0];
  BgradQ[2] = 0;

  BgradQ[3]= -gradQ[5];
  BgradQ[4]= 0;
  BgradQ[5]= -gradQ[3];

}


void Linear::MyLinearSolver::coefficientMatrix(const double* const Q,const int d,double* Bn){
  static tarch::logging::Log _log("MyLinearSolver::coefficientMatrix");

  double nv[2] = {0.0};

  nv[d] = 1.0;
  
  double B1[3][3];
  double B2[3][3];
   
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      B1[i][j] = 0.0;
      B2[i][j] = 0.0;
    }
  }
  
  B1[0][1] = -1.0; 
  B1[1][0] = -1.0;
  
  B2[0][2] = -1.0; 
  B2[2][0] = -1.0;
  
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      Bn[i*3+ j] = nv[0]*B1[i][j] + nv[1]*B2[i][j];
    }
  }

}


void Linear::MyLinearSolver::riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,double* tempFaceUnknownsArray,double** tempStateSizedVectors,double** tempStateSizedSquareMatrices,const double dt,const int normalNonZeroIndex){

  constexpr int numberOfVariables  = MyLinearSolver::NumberOfVariables;
  constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
  constexpr int numberOfParameters = MyLinearSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int basisSize          = MyLinearSolver::Order+1;
  constexpr int order              = basisSize - 1;
  
  kernels::idx2 idx_QLR(basisSize, numberOfVariables);

  kernels::idx2 idx_FLR(basisSize, numberOfVariables);

  double n[2]={0,0};
  n[normalNonZeroIndex]=1;


    for (int i = 0; i < basisSize; i++) {
      double v_p=QL[idx_QLR(i,1)]*n[0]+QL[idx_QLR(i,2)]*n[1];

      double v_m=QR[idx_QLR(i,1)]*n[0]+QR[idx_QLR(i,2)]*n[1];

      // printf("QR1: %f \n",QR[idx_QLR(i,1)]);
      // printf("QR2: %f \n",QR[idx_QLR(i,2)]);
      // printf("n1: %f \n",n[0]);
      // printf("n2: %f \n",n[1]);
      // printf("nnZI: %d \n",normalNonZeroIndex);
      
      //      printf("v_p: %f \n",v_p);
      //      printf("v_m: %f \n",v_m);

      double sigma_p = QL[idx_QLR(i,0)];
      double sigma_m = QR[idx_QLR(i,0)];

      //      printf("sigma_p: %f \n",sigma_p);
      //      printf("sigma_m: %f \n",sigma_m);

      double z_p=1.0;
      double z_m=1.0;

      double v_hat_p=0;
      double v_hat_m=0;
      double sigma_hat_p=0;
      double sigma_hat_m=0;

      // riemannSolver_Nodal(v_p,v_m, sigma_p, sigma_m, z_p , z_m, v_hat_p , v_hat_m, sigma_hat_p, sigma_hat_m);

      v_hat_p= v_m;
      v_hat_m= v_p;
      sigma_hat_m=sigma_p;
      sigma_hat_p=sigma_m;

      FR[idx_FLR(i, 0)] = -0.5*(z_m *(v_m-v_hat_m) - (sigma_m-sigma_hat_m));
      FL[idx_FLR(i, 0)] = 0.5*(z_p *(v_p-v_hat_p) + (sigma_p-sigma_hat_p));

      // printf("FR: %f \n",FR[idx_FLR(i, 0)]);
      // printf("FL: %f \n",FL[idx_FLR(i, 0)]);


      FR[idx_FLR(i, 1)] = 0.5*((v_m-v_hat_m) - (sigma_m-sigma_hat_m)/z_m)*n[0];
      FL[idx_FLR(i, 1)] = 0.5*((v_p-v_hat_p) + (sigma_p-sigma_hat_p)/z_p)*n[0];

      FR[idx_FLR(i, 2)] = 0.5*((v_m-v_hat_m) - (sigma_m-sigma_hat_m)/z_m)*n[1];
      FL[idx_FLR(i, 2)] = 0.5*((v_p-v_hat_p) + (sigma_p-sigma_hat_p)/z_p)*n[1];

    }

}



void Linear::MyLinearSolver::pointSource(const double* const x,const double t,const double dt, double* forceVector, double* x0){
     double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.7;
  //double t0 = 0.1;
  double f = 0.0;
  double M0 = 1000.0;
  

  //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-pow((t-t0)/(std::sqrt(2)*sigma),2.0)));

  f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));


  //f = M0*t/(t0*t0)*std::exp(-t/t0);

  x0[0] = 0.5;
  x0[1] = 0.5;

  forceVector[0] = 1.*f;
  forceVector[1] = 0.0;
  forceVector[2] = 0.0;
}



 void Linear::MyLinearSolver::riemannSolver_Nodal(double v_p,double v_m, double sigma_p, double sigma_m, double z_p , double z_m, double v_hat_p , double v_hat_m, double sigma_hat_p, double sigma_hat_m){
   double p=0;
   double q=0;
   double phi=0;
   double v_hat=0;
   double eta=0;

   p=z_p*v_m + sigma_m;
   q=z_m*v_p - sigma_p;

   eta=(z_p*z_m)/(z_p+z_m);

   phi= p/z_p - q/z_m;

   sigma_hat_p=phi;
   sigma_hat_m=phi;

   v_hat_m=(p-phi)/z_p;
   v_hat_p=(q+phi)/z_m;
 }
