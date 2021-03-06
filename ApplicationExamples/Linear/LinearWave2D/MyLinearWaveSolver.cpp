#include "MyLinearWaveSolver.h"

#include "MyLinearWaveSolver_Variables.h"

#include "kernels/KernelUtils.h"

#include "kernels/aderdg/generic/Kernels.h"



tarch::logging::Log Linear::MyLinearWaveSolver::_log( "Linear::MyLinearWaveSolver" );


void Linear::MyLinearWaveSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void Linear::MyLinearWaveSolver::adjustSolution(double* const luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 3 + 0
  // @todo Please implement/augment if required
  
  if (tarch::la::equals(t,0.0)) {
    //initialise luh
    
    constexpr int basisSize = MyLinearWaveSolver::Order+1;
    int numberOfData=MyLinearWaveSolver::NumberOfParameters+MyLinearWaveSolver::NumberOfVariables;


     kernels::idx3 id_3(basisSize,basisSize,numberOfData);

     double offset_x=center[0]-0.5*dx[0];
     double offset_y=center[1]-0.5*dx[1];
    

     double width_x=dx[0];
     double width_y=dx[1];
   

    for (int j=0; j< basisSize; j++){
      for (int i=0; i< basisSize; i++){
     
	double x  =  (offset_x+width_x*kernels::gaussLegendreNodes[basisSize-1][i]);
	double y  =  (offset_y+width_y*kernels::gaussLegendreNodes[basisSize-1][j]);


	//Pressure
	luh[id_3(j,i,0)]  = 0*std::exp(-100*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));

	// //Velocity
	luh[id_3(j,i,1)]  = 0;
	luh[id_3(j,i,2)]  = 0;

	// Parameters
	//	luh[id_3(j,i,3)]  = 2.7;     // density
	//	luh[id_3(j,i,4)]  = 6.0;   // wavespeed
	
	// Parameters
	luh[id_3(j,i,3)]  = 1.0;     // density
	luh[id_3(j,i,4)]  = 1.0;   // wavespeed
	
      }
    }

    
  }
  
}

void Linear::MyLinearWaveSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double* const fluxOut,double* const stateOut) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 3 + 0

  // @todo Please implement/augment if required

  constexpr int numberOfVariables  = MyLinearWaveSolver::NumberOfVariables;
  constexpr int numberOfParameters = MyLinearWaveSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
 

  for (int i = 0; i<numberOfData; i++){
    stateOut[i] = stateIn[i];
  }
 
  for (int i = 0; i< numberOfVariables; i++){
  fluxOut[i] =  fluxIn[i];
 }
 
}

exahype::solvers::Solver::RefinementControl Linear::MyLinearWaveSolver::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void Linear::MyLinearWaveSolver::eigenvalues(const double* const Q,const int d,double* const lambda) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 3 + 0
  
  // @todo Please implement/augment if required
  double cp = Q[4];
  //double cp = 1.0;
  
 
  lambda[0] =  cp;
  lambda[1] = -cp;
  lambda[2] =  0.0;
}


void Linear::MyLinearWaveSolver::flux(const double* const Q,double** const F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 3 + 0
  
  // @todo Please implement/augment if required
  double p = Q[0];
  double u = Q[1];
  double v = Q[2];
  
  F[0][0] = -u;
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  
  F[1][0] = -v;
  F[1][1] = 0.0;
  F[1][2] = 0.0;
  
}



void  Linear::MyLinearWaveSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // @todo Please implement/augment if required
  static tarch::logging::Log _log("MyLinearSolver::nonConservativeProduct");

  double p_x = gradQ[0];
  double u_x = gradQ[1];
  double v_x = gradQ[2];

  double p_y = gradQ[3];
  double u_y = gradQ[4];
  double v_y = gradQ[5];
  
  BgradQ[0] = 0.0; 
  BgradQ[1] = -p_x;
  BgradQ[2] = 0.0;

  BgradQ[3]= 0.0; 
  BgradQ[4]= 0.0;
  BgradQ[5]= -p_y;
  
}


    /**
     * @TODO LR : document
     */
void Linear::MyLinearWaveSolver::multiplyMaterialParameterMatrix(const double* const Q, double* const rhs) {
  // @todo Please implement/augment if required

  //double cp = 1.0;
  //double rho = 1.0;
 
  double rho = Q[3];
  double cp  = Q[4];
  double lam = rho*cp*cp;
  
  rhs[0] = lam*rhs[0];
  rhs[1] = 1.0/rho*rhs[1];
  rhs[2] = 1.0/rho*rhs[2];;

  rhs[3]= lam*rhs[3];
  rhs[4]= 1.0/rho*rhs[4];
  rhs[5]= 1.0/rho*rhs[5];
  
}


void Linear::MyLinearWaveSolver::coefficientMatrix(const double* const Q,const int d,double* const Bn){
  static tarch::logging::Log _log("MyLinearWaveSolver::coefficientMatrix");

  double cp = Q[4];
  double rho = Q[3];

  //double cp = 1.0;
  //double rho = 1.0;
  double lam = rho*cp*cp;

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
  
  B1[0][1] = -lam; 
  B1[1][0] = -1/rho;
  
  B2[0][2] = -lam; 
  B2[2][0] = -1/rho;
  
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      Bn[i*3+ j] = nv[0]*B1[i][j] + nv[1]*B2[i][j];
    }
  }

}


void Linear::MyLinearWaveSolver::pointSource(const double* const x,const double t,const double dt, double* const forceVector, double* const x0, int n){
  static tarch::logging::Log _log("MyLinearWaveSolver::pointSource");
  double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.7;
  //double t0 = 0.1;
  double f = 0.0;
  double M0 = 1000.0;
  

  //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-pow((t-t0)/(std::sqrt(2)*sigma),2.0)));

  f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));


  //f = M0*t/(t0*t0)*std::exp(-t/t0);

  if(n == 0){
  
    x0[0] = 5.0;
    x0[1] = 5.0;
    
    forceVector[0] = 1.*f;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
  }else if(n == 1){
  
    x0[0] = 5.0;
    x0[1] = 2.5;
    
    forceVector[0] = 0.*f;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;

  }    
}

void Linear::MyLinearWaveSolver::riemannSolver_BC0(double v, double sigma, double z,  double r, double& v_hat, double& sigma_hat){
  
   double p = 0.5*(z*v + sigma);

   v_hat = (1+r)/z*p;
   sigma_hat = (1-r)*p;

}


void Linear::MyLinearWaveSolver::riemannSolver_BCn(double v,double sigma, double z, double r, double& v_hat, double& sigma_hat){
  
   double q = 0.5*(z*v - sigma);

   v_hat = (1+r)/z*q;
   sigma_hat = -(1-r)*q;

}


void Linear::MyLinearWaveSolver::riemannSolver_Nodal(double v_p,double v_m, double sigma_p, double sigma_m, double z_p , double z_m, double& v_hat_p , double& v_hat_m, double& sigma_hat_p, double& sigma_hat_m){
   double p=0;
   double q=0;
   double phi=0;
   double v_hat=0;
   double eta=0;

   p = z_m*v_p + sigma_p;
   q = z_p*v_m - sigma_m;

   eta=(z_p*z_m)/(z_p+z_m);

   phi= eta*(p/z_p - q/z_m);

   sigma_hat_p=phi;
   sigma_hat_m=phi;

   v_hat_m=(p-phi)/z_p;
   v_hat_p=(q+phi)/z_m;

 }



// void Linear::MyLinearWaveSolver::riemannSolver(double* const FL,double* const FR,const double* const QL,const double* const QR,const double dt,const int normalNonZeroIndex, bool isBoundaryFace, int faceIndex){

//   constexpr int numberOfVariables  = MyLinearWaveSolver::NumberOfVariables;
//   constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
//   constexpr int numberOfParameters = MyLinearWaveSolver::NumberOfParameters;
//   constexpr int numberOfData       = numberOfVariables+numberOfParameters;
//   constexpr int basisSize          = MyLinearWaveSolver::Order+1;
//   constexpr int order              = basisSize - 1;
  
//   kernels::idx2 idx_QLR(basisSize, numberOfData);
  
//   kernels::idx2 idx_FLR(basisSize, numberOfVariables);
  
//   double n[2]={0,0};
//   n[normalNonZeroIndex]=1;
  
//   double cp_p;
//   double rho_p;
//   double lam_p;

//   double cp_m;
//   double rho_m;
//   double lam_m;
  
//   double n_p[2]={0,0};
//   double n_m[2]={0,0};

//   double m_p[2]={0,0};
//   double m_m[2]={0,0};
  
//   double norm_p_qr;
//   double norm_m_qr;
  
//   //std::cout<<" riemann solver called"<<std::endl;
//   //std::exit(-1);
  

//   for (int k = 0; k < 2; k++){
    
//     n_m[k] = n[k];
//     n_p[k] = n[k];
    
//   }
  
//   norm_m_qr = 1.0;
//   norm_p_qr = 1.0;   
  
  
//   for (int i = 0; i < basisSize; i++) {
    
//     rho_m = QL[idx_QLR(i,3)];
//     cp_m  = QL[idx_QLR(i,4)];
//     lam_m = cp_m*cp_m*rho_m;

//     rho_p = QR[idx_QLR(i,3)];
//     cp_p  = QR[idx_QLR(i,4)];
//     lam_p = cp_p*cp_p*rho_p;
    
    
    
//     double v_m=QL[idx_QLR(i,1)]*n_m[0]+QL[idx_QLR(i,2)]*n_m[1];
//     double v_p=QR[idx_QLR(i,1)]*n_p[0]+QR[idx_QLR(i,2)]*n_p[1];
    
//     double sigma_m = QL[idx_QLR(i,0)];
//     double sigma_p = QR[idx_QLR(i,0)];
    
    
//     double z_p=rho_p*cp_p;
//     double z_m=rho_m*cp_m;
    
//     double v_hat_p=0;
//     double v_hat_m=0;
//     double sigma_hat_p=0;
//     double sigma_hat_m=0;

//     double r = 0.;

//     if (faceIndex  == 2) { r = 1.0;} 
    
//     if (isBoundaryFace) { // external boundaries
//       if (faceIndex % 2  == 0) {
	
// 	riemannSolver_BC0(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
// 	riemannSolver_BC0(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
//       }
      
      
//       if  (faceIndex % 2  == 1) {
	
	
// 	riemannSolver_BCn(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
// 	riemannSolver_BCn(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
//       }
      
//     }
//     else {// inter-element boundaries
//       riemannSolver_Nodal(v_p,v_m, sigma_p, sigma_m, z_p , z_m, v_hat_p , v_hat_m, sigma_hat_p, sigma_hat_m);
//     }
        
//     FR[idx_FLR(i,0)] = -norm_p_qr*0.5*lam_p*((v_p-v_hat_p) - (sigma_p-sigma_hat_p)/z_p);
//     FL[idx_FLR(i,0)] =  norm_m_qr*0.5*lam_m*((v_m-v_hat_m) + (sigma_m-sigma_hat_m)/z_m);
    
//     FR[idx_FLR(i,1)] = norm_p_qr*0.5/rho_p*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n_p[0];
//     FL[idx_FLR(i,1)] = norm_m_qr*0.5/rho_m*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n_m[0];
    
//     FR[idx_FLR(i,2)] = norm_p_qr*0.5/rho_p*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n_p[1];
//     FL[idx_FLR(i,2)] = norm_m_qr*0.5/rho_m*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n_m[1];
        
//   }
  
// }

