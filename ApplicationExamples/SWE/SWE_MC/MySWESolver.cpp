// This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
// 
// ========================
//   www.exahype.eu
// ========================

#include "MySWESolver.h"
#include "bathymetry.h"
#include "MySWESolver_Variables.h"
#include "kernels/KernelUtils.h"
#include "kernels/GaussLegendreQuadrature.h"
#include "readCsv.h"

double grav;
tarch::logging::Log SWE::MySWESolver::_log( "SWE::MySWESolver" );


void SWE::MySWESolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
    grav = 9.81;

    //for testing csv writer/reader
    const int nelem = 9*3;
    std::vector<double> a((nelem+1)*(nelem+1));
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,0.4);
    //create uniformly distributed pertubation for each nodal point
    for(int i = 0; i < nelem+1; i++){
        for(int j = 0; j < nelem+1; j++){
            a[(nelem+1)*i+j] = distribution(generator);
        }
    }
    writeCsv("Input/parameters.csv", a);

    std::vector<double> a_bath = {3.0,6.0,4.0,-1.0,-3.0,-6.0,18.0,1.0,-4.0,15.0,9.0,12.0,17.0,-18.0,24.0,6.0};
    writeCsv("Input/bathymetry.csv", a_bath);

    std::vector<std::vector<double>> a_measurements(4);
    a_measurements[0] = {0.3,0.6};
    a_measurements[1] = {0.2,0.6};
    a_measurements[2] = {0.1,0.6};
    a_measurements[3] = {0.1,0.5};
    writeCsv("Input/measurements.csv", a_measurements);
}

void SWE::MySWESolver::adjustSolution(double *luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
    // Dimensions                        = 2
    // Number of variables + parameters  = 4 + 0
    if (tarch::la::equals(t,0.0)) {
        constexpr int basisSize = MySWESolver::Order+1;
        constexpr int numberOfData=MySWESolver::NumberOfParameters+MySWESolver::NumberOfVariables;

        kernels::idx3 id_xyf(basisSize,basisSize,numberOfData);
        kernels::idx2 id_xy(basisSize,basisSize);

        int num_nodes = basisSize;

        double offset_x=center[0]-0.5*dx[0];
        double offset_y=center[1]-0.5*dx[1];

        for (int i=0; i< num_nodes; i++){
            for (int j=0; j< num_nodes; j++){

                double x  =  (offset_x+dx[0]*kernels::gaussLegendreNodes[basisSize-1][i]);
                double y  =  (offset_y+dx[1]*kernels::gaussLegendreNodes[basisSize-1][j]);

                double b =  SWE::bathymetry(x,y) + linearInterpolation(x,y);

                if(x < 0.5) {
                    luh[id_xyf(i,j,0)] = 2.0 - b; //h
                }
                else{
                    luh[id_xyf(i,j,0)] = 1.5 - b; //h
                }
                luh[id_xyf(i,j,1)] = 0; //hu
                luh[id_xyf(i,j,2)] = 0; //hv
                luh[id_xyf(i,j,3)] = b; //b
            }
        }
    }
}

void SWE::MySWESolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0
  // Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0
    //Wall
    std::copy_n(stateIn, NumberOfVariables, stateOut);
    stateOut[1+normalNonZero] =  -stateOut[1+normalNonZero];
    double _F[2][NumberOfVariables]={0.0};
    double* F[2] = {_F[0], _F[1]};
    flux(stateOut,F);
    std::copy_n(F[normalNonZero], NumberOfVariables, fluxOut);
}

exahype::solvers::Solver::RefinementControl SWE::MySWESolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void SWE::MySWESolver::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0
    ReadOnlyVariables vars(Q);
    Variables eigs(lambda);

    const double c = std::sqrt(grav*vars.h());
    const double ih = 1./vars.h();
    double u_n = Q[d + 1] * ih;

    eigs.h() = u_n + c;
    eigs.hu() = u_n - c;
    eigs.hv() = u_n;
    eigs.b() = 0.0;
}





void SWE::MySWESolver::flux(const double* const Q,double** F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0
    ReadOnlyVariables vars(Q);

    const double ih = 1./vars.h();

    double* f = F[0];
    double* g = F[1];

    f[0] = vars.hu();
    f[1] = vars.hu()*vars.hu()*ih + 0.5*grav*vars.h()*vars.h();
    f[2] = vars.hu()*vars.hv()*ih;
    f[3] = 0.0;

    g[0] = vars.hv();
    g[1] = vars.hu()*vars.hv()*ih;
    g[2] = vars.hv()*vars.hv()*ih + 0.5*grav*vars.h()*vars.h();
    g[3] = 0.0;
}



void  SWE::MySWESolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  // @todo Please implement/augment if required
    kernels::idx2 idx_gradQ(DIMENSIONS,NumberOfVariables);

    BgradQ[0] = 0.0;
    BgradQ[1] = grav*Q[0]*gradQ[idx_gradQ(0,3)];
    BgradQ[2] = grav*Q[0]*gradQ[idx_gradQ(1,3)];
    BgradQ[3] = 0.0;
}

