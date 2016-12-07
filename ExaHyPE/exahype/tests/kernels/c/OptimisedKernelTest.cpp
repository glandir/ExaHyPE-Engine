/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#ifdef TEST_OPT_KERNEL

#include <sstream>
#include <iomanip> 
#include <random>
#include <cstring>

#include "exahype/tests/kernels/c/OptimisedKernelTest.h"


registerTest(exahype::tests::c::OptimisedKernelTest)

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", off)
#endif

tarch::logging::Log exahype::tests::c::OptimisedKernelTest::_log( "exahype::tests::c::OptimisedKernelTest" );

namespace exahype {
namespace tests {
namespace c {

int OptimisedKernelTest::_numberOfVariables;
int OptimisedKernelTest::_basisSize;
int OptimisedKernelTest::_order; 
int OptimisedKernelTest::_dim; 
bool OptimisedKernelTest::_isLinear;
 
const double OptimisedKernelTest::eps  = 1.0e-14;
const double OptimisedKernelTest::eps2 = 1.0e-12; //for known reordered operations
const double OptimisedKernelTest::eps3 = 1.0e-5;  //for known unprecise operations (to fix)
#ifdef Dim2
const std::string OptimisedKernelTest::dim = "2";
#endif
#ifdef Dim3
const std::string OptimisedKernelTest::dim = "3";
#endif

OptimisedKernelTest::OptimisedKernelTest()
    : tarch::tests::TestCase("exahype::tests::c::OptimisedKernelTest") {
      
  _numberOfVariables = kernels::aderdg::optimised::getNumberOfVariable();
  _basisSize         = kernels::aderdg::optimised::getBasisSize();
  _order             = kernels::aderdg::optimised::getBasisSize() - 1;
  _dim               = kernels::aderdg::optimised::getDimension();
  _isLinear          = kernels::aderdg::optimised::isLinear();
}

OptimisedKernelTest::~OptimisedKernelTest() {}


void OptimisedKernelTest::adjustedSolutionValues(const double* const x, const double w, const double t, const double dt, double* Q) {
  adjustedSolutionValues_Euler(x,w,t,dt,Q);
}

void OptimisedKernelTest::flux(const double* const Q, double** F) {
  flux_Euler(Q, F);
}

void OptimisedKernelTest::source(const double* const Q, double* S) {
  source_Euler(Q,S);
}

void OptimisedKernelTest::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  eigenvalues_Euler(Q, normalNonZeroIndex, lambda);
}

void OptimisedKernelTest::ncp(const double* const Q, const double* const gradQ, double* BgradQ) {
	std::memset(BgradQ, 0, _numberOfVariables * sizeof(double));
}

void OptimisedKernelTest::matrixb(const double* const Q, const int normalNonZero, double* Bn) {
	std::memset(Bn, 0, _numberOfVariables * sizeof(double));
}


void OptimisedKernelTest::fluxSplitted(const double* const Q, double* f, double* g
#ifdef Dim3
            , double* h
#endif    
) {
  
#ifdef Dim3
  double* F[3];
  F[2] = h;
#else
  double* F[2];
#endif    
  F[0] = f;
  F[1] = g;
   
  flux(Q, F);
}


void OptimisedKernelTest::adjustedSolutionValues_Euler(const double* const x,
                                                  const double w,
                                                  const double t,
                                                  const double dt, double* Q) {

  double GAMMA = 1.4;
  int repetition = _numberOfVariables / 5;
  int rest = _numberOfVariables % 5;
  
  for(int i=0; i<repetition; i++) {
    Q[5*i+0] = 1.;
    Q[5*i+1] = 0.;
    Q[5*i+2] = 0.;
    Q[5*i+3] = 0.;  
    Q[5*i+4] = 1. / (GAMMA -1) +
          std::exp(-((x[0] -0.5) *(x[0] -0.5) + (x[1] -0.5) *(x[1] -0.5) 
#ifdef Dim3
            + (x[2] -0.5) *(x[2] -0.5)
#endif        
          ) /
          (0.5 *0.5
#ifdef Dim3
            *0.5
#endif        
          )) *
          1.0e-1;
  }
  
  //fill the rest with 1
  for(int i=0; i<rest; i++) {
    Q[repetition*5+i] = 1.0;
  }
}


//TODO JMG adapt for nVar != 5
void OptimisedKernelTest::flux_Euler(const double* const Q, double** F) {
  // Dimensions             = 2/3
  // Number of variables    = 5 (#unknowns + #parameters)

  const double GAMMA = 1.4;

  const double irho = 1.0 / Q[0];
#ifdef Dim2
  const double p =
      (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);
#elif Dim3
  const double p =
      (GAMMA - 1) *
      (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * irho);
#else
#error Dim2 or Dim3 must be defined
#endif

  double* f = F[0];
  double* g = F[1];

  // f
  f[0] = Q[1];
  f[1] = irho * Q[1] * Q[1] + p;
  f[2] = irho * Q[1] * Q[2];
  f[3] = irho * Q[1] * Q[3];
  f[4] = irho * Q[1] * (Q[4] + p);
  // g
  g[0] = Q[2];
  g[1] = irho * Q[2] * Q[1];
  g[2] = irho * Q[2] * Q[2] + p;
  g[3] = irho * Q[2] * Q[3];
  g[4] = irho * Q[2] * (Q[4] + p);

#ifdef Dim3
  double* h = F[2];
  // h
  h[0] = Q[3];
  h[1] = irho * Q[3] * Q[1];
  h[2] = irho * Q[3] * Q[2];
  h[3] = irho * Q[3] * Q[3] + p;
  h[4] = irho * Q[3] * (Q[4] + p);
#endif
}


//TODO JMG adapt for nVar != 5
void OptimisedKernelTest::eigenvalues_Euler(const double* const Q,
                                       const int normalNonZeroIndex,
                                       double* lambda) {
  // Dimensions             = 2/3
  // Number of variables    = 5 (#unknowns + #parameters)
  const double GAMMA = 1.4;

  double irho = 1.0 / Q[0];

#ifdef Dim2
  double p = (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);
#elif Dim3
  double p = (GAMMA - 1) *
             (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * irho);
#else
#error Dim2 or Dim3 must be defined
#endif

  double u_n = Q[normalNonZeroIndex + 1] * irho;
  double c = std::sqrt(GAMMA * p * irho);

  lambda[0] = u_n - c;
  lambda[1] = u_n;
  lambda[2] = u_n;
  lambda[3] = u_n;
  lambda[4] = u_n + c;
}

//TODO JMG adapt to nVar != 5
void OptimisedKernelTest::source_Euler(const double* const Q, double* S) {
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
}

int OptimisedKernelTest::getNumberOfVariables() {
  return _numberOfVariables;
}

int OptimisedKernelTest::getNodesPerCoordinateAxis() {
  return _basisSize; //basisSize
}



void OptimisedKernelTest::run() {
  _log.info("OptimisedKernelTest::run()", "OptimisedKernelTest is active");

  _dt = 0.01;
  
#ifdef Dim2    
  _luhSize = _numberOfVariables*_basisSize*_basisSize;
  const int basisSizePowDim = _basisSize * _basisSize;
#else
  _luhSize = _numberOfVariables*_basisSize*_basisSize*_basisSize;
  const int basisSizePowDim = _basisSize * _basisSize * _basisSize;
#endif
  _luh  = new double[_luhSize]();
  _lduh = new double[_luhSize]();
  _lFhi_gen = new double[basisSizePowDim*_numberOfVariables*(_dim+1)]();
  _lFhi_opt = new double[basisSizePowDim*_numberOfVariables*(_dim+1)]();
  _lQhi_gen = new double[basisSizePowDim*_numberOfVariables]();
  _lQhi_opt = new double[basisSizePowDim*_numberOfVariables]();

  testMethod(testSolutionAdjustment); //initialize _luh
  testMethod(testStableTimeStep); //initialize dt
  if(!_isLinear) {
    testMethod(testSpaceTimePredictorNonLinear);
  }
  testMethod(testVolumeIntegral);
  testMethod(testSolutionUpdate);


  delete[] _luh;
  delete[] _lduh;
  delete[] _lFhi_gen;
  delete[] _lFhi_opt;
  delete[] _lQhi_gen;
  delete[] _lQhi_opt;
}


void OptimisedKernelTest::testSolutionAdjustment() {
  std::ostringstream out;
  out << "Test solutionAdjustment with gaussian pulse on Q[0], ORDER="<< _order <<", NVAR=" << _numberOfVariables;
  logInfo("OptimisedKernelTest::testSolutionAdjustment()", out.str());
  
#ifdef Dim2    
  const double dx[2] = {1.0, 1.0};
  const double center[2] = {0.5, 0.5};
#else
  const double dx[3] = {1.0, 1.0, 1.0};
  const double center[3] = {0.5, 0.5, 0.5};
#endif 
  
  double* luh_optimised = new double[_luhSize]();

  double t = 0.0;
  double dt = 0.0;
  
  kernels::aderdg::generic::c::solutionAdjustment<OptimisedKernelTest>( *this, _luh, center[0], dx[0], t, dt );
  kernels::aderdg::optimised::solutionAdjustment<OptimisedKernelTest::adjustedSolutionValues>( luh_optimised, center[0], dx[0], t, dt );
   
  for(int i=0; i<_luhSize; i++) {
    validateNumericalEqualsWithEps(luh_optimised[i], _luh[i], eps);
  }
  
  delete[] luh_optimised;
}


void OptimisedKernelTest::testStableTimeStep() {
  logInfo("OptimisedKernelTest::testStableTimeStep()", "Test stableTimeStep");
#ifdef Dim2    
  const double dx[2] = {1.0, 1.0};
#else
  const double dx[3] = {1.0, 1.0, 1.0};
#endif   
  double* tempEigenvalues = new double[_numberOfVariables];
  
  double dt_opt = kernels::aderdg::optimised::stableTimeStepSize<eigenvalues>( _luh, dx[0] );
  _dt = kernels::aderdg::generic::c::stableTimeStepSize<OptimisedKernelTest>( *this, _luh, tempEigenvalues, dx[0] );
  
  validateNumericalEqualsWithEps(dt_opt, _dt, eps);
  
  delete[] tempEigenvalues;
}
 

void OptimisedKernelTest::testSpaceTimePredictorNonLinear() {
    std::ostringstream out;
  out << "Test spaceTimePredictorNonLinear with gaussian pulse on Q[0], ORDER="<< _order <<", NVAR=" << _numberOfVariables;
  logInfo("OptimisedKernelTest::testSpaceTimePredictorNonLinear()", out.str());

#ifdef Dim2    
  const double dx[2] = {1.0, 1.0};
  const double center[2] = {0.5, 0.5};
  const int basisSizePowDim = _basisSize * _basisSize;
  const int sizeLFi = basisSizePowDim * _basisSize * (_dim+1) * _numberOfVariables; // idx_lFi(t, y, x, nDim + 1 for Source, nVar)
  const int sizeFhbnd = 2 * 2 * _basisSize * _numberOfVariables;
#else
  const double dx[3] = {1.0, 1.0, 1.0};
  const double center[3] = {0.5, 0.5, 0.5};
  const int basisSizePowDim = _basisSize * _basisSize * _basisSize;
  const int sizeLFi = basisSizePowDim * _basisSize * (_dim+1) * _numberOfVariables; // idx_lFi(t, z, y, x, nDim + 1 for Source, nVar)
  const int sizeFhbnd = 2 * 3 * _basisSize * _basisSize * _numberOfVariables;
#endif 
  const int sizeLQi = basisSizePowDim * _basisSize * _numberOfVariables; // idx_lQi(z?,y,x,t,nVar)
  
  //tmp storage
  double** tempSpaceTimeUnknowns = new double*[4];
  for (int i=0; i<4; ++i) { // 0: lQi, 1: lQi_old, 2: rhs, 3: rhs_0 (spaceTimePredictorNonLinear.cpph)
    tempSpaceTimeUnknowns[i] = new double[sizeLQi]();
  }
  double** tempSpaceTimeFluxUnknowns = new double*[2];
  for (int i=0; i<2; ++i) { //0: lFi, 1: gradQ
    tempSpaceTimeFluxUnknowns[i] = new double[sizeLFi]();
  }

  double* tempStateSizedVectors = new double[_numberOfVariables]();
  
  //out
  double* lQhbnd_generic = new double[sizeFhbnd]();
  double* lFhbnd_generic = new double[sizeFhbnd]();
  double* lQhbnd_optimised = new double[sizeFhbnd]();
  double* lFhbnd_optimised = new double[sizeFhbnd]();
  
  
  //compute
  kernels::aderdg::optimised::picardLoop<fluxSplitted>( tempSpaceTimeUnknowns[0], tempSpaceTimeFluxUnknowns[0], _luh, dx[0], _dt );
  kernels::aderdg::optimised::predictor( _lQhi_opt, _lFhi_opt, tempSpaceTimeUnknowns[0], tempSpaceTimeFluxUnknowns[0] );
  kernels::aderdg::optimised::extrapolator( lQhbnd_optimised, lFhbnd_optimised, _lQhi_opt, _lFhi_opt );
  
  kernels::aderdg::generic::c::spaceTimePredictorNonlinear<OptimisedKernelTest>( *this, lQhbnd_generic, lFhbnd_generic, tempSpaceTimeUnknowns, tempSpaceTimeFluxUnknowns, _lQhi_gen, _lFhi_gen, tempStateSizedVectors, _luh, dx[0], _dt );
  
  //compare
#ifdef Dim2  
  kernels::idx3 idx_bnd_gen(2 * _dim, _basisSize, _numberOfVariables);
  kernels::idx3 idx_bnd_opt(2 * _dim, _numberOfVariables, _basisSize);
  for(int i=0; i<2*_dim; i++) {
    for(int j=0; j<_basisSize; j++) {
      for(int k=0; k<_numberOfVariables; k++) {
        validateNumericalEqualsWithEps(lQhbnd_optimised[idx_bnd_opt(i,k,j)], lQhbnd_generic[idx_bnd_gen(i,j,k)], eps3);
        validateNumericalEqualsWithEps(lFhbnd_optimised[idx_bnd_opt(i,k,j)], lFhbnd_generic[idx_bnd_gen(i,j,k)], eps3);
      }
    }    
  }
  
  kernels::idx3 idx_lQhi(_basisSize, _basisSize, _numberOfVariables);
  kernels::idx3 idx_lFhi(_basisSize, _basisSize, _numberOfVariables);
  const int lFhi_shift = _basisSize*_basisSize*_numberOfVariables;
  for(int i=0; i<_basisSize; i++) {
    for(int j=0; j<_basisSize; j++) {
      for(int k=0; k<_numberOfVariables; k++) {
        validateNumericalEqualsWithEps(_lQhi_opt[idx_lQhi(i,j,k)], _lQhi_gen[idx_lQhi(i,j,k)], eps3);
        validateNumericalEqualsWithEps(_lFhi_opt[idx_lFhi(i,j,k)           ], _lFhi_gen[idx_lFhi(i,j,k)           ], eps3); //lFhi_x
        validateNumericalEqualsWithEps(_lFhi_opt[idx_lFhi(i,j,k)+lFhi_shift], _lFhi_gen[idx_lFhi(i,j,k)+lFhi_shift], eps3); //lFhi_y
      }
    }    
  }
  
#else
  kernels::idx4 idx_bnd_gen(2 * _dim, _basisSize, _basisSize, _numberOfVariables);
  kernels::idx4 idx_bnd_opt(2 * _dim, _numberOfVariables, _basisSize, _basisSize);  
  for(int i=0; i<2*_dim; i++) {
    for(int j=0; j<_basisSize; j++) {
      for(int k=0; k<_basisSize; k++) {
        for(int l=0; l<_numberOfVariables; l++) {
          validateNumericalEqualsWithEps(lQhbnd_optimised[idx_bnd_opt(i,l,j,k)], lQhbnd_generic[idx_bnd_gen(i,j,k,l)], eps3);
          validateNumericalEqualsWithEps(lFhbnd_optimised[idx_bnd_opt(i,l,j,k)], lFhbnd_generic[idx_bnd_gen(i,j,k,l)], eps3);
        }
      }
    }    
  }
  
  kernels::idx4 idx_lQhi(_basisSize, _basisSize, _basisSize, _numberOfVariables);
  kernels::idx4 idx_lFhi(_basisSize, _basisSize, _basisSize, _numberOfVariables);
  const int lFhi_shift = _basisSize*_basisSize*_basisSize*_numberOfVariables;
  for(int i=0; i<_basisSize; i++) {
    for(int j=0; j<_basisSize; j++) {
      for(int k=0; k<_basisSize; k++) {
        for(int l=0; l<_numberOfVariables; l++) {
          validateNumericalEqualsWithEps(_lQhi_opt[idx_lQhi(i,j,k,l)], _lQhi_gen[idx_lQhi(i,j,k,l)], eps3);
          validateNumericalEqualsWithEps(_lFhi_opt[idx_lFhi(i,j,k,l)+0*lFhi_shift], _lFhi_gen[idx_lFhi(i,j,k,l)+0*lFhi_shift], eps3); //lFhi_x
          validateNumericalEqualsWithEps(_lFhi_opt[idx_lFhi(i,j,k,l)+1*lFhi_shift], _lFhi_gen[idx_lFhi(i,j,k,l)+1*lFhi_shift], eps3); //lFhi_y
          validateNumericalEqualsWithEps(_lFhi_opt[idx_lFhi(i,j,k,l)+2*lFhi_shift], _lFhi_gen[idx_lFhi(i,j,k,l)+2*lFhi_shift], eps3); //lFhi_z
        }
      }
    }    
  }
#endif
  
  
  //clean
  for (int i=0; i<4; ++i) { 
    delete[] tempSpaceTimeUnknowns[i];
  }
  delete[] tempSpaceTimeUnknowns;
  for (int i=0; i<2; ++i) { 
    delete[] tempSpaceTimeFluxUnknowns[i];
  }
  delete[] tempSpaceTimeFluxUnknowns;
  delete[] tempStateSizedVectors;
  delete[] lQhbnd_generic;
  delete[] lFhbnd_generic;
  delete[] lQhbnd_optimised;
  delete[] lFhbnd_optimised;
}


void OptimisedKernelTest::testVolumeIntegral() {
  std::ostringstream out;
  out << "Test testVolumeIntegral with random values, ORDER="<< _order <<", NVAR=" << _numberOfVariables;
  logInfo("OptimisedKernelTest::testVolumeIntegral()", out.str());
  
#ifdef Dim2    
  const double dx[2] = {1.0, 1.0};
#else
  const double dx[3] = {1.0, 1.0, 1.0};
#endif

  double* lduh_opt = new double[_luhSize];
  
  kernels::aderdg::generic::c::volumeIntegralNonlinear( _lduh, _lFhi_gen, dx[0], _numberOfVariables, 0, _basisSize );
  kernels::aderdg::optimised::volumeIntegral( lduh_opt, _lFhi_gen, dx[0] );
  
  for(int i=0; i<_luhSize; i++) {
    validateNumericalEqualsWithEps(lduh_opt[i], _lduh[i], eps2);
  }
  
  delete[] lduh_opt;
}


void OptimisedKernelTest::testSolutionUpdate() {
  std::ostringstream out;
  out << "Test testSolutionUpdate with random values, ORDER="<< _order <<", NVAR=" << _numberOfVariables;
  logInfo("OptimisedKernelTest::testSolutionUpdate()", out.str());

  const double dt = 0.05;
  double* luh_generic = new double[_luhSize];
  double* luh_optimised = new double[_luhSize];
  double* lduh_generic = new double[_luhSize];
  double* lduh_optimised = new double[_luhSize];
  
  std::random_device rd; //to generate a randome seed
  std::mt19937 mt(rd()); //mersenne twister random number generator with random seed
  std::uniform_real_distribution<double> dist(-1.0, 1.0); // [-1.0,1.0)
  
  for(int i=0; i<_luhSize; i++) {
    luh_generic[i] = dist(mt);
    lduh_generic[i] = dist(mt);
  }

  std::memcpy(luh_optimised, luh_generic, _luhSize*sizeof(double));
  std::memcpy(lduh_optimised, lduh_generic, _luhSize*sizeof(double));

  kernels::aderdg::generic::c::solutionUpdate( luh_generic, lduh_generic, dt, _numberOfVariables, 0, _basisSize );
  kernels::aderdg::optimised::solutionUpdate( luh_optimised, lduh_optimised, dt );
  
  for(int i=0; i<_luhSize; i++) {
    validateNumericalEqualsWithEps(luh_optimised[i], luh_generic[i], eps2);
  }
  
  delete[] luh_generic;
  delete[] luh_optimised;
  delete[] lduh_generic;
  delete[] lduh_optimised;

}





}  // namespace c
}  // namespace tests
}  // namespace exahype

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", on)
#endif

#endif //TEST_OPT_KERNEL