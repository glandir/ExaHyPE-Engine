// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ErrorWriter.h"
#include "PDE.h"
#include "kernels/GaussLegendreQuadrature.h"
#include <cmath>

#include "kernels/aderdg/generic/c/sizes.cpph"
#include "kernels/KernelUtils.h" // matrix indexing

GPRDIM::ErrorWriter::ErrorWriter() : exahype::plotters::ADERDG2UserDefined::ADERDG2UserDefined(),
	errors("output/error-"){
  // @TODO Please insert your code here.
  	char name[10];
	for(int m=0; m < nVar; m++) {
		pdevarname_(name,&m);
		//printf("***********************************************************");
		//printf(name);
		//printf("\n");
		errors.add(m, name);
    }
}


void GPRDIM::ErrorWriter::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
    double timeStamp) {
 // @TODO Please insert your code here.
  double dV=sizeOfPatch[0]*sizeOfPatch[1];
#if DIMENSIONS==3
	dV = dV * sizeOfPatch[2];
	kernels::idx4 id_xyz_dof(basisSize,basisSize,basisSize,nVar);
#else
	kernels::idx3 id_xy_dof(basisSize,basisSize,nVar);
#endif

	
	double localError[nVar]={0.};
	double w_x, w_y;
	double w_z = 1.0; //constant in two dimensions
	
	double pos[DIMENSIONS] = {0.};
	
	for(int i = 0; i < basisSize; i++){
		w_x= kernels::gaussLegendreWeights[order][i];
		pos[0] = kernels::gaussLegendreNodes[order][i]*sizeOfPatch[0]+offsetOfPatch[0];
		for(int j = 0; j <basisSize ; j++){
			w_y = kernels::gaussLegendreWeights[order][j];
			pos[1] = kernels::gaussLegendreNodes[order][j]*sizeOfPatch[1]+offsetOfPatch[1];
#if DIMENSIONS==3	
			for(int k = 0; k <basisSize ; k++){
				w_z = kernels::gaussLegendreWeights[order][k];
				pos[2] = kernels::gaussLegendreNodes[order][k]*sizeOfPatch[2]+offsetOfPatch[2];
#endif	
				double numerical[nVar];
#if DIMENSIONS==3
				getNumericalSolution(numerical,&u[id_xyz_dof(i,j,k,0)]);
#else
				getNumericalSolution(numerical,&u[id_xy_dof(i,j,0)]);
#endif	
				double exact[nVar];
				getExactSolution(exact,pos,timeStamp);
				
						
				for(int m = 0; m <nVar ; m++){
					localError[m] += std::abs(numerical[m]-exact[m]) * w_x * w_y * w_z;
				}
#if DIMENSIONS==3				
			}
#endif			
		}
	}
  	errors.addValue(localError, dV);
}


void GPRDIM::ErrorWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
  errors.startRow(time);
}


void GPRDIM::ErrorWriter::finishPlotting() {
  // @TODO Please insert your code here.
  errors.finishRow();
}

void GPRDIM::ErrorWriter::getNumericalSolution(double* numerical, double* u){
	getnumericalsolution_(numerical,u);
}

void GPRDIM::ErrorWriter::getExactSolution(double* exact, double* pos, double timeStamp){
	getexactsolution_(exact,pos,&timeStamp);
}
