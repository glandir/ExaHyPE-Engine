// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ConstraintsWriter.h"

#include "kernels/GaussLegendreBasis.h"
#include "kernels/KernelUtils.h"

#include "PDE.h" // ADMConstraints()
#include "kernels/aderdg/generic/c/computeGradients.cpph"


using namespace kernels;
using namespace exahype::plotters::ascii;

FOCCZ4::ConstraintsWriter::ConstraintsWriter() :
	exahype::plotters::ADERDG2UserDefined::ADERDG2UserDefined(),
	constraintReductions("output/constraint-", ".asc", exahype::plotters::ascii::parallel::postprocess)
{
	constraintReductions.add(0, "hamiltonian");
	constraintReductions.add(1, "momentum1");
	constraintReductions.add(2, "momentum2");
	constraintReductions.add(3, "momentum3");
	constraintReductions.add(4, "det1");
	constraintReductions.add(5, "traceA");
}

void FOCCZ4::ConstraintsWriter::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& dx, // sizeOfPatch
    double* u,
    double timeStamp) {

	double gradQ[basisSize3 * DIMENSIONS * numberOfVariables];
	
	kernels::aderdg::generic::c::computeGradQ<FOCCZ4::AbstractFOCCZ4Solver_ADERDG>(gradQ, u, dx);

	// volume form for integration
	double scaling = tarch::la::volume(dx);
	double wx, wy, wz, constraints[6]; // 6 constraints!

	idx5 idx_gradQ(basisSize, basisSize, basisSize, DIMENSIONS, numberOfVariables);
	idx4 idx_u(basisSize, basisSize, basisSize, numberOfVariables);

	for (int iz = 0; iz < basisSize; iz++) {
	for (int iy = 0; iy < basisSize; iy++) {
	for (int ix = 0; ix < basisSize; ix++) {
		// Gauss-Legendre weights from pos argument
		wx = kernels::legendre::weights[FOCCZ4::AbstractFOCCZ4Solver_ADERDG::Order][ix];
		wy = kernels::legendre::weights[FOCCZ4::AbstractFOCCZ4Solver_ADERDG::Order][iy];
		wz = kernels::legendre::weights[FOCCZ4::AbstractFOCCZ4Solver_ADERDG::Order][iz];

		admconstraints_(constraints, &u[idx_u(iz,iy,ix,0)], &gradQ[idx_gradQ(iz,iy,ix,0,0)]);
		constraintReductions.addValue(constraints, scaling*wx*wy*wz);
	} // x
	} // y
	} // z

}


void FOCCZ4::ConstraintsWriter::startPlotting( double time) {
	constraintReductions.startRow(time);
}


void FOCCZ4::ConstraintsWriter::finishPlotting() {
	constraintReductions.finishRow();
}
