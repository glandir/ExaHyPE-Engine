#!/bin/env python
##
# @file This file is part of the ExaHyPE project.
# @author ExaHyPE Group (exahype@lists.lrz.de)
#
# @section LICENSE
#
# Copyright (c) 2016  http://exahype.eu
# All rights reserved.
#
# The project has received funding from the European Union's Horizon 
# 2020 research and innovation programme under grant agreement
# No 671698. For copyrights and licensing, please consult the webpage.
#
# Released under the BSD 3 Open Source License.
# For the full license text, see LICENSE.txt
#
#
# @section DESCRIPTION
#
# Generates the data structures for the system matrices
# for a specific configuration. In particular, it decides
# on padding.
#
# @note Utilises the same procedures as the generic code base
# to generate the matrices, but internally transposes the
# data or introduces padding.
#

import numpy as np
import Backend
import re
import sys
import os
dir = os.path.dirname(__file__)
sys.path.insert(1, dir+'/../Miscellaneous/aderdg')
import aderdg

class DGMatrixGenerator:
    m_config = {}

    # order of the approximation polynomial
    m_order      = -1

    # number of degrees of freedom
    m_nDof       = -1

    # number of dimensions we simulate
    m_nDim       = -1

    # quadrature nodes and weights mapped onto [0,1]
    m_xGPN       = []
    m_wGPN       = []

    # linear/nonlinear
    m_numerics  = ''

    # name of generated output files
    m_sourceName = "DGMatrices.cpp"
    m_headerName = "DGMatrices.h"


    def __init__(self, i_config, i_numerics):
        self.m_order     = i_config['nDof']-1
        self.m_config    = i_config
        self.m_nDim      = i_config['nDim']
        self.m_numerics  = i_numerics

        # compute the Gauss-Legendre weights
        x, w = np.polynomial.legendre.leggauss(self.m_order+1)
        # map onto [0,1]
        self.m_xGPN = 0.5*(x+1)
        self.m_wGPN = 0.5*w


    def generateCode(self):
        # padding of the system matrices is the same
        l_padHeight = Backend.getPadWidth(self.m_config['nDof'])
        l_padMatrix = np.zeros((l_padHeight, self.m_config['nDof']))

        l_matrices = {}

        # Kxi
        Kxi = aderdg.assembleStiffnessMatrix(self.m_xGPN, self.m_wGPN, self.m_order)
        l_linearisedKxi = np.concatenate((Kxi,l_padMatrix),axis=0).flatten('F')
        # Kxi is now
        # [ Kxi ]
        # [0...0]
        # [0...0]
        # stored in column major order
        l_matrices['Kxi'] = l_linearisedKxi

        # Kxi transposed
        # for the time being we need both Kxi and its tranpose
        Kxi_T = np.transpose(Kxi)
        l_linearised_Kxi_T = np.concatenate((Kxi_T,l_padMatrix),axis=0).flatten('F')
        l_matrices['Kxi_T'] = l_linearised_Kxi_T

        # iK1
        iK1 = np.transpose(np.linalg.inv(aderdg.assembleK1(Kxi, self.m_xGPN, self.m_order)))
        l_linearisediK1 = np.concatenate((iK1,l_padMatrix),axis=0).flatten('F')
        l_matrices['iK1'] = l_linearisediK1

        # dudx
        MM   = aderdg.assembleMassMatrix(self.m_xGPN, self.m_wGPN, self.m_order)
        dudx = aderdg.assembleDiscreteDerivativeOperator(MM,Kxi)
        # we multiply always from the left -> no transpose
        l_lineariseddudx  = np.concatenate((dudx,l_padMatrix),axis=0).flatten('F')
        l_matrices['dudx'] = l_lineariseddudx
        
        # dudx_t
        dudx_T = np.transpose(dudx)
        l_lineariseddudx_T  = np.concatenate((dudx_T,l_padMatrix),axis=0).flatten('F')
        l_matrices['dudx_T'] = l_lineariseddudx_T

        # tmp memory for *s*caled *m*atrices
        s_m = np.zeros(np.concatenate((Kxi,l_padMatrix),axis=0).shape)
        l_matrices['s_m'] = s_m.flatten('F')

        # tmp memory for *s*caled *v*ector
        s_v = np.zeros((1,Backend.getSizeWithPadding(self.m_config['nVar'])))
        l_matrices['s_v'] = s_v.flatten('F')

        # tmp memory for simd vectorisation on boundaries, sized as e.g. lQbndL
        # 2D: nVar * getSizeWithPadding(nDOF)
        # 3D: nVar * getSizeWithPadding(nDOF*nDOF)
        l_nTotalDof = self.m_config['nDof']**(self.m_nDim-1)
        l_nEntries  = Backend.getSizeWithPadding(self.m_config['nVar']) * Backend.getSizeWithPadding(l_nTotalDof)
        tmp_bnd = np.zeros((1,l_nEntries))
        l_matrices['tmp_bnd'] = tmp_bnd.flatten('F')


        # [FLCoeff 0...0]; [FRCoeff 0...0];
        # right now FLCoeff, FRCoeff no pad (gives no benefit w.r.t libxsmm)
        FLCoeff, _ = np.array(aderdg.BaseFunc1d(0.0, self.m_xGPN, self.m_order))
        FRCoeff, _ = np.array(aderdg.BaseFunc1d(1.0, self.m_xGPN, self.m_order))
        #l_paddedFLCoeff = np.pad(FLCoeff, [0, Backend.getPadWidth(self.m_config['nDof'])], mode='constant')
        #l_paddedFRCoeff = np.pad(FRCoeff, [0, Backend.getPadWidth(self.m_config['nDof'])], mode='constant')
        #l_matrices['FLCoeff'] = l_paddedFLCoeff
        #l_matrices['FRCoeff'] = l_paddedFRCoeff
        l_matrices['FLCoeff'] = FLCoeff
        l_matrices['FRCoeff'] = FRCoeff

        # F0 with padding for DSCAL
        F0, _ = np.array(aderdg.BaseFunc1d(0.0, self.m_xGPN, self.m_order))
        l_paddedF0 = np.pad(F0, [0, Backend.getPadWidth(self.m_config['nDof'])], mode='constant')
        l_matrices['F0'] = l_paddedF0

        self.__generateHeaderFile()
        self.__writeToFile(l_matrices)


    def __generateHeaderFile(self):
        l_sourceFile = open(self.m_headerName, 'a')
        # copied from Dominic's DGMatrices.h; matrices are linearised
        l_includeGuard = '#ifndef _EXAHYPE_KERNELS_ADERDG_OPTIMISED_DGMATRICES_H_\n'   \
                         '#define _EXAHYPE_KERNELS_ADERDG_OPTIMISED_DGMATRICES_H_\n\n'
        l_sourceFile.write(l_includeGuard)
        l_sourceFile.write('#include <set>\n\n')
        l_sourceFile.write('namespace kernels { \n'\
                           'namespace aderdg {\n'\
                           'namespace optimised {\n\n'\
                           'void initDGMatrices(const std::set<int>& orders);\n' \
                           'void freeDGMatrices(const std::set<int>& orders);\n\n' \
                           'extern double *Kxi;\n'      \
                           'extern double *Kxi_T;\n'    \
                           'extern double *iK1;\n'      \
                           'extern double *dudx;\n'     \
                           'extern double *dudx_T;\n'     \
                           'extern double *s_m;\n'      \
                           'extern double *s_v;\n'      \
                           'extern double *tmp_bnd;\n'  \
                           'extern double *F0; \n'      \
                           'extern double *FLCoeff;\n'  \
                           'extern double *FRCoeff;\n'  \
                           'extern double ***equidistantGridProjector1d;\n' \
                           'extern double **** fineGridProjector1d;\n\n}\n}\n}\n')
        # close include guard
        l_sourceFile.write('#endif /* _EXAHYPE_KERNELS_ADERDG_OPTIMISED_DGMATRICES_H_ */')
        l_sourceFile.close()


    def __writeToFile(self, i_matrices):
        l_sourceFile = open(self.m_sourceName, 'a')
        l_sourceFile.write('#include "kernels/aderdg/optimised/'+self.m_headerName+'"\n' \
                           '#include <mm_malloc.h> //g++\n\n')
        l_sourceFile.write('double* kernels::aderdg::optimised::Kxi;\n'     \
                           'double* kernels::aderdg::optimised::Kxi_T;\n'   \
                           'double* kernels::aderdg::optimised::iK1;\n'     \
                           'double* kernels::aderdg::optimised::dudx;\n'    \
                           'double* kernels::aderdg::optimised::dudx_T;\n'    \
                           'double* kernels::aderdg::optimised::s_m;\n'     \
                           'double* kernels::aderdg::optimised::s_v;\n'     \
                           'double* kernels::aderdg::optimised::tmp_bnd;\n' \
                           'double* kernels::aderdg::optimised::F0;\n'      \
                           'double* kernels::aderdg::optimised::FLCoeff;\n' \
                           'double* kernels::aderdg::optimised::FRCoeff;\n' \
                           'double*** kernels::aderdg::optimised::equidistantGridProjector1d;\n' \
                           'double**** kernels::aderdg::optimised::fineGridProjector1d;\n\n')


        l_sourceFile.write('void kernels::aderdg::optimised::freeDGMatrices(const std::set<int>& orders) {\n')
        for matrix in i_matrices:
            l_sourceFile.write('_mm_free('+str(matrix)+');\n')
        l_sourceFile.write('  // copied from generic DGMatrices.cpp\n'\
                           '  const int MAX_ORDER=9;\n'\
                           '  for (int ii = 0; ii < MAX_ORDER + 1; ii++) {\n' \
                           '    for (int jj = 0; jj < ii + 1; jj++) {\n'\
                           '      delete [] equidistantGridProjector1d[ii][jj];\n'\
                           '      delete [] fineGridProjector1d[ii][0][jj];\n'\
                           '      delete [] fineGridProjector1d[ii][1][jj];\n'\
                           '      delete [] fineGridProjector1d[ii][2][jj];\n'\
                           '    }\n'\
                           '    delete [] fineGridProjector1d[ii][0];\n'\
                           '    delete [] fineGridProjector1d[ii][1];\n'\
                           '    delete [] fineGridProjector1d[ii][2];\n'\
                           '    delete [] fineGridProjector1d[ii];\n'\
                           '    delete [] equidistantGridProjector1d[ii];\n'\
                           '  }\n'\
                           '  delete [] equidistantGridProjector1d;\n'\
                           '  delete [] fineGridProjector1d;\n')
        l_sourceFile.write('}\n\n')

        l_sourceFile.write('void kernels::aderdg::optimised::initDGMatrices(const std::set<int>& orders) {\n')

        for matrix in i_matrices:
            l_elemCount = i_matrices[matrix].size
            l_sourceFile.write(str(matrix)+' = (double *) _mm_malloc(sizeof(double)*'+str(l_elemCount)+', ALIGNMENT);\n')

        for matrix in i_matrices:
            for idx in range(0, i_matrices[matrix].size):
                l_sourceFile.write(str(matrix) + '['+str(idx)+'] = %.12e' % i_matrices[matrix].item(idx)+';\n')


        l_sourceFile.write('const int MAX_ORDER = 9;\n')
        l_sourceFile.write(
            '   equidistantGridProjector1d = new double**  [MAX_ORDER + 1];  // ***\n '\
            '   fineGridProjector1d        = new double*** [MAX_ORDER + 1];  // ****\n '\
            ' \n '\
            '   for (int ii = 0; ii < MAX_ORDER + 1; ii++) {\n '\
            '     equidistantGridProjector1d[ii] = new double* [ii + 1];\n '\
            '     fineGridProjector1d[ii] = new double** [3];\n '\
            '     fineGridProjector1d[ii][0] = new double* [ii+1];\n '\
            '     fineGridProjector1d[ii][1] = new double* [ii+1];\n '\
            '     fineGridProjector1d[ii][2] = new double* [ii+1];\n '\
            ' \n '\
            '     for (int jj = 0; jj < ii + 1; jj++) {\n '\
            '       equidistantGridProjector1d[ii][jj] = new double[ii + 1];\n '\
            '       fineGridProjector1d[ii][0][jj] = new double [ii+1];\n '\
            '       fineGridProjector1d[ii][1][jj] = new double [ii+1];\n '\
            '       fineGridProjector1d[ii][2][jj] = new double [ii+1];\n '\
            '     }\n '\
            '   }\n '\
            '   equidistantGridProjector1d[0][0][0] = 1.000000000000e+00;\n '\
            '   fineGridProjector1d[0][0][0][0] = 1.000000000000e+00;\n '\
            '   fineGridProjector1d[0][1][0][0] = 1.000000000000e+00;\n '\
            '   fineGridProjector1d[0][2][0][0] = 1.000000000000e+00;\n '\
            '   \n '\
            '   equidistantGridProjector1d[1][0][0] = 1.366025403784e+00;\n'\
            '   equidistantGridProjector1d[1][0][1] = -3.660254037844e-01;\n'\
            '   equidistantGridProjector1d[1][1][0] = -3.660254037844e-01;\n'\
            '   equidistantGridProjector1d[1][1][1] = 1.366025403784e+00;\n'\
            '   fineGridProjector1d[1][0][0][0] = 1.244016935856e+00;\n '\
            '   fineGridProjector1d[1][0][0][1] = 9.106836025230e-01;\n '\
            '   fineGridProjector1d[1][0][1][0] = -2.440169358563e-01;\n '\
            '   fineGridProjector1d[1][0][1][1] = 8.931639747704e-02;\n '\
            '   fineGridProjector1d[1][1][0][0] = 6.666666666667e-01;\n '\
            '   fineGridProjector1d[1][1][0][1] = 3.333333333333e-01;\n '\
            '   fineGridProjector1d[1][1][1][0] = 3.333333333333e-01;\n '\
            '   fineGridProjector1d[1][1][1][1] = 6.666666666667e-01;\n '\
            '   fineGridProjector1d[1][2][0][0] = 8.931639747704e-02;\n '\
            '   fineGridProjector1d[1][2][0][1] = -2.440169358563e-01;\n '\
            '   fineGridProjector1d[1][2][1][0] = 9.106836025230e-01;\n '\
            '   fineGridProjector1d[1][2][1][1] = 1.244016935856e+00;\n '\
            '   \n '\
			'   equidistantGridProjector1d[2][0][0] = 1.478830557701e+00;\n'\
			'   equidistantGridProjector1d[2][0][1] = -0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[2][0][2] = 1.878361089654e-01;\n'\
			'   equidistantGridProjector1d[2][1][0] = -6.666666666667e-01;\n'\
			'   equidistantGridProjector1d[2][1][1] = 1.000000000000e+00;\n'\
			'   equidistantGridProjector1d[2][1][2] = -6.666666666667e-01;\n'\
			'   equidistantGridProjector1d[2][2][0] = 1.878361089654e-01;\n'\
			'   equidistantGridProjector1d[2][2][1] = 0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[2][2][2] = 1.478830557701e+00;\n'\
            '   fineGridProjector1d[2][0][0][0] = 1.309811730779e+00;\n '\
            '   fineGridProjector1d[2][0][0][1] = 8.007018532823e-01;\n '\
            '   fineGridProjector1d[2][0][0][2] = 4.027030868966e-01;\n '\
            '   fineGridProjector1d[2][0][1][0] = -4.256271624011e-01;\n '\
            '   fineGridProjector1d[2][0][1][1] = 2.592592592593e-01;\n '\
            '   fineGridProjector1d[2][0][1][2] = 7.219234586974e-01;\n '\
            '   fineGridProjector1d[2][0][2][0] = 1.158154316219e-01;\n '\
            '   fineGridProjector1d[2][0][2][1] = -5.996111254156e-02;\n '\
            '   fineGridProjector1d[2][0][2][2] = -1.246265455940e-01;\n '\
            '   fineGridProjector1d[2][1][0][0] = 2.222222222222e-01;\n '\
            '   fineGridProjector1d[2][1][0][1] = -0.000000000000e+00;\n '\
            '   fineGridProjector1d[2][1][0][2] = -1.111111111111e-01;\n '\
            '   fineGridProjector1d[2][1][1][0] = 8.888888888889e-01;\n '\
            '   fineGridProjector1d[2][1][1][1] = 1.000000000000e+00;\n '\
            '   fineGridProjector1d[2][1][1][2] = 8.888888888889e-01;\n '\
            '   fineGridProjector1d[2][1][2][0] = -1.111111111111e-01;\n '\
            '   fineGridProjector1d[2][1][2][1] = 0.000000000000e+00;\n '\
            '   fineGridProjector1d[2][1][2][2] = 2.222222222222e-01;\n '\
            '   fineGridProjector1d[2][2][0][0] = -1.246265455940e-01;\n '\
            '   fineGridProjector1d[2][2][0][1] = -5.996111254156e-02;\n '\
            '   fineGridProjector1d[2][2][0][2] = 1.158154316219e-01;\n '\
            '   fineGridProjector1d[2][2][1][0] = 7.219234586974e-01;\n '\
            '   fineGridProjector1d[2][2][1][1] = 2.592592592593e-01;\n '\
            '   fineGridProjector1d[2][2][1][2] = -4.256271624011e-01;\n '\
            '   fineGridProjector1d[2][2][2][0] = 4.027030868966e-01;\n '\
            '   fineGridProjector1d[2][2][2][1] = 8.007018532823e-01;\n '\
            '   fineGridProjector1d[2][2][2][2] = 1.309811730779e+00;\n '\
            '   \n '\
			'   equidistantGridProjector1d[3][0][0] = 1.526788125457e+00;\n'\
			'   equidistantGridProjector1d[3][0][1] = -4.959188482477e-03;\n'\
			'   equidistantGridProjector1d[3][0][2] = -2.191327726011e-03;\n'\
			'   equidistantGridProjector1d[3][0][3] = -1.139171962820e-01;\n'\
			'   equidistantGridProjector1d[3][1][0] = -8.136324494869e-01;\n'\
			'   equidistantGridProjector1d[3][1][1] = 9.973040189732e-01;\n'\
			'   equidistantGridProjector1d[3][1][2] = 9.846497235301e-03;\n'\
			'   equidistantGridProjector1d[3][1][3] = 4.007615203117e-01;\n'\
			'   equidistantGridProjector1d[3][2][0] = 4.007615203117e-01;\n'\
			'   equidistantGridProjector1d[3][2][1] = 9.846497235300e-03;\n'\
			'   equidistantGridProjector1d[3][2][2] = 9.973040189732e-01;\n'\
			'   equidistantGridProjector1d[3][2][3] = -8.136324494869e-01;\n'\
			'   equidistantGridProjector1d[3][3][0] = -1.139171962820e-01;\n'\
			'   equidistantGridProjector1d[3][3][1] = -2.191327726011e-03;\n'\
			'   equidistantGridProjector1d[3][3][2] = -4.959188482477e-03;\n'\
			'   equidistantGridProjector1d[3][3][3] = 1.526788125457e+00;\n'\
            '   fineGridProjector1d[3][0][0][0] = 1.336580943538e+00;\n '\
            '   fineGridProjector1d[3][0][0][1] = 7.501737820555e-01;\n '\
            '   fineGridProjector1d[3][0][0][2] = 2.500683135133e-01;\n '\
            '   fineGridProjector1d[3][0][0][3] = 3.282922721363e-02;\n '\
            '   fineGridProjector1d[3][0][1][0] = -5.106599509271e-01;\n '\
            '   fineGridProjector1d[3][0][1][1] = 3.503991256637e-01;\n '\
            '   fineGridProjector1d[3][0][1][2] = 9.137546431442e-01;\n '\
            '   fineGridProjector1d[3][0][1][3] = 1.010071394721e+00;\n '\
            '   fineGridProjector1d[3][0][2][0] = 2.422582771987e-01;\n '\
            '   fineGridProjector1d[3][0][2][1] = -1.376638598071e-01;\n '\
            '   fineGridProjector1d[3][0][2][2] = -2.182390043732e-01;\n '\
            '   fineGridProjector1d[3][0][2][3] = -5.564103858755e-02;\n '\
            '   fineGridProjector1d[3][0][3][0] = -6.817926981005e-02;\n '\
            '   fineGridProjector1d[3][0][3][1] = 3.709095208782e-02;\n '\
            '   fineGridProjector1d[3][0][3][2] = 5.441604771565e-02;\n '\
            '   fineGridProjector1d[3][0][3][3] = 1.274041665251e-02;\n '\
            '   fineGridProjector1d[3][1][0][0] = -3.535004259642e-02;\n '\
            '   fineGridProjector1d[3][1][0][1] = -9.286838847674e-02;\n '\
            '   fineGridProjector1d[3][1][0][2] = -7.126778652900e-02;\n '\
            '   fineGridProjector1d[3][1][0][3] = -1.767502129821e-02;\n '\
            '   fineGridProjector1d[3][1][1][0] = 9.710461986761e-01;\n '\
            '   fineGridProjector1d[3][1][1][1] = 7.760907833372e-01;\n '\
            '   fineGridProjector1d[3][1][1][2] = 3.880453916686e-01;\n '\
            '   fineGridProjector1d[3][1][1][3] = 8.197886521854e-02;\n '\
            '   fineGridProjector1d[3][1][2][0] = 8.197886521854e-02;\n '\
            '   fineGridProjector1d[3][1][2][1] = 3.880453916686e-01;\n '\
            '   fineGridProjector1d[3][1][2][2] = 7.760907833372e-01;\n '\
            '   fineGridProjector1d[3][1][2][3] = 9.710461986761e-01;\n '\
            '   fineGridProjector1d[3][1][3][0] = -1.767502129821e-02;\n '\
            '   fineGridProjector1d[3][1][3][1] = -7.126778652900e-02;\n '\
            '   fineGridProjector1d[3][1][3][2] = -9.286838847674e-02;\n '\
            '   fineGridProjector1d[3][1][3][3] = -3.535004259642e-02;\n '\
            '   fineGridProjector1d[3][2][0][0] = 1.274041665251e-02;\n '\
            '   fineGridProjector1d[3][2][0][1] = 5.441604771565e-02;\n '\
            '   fineGridProjector1d[3][2][0][2] = 3.709095208782e-02;\n '\
            '   fineGridProjector1d[3][2][0][3] = -6.817926981005e-02;\n '\
            '   fineGridProjector1d[3][2][1][0] = -5.564103858755e-02;\n '\
            '   fineGridProjector1d[3][2][1][1] = -2.182390043732e-01;\n '\
            '   fineGridProjector1d[3][2][1][2] = -1.376638598071e-01;\n '\
            '   fineGridProjector1d[3][2][1][3] = 2.422582771987e-01;\n '\
            '   fineGridProjector1d[3][2][2][0] = 1.010071394721e+00;\n '\
            '   fineGridProjector1d[3][2][2][1] = 9.137546431442e-01;\n '\
            '   fineGridProjector1d[3][2][2][2] = 3.503991256637e-01;\n '\
            '   fineGridProjector1d[3][2][2][3] = -5.106599509271e-01;\n '\
            '   fineGridProjector1d[3][2][3][0] = 3.282922721363e-02;\n '\
            '   fineGridProjector1d[3][2][3][1] = 2.500683135133e-01;\n '\
            '   fineGridProjector1d[3][2][3][2] = 7.501737820555e-01;\n '\
            '   fineGridProjector1d[3][2][3][3] = 1.336580943538e+00;\n '\
            '   \n '\
			'   equidistantGridProjector1d[4][0][0] = 1.551408049094e+00;\n'\
			'   equidistantGridProjector1d[4][0][1] = -3.219524114026e-02;\n'\
			'   equidistantGridProjector1d[4][0][2] = 0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[4][0][3] = 9.299705243306e-03;\n'\
			'   equidistantGridProjector1d[4][0][4] = 7.635866179581e-02;\n'\
			'   equidistantGridProjector1d[4][1][0] = -8.931583920001e-01;\n'\
			'   equidistantGridProjector1d[4][1][1] = 9.627256459213e-01;\n'\
			'   equidistantGridProjector1d[4][1][2] = -0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[4][1][3] = -3.566344335768e-02;\n'\
			'   equidistantGridProjector1d[4][1][4] = -2.679416522234e-01;\n'\
			'   equidistantGridProjector1d[4][2][0] = 5.333333333333e-01;\n'\
			'   equidistantGridProjector1d[4][2][1] = 9.583333333333e-02;\n'\
			'   equidistantGridProjector1d[4][2][2] = 1.000000000000e+00;\n'\
			'   equidistantGridProjector1d[4][2][3] = 9.583333333333e-02;\n'\
			'   equidistantGridProjector1d[4][2][4] = 5.333333333333e-01;\n'\
			'   equidistantGridProjector1d[4][3][0] = -2.679416522234e-01;\n'\
			'   equidistantGridProjector1d[4][3][1] = -3.566344335768e-02;\n'\
			'   equidistantGridProjector1d[4][3][2] = 0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[4][3][3] = 9.627256459213e-01;\n'\
			'   equidistantGridProjector1d[4][3][4] = -8.931583920001e-01;\n'\
			'   equidistantGridProjector1d[4][4][0] = 7.635866179581e-02;\n'\
			'   equidistantGridProjector1d[4][4][1] = 9.299705243306e-03;\n'\
			'   equidistantGridProjector1d[4][4][2] = -0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[4][4][3] = -3.219524114025e-02;\n'\
			'   equidistantGridProjector1d[4][4][4] = 1.551408049094e+00;\n'\
			'   fineGridProjector1d[4][0][0][0] = 1.350055263641e+00;\n '\
            '   fineGridProjector1d[4][0][0][1] = 7.240733933961e-01;\n '\
            '   fineGridProjector1d[4][0][0][2] = 1.856876204782e-01;\n '\
            '   fineGridProjector1d[4][0][0][3] = -4.093289536960e-02;\n '\
            '   fineGridProjector1d[4][0][0][4] = -8.338741092706e-02;\n '\
            '   fineGridProjector1d[4][0][1][0] = -5.558211423732e-01;\n '\
            '   fineGridProjector1d[4][0][1][1] = 4.000375957323e-01;\n '\
            '   fineGridProjector1d[4][0][1][2] = 9.825172467869e-01;\n '\
            '   fineGridProjector1d[4][0][1][3] = 9.469854295406e-01;\n '\
            '   fineGridProjector1d[4][0][1][4] = 7.356281382303e-01;\n '\
            '   fineGridProjector1d[4][0][2][0] = 3.193976695423e-01;\n '\
            '   fineGridProjector1d[4][0][2][1] = -1.882041223073e-01;\n '\
            '   fineGridProjector1d[4][0][2][2] = -2.444444444444e-01;\n '\
            '   fineGridProjector1d[4][0][2][3] = 1.289969394682e-01;\n '\
            '   fineGridProjector1d[4][0][2][4] = 4.538470075812e-01;\n '\
            '   fineGridProjector1d[4][0][3][0] = -1.586695550735e-01;\n '\
            '   fineGridProjector1d[4][0][3][1] = 8.889508132080e-02;\n '\
            '   fineGridProjector1d[4][0][3][2] = 1.045161012876e-01;\n '\
            '   fineGridProjector1d[4][0][3][3] = -4.735859976571e-02;\n '\
            '   fineGridProjector1d[4][0][3][4] = -1.416250802149e-01;\n '\
            '   fineGridProjector1d[4][0][4][0] = 4.503776426321e-02;\n '\
            '   fineGridProjector1d[4][0][4][1] = -2.480194814193e-02;\n '\
            '   fineGridProjector1d[4][0][4][2] = -2.827652410819e-02;\n '\
            '   fineGridProjector1d[4][0][4][3] = 1.230912612653e-02;\n '\
            '   fineGridProjector1d[4][0][4][4] = 3.553734533050e-02;\n '\
            '   fineGridProjector1d[4][1][0][0] = -8.312593247755e-02;\n '\
            '   fineGridProjector1d[4][1][0][1] = -5.756778485668e-02;\n '\
            '   fineGridProjector1d[4][1][0][2] = 0.000000000000e+00;\n '\
            '   fineGridProjector1d[4][1][0][3] = 3.853284399929e-02;\n '\
            '   fineGridProjector1d[4][1][0][4] = 4.156296623877e-02;\n '\
            '   fineGridProjector1d[4][1][1][0] = 6.015917698357e-01;\n '\
            '   fineGridProjector1d[4][1][1][1] = 3.300395127245e-01;\n '\
            '   fineGridProjector1d[4][1][1][2] = -0.000000000000e+00;\n '\
            '   fineGridProjector1d[4][1][1][3] = -1.650197563622e-01;\n '\
            '   fineGridProjector1d[4][1][1][4] = -1.692052057501e-01;\n '\
            '   fineGridProjector1d[4][1][2][0] = 6.091764021532e-01;\n '\
            '   fineGridProjector1d[4][1][2][1] = 8.540151844952e-01;\n '\
            '   fineGridProjector1d[4][1][2][2] = 1.000000000000e+00;\n '\
            '   fineGridProjector1d[4][1][2][3] = 8.540151844952e-01;\n '\
            '   fineGridProjector1d[4][1][2][4] = 6.091764021532e-01;\n '\
            '   fineGridProjector1d[4][1][3][0] = -1.692052057501e-01;\n '\
            '   fineGridProjector1d[4][1][3][1] = -1.650197563622e-01;\n '\
            '   fineGridProjector1d[4][1][3][2] = 0.000000000000e+00;\n '\
            '   fineGridProjector1d[4][1][3][3] = 3.300395127245e-01;\n '\
            '   fineGridProjector1d[4][1][3][4] = 6.015917698357e-01;\n '\
            '   fineGridProjector1d[4][1][4][0] = 4.156296623877e-02;\n '\
            '   fineGridProjector1d[4][1][4][1] = 3.853284399929e-02;\n '\
            '   fineGridProjector1d[4][1][4][2] = -0.000000000000e+00;\n '\
            '   fineGridProjector1d[4][1][4][3] = -5.756778485668e-02;\n '\
            '   fineGridProjector1d[4][1][4][4] = -8.312593247755e-02;\n '\
            '   fineGridProjector1d[4][2][0][0] = 3.553734533050e-02;\n '\
            '   fineGridProjector1d[4][2][0][1] = 1.230912612653e-02;\n '\
            '   fineGridProjector1d[4][2][0][2] = -2.827652410819e-02;\n '\
            '   fineGridProjector1d[4][2][0][3] = -2.480194814193e-02;\n '\
            '   fineGridProjector1d[4][2][0][4] = 4.503776426321e-02;\n '\
            '   fineGridProjector1d[4][2][1][0] = -1.416250802149e-01;\n '\
            '   fineGridProjector1d[4][2][1][1] = -4.735859976571e-02;\n '\
            '   fineGridProjector1d[4][2][1][2] = 1.045161012876e-01;\n '\
            '   fineGridProjector1d[4][2][1][3] = 8.889508132080e-02;\n '\
            '   fineGridProjector1d[4][2][1][4] = -1.586695550735e-01;\n '\
            '   fineGridProjector1d[4][2][2][0] = 4.538470075812e-01;\n '\
            '   fineGridProjector1d[4][2][2][1] = 1.289969394682e-01;\n '\
            '   fineGridProjector1d[4][2][2][2] = -2.444444444444e-01;\n '\
            '   fineGridProjector1d[4][2][2][3] = -1.882041223073e-01;\n '\
            '   fineGridProjector1d[4][2][2][4] = 3.193976695423e-01;\n '\
            '   fineGridProjector1d[4][2][3][0] = 7.356281382303e-01;\n '\
            '   fineGridProjector1d[4][2][3][1] = 9.469854295406e-01;\n '\
            '   fineGridProjector1d[4][2][3][2] = 9.825172467869e-01;\n '\
            '   fineGridProjector1d[4][2][3][3] = 4.000375957323e-01;\n '\
            '   fineGridProjector1d[4][2][3][4] = -5.558211423732e-01;\n '\
            '   fineGridProjector1d[4][2][4][0] = -8.338741092706e-02;\n '\
            '   fineGridProjector1d[4][2][4][1] = -4.093289536960e-02;\n '\
            '   fineGridProjector1d[4][2][4][2] = 1.856876204782e-01;\n '\
            '   fineGridProjector1d[4][2][4][3] = 7.240733933961e-01;\n '\
            '   fineGridProjector1d[4][2][4][4] = 1.350055263641e+00;\n '\
            '   \n '\
			'   equidistantGridProjector1d[5][0][0] = 1.565673200151e+00;\n'\
			'   equidistantGridProjector1d[5][0][1] = -5.472930574888e-02;\n'\
			'   equidistantGridProjector1d[5][0][2] = 1.163100230957e-02;\n'\
			'   equidistantGridProjector1d[5][0][3] = 7.522811435131e-03;\n'\
			'   equidistantGridProjector1d[5][0][4] = -1.187353192110e-02;\n'\
			'   equidistantGridProjector1d[5][0][5] = -5.471272432927e-02;\n'\
			'   equidistantGridProjector1d[5][1][0] = -9.404628431763e-01;\n'\
			'   equidistantGridProjector1d[5][1][1] = 8.958303271943e-01;\n'\
			'   equidistantGridProjector1d[5][1][2] = -5.566479034374e-02;\n'\
			'   equidistantGridProjector1d[5][1][3] = -2.981054805679e-02;\n'\
			'   equidistantGridProjector1d[5][1][4] = 4.347670203999e-02;\n'\
			'   equidistantGridProjector1d[5][1][5] = 1.918000140387e-01;\n'\
			'   equidistantGridProjector1d[5][2][0] = 6.169300554305e-01;\n'\
			'   equidistantGridProjector1d[5][2][1] = 2.236884405952e-01;\n'\
			'   equidistantGridProjector1d[5][2][2] = 9.800324251472e-01;\n'\
			'   equidistantGridProjector1d[5][2][3] = 8.628909950863e-02;\n'\
			'   equidistantGridProjector1d[5][2][4] = -9.639263215958e-02;\n'\
			'   equidistantGridProjector1d[5][2][5] = -3.792277021146e-01;\n'\
			'   equidistantGridProjector1d[5][3][0] = -3.792277021146e-01;\n'\
			'   equidistantGridProjector1d[5][3][1] = -9.639263215958e-02;\n'\
			'   equidistantGridProjector1d[5][3][2] = 8.628909950863e-02;\n'\
			'   equidistantGridProjector1d[5][3][3] = 9.800324251472e-01;\n'\
			'   equidistantGridProjector1d[5][3][4] = 2.236884405952e-01;\n'\
			'   equidistantGridProjector1d[5][3][5] = 6.169300554305e-01;\n'\
			'   equidistantGridProjector1d[5][4][0] = 1.918000140387e-01;\n'\
			'   equidistantGridProjector1d[5][4][1] = 4.347670203999e-02;\n'\
			'   equidistantGridProjector1d[5][4][2] = -2.981054805679e-02;\n'\
			'   equidistantGridProjector1d[5][4][3] = -5.566479034374e-02;\n'\
			'   equidistantGridProjector1d[5][4][4] = 8.958303271943e-01;\n'\
			'   equidistantGridProjector1d[5][4][5] = -9.404628431763e-01;\n'\
			'   equidistantGridProjector1d[5][5][0] = -5.471272432927e-02;\n'\
			'   equidistantGridProjector1d[5][5][1] = -1.187353192110e-02;\n'\
			'   equidistantGridProjector1d[5][5][2] = 7.522811435131e-03;\n'\
			'   equidistantGridProjector1d[5][5][3] = 1.163100230957e-02;\n'\
			'   equidistantGridProjector1d[5][5][4] = -5.472930574888e-02;\n'\
			'   equidistantGridProjector1d[5][5][5] = 1.565673200151e+00;\n'\
            '   fineGridProjector1d[5][0][0][0] = 1.357780727356e+00;\n '\
            '   fineGridProjector1d[5][0][0][1] = 7.089884072874e-01;\n '\
            '   fineGridProjector1d[5][0][0][2] = 1.532342677131e-01;\n '\
            '   fineGridProjector1d[5][0][0][3] = -6.173439294873e-02;\n '\
            '   fineGridProjector1d[5][0][0][4] = -7.124783768985e-02;\n '\
            '   fineGridProjector1d[5][0][0][5] = -4.256187183674e-02;\n '\
            '   fineGridProjector1d[5][0][1][0] = -5.824222045238e-01;\n '\
            '   fineGridProjector1d[5][0][1][1] = 4.294614688616e-01;\n '\
            '   fineGridProjector1d[5][0][1][2] = 1.011931568114e+00;\n '\
            '   fineGridProjector1d[5][0][1][3] = 8.672294110247e-01;\n '\
            '   fineGridProjector1d[5][0][1][4] = 4.856629111778e-01;\n '\
            '   fineGridProjector1d[5][0][1][5] = 2.421957112211e-01;\n '\
            '   fineGridProjector1d[5][0][2][0] = 3.675417711406e-01;\n '\
            '   fineGridProjector1d[5][0][2][1] = -2.205221708713e-01;\n '\
            '   fineGridProjector1d[5][0][2][2] = -2.498095078858e-01;\n '\
            '   fineGridProjector1d[5][0][2][3] = 2.717701346628e-01;\n '\
            '   fineGridProjector1d[5][0][2][4] = 7.411544572976e-01;\n '\
            '   fineGridProjector1d[5][0][2][5] = 9.301096244275e-01;\n '\
            '   fineGridProjector1d[5][0][3][0] = -2.233071399694e-01;\n '\
            '   fineGridProjector1d[5][0][3][1] = 1.270313012119e-01;\n '\
            '   fineGridProjector1d[5][0][3][2] = 1.287538761665e-01;\n '\
            '   fineGridProjector1d[5][0][3][3] = -1.147011120037e-01;\n '\
            '   fineGridProjector1d[5][0][3][4] = -2.247049479625e-01;\n '\
            '   fineGridProjector1d[5][0][3][5] = -1.834117794241e-01;\n '\
            '   fineGridProjector1d[5][0][4][0] = 1.124115732232e-01;\n '\
            '   fineGridProjector1d[5][0][4][1] = -6.264912952844e-02;\n '\
            '   fineGridProjector1d[5][0][4][2] = -6.111282650553e-02;\n '\
            '   fineGridProjector1d[5][0][4][3] = 5.146568029542e-02;\n '\
            '   fineGridProjector1d[5][0][4][4] = 9.426074831775e-02;\n '\
            '   fineGridProjector1d[5][0][4][5] = 7.271825342658e-02;\n '\
            '   fineGridProjector1d[5][0][5][0] = -3.200472722623e-02;\n '\
            '   fineGridProjector1d[5][0][5][1] = 1.769012303893e-02;\n '\
            '   fineGridProjector1d[5][0][5][2] = 1.700262239728e-02;\n '\
            '   fineGridProjector1d[5][0][5][3] = -1.402972103050e-02;\n '\
            '   fineGridProjector1d[5][0][5][4] = -2.512533114078e-02;\n '\
            '   fineGridProjector1d[5][0][5][5] = -1.904993781431e-02;\n '\
            '   fineGridProjector1d[5][1][0][0] = -2.564351162208e-02;\n '\
            '   fineGridProjector1d[5][1][0][1] = 5.718577381192e-03;\n '\
            '   fineGridProjector1d[5][1][0][2] = 3.368861965815e-02;\n '\
            '   fineGridProjector1d[5][1][0][3] = 2.839304812755e-02;\n '\
            '   fineGridProjector1d[5][1][0][4] = 3.532052521137e-03;\n '\
            '   fineGridProjector1d[5][1][0][5] = -1.282175581104e-02;\n '\
            '   fineGridProjector1d[5][1][1][0] = 1.371026621756e-01;\n '\
            '   fineGridProjector1d[5][1][1][1] = -2.783766183457e-02;\n '\
            '   fineGridProjector1d[5][1][1][2] = -1.488648304361e-01;\n '\
            '   fineGridProjector1d[5][1][1][3] = -1.168953876580e-01;\n '\
            '   fineGridProjector1d[5][1][1][4] = -1.391883091728e-02;\n '\
            '   fineGridProjector1d[5][1][1][5] = 4.942106361179e-02;\n '\
            '   fineGridProjector1d[5][1][2][0] = 9.808364045743e-01;\n '\
            '   fineGridProjector1d[5][1][2][1] = 9.930953983583e-01;\n '\
            '   fineGridProjector1d[5][1][2][2] = 8.024523668723e-01;\n '\
            '   fineGridProjector1d[5][1][2][3] = 4.012261834361e-01;\n '\
            '   fineGridProjector1d[5][1][2][4] = 3.941046449119e-02;\n '\
            '   fineGridProjector1d[5][1][2][5] = -1.288948629286e-01;\n '\
            '   fineGridProjector1d[5][1][3][0] = -1.288948629286e-01;\n '\
            '   fineGridProjector1d[5][1][3][1] = 3.941046449119e-02;\n '\
            '   fineGridProjector1d[5][1][3][2] = 4.012261834361e-01;\n '\
            '   fineGridProjector1d[5][1][3][3] = 8.024523668723e-01;\n '\
            '   fineGridProjector1d[5][1][3][4] = 9.930953983583e-01;\n '\
            '   fineGridProjector1d[5][1][3][5] = 9.808364045743e-01;\n '\
            '   fineGridProjector1d[5][1][4][0] = 4.942106361179e-02;\n '\
            '   fineGridProjector1d[5][1][4][1] = -1.391883091728e-02;\n '\
            '   fineGridProjector1d[5][1][4][2] = -1.168953876580e-01;\n '\
            '   fineGridProjector1d[5][1][4][3] = -1.488648304361e-01;\n '\
            '   fineGridProjector1d[5][1][4][4] = -2.783766183457e-02;\n '\
            '   fineGridProjector1d[5][1][4][5] = 1.371026621756e-01;\n '\
            '   fineGridProjector1d[5][1][5][0] = -1.282175581104e-02;\n '\
            '   fineGridProjector1d[5][1][5][1] = 3.532052521137e-03;\n '\
            '   fineGridProjector1d[5][1][5][2] = 2.839304812755e-02;\n '\
            '   fineGridProjector1d[5][1][5][3] = 3.368861965815e-02;\n '\
            '   fineGridProjector1d[5][1][5][4] = 5.718577381192e-03;\n '\
            '   fineGridProjector1d[5][1][5][5] = -2.564351162208e-02;\n '\
            '   fineGridProjector1d[5][2][0][0] = -1.904993781431e-02;\n '\
            '   fineGridProjector1d[5][2][0][1] = -2.512533114078e-02;\n '\
            '   fineGridProjector1d[5][2][0][2] = -1.402972103050e-02;\n '\
            '   fineGridProjector1d[5][2][0][3] = 1.700262239728e-02;\n '\
            '   fineGridProjector1d[5][2][0][4] = 1.769012303893e-02;\n '\
            '   fineGridProjector1d[5][2][0][5] = -3.200472722623e-02;\n '\
            '   fineGridProjector1d[5][2][1][0] = 7.271825342658e-02;\n '\
            '   fineGridProjector1d[5][2][1][1] = 9.426074831775e-02;\n '\
            '   fineGridProjector1d[5][2][1][2] = 5.146568029542e-02;\n '\
            '   fineGridProjector1d[5][2][1][3] = -6.111282650553e-02;\n '\
            '   fineGridProjector1d[5][2][1][4] = -6.264912952844e-02;\n '\
            '   fineGridProjector1d[5][2][1][5] = 1.124115732232e-01;\n '\
            '   fineGridProjector1d[5][2][2][0] = -1.834117794241e-01;\n '\
            '   fineGridProjector1d[5][2][2][1] = -2.247049479625e-01;\n '\
            '   fineGridProjector1d[5][2][2][2] = -1.147011120037e-01;\n '\
            '   fineGridProjector1d[5][2][2][3] = 1.287538761665e-01;\n '\
            '   fineGridProjector1d[5][2][2][4] = 1.270313012119e-01;\n '\
            '   fineGridProjector1d[5][2][2][5] = -2.233071399694e-01;\n '\
            '   fineGridProjector1d[5][2][3][0] = 9.301096244275e-01;\n '\
            '   fineGridProjector1d[5][2][3][1] = 7.411544572976e-01;\n '\
            '   fineGridProjector1d[5][2][3][2] = 2.717701346628e-01;\n '\
            '   fineGridProjector1d[5][2][3][3] = -2.498095078858e-01;\n '\
            '   fineGridProjector1d[5][2][3][4] = -2.205221708713e-01;\n '\
            '   fineGridProjector1d[5][2][3][5] = 3.675417711406e-01;\n '\
            '   fineGridProjector1d[5][2][4][0] = 2.421957112211e-01;\n '\
            '   fineGridProjector1d[5][2][4][1] = 4.856629111778e-01;\n '\
            '   fineGridProjector1d[5][2][4][2] = 8.672294110247e-01;\n '\
            '   fineGridProjector1d[5][2][4][3] = 1.011931568114e+00;\n '\
            '   fineGridProjector1d[5][2][4][4] = 4.294614688616e-01;\n '\
            '   fineGridProjector1d[5][2][4][5] = -5.824222045238e-01;\n '\
            '   fineGridProjector1d[5][2][5][0] = -4.256187183674e-02;\n '\
            '   fineGridProjector1d[5][2][5][1] = -7.124783768985e-02;\n '\
            '   fineGridProjector1d[5][2][5][2] = -6.173439294873e-02;\n '\
            '   fineGridProjector1d[5][2][5][3] = 1.532342677131e-01;\n '\
            '   fineGridProjector1d[5][2][5][4] = 7.089884072874e-01;\n '\
            '   fineGridProjector1d[5][2][5][5] = 1.357780727356e+00;\n '\
            '   \n '\
			'   equidistantGridProjector1d[6][0][0] = 1.574662499711e+00;\n'\
			'   equidistantGridProjector1d[6][0][1] = -6.825740746497e-02;\n'\
			'   equidistantGridProjector1d[6][0][2] = 2.160099903605e-02;\n'\
			'   equidistantGridProjector1d[6][0][3] = -0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[6][0][4] = -1.037189511212e-02;\n'\
			'   equidistantGridProjector1d[6][0][5] = 1.193155743473e-02;\n'\
			'   equidistantGridProjector1d[6][0][6] = 4.111514886291e-02;\n'\
			'   equidistantGridProjector1d[6][1][0] = -9.707266965061e-01;\n'\
			'   equidistantGridProjector1d[6][1][1] = 8.062490980084e-01;\n'\
			'   equidistantGridProjector1d[6][1][2] = -1.020216283465e-01;\n'\
			'   equidistantGridProjector1d[6][1][3] = 0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[6][1][4] = 3.874442670885e-02;\n'\
			'   equidistantGridProjector1d[6][1][5] = -4.286290507063e-02;\n'\
			'   equidistantGridProjector1d[6][1][6] = -1.440701036121e-01;\n'\
			'   equidistantGridProjector1d[6][2][0] = 6.721078619224e-01;\n'\
			'   equidistantGridProjector1d[6][2][1] = 3.683285328645e-01;\n'\
			'   equidistantGridProjector1d[6][2][2] = 9.140872873432e-01;\n'\
			'   equidistantGridProjector1d[6][2][3] = -0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[6][2][4] = -8.966999503802e-02;\n'\
			'   equidistantGridProjector1d[6][2][5] = 8.957291141200e-02;\n'\
			'   equidistantGridProjector1d[6][2][6] = 2.840541467652e-01;\n'\
			'   equidistantGridProjector1d[6][3][0] = -4.571428571429e-01;\n'\
			'   equidistantGridProjector1d[6][3][1] = -1.649617871840e-01;\n'\
			'   equidistantGridProjector1d[6][3][2] = 2.276308054086e-01;\n'\
			'   equidistantGridProjector1d[6][3][3] = 1.000000000000e+00;\n'\
			'   equidistantGridProjector1d[6][3][4] = 2.276308054086e-01;\n'\
			'   equidistantGridProjector1d[6][3][5] = -1.649617871840e-01;\n'\
			'   equidistantGridProjector1d[6][3][6] = -4.571428571429e-01;\n'\
			'   equidistantGridProjector1d[6][4][0] = 2.840541467652e-01;\n'\
			'   equidistantGridProjector1d[6][4][1] = 8.957291141200e-02;\n'\
			'   equidistantGridProjector1d[6][4][2] = -8.966999503802e-02;\n'\
			'   equidistantGridProjector1d[6][4][3] = 0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[6][4][4] = 9.140872873432e-01;\n'\
			'   equidistantGridProjector1d[6][4][5] = 3.683285328645e-01;\n'\
			'   equidistantGridProjector1d[6][4][6] = 6.721078619224e-01;\n'\
			'   equidistantGridProjector1d[6][5][0] = -1.440701036121e-01;\n'\
			'   equidistantGridProjector1d[6][5][1] = -4.286290507063e-02;\n'\
			'   equidistantGridProjector1d[6][5][2] = 3.874442670885e-02;\n'\
			'   equidistantGridProjector1d[6][5][3] = -0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[6][5][4] = -1.020216283465e-01;\n'\
			'   equidistantGridProjector1d[6][5][5] = 8.062490980084e-01;\n'\
			'   equidistantGridProjector1d[6][5][6] = -9.707266965061e-01;\n'\
			'   equidistantGridProjector1d[6][6][0] = 4.111514886291e-02;\n'\
			'   equidistantGridProjector1d[6][6][1] = 1.193155743473e-02;\n'\
			'   equidistantGridProjector1d[6][6][2] = -1.037189511212e-02;\n'\
			'   equidistantGridProjector1d[6][6][3] = 0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[6][6][4] = 2.160099903605e-02;\n'\
			'   equidistantGridProjector1d[6][6][5] = -6.825740746497e-02;\n'\
			'   equidistantGridProjector1d[6][6][6] = 1.574662499711e+00;\n'\
            '   fineGridProjector1d[6][0][0][0] = 1.362618484193e+00;\n '\
            '   fineGridProjector1d[6][0][0][1] = 6.995131651884e-01;\n '\
            '   fineGridProjector1d[6][0][0][2] = 1.346121125604e-01;\n '\
            '   fineGridProjector1d[6][0][0][3] = -6.825740746497e-02;\n '\
            '   fineGridProjector1d[6][0][0][4] = -5.320343713617e-02;\n '\
            '   fineGridProjector1d[6][0][0][5] = -5.194099155088e-03;\n '\
            '   fineGridProjector1d[6][0][0][6] = 1.754760402688e-02;\n '\
            '   fineGridProjector1d[6][0][1][0] = -5.993422803481e-01;\n '\
            '   fineGridProjector1d[6][0][1][1] = 4.482086256914e-01;\n '\
            '   fineGridProjector1d[6][0][1][2] = 1.026548137308e+00;\n '\
            '   fineGridProjector1d[6][0][1][3] = 8.062490980084e-01;\n '\
            '   fineGridProjector1d[6][0][1][4] = 3.311111433515e-01;\n '\
            '   fineGridProjector1d[6][0][1][5] = 2.674414068454e-02;\n '\
            '   fineGridProjector1d[6][0][1][6] = -8.408881805116e-02;\n '\
            '   fineGridProjector1d[6][0][2][0] = 3.991303365277e-01;\n '\
            '   fineGridProjector1d[6][0][2][1] = -2.419736651777e-01;\n '\
            '   fineGridProjector1d[6][0][2][2] = -2.492097393249e-01;\n '\
            '   fineGridProjector1d[6][0][2][3] = 3.683285328645e-01;\n '\
            '   fineGridProjector1d[6][0][2][4] = 8.821602474804e-01;\n '\
            '   fineGridProjector1d[6][0][2][5] = 1.004657317050e+00;\n '\
            '   fineGridProjector1d[6][0][2][6] = 9.426302282836e-01;\n '\
            '   fineGridProjector1d[6][0][3][0] = -2.682736762509e-01;\n '\
            '   fineGridProjector1d[6][0][3][1] = 1.539827176121e-01;\n '\
            '   fineGridProjector1d[6][0][3][2] = 1.409095078758e-01;\n '\
            '   fineGridProjector1d[6][0][3][3] = -1.649617871840e-01;\n '\
            '   fineGridProjector1d[6][0][3][4] = -2.385794191977e-01;\n '\
            '   fineGridProjector1d[6][0][3][5] = -3.740814808829e-02;\n '\
            '   fineGridProjector1d[6][0][3][6] = 1.711136235543e-01;\n '\
            '   fineGridProjector1d[6][0][4][0] = 1.658704710675e-01;\n '\
            '   fineGridProjector1d[6][0][4][1] = -9.314489156547e-02;\n '\
            '   fineGridProjector1d[6][0][4][2] = -8.172984903586e-02;\n '\
            '   fineGridProjector1d[6][0][4][3] = 8.957291141200e-02;\n '\
            '   fineGridProjector1d[6][0][4][4] = 1.181633190270e-01;\n '\
            '   fineGridProjector1d[6][0][4][5] = 1.660901817810e-02;\n '\
            '   fineGridProjector1d[6][0][4][6] = -6.924756100357e-02;\n '\
            '   fineGridProjector1d[6][0][5][0] = -8.393064203396e-02;\n '\
            '   fineGridProjector1d[6][0][5][1] = 4.665527271095e-02;\n '\
            '   fineGridProjector1d[6][0][5][2] = 4.018269078282e-02;\n '\
            '   fineGridProjector1d[6][0][5][3] = -4.286290507063e-02;\n '\
            '   fineGridProjector1d[6][0][5][4] = -5.466327165748e-02;\n '\
            '   fineGridProjector1d[6][0][5][5] = -7.418235241249e-03;\n '\
            '   fineGridProjector1d[6][0][5][6] = 3.013146140696e-02;\n '\
            '   fineGridProjector1d[6][0][6][0] = 2.392730684521e-02;\n '\
            '   fineGridProjector1d[6][0][6][1] = -1.324122445973e-02;\n '\
            '   fineGridProjector1d[6][0][6][2] = -1.131286016633e-02;\n '\
            '   fineGridProjector1d[6][0][6][3] = 1.193155743473e-02;\n '\
            '   fineGridProjector1d[6][0][6][4] = 1.501141813240e-02;\n '\
            '   fineGridProjector1d[6][0][6][5] = 2.010006571887e-03;\n '\
            '   fineGridProjector1d[6][0][6][6] = -8.086538217023e-03;\n '\
            '   fineGridProjector1d[6][1][0][0] = 2.500315801791e-02;\n '\
            '   fineGridProjector1d[6][1][0][1] = 3.217671765537e-02;\n '\
            '   fineGridProjector1d[6][1][0][2] = 2.453160864937e-02;\n '\
            '   fineGridProjector1d[6][1][0][3] = -0.000000000000e+00;\n '\
            '   fineGridProjector1d[6][1][0][4] = -1.841078629593e-02;\n '\
            '   fineGridProjector1d[6][1][0][5] = -1.887997532051e-02;\n '\
            '   fineGridProjector1d[6][1][0][6] = -1.250157900895e-02;\n '\
            '   fineGridProjector1d[6][1][1][0] = -1.165016974039e-01;\n '\
            '   fineGridProjector1d[6][1][1][1] = -1.430425625246e-01;\n '\
            '   fineGridProjector1d[6][1][1][2] = -1.031035679286e-01;\n '\
            '   fineGridProjector1d[6][1][1][3] = 0.000000000000e+00;\n '\
            '   fineGridProjector1d[6][1][1][4] = 7.128827947233e-02;\n '\
            '   fineGridProjector1d[6][1][1][5] = 7.152128126230e-02;\n '\
            '   fineGridProjector1d[6][1][1][6] = 4.682111535637e-02;\n '\
            '   fineGridProjector1d[6][1][2][0] = 8.810777835154e-01;\n '\
            '   fineGridProjector1d[6][1][2][1] = 7.093283979959e-01;\n '\
            '   fineGridProjector1d[6][1][2][2] = 3.676964689786e-01;\n '\
            '   fineGridProjector1d[6][1][2][3] = -0.000000000000e+00;\n '\
            '   fineGridProjector1d[6][1][2][4] = -1.838482344893e-01;\n '\
            '   fineGridProjector1d[6][1][2][5] = -1.723490862383e-01;\n '\
            '   fineGridProjector1d[6][1][2][6] = -1.091575785908e-01;\n '\
            '   fineGridProjector1d[6][1][3][0] = 2.852587981140e-01;\n '\
            '   fineGridProjector1d[6][1][3][1] = 5.212452271698e-01;\n '\
            '   fineGridProjector1d[6][1][3][2] = 8.418462316135e-01;\n '\
            '   fineGridProjector1d[6][1][3][3] = 1.000000000000e+00;\n '\
            '   fineGridProjector1d[6][1][3][4] = 8.418462316135e-01;\n '\
            '   fineGridProjector1d[6][1][3][5] = 5.212452271698e-01;\n '\
            '   fineGridProjector1d[6][1][3][6] = 2.852587981140e-01;\n '\
            '   fineGridProjector1d[6][1][4][0] = -1.091575785908e-01;\n '\
            '   fineGridProjector1d[6][1][4][1] = -1.723490862383e-01;\n '\
            '   fineGridProjector1d[6][1][4][2] = -1.838482344893e-01;\n '\
            '   fineGridProjector1d[6][1][4][3] = 0.000000000000e+00;\n '\
            '   fineGridProjector1d[6][1][4][4] = 3.676964689786e-01;\n '\
            '   fineGridProjector1d[6][1][4][5] = 7.093283979959e-01;\n '\
            '   fineGridProjector1d[6][1][4][6] = 8.810777835154e-01;\n '\
            '   fineGridProjector1d[6][1][5][0] = 4.682111535637e-02;\n '\
            '   fineGridProjector1d[6][1][5][1] = 7.152128126230e-02;\n '\
            '   fineGridProjector1d[6][1][5][2] = 7.128827947233e-02;\n '\
            '   fineGridProjector1d[6][1][5][3] = -0.000000000000e+00;\n '\
            '   fineGridProjector1d[6][1][5][4] = -1.031035679286e-01;\n '\
            '   fineGridProjector1d[6][1][5][5] = -1.430425625246e-01;\n '\
            '   fineGridProjector1d[6][1][5][6] = -1.165016974039e-01;\n '\
            '   fineGridProjector1d[6][1][6][0] = -1.250157900895e-02;\n '\
            '   fineGridProjector1d[6][1][6][1] = -1.887997532051e-02;\n '\
            '   fineGridProjector1d[6][1][6][2] = -1.841078629593e-02;\n '\
            '   fineGridProjector1d[6][1][6][3] = 0.000000000000e+00;\n '\
            '   fineGridProjector1d[6][1][6][4] = 2.453160864937e-02;\n '\
            '   fineGridProjector1d[6][1][6][5] = 3.217671765537e-02;\n '\
            '   fineGridProjector1d[6][1][6][6] = 2.500315801791e-02;\n '\
            '   fineGridProjector1d[6][2][0][0] = -8.086538217023e-03;\n '\
            '   fineGridProjector1d[6][2][0][1] = 2.010006571887e-03;\n '\
            '   fineGridProjector1d[6][2][0][2] = 1.501141813240e-02;\n '\
            '   fineGridProjector1d[6][2][0][3] = 1.193155743473e-02;\n '\
            '   fineGridProjector1d[6][2][0][4] = -1.131286016633e-02;\n '\
            '   fineGridProjector1d[6][2][0][5] = -1.324122445973e-02;\n '\
            '   fineGridProjector1d[6][2][0][6] = 2.392730684521e-02;\n '\
            '   fineGridProjector1d[6][2][1][0] = 3.013146140696e-02;\n '\
            '   fineGridProjector1d[6][2][1][1] = -7.418235241249e-03;\n '\
            '   fineGridProjector1d[6][2][1][2] = -5.466327165748e-02;\n '\
            '   fineGridProjector1d[6][2][1][3] = -4.286290507063e-02;\n '\
            '   fineGridProjector1d[6][2][1][4] = 4.018269078282e-02;\n '\
            '   fineGridProjector1d[6][2][1][5] = 4.665527271095e-02;\n '\
            '   fineGridProjector1d[6][2][1][6] = -8.393064203396e-02;\n '\
            '   fineGridProjector1d[6][2][2][0] = -6.924756100357e-02;\n '\
            '   fineGridProjector1d[6][2][2][1] = 1.660901817810e-02;\n '\
            '   fineGridProjector1d[6][2][2][2] = 1.181633190270e-01;\n '\
            '   fineGridProjector1d[6][2][2][3] = 8.957291141200e-02;\n '\
            '   fineGridProjector1d[6][2][2][4] = -8.172984903586e-02;\n '\
            '   fineGridProjector1d[6][2][2][5] = -9.314489156547e-02;\n '\
            '   fineGridProjector1d[6][2][2][6] = 1.658704710675e-01;\n '\
            '   fineGridProjector1d[6][2][3][0] = 1.711136235543e-01;\n '\
            '   fineGridProjector1d[6][2][3][1] = -3.740814808829e-02;\n '\
            '   fineGridProjector1d[6][2][3][2] = -2.385794191977e-01;\n '\
            '   fineGridProjector1d[6][2][3][3] = -1.649617871840e-01;\n '\
            '   fineGridProjector1d[6][2][3][4] = 1.409095078758e-01;\n '\
            '   fineGridProjector1d[6][2][3][5] = 1.539827176121e-01;\n '\
            '   fineGridProjector1d[6][2][3][6] = -2.682736762509e-01;\n '\
            '   fineGridProjector1d[6][2][4][0] = 9.426302282836e-01;\n '\
            '   fineGridProjector1d[6][2][4][1] = 1.004657317050e+00;\n '\
            '   fineGridProjector1d[6][2][4][2] = 8.821602474804e-01;\n '\
            '   fineGridProjector1d[6][2][4][3] = 3.683285328645e-01;\n '\
            '   fineGridProjector1d[6][2][4][4] = -2.492097393249e-01;\n '\
            '   fineGridProjector1d[6][2][4][5] = -2.419736651777e-01;\n '\
            '   fineGridProjector1d[6][2][4][6] = 3.991303365277e-01;\n '\
            '   fineGridProjector1d[6][2][5][0] = -8.408881805116e-02;\n '\
            '   fineGridProjector1d[6][2][5][1] = 2.674414068454e-02;\n '\
            '   fineGridProjector1d[6][2][5][2] = 3.311111433515e-01;\n '\
            '   fineGridProjector1d[6][2][5][3] = 8.062490980084e-01;\n '\
            '   fineGridProjector1d[6][2][5][4] = 1.026548137308e+00;\n '\
            '   fineGridProjector1d[6][2][5][5] = 4.482086256914e-01;\n '\
            '   fineGridProjector1d[6][2][5][6] = -5.993422803481e-01;\n '\
            '   fineGridProjector1d[6][2][6][0] = 1.754760402688e-02;\n '\
            '   fineGridProjector1d[6][2][6][1] = -5.194099155088e-03;\n '\
            '   fineGridProjector1d[6][2][6][2] = -5.320343713617e-02;\n '\
            '   fineGridProjector1d[6][2][6][3] = -6.825740746497e-02;\n '\
            '   fineGridProjector1d[6][2][6][4] = 1.346121125604e-01;\n '\
            '   fineGridProjector1d[6][2][6][5] = 6.995131651884e-01;\n '\
            '   fineGridProjector1d[6][2][6][6] = 1.362618484193e+00;\n '\
            '   \n '\
			'   equidistantGridProjector1d[7][0][0] = 1.580687063031e+00;\n'\
			'   equidistantGridProjector1d[7][0][1] = -7.335564462275e-02;\n'\
			'   equidistantGridProjector1d[7][0][2] = 2.743355394046e-02;\n'\
			'   equidistantGridProjector1d[7][0][3] = -7.241518001033e-03;\n'\
			'   equidistantGridProjector1d[7][0][4] = -5.365970005741e-03;\n'\
			'   equidistantGridProjector1d[7][0][5] = 1.050279558546e-02;\n'\
			'   equidistantGridProjector1d[7][0][6] = -1.077633804450e-02;\n'\
			'   equidistantGridProjector1d[7][0][7] = -3.202042284583e-02;\n'\
			'   equidistantGridProjector1d[7][1][0] = -9.912041583117e-01;\n'\
			'   equidistantGridProjector1d[7][1][1] = 7.033555058504e-01;\n'\
			'   equidistantGridProjector1d[7][1][2] = -1.272414335879e-01;\n'\
			'   equidistantGridProjector1d[7][1][3] = 2.907068440897e-02;\n'\
			'   equidistantGridProjector1d[7][1][4] = 2.023012984288e-02;\n'\
			'   equidistantGridProjector1d[7][1][5] = -3.822681414156e-02;\n'\
			'   equidistantGridProjector1d[7][1][6] = 3.834864110228e-02;\n'\
			'   equidistantGridProjector1d[7][1][7] = 1.121772102087e-01;\n'\
			'   equidistantGridProjector1d[7][2][0] = 7.101568903172e-01;\n'\
			'   equidistantGridProjector1d[7][2][1] = 5.132097472522e-01;\n'\
			'   equidistantGridProjector1d[7][2][2] = 8.075701993243e-01;\n'\
			'   equidistantGridProjector1d[7][2][3] = -8.303567809735e-02;\n'\
			'   equidistantGridProjector1d[7][2][4] = -4.754068964647e-02;\n'\
			'   equidistantGridProjector1d[7][2][5] = 8.206947281072e-02;\n'\
			'   equidistantGridProjector1d[7][2][6] = -7.813245647685e-02;\n'\
			'   equidistantGridProjector1d[7][2][7] = -2.208713667044e-01;\n'\
			'   equidistantGridProjector1d[7][3][0] = -5.126556338014e-01;\n'\
			'   equidistantGridProjector1d[7][3][1] = -2.267109983185e-01;\n'\
			'   equidistantGridProjector1d[7][3][2] = 3.968483941550e-01;\n'\
			'   equidistantGridProjector1d[7][3][3] = 9.728943391547e-01;\n'\
			'   equidistantGridProjector1d[7][3][4] = 1.209887023440e-01;\n'\
			'   equidistantGridProjector1d[7][3][5] = -1.589561680865e-01;\n'\
			'   equidistantGridProjector1d[7][3][6] = 1.340615432577e-01;\n'\
			'   equidistantGridProjector1d[7][3][7] = 3.537304181064e-01;\n'\
			'   equidistantGridProjector1d[7][4][0] = 3.537304181064e-01;\n'\
			'   equidistantGridProjector1d[7][4][1] = 1.340615432577e-01;\n'\
			'   equidistantGridProjector1d[7][4][2] = -1.589561680865e-01;\n'\
			'   equidistantGridProjector1d[7][4][3] = 1.209887023440e-01;\n'\
			'   equidistantGridProjector1d[7][4][4] = 9.728943391547e-01;\n'\
			'   equidistantGridProjector1d[7][4][5] = 3.968483941550e-01;\n'\
			'   equidistantGridProjector1d[7][4][6] = -2.267109983185e-01;\n'\
			'   equidistantGridProjector1d[7][4][7] = -5.126556338014e-01;\n'\
			'   equidistantGridProjector1d[7][5][0] = -2.208713667044e-01;\n'\
			'   equidistantGridProjector1d[7][5][1] = -7.813245647685e-02;\n'\
			'   equidistantGridProjector1d[7][5][2] = 8.206947281072e-02;\n'\
			'   equidistantGridProjector1d[7][5][3] = -4.754068964647e-02;\n'\
			'   equidistantGridProjector1d[7][5][4] = -8.303567809735e-02;\n'\
			'   equidistantGridProjector1d[7][5][5] = 8.075701993243e-01;\n'\
			'   equidistantGridProjector1d[7][5][6] = 5.132097472522e-01;\n'\
			'   equidistantGridProjector1d[7][5][7] = 7.101568903172e-01;\n'\
			'   equidistantGridProjector1d[7][6][0] = 1.121772102087e-01;\n'\
			'   equidistantGridProjector1d[7][6][1] = 3.834864110228e-02;\n'\
			'   equidistantGridProjector1d[7][6][2] = -3.822681414156e-02;\n'\
			'   equidistantGridProjector1d[7][6][3] = 2.023012984288e-02;\n'\
			'   equidistantGridProjector1d[7][6][4] = 2.907068440897e-02;\n'\
			'   equidistantGridProjector1d[7][6][5] = -1.272414335879e-01;\n'\
			'   equidistantGridProjector1d[7][6][6] = 7.033555058504e-01;\n'\
			'   equidistantGridProjector1d[7][6][7] = -9.912041583117e-01;\n'\
			'   equidistantGridProjector1d[7][7][0] = -3.202042284583e-02;\n'\
			'   equidistantGridProjector1d[7][7][1] = -1.077633804450e-02;\n'\
			'   equidistantGridProjector1d[7][7][2] = 1.050279558546e-02;\n'\
			'   equidistantGridProjector1d[7][7][3] = -5.365970005741e-03;\n'\
			'   equidistantGridProjector1d[7][7][4] = -7.241518001033e-03;\n'\
			'   equidistantGridProjector1d[7][7][5] = 2.743355394046e-02;\n'\
			'   equidistantGridProjector1d[7][7][6] = -7.335564462275e-02;\n'\
			'   equidistantGridProjector1d[7][7][7] = 1.580687063031e+00;\n'\
            '   fineGridProjector1d[7][0][0][0] = 1.365847593693e+00;\n '\
            '   fineGridProjector1d[7][0][0][1] = 6.931802462555e-01;\n '\
            '   fineGridProjector1d[7][0][0][2] = 1.229100718913e-01;\n '\
            '   fineGridProjector1d[7][0][0][3] = -7.032301133937e-02;\n '\
            '   fineGridProjector1d[7][0][0][4] = -3.989926453250e-02;\n '\
            '   fineGridProjector1d[7][0][0][5] = 1.296433184770e-02;\n '\
            '   fineGridProjector1d[7][0][0][6] = 2.986020799602e-02;\n '\
            '   fineGridProjector1d[7][0][0][7] = 2.863609422315e-02;\n '\
            '   fineGridProjector1d[7][0][1][0] = -6.107484297948e-01;\n '\
            '   fineGridProjector1d[7][0][1][1] = 4.608504141649e-01;\n '\
            '   fineGridProjector1d[7][0][1][2] = 1.034684310493e+00;\n '\
            '   fineGridProjector1d[7][0][1][3] = 7.623766767031e-01;\n '\
            '   fineGridProjector1d[7][0][1][4] = 2.377781367632e-01;\n '\
            '   fineGridProjector1d[7][0][1][5] = -6.394563419373e-02;\n '\
            '   fineGridProjector1d[7][0][1][6] = -1.355379672259e-01;\n '\
            '   fineGridProjector1d[7][0][1][7] = -1.253727350632e-01;\n '\
            '   fineGridProjector1d[7][0][2][0] = 4.208308498076e-01;\n '\
            '   fineGridProjector1d[7][0][2][1] = -2.568048602856e-01;\n '\
            '   fineGridProjector1d[7][0][2][2] = -2.470614882914e-01;\n '\
            '   fineGridProjector1d[7][0][2][3] = 4.338527015142e-01;\n '\
            '   fineGridProjector1d[7][0][2][4] = 9.499278267926e-01;\n '\
            '   fineGridProjector1d[7][0][2][5] = 9.583437398281e-01;\n '\
            '   fineGridProjector1d[7][0][2][6] = 7.203812928376e-01;\n '\
            '   fineGridProjector1d[7][0][2][7] = 5.271518345152e-01;\n '\
            '   fineGridProjector1d[7][0][3][0] = -3.001846752984e-01;\n '\
            '   fineGridProjector1d[7][0][3][1] = 1.732860881489e-01;\n '\
            '   fineGridProjector1d[7][0][3][2] = 1.474619898397e-01;\n '\
            '   fineGridProjector1d[7][0][3][3] = -2.002855123277e-01;\n '\
            '   fineGridProjector1d[7][0][3][4] = -2.236539758146e-01;\n '\
            '   fineGridProjector1d[7][0][3][5] = 1.315770211480e-01;\n '\
            '   fineGridProjector1d[7][0][3][6] = 5.115653916979e-01;\n '\
            '   fineGridProjector1d[7][0][3][7] = 7.184640259508e-01;\n '\
            '   fineGridProjector1d[7][0][4][0] = 2.060736411242e-01;\n '\
            '   fineGridProjector1d[7][0][4][1] = -1.163032031210e-01;\n '\
            '   fineGridProjector1d[7][0][4][2] = -9.469656069316e-02;\n '\
            '   fineGridProjector1d[7][0][4][3] = 1.196502202314e-01;\n '\
            '   fineGridProjector1d[7][0][4][4] = 1.196535884089e-01;\n '\
            '   fineGridProjector1d[7][0][4][5] = -6.005551773248e-02;\n '\
            '   fineGridProjector1d[7][0][4][6] = -1.904996614855e-01;\n '\
            '   fineGridProjector1d[7][0][4][7] = -2.211431752079e-01;\n '\
            '   fineGridProjector1d[7][0][5][0] = -1.283480353197e-01;\n '\
            '   fineGridProjector1d[7][0][5][1] = 7.164436760859e-02;\n '\
            '   fineGridProjector1d[7][0][5][2] = 5.715210354087e-02;\n '\
            '   fineGridProjector1d[7][0][5][3] = -7.002016468091e-02;\n '\
            '   fineGridProjector1d[7][0][5][4] = -6.717992322398e-02;\n '\
            '   fineGridProjector1d[7][0][5][5] = 3.207906308248e-02;\n '\
            '   fineGridProjector1d[7][0][5][6] = 9.672621942444e-02;\n '\
            '   fineGridProjector1d[7][0][5][7] = 1.081757500422e-01;\n '\
            '   fineGridProjector1d[7][0][6][0] = 6.510002774521e-02;\n '\
            '   fineGridProjector1d[7][0][6][1] = -3.613355466876e-02;\n '\
            '   fineGridProjector1d[7][0][6][2] = -2.852872284388e-02;\n '\
            '   fineGridProjector1d[7][0][6][3] = 3.443368273302e-02;\n '\
            '   fineGridProjector1d[7][0][6][4] = 3.241364216707e-02;\n '\
            '   fineGridProjector1d[7][0][6][5] = -1.514937608035e-02;\n '\
            '   fineGridProjector1d[7][0][6][6] = -4.476019539319e-02;\n '\
            '   fineGridProjector1d[7][0][6][7] = -4.935970743790e-02;\n '\
            '   fineGridProjector1d[7][0][7][0] = -1.857097195694e-02;\n '\
            '   fineGridProjector1d[7][0][7][1] = 1.028050189740e-02;\n '\
            '   fineGridProjector1d[7][0][7][2] = 8.078296063724e-03;\n '\
            '   fineGridProjector1d[7][0][7][3] = -9.684592833815e-03;\n '\
            '   fineGridProjector1d[7][0][7][4] = -9.040030560746e-03;\n '\
            '   fineGridProjector1d[7][0][7][5] = 4.186372100316e-03;\n '\
            '   fineGridProjector1d[7][0][7][6] = 1.226471214876e-02;\n '\
            '   fineGridProjector1d[7][0][7][7] = 1.344791297770e-02;\n '\
            '   fineGridProjector1d[7][1][0][0] = 2.570660756470e-02;\n '\
            '   fineGridProjector1d[7][1][0][1] = 1.654900747731e-02;\n '\
            '   fineGridProjector1d[7][1][0][2] = -1.573234145260e-03;\n '\
            '   fineGridProjector1d[7][1][0][3] = -1.657195712217e-02;\n '\
            '   fineGridProjector1d[7][1][0][4] = -1.458790409097e-02;\n '\
            '   fineGridProjector1d[7][1][0][5] = -1.087803748778e-03;\n '\
            '   fineGridProjector1d[7][1][0][6] = 9.378967905442e-03;\n '\
            '   fineGridProjector1d[7][1][0][7] = 1.285330378235e-02;\n '\
            '   fineGridProjector1d[7][1][1][0] = -1.108802057964e-01;\n '\
            '   fineGridProjector1d[7][1][1][1] = -6.950726400108e-02;\n '\
            '   fineGridProjector1d[7][1][1][2] = 6.381407609524e-03;\n '\
            '   fineGridProjector1d[7][1][1][3] = 6.504783375660e-02;\n '\
            '   fineGridProjector1d[7][1][1][4] = 5.577460134906e-02;\n '\
            '   fineGridProjector1d[7][1][1][5] = 4.080876041870e-03;\n '\
            '   fineGridProjector1d[7][1][1][6] = -3.475363200054e-02;\n '\
            '   fineGridProjector1d[7][1][1][7] = -4.731725481747e-02;\n '\
            '   fineGridProjector1d[7][1][2][0] = 4.300245431307e-01;\n '\
            '   fineGridProjector1d[7][1][2][1] = 2.373938397973e-01;\n '\
            '   fineGridProjector1d[7][1][2][2] = -1.892476299597e-02;\n '\
            '   fineGridProjector1d[7][1][2][3] = -1.722409753590e-01;\n '\
            '   fineGridProjector1d[7][1][2][4] = -1.363382599180e-01;\n '\
            '   fineGridProjector1d[7][1][2][5] = -9.462381497985e-03;\n '\
            '   fineGridProjector1d[7][1][2][6] = 7.801524774262e-02;\n '\
            '   fineGridProjector1d[7][1][2][7] = 1.044694886293e-01;\n '\
            '   fineGridProjector1d[7][1][3][0] = 8.031145301117e-01;\n '\
            '   fineGridProjector1d[7][1][3][1] = 9.336984282881e-01;\n '\
            '   fineGridProjector1d[7][1][3][2] = 9.976154452652e-01;\n '\
            '   fineGridProjector1d[7][1][3][3] = 8.126111075896e-01;\n '\
            '   fineGridProjector1d[7][1][3][4] = 4.063055537948e-01;\n '\
            '   fineGridProjector1d[7][1][3][5] = 2.297045347135e-02;\n '\
            '   fineGridProjector1d[7][1][3][6] = -1.707745952091e-01;\n '\
            '   fineGridProjector1d[7][1][3][7] = -2.179710126049e-01;\n '\
            '   fineGridProjector1d[7][1][4][0] = -2.179710126049e-01;\n '\
            '   fineGridProjector1d[7][1][4][1] = -1.707745952091e-01;\n '\
            '   fineGridProjector1d[7][1][4][2] = 2.297045347135e-02;\n '\
            '   fineGridProjector1d[7][1][4][3] = 4.063055537948e-01;\n '\
            '   fineGridProjector1d[7][1][4][4] = 8.126111075896e-01;\n '\
            '   fineGridProjector1d[7][1][4][5] = 9.976154452652e-01;\n '\
            '   fineGridProjector1d[7][1][4][6] = 9.336984282881e-01;\n '\
            '   fineGridProjector1d[7][1][4][7] = 8.031145301117e-01;\n '\
            '   fineGridProjector1d[7][1][5][0] = 1.044694886293e-01;\n '\
            '   fineGridProjector1d[7][1][5][1] = 7.801524774262e-02;\n '\
            '   fineGridProjector1d[7][1][5][2] = -9.462381497984e-03;\n '\
            '   fineGridProjector1d[7][1][5][3] = -1.363382599180e-01;\n '\
            '   fineGridProjector1d[7][1][5][4] = -1.722409753590e-01;\n '\
            '   fineGridProjector1d[7][1][5][5] = -1.892476299597e-02;\n '\
            '   fineGridProjector1d[7][1][5][6] = 2.373938397973e-01;\n '\
            '   fineGridProjector1d[7][1][5][7] = 4.300245431307e-01;\n '\
            '   fineGridProjector1d[7][1][6][0] = -4.731725481747e-02;\n '\
            '   fineGridProjector1d[7][1][6][1] = -3.475363200054e-02;\n '\
            '   fineGridProjector1d[7][1][6][2] = 4.080876041870e-03;\n '\
            '   fineGridProjector1d[7][1][6][3] = 5.577460134906e-02;\n '\
            '   fineGridProjector1d[7][1][6][4] = 6.504783375660e-02;\n '\
            '   fineGridProjector1d[7][1][6][5] = 6.381407609525e-03;\n '\
            '   fineGridProjector1d[7][1][6][6] = -6.950726400108e-02;\n '\
            '   fineGridProjector1d[7][1][6][7] = -1.108802057964e-01;\n '\
            '   fineGridProjector1d[7][1][7][0] = 1.285330378235e-02;\n '\
            '   fineGridProjector1d[7][1][7][1] = 9.378967905442e-03;\n '\
            '   fineGridProjector1d[7][1][7][2] = -1.087803748778e-03;\n '\
            '   fineGridProjector1d[7][1][7][3] = -1.458790409097e-02;\n '\
            '   fineGridProjector1d[7][1][7][4] = -1.657195712217e-02;\n '\
            '   fineGridProjector1d[7][1][7][5] = -1.573234145260e-03;\n '\
            '   fineGridProjector1d[7][1][7][6] = 1.654900747731e-02;\n '\
            '   fineGridProjector1d[7][1][7][7] = 2.570660756470e-02;\n '\
            '   fineGridProjector1d[7][2][0][0] = 1.344791297770e-02;\n '\
            '   fineGridProjector1d[7][2][0][1] = 1.226471214876e-02;\n '\
            '   fineGridProjector1d[7][2][0][2] = 4.186372100316e-03;\n '\
            '   fineGridProjector1d[7][2][0][3] = -9.040030560746e-03;\n '\
            '   fineGridProjector1d[7][2][0][4] = -9.684592833815e-03;\n '\
            '   fineGridProjector1d[7][2][0][5] = 8.078296063724e-03;\n '\
            '   fineGridProjector1d[7][2][0][6] = 1.028050189740e-02;\n '\
            '   fineGridProjector1d[7][2][0][7] = -1.857097195694e-02;\n '\
            '   fineGridProjector1d[7][2][1][0] = -4.935970743790e-02;\n '\
            '   fineGridProjector1d[7][2][1][1] = -4.476019539319e-02;\n '\
            '   fineGridProjector1d[7][2][1][2] = -1.514937608035e-02;\n '\
            '   fineGridProjector1d[7][2][1][3] = 3.241364216707e-02;\n '\
            '   fineGridProjector1d[7][2][1][4] = 3.443368273302e-02;\n '\
            '   fineGridProjector1d[7][2][1][5] = -2.852872284388e-02;\n '\
            '   fineGridProjector1d[7][2][1][6] = -3.613355466876e-02;\n '\
            '   fineGridProjector1d[7][2][1][7] = 6.510002774521e-02;\n '\
            '   fineGridProjector1d[7][2][2][0] = 1.081757500422e-01;\n '\
            '   fineGridProjector1d[7][2][2][1] = 9.672621942444e-02;\n '\
            '   fineGridProjector1d[7][2][2][2] = 3.207906308248e-02;\n '\
            '   fineGridProjector1d[7][2][2][3] = -6.717992322398e-02;\n '\
            '   fineGridProjector1d[7][2][2][4] = -7.002016468091e-02;\n '\
            '   fineGridProjector1d[7][2][2][5] = 5.715210354087e-02;\n '\
            '   fineGridProjector1d[7][2][2][6] = 7.164436760859e-02;\n '\
            '   fineGridProjector1d[7][2][2][7] = -1.283480353197e-01;\n '\
            '   fineGridProjector1d[7][2][3][0] = -2.211431752079e-01;\n '\
            '   fineGridProjector1d[7][2][3][1] = -1.904996614855e-01;\n '\
            '   fineGridProjector1d[7][2][3][2] = -6.005551773248e-02;\n '\
            '   fineGridProjector1d[7][2][3][3] = 1.196535884089e-01;\n '\
            '   fineGridProjector1d[7][2][3][4] = 1.196502202314e-01;\n '\
            '   fineGridProjector1d[7][2][3][5] = -9.469656069315e-02;\n '\
            '   fineGridProjector1d[7][2][3][6] = -1.163032031210e-01;\n '\
            '   fineGridProjector1d[7][2][3][7] = 2.060736411242e-01;\n '\
            '   fineGridProjector1d[7][2][4][0] = 7.184640259508e-01;\n '\
            '   fineGridProjector1d[7][2][4][1] = 5.115653916979e-01;\n '\
            '   fineGridProjector1d[7][2][4][2] = 1.315770211480e-01;\n '\
            '   fineGridProjector1d[7][2][4][3] = -2.236539758146e-01;\n '\
            '   fineGridProjector1d[7][2][4][4] = -2.002855123277e-01;\n '\
            '   fineGridProjector1d[7][2][4][5] = 1.474619898397e-01;\n '\
            '   fineGridProjector1d[7][2][4][6] = 1.732860881489e-01;\n '\
            '   fineGridProjector1d[7][2][4][7] = -3.001846752984e-01;\n '\
            '   fineGridProjector1d[7][2][5][0] = 5.271518345152e-01;\n '\
            '   fineGridProjector1d[7][2][5][1] = 7.203812928376e-01;\n '\
            '   fineGridProjector1d[7][2][5][2] = 9.583437398281e-01;\n '\
            '   fineGridProjector1d[7][2][5][3] = 9.499278267926e-01;\n '\
            '   fineGridProjector1d[7][2][5][4] = 4.338527015142e-01;\n '\
            '   fineGridProjector1d[7][2][5][5] = -2.470614882914e-01;\n '\
            '   fineGridProjector1d[7][2][5][6] = -2.568048602856e-01;\n '\
            '   fineGridProjector1d[7][2][5][7] = 4.208308498076e-01;\n '\
            '   fineGridProjector1d[7][2][6][0] = -1.253727350632e-01;\n '\
            '   fineGridProjector1d[7][2][6][1] = -1.355379672259e-01;\n '\
            '   fineGridProjector1d[7][2][6][2] = -6.394563419373e-02;\n '\
            '   fineGridProjector1d[7][2][6][3] = 2.377781367632e-01;\n '\
            '   fineGridProjector1d[7][2][6][4] = 7.623766767031e-01;\n '\
            '   fineGridProjector1d[7][2][6][5] = 1.034684310493e+00;\n '\
            '   fineGridProjector1d[7][2][6][6] = 4.608504141649e-01;\n '\
            '   fineGridProjector1d[7][2][6][7] = -6.107484297948e-01;\n '\
            '   fineGridProjector1d[7][2][7][0] = 2.863609422315e-02;\n '\
            '   fineGridProjector1d[7][2][7][1] = 2.986020799602e-02;\n '\
            '   fineGridProjector1d[7][2][7][2] = 1.296433184770e-02;\n '\
            '   fineGridProjector1d[7][2][7][3] = -3.989926453250e-02;\n '\
            '   fineGridProjector1d[7][2][7][4] = -7.032301133937e-02;\n '\
            '   fineGridProjector1d[7][2][7][5] = 1.229100718913e-01;\n '\
            '   fineGridProjector1d[7][2][7][6] = 6.931802462555e-01;\n '\
            '   fineGridProjector1d[7][2][7][7] = 1.365847593693e+00;\n '\
            '   \n '\
			'   equidistantGridProjector1d[8][0][0] = 1.584919424220e+00;\n'\
			'   equidistantGridProjector1d[8][0][1] = -7.178401254211e-02;\n'\
			'   equidistantGridProjector1d[8][0][2] = 2.887705048168e-02;\n'\
			'   equidistantGridProjector1d[8][0][3] = -1.242504204771e-02;\n'\
			'   equidistantGridProjector1d[8][0][4] = 0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[8][0][5] = 7.325121017319e-03;\n'\
			'   equidistantGridProjector1d[8][0][6] = -9.208182122077e-03;\n'\
			'   equidistantGridProjector1d[8][0][7] = 9.114643098419e-03;\n'\
			'   equidistantGridProjector1d[8][0][8] = 2.563991175815e-02;\n'\
			'   equidistantGridProjector1d[8][1][0] = -1.005682886391e+00;\n'\
			'   equidistantGridProjector1d[8][1][1] = 5.948306835392e-01;\n'\
			'   equidistantGridProjector1d[8][1][2] = -1.314659908569e-01;\n'\
			'   equidistantGridProjector1d[8][1][3] = 4.975584405363e-02;\n'\
			'   equidistantGridProjector1d[8][1][4] = -0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[8][1][5] = -2.684865303582e-02;\n'\
			'   equidistantGridProjector1d[8][1][6] = 3.306559423741e-02;\n'\
			'   equidistantGridProjector1d[8][1][7] = -3.226540900764e-02;\n'\
			'   equidistantGridProjector1d[8][1][8] = -8.981367941100e-02;\n'\
			'   equidistantGridProjector1d[8][2][0] = 7.373969413835e-01;\n'\
			'   equidistantGridProjector1d[8][2][1] = 6.475601423470e-01;\n'\
			'   equidistantGridProjector1d[8][2][2] = 6.736929701053e-01;\n'\
			'   equidistantGridProjector1d[8][2][3] = -1.387352224019e-01;\n'\
			'   equidistantGridProjector1d[8][2][4] = 0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[8][2][5] = 5.839018366929e-02;\n'\
			'   equidistantGridProjector1d[8][2][6] = -6.860023077464e-02;\n'\
			'   equidistantGridProjector1d[8][2][7] = 6.489443182319e-02;\n'\
			'   equidistantGridProjector1d[8][2][8] = 1.767099114311e-01;\n'\
			'   equidistantGridProjector1d[8][3][0] = -5.532193350361e-01;\n'\
			'   equidistantGridProjector1d[8][3][1] = -2.724936537078e-01;\n'\
			'   equidistantGridProjector1d[8][3][2] = 5.698554536819e-01;\n'\
			'   equidistantGridProjector1d[8][3][3] = 8.902397194841e-01;\n'\
			'   equidistantGridProjector1d[8][3][4] = -0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[8][3][5] = -1.151118027820e-01;\n'\
			'   equidistantGridProjector1d[8][3][6] = 1.215040693743e-01;\n'\
			'   equidistantGridProjector1d[8][3][7] = -1.079942941608e-01;\n'\
			'   equidistantGridProjector1d[8][3][8] = -2.822994943042e-01;\n'\
			'   equidistantGridProjector1d[8][4][0] = 4.063492063492e-01;\n'\
			'   equidistantGridProjector1d[8][4][1] = 1.681374686105e-01;\n'\
			'   equidistantGridProjector1d[8][4][2] = -2.177207341270e-01;\n'\
			'   equidistantGridProjector1d[8][4][3] = 2.874098520430e-01;\n'\
			'   equidistantGridProjector1d[8][4][4] = 1.000000000000e+00;\n'\
			'   equidistantGridProjector1d[8][4][5] = 2.874098520430e-01;\n'\
			'   equidistantGridProjector1d[8][4][6] = -2.177207341270e-01;\n'\
			'   equidistantGridProjector1d[8][4][7] = 1.681374686105e-01;\n'\
			'   equidistantGridProjector1d[8][4][8] = 4.063492063492e-01;\n'\
			'   equidistantGridProjector1d[8][5][0] = -2.822994943042e-01;\n'\
			'   equidistantGridProjector1d[8][5][1] = -1.079942941608e-01;\n'\
			'   equidistantGridProjector1d[8][5][2] = 1.215040693743e-01;\n'\
			'   equidistantGridProjector1d[8][5][3] = -1.151118027820e-01;\n'\
			'   equidistantGridProjector1d[8][5][4] = 0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[8][5][5] = 8.902397194841e-01;\n'\
			'   equidistantGridProjector1d[8][5][6] = 5.698554536819e-01;\n'\
			'   equidistantGridProjector1d[8][5][7] = -2.724936537078e-01;\n'\
			'   equidistantGridProjector1d[8][5][8] = -5.532193350361e-01;\n'\
			'   equidistantGridProjector1d[8][6][0] = 1.767099114311e-01;\n'\
			'   equidistantGridProjector1d[8][6][1] = 6.489443182319e-02;\n'\
			'   equidistantGridProjector1d[8][6][2] = -6.860023077464e-02;\n'\
			'   equidistantGridProjector1d[8][6][3] = 5.839018366929e-02;\n'\
			'   equidistantGridProjector1d[8][6][4] = -0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[8][6][5] = -1.387352224019e-01;\n'\
			'   equidistantGridProjector1d[8][6][6] = 6.736929701053e-01;\n'\
			'   equidistantGridProjector1d[8][6][7] = 6.475601423470e-01;\n'\
			'   equidistantGridProjector1d[8][6][8] = 7.373969413835e-01;\n'\
			'   equidistantGridProjector1d[8][7][0] = -8.981367941100e-02;\n'\
			'   equidistantGridProjector1d[8][7][1] = -3.226540900764e-02;\n'\
			'   equidistantGridProjector1d[8][7][2] = 3.306559423741e-02;\n'\
			'   equidistantGridProjector1d[8][7][3] = -2.684865303582e-02;\n'\
			'   equidistantGridProjector1d[8][7][4] = 0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[8][7][5] = 4.975584405363e-02;\n'\
			'   equidistantGridProjector1d[8][7][6] = -1.314659908569e-01;\n'\
			'   equidistantGridProjector1d[8][7][7] = 5.948306835392e-01;\n'\
			'   equidistantGridProjector1d[8][7][8] = -1.005682886391e+00;\n'\
			'   equidistantGridProjector1d[8][8][0] = 2.563991175815e-02;\n'\
			'   equidistantGridProjector1d[8][8][1] = 9.114643098419e-03;\n'\
			'   equidistantGridProjector1d[8][8][2] = -9.208182122077e-03;\n'\
			'   equidistantGridProjector1d[8][8][3] = 7.325121017319e-03;\n'\
			'   equidistantGridProjector1d[8][8][4] = -0.000000000000e+00;\n'\
			'   equidistantGridProjector1d[8][8][5] = -1.242504204771e-02;\n'\
			'   equidistantGridProjector1d[8][8][6] = 2.887705048168e-02;\n'\
			'   equidistantGridProjector1d[8][8][7] = -7.178401254211e-02;\n'\
			'   equidistantGridProjector1d[8][8][8] = 1.584919424220e+00;\n'\
            '   fineGridProjector1d[8][0][0][0] = 1.368109827999e+00;\n '\
            '   fineGridProjector1d[8][0][0][1] = 6.887407902660e-01;\n '\
            '   fineGridProjector1d[8][0][0][2] = 1.150565385411e-01;\n '\
            '   fineGridProjector1d[8][0][0][3] = -7.084083809576e-02;\n '\
            '   fineGridProjector1d[8][0][0][4] = -3.076782684190e-02;\n '\
            '   fineGridProjector1d[8][0][0][5] = 2.100500706413e-02;\n '\
            '   fineGridProjector1d[8][0][0][6] = 2.712272997683e-02;\n '\
            '   fineGridProjector1d[8][0][0][7] = 1.416709359521e-02;\n '\
            '   fineGridProjector1d[8][0][0][8] = 4.261981043666e-03;\n '\
            '   fineGridProjector1d[8][0][1][0] = -6.187928527741e-01;\n '\
            '   fineGridProjector1d[8][0][1][1] = 4.697652355865e-01;\n '\
            '   fineGridProjector1d[8][0][1][2] = 1.039620473952e+00;\n '\
            '   fineGridProjector1d[8][0][1][3] = 7.306133731926e-01;\n '\
            '   fineGridProjector1d[8][0][1][4] = 1.789773990040e-01;\n '\
            '   fineGridProjector1d[8][0][1][5] = -1.013262605992e-01;\n '\
            '   fineGridProjector1d[8][0][1][6] = -1.199562235996e-01;\n '\
            '   fineGridProjector1d[8][0][1][7] = -5.994649707373e-02;\n '\
            '   fineGridProjector1d[8][0][1][8] = -1.766650402728e-02;\n '\
            '   fineGridProjector1d[8][0][2][0] = 4.363271678107e-01;\n '\
            '   fineGridProjector1d[8][0][2][1] = -2.674373249127e-01;\n '\
            '   fineGridProjector1d[8][0][2][2] = -2.447156358544e-01;\n '\
            '   fineGridProjector1d[8][0][2][3] = 4.796560826380e-01;\n '\
            '   fineGridProjector1d[8][0][2][4] = 9.833424775980e-01;\n '\
            '   fineGridProjector1d[8][0][2][5] = 8.871202945778e-01;\n '\
            '   fineGridProjector1d[8][0][2][6] = 5.128820931979e-01;\n '\
            '   fineGridProjector1d[8][0][2][7] = 2.060325947833e-01;\n '\
            '   fineGridProjector1d[8][0][2][8] = 5.578602610345e-02;\n '\
            '   fineGridProjector1d[8][0][3][0] = -3.234409572727e-01;\n '\
            '   fineGridProjector1d[8][0][3][1] = 1.874368719192e-01;\n '\
            '   fineGridProjector1d[8][0][3][2] = 1.512397410921e-01;\n '\
            '   fineGridProjector1d[8][0][3][3] = -2.253057212804e-01;\n '\
            '   fineGridProjector1d[8][0][3][4] = -2.006913028073e-01;\n '\
            '   fineGridProjector1d[8][0][3][5] = 2.719807184868e-01;\n '\
            '   fineGridProjector1d[8][0][3][6] = 7.368938551619e-01;\n '\
            '   fineGridProjector1d[8][0][3][7] = 9.553323694180e-01;\n '\
            '   fineGridProjector1d[8][0][3][8] = 1.000768244617e+00;\n '\
            '   fineGridProjector1d[8][0][4][0] = 2.363501286965e-01;\n '\
            '   fineGridProjector1d[8][0][4][1] = -1.338561210733e-01;\n '\
            '   fineGridProjector1d[8][0][4][2] = -1.032021296234e-01;\n '\
            '   fineGridProjector1d[8][0][4][3] = 1.424035453955e-01;\n '\
            '   fineGridProjector1d[8][0][4][4] = 1.120439925195e-01;\n '\
            '   fineGridProjector1d[8][0][4][5] = -1.240208925713e-01;\n '\
            '   fineGridProjector1d[8][0][4][6] = -2.390700449409e-01;\n '\
            '   fineGridProjector1d[8][0][4][7] = -1.705842000953e-01;\n '\
            '   fineGridProjector1d[8][0][4][8] = -6.228395877434e-02;\n '\
            '   fineGridProjector1d[8][0][5][0] = -1.637673418227e-01;\n '\
            '   fineGridProjector1d[8][0][5][1] = 9.169461454539e-02;\n '\
            '   fineGridProjector1d[8][0][5][2] = 6.919032205266e-02;\n '\
            '   fineGridProjector1d[8][0][5][3] = -9.235604337196e-02;\n '\
            '   fineGridProjector1d[8][0][5][4] = -6.934904295262e-02;\n '\
            '   fineGridProjector1d[8][0][5][5] = 7.219113702779e-02;\n '\
            '   fineGridProjector1d[8][0][5][6] = 1.292610188375e-01;\n '\
            '   fineGridProjector1d[8][0][5][7] = 8.548944317449e-02;\n '\
            '   fineGridProjector1d[8][0][5][8] = 2.949460085591e-02;\n '\
            '   fineGridProjector1d[8][0][6][0] = 1.023645981368e-01;\n '\
            '   fineGridProjector1d[8][0][6][1] = -5.695830736689e-02;\n '\
            '   fineGridProjector1d[8][0][6][2] = -4.248988832803e-02;\n '\
            '   fineGridProjector1d[8][0][6][3] = 5.576363371921e-02;\n '\
            '   fineGridProjector1d[8][0][6][4] = 4.094211526850e-02;\n '\
            '   fineGridProjector1d[8][0][6][5] = -4.147318301770e-02;\n '\
            '   fineGridProjector1d[8][0][6][6] = -7.208033628684e-02;\n '\
            '   fineGridProjector1d[8][0][6][7] = -4.637288100377e-02;\n '\
            '   fineGridProjector1d[8][0][6][8] = -1.570025966080e-02;\n '\
            '   fineGridProjector1d[8][0][7][0] = -5.198556224023e-02;\n '\
            '   fineGridProjector1d[8][0][7][1] = 2.882675600633e-02;\n '\
            '   fineGridProjector1d[8][0][7][2] = 2.137069854170e-02;\n '\
            '   fineGridProjector1d[8][0][7][3] = -2.779516775248e-02;\n '\
            '   fineGridProjector1d[8][0][7][4] = -2.017199173226e-02;\n '\
            '   fineGridProjector1d[8][0][7][5] = 2.015817652592e-02;\n '\
            '   fineGridProjector1d[8][0][7][6] = 3.454076586106e-02;\n '\
            '   fineGridProjector1d[8][0][7][7] = 2.194286111129e-02;\n '\
            '   fineGridProjector1d[8][0][7][8] = 7.367438559345e-03;\n '\
            '   fineGridProjector1d[8][0][8][0] = 1.483499146706e-02;\n '\
            '   fineGridProjector1d[8][0][8][1] = -8.212514970601e-03;\n '\
            '   fineGridProjector1d[8][0][8][2] = -6.070120374007e-03;\n '\
            '   fineGridProjector1d[8][0][8][3] = 7.861135555292e-03;\n '\
            '   fineGridProjector1d[8][0][8][4] = 5.674179944117e-03;\n '\
            '   fineGridProjector1d[8][0][8][5] = -5.634997494254e-03;\n '\
            '   fineGridProjector1d[8][0][8][6] = -9.593858207808e-03;\n '\
            '   fineGridProjector1d[8][0][8][7] = -6.060783909557e-03;\n '\
            '   fineGridProjector1d[8][0][8][8] = -2.027568717351e-03;\n '\
            '   fineGridProjector1d[8][1][0][0] = -3.204950709587e-04;\n '\
            '   fineGridProjector1d[8][1][0][1] = -8.487295119646e-03;\n '\
            '   fineGridProjector1d[8][1][0][2] = -1.599181949454e-02;\n '\
            '   fineGridProjector1d[8][1][0][3] = -1.305225432363e-02;\n '\
            '   fineGridProjector1d[8][1][0][4] = 0.000000000000e+00;\n '\
            '   fineGridProjector1d[8][1][0][5] = 1.043064600190e-02;\n '\
            '   fineGridProjector1d[8][1][0][6] = 1.041516373233e-02;\n '\
            '   fineGridProjector1d[8][1][0][7] = 4.693352970847e-03;\n '\
            '   fineGridProjector1d[8][1][0][8] = 1.602475354794e-04;\n '\
            '   fineGridProjector1d[8][1][1][0] = 1.316868131421e-03;\n '\
            '   fineGridProjector1d[8][1][1][1] = 3.430894470685e-02;\n '\
            '   fineGridProjector1d[8][1][1][2] = 6.318933073355e-02;\n '\
            '   fineGridProjector1d[8][1][1][3] = 5.039277096846e-02;\n '\
            '   fineGridProjector1d[8][1][1][4] = -0.000000000000e+00;\n '\
            '   fineGridProjector1d[8][1][1][5] = -3.885461553042e-02;\n '\
            '   fineGridProjector1d[8][1][1][6] = -3.835577766183e-02;\n '\
            '   fineGridProjector1d[8][1][1][7] = -1.715447235343e-02;\n '\
            '   fineGridProjector1d[8][1][1][8] = -5.833546848525e-04;\n '\
            '   fineGridProjector1d[8][1][2][0] = -4.020900087080e-03;\n '\
            '   fineGridProjector1d[8][1][2][1] = -9.877859498154e-02;\n '\
            '   fineGridProjector1d[8][1][2][2] = -1.687362665017e-01;\n '\
            '   fineGridProjector1d[8][1][2][3] = -1.255169682468e-01;\n '\
            '   fineGridProjector1d[8][1][2][4] = 0.000000000000e+00;\n '\
            '   fineGridProjector1d[8][1][2][5] = 8.790847293900e-02;\n '\
            '   fineGridProjector1d[8][1][2][6] = 8.436813325085e-02;\n '\
            '   fineGridProjector1d[8][1][2][7] = 3.706148650524e-02;\n '\
            '   fineGridProjector1d[8][1][2][8] = 1.248467734610e-03;\n '\
            '   fineGridProjector1d[8][1][3][0] = 9.994063101772e-01;\n '\
            '   fineGridProjector1d[8][1][3][1] = 9.511681923937e-01;\n '\
            '   fineGridProjector1d[8][1][3][2] = 7.552371157029e-01;\n '\
            '   fineGridProjector1d[8][1][3][3] = 3.847097129182e-01;\n '\
            '   fineGridProjector1d[8][1][3][4] = -0.000000000000e+00;\n '\
            '   fineGridProjector1d[8][1][3][5] = -1.923548564591e-01;\n '\
            '   fineGridProjector1d[8][1][3][6] = -1.711231094865e-01;\n '\
            '   fineGridProjector1d[8][1][3][7] = -7.190018280917e-02;\n '\
            '   fineGridProjector1d[8][1][3][8] = -2.368618235197e-03;\n '\
            '   fineGridProjector1d[8][1][4][0] = 5.161474499376e-03;\n '\
            '   fineGridProjector1d[8][1][4][1] = 1.690885686871e-01;\n '\
            '   fineGridProjector1d[8][1][4][2] = 4.809972297250e-01;\n '\
            '   fineGridProjector1d[8][1][4][3] = 8.363370917324e-01;\n '\
            '   fineGridProjector1d[8][1][4][4] = 1.000000000000e+00;\n '\
            '   fineGridProjector1d[8][1][4][5] = 8.363370917324e-01;\n '\
            '   fineGridProjector1d[8][1][4][6] = 4.809972297250e-01;\n '\
            '   fineGridProjector1d[8][1][4][7] = 1.690885686871e-01;\n '\
            '   fineGridProjector1d[8][1][4][8] = 5.161474499377e-03;\n '\
            '   fineGridProjector1d[8][1][5][0] = -2.368618235196e-03;\n '\
            '   fineGridProjector1d[8][1][5][1] = -7.190018280917e-02;\n '\
            '   fineGridProjector1d[8][1][5][2] = -1.711231094865e-01;\n '\
            '   fineGridProjector1d[8][1][5][3] = -1.923548564591e-01;\n '\
            '   fineGridProjector1d[8][1][5][4] = 0.000000000000e+00;\n '\
            '   fineGridProjector1d[8][1][5][5] = 3.847097129182e-01;\n '\
            '   fineGridProjector1d[8][1][5][6] = 7.552371157029e-01;\n '\
            '   fineGridProjector1d[8][1][5][7] = 9.511681923937e-01;\n '\
            '   fineGridProjector1d[8][1][5][8] = 9.994063101772e-01;\n '\
            '   fineGridProjector1d[8][1][6][0] = 1.248467734609e-03;\n '\
            '   fineGridProjector1d[8][1][6][1] = 3.706148650524e-02;\n '\
            '   fineGridProjector1d[8][1][6][2] = 8.436813325085e-02;\n '\
            '   fineGridProjector1d[8][1][6][3] = 8.790847293900e-02;\n '\
            '   fineGridProjector1d[8][1][6][4] = -0.000000000000e+00;\n '\
            '   fineGridProjector1d[8][1][6][5] = -1.255169682468e-01;\n '\
            '   fineGridProjector1d[8][1][6][6] = -1.687362665017e-01;\n '\
            '   fineGridProjector1d[8][1][6][7] = -9.877859498154e-02;\n '\
            '   fineGridProjector1d[8][1][6][8] = -4.020900087080e-03;\n '\
            '   fineGridProjector1d[8][1][7][0] = -5.833546848524e-04;\n '\
            '   fineGridProjector1d[8][1][7][1] = -1.715447235343e-02;\n '\
            '   fineGridProjector1d[8][1][7][2] = -3.835577766183e-02;\n '\
            '   fineGridProjector1d[8][1][7][3] = -3.885461553042e-02;\n '\
            '   fineGridProjector1d[8][1][7][4] = 0.000000000000e+00;\n '\
            '   fineGridProjector1d[8][1][7][5] = 5.039277096846e-02;\n '\
            '   fineGridProjector1d[8][1][7][6] = 6.318933073355e-02;\n '\
            '   fineGridProjector1d[8][1][7][7] = 3.430894470685e-02;\n '\
            '   fineGridProjector1d[8][1][7][8] = 1.316868131421e-03;\n '\
            '   fineGridProjector1d[8][1][8][0] = 1.602475354794e-04;\n '\
            '   fineGridProjector1d[8][1][8][1] = 4.693352970847e-03;\n '\
            '   fineGridProjector1d[8][1][8][2] = 1.041516373233e-02;\n '\
            '   fineGridProjector1d[8][1][8][3] = 1.043064600190e-02;\n '\
            '   fineGridProjector1d[8][1][8][4] = -0.000000000000e+00;\n '\
            '   fineGridProjector1d[8][1][8][5] = -1.305225432363e-02;\n '\
            '   fineGridProjector1d[8][1][8][6] = -1.599181949454e-02;\n '\
            '   fineGridProjector1d[8][1][8][7] = -8.487295119646e-03;\n '\
            '   fineGridProjector1d[8][1][8][8] = -3.204950709588e-04;\n '\
            '   fineGridProjector1d[8][2][0][0] = -2.027568717351e-03;\n '\
            '   fineGridProjector1d[8][2][0][1] = -6.060783909557e-03;\n '\
            '   fineGridProjector1d[8][2][0][2] = -9.593858207808e-03;\n '\
            '   fineGridProjector1d[8][2][0][3] = -5.634997494254e-03;\n '\
            '   fineGridProjector1d[8][2][0][4] = 5.674179944117e-03;\n '\
            '   fineGridProjector1d[8][2][0][5] = 7.861135555292e-03;\n '\
            '   fineGridProjector1d[8][2][0][6] = -6.070120374007e-03;\n '\
            '   fineGridProjector1d[8][2][0][7] = -8.212514970601e-03;\n '\
            '   fineGridProjector1d[8][2][0][8] = 1.483499146706e-02;\n '\
            '   fineGridProjector1d[8][2][1][0] = 7.367438559345e-03;\n '\
            '   fineGridProjector1d[8][2][1][1] = 2.194286111129e-02;\n '\
            '   fineGridProjector1d[8][2][1][2] = 3.454076586106e-02;\n '\
            '   fineGridProjector1d[8][2][1][3] = 2.015817652592e-02;\n '\
            '   fineGridProjector1d[8][2][1][4] = -2.017199173226e-02;\n '\
            '   fineGridProjector1d[8][2][1][5] = -2.779516775248e-02;\n '\
            '   fineGridProjector1d[8][2][1][6] = 2.137069854170e-02;\n '\
            '   fineGridProjector1d[8][2][1][7] = 2.882675600633e-02;\n '\
            '   fineGridProjector1d[8][2][1][8] = -5.198556224023e-02;\n '\
            '   fineGridProjector1d[8][2][2][0] = -1.570025966080e-02;\n '\
            '   fineGridProjector1d[8][2][2][1] = -4.637288100377e-02;\n '\
            '   fineGridProjector1d[8][2][2][2] = -7.208033628684e-02;\n '\
            '   fineGridProjector1d[8][2][2][3] = -4.147318301770e-02;\n '\
            '   fineGridProjector1d[8][2][2][4] = 4.094211526850e-02;\n '\
            '   fineGridProjector1d[8][2][2][5] = 5.576363371921e-02;\n '\
            '   fineGridProjector1d[8][2][2][6] = -4.248988832803e-02;\n '\
            '   fineGridProjector1d[8][2][2][7] = -5.695830736689e-02;\n '\
            '   fineGridProjector1d[8][2][2][8] = 1.023645981368e-01;\n '\
            '   fineGridProjector1d[8][2][3][0] = 2.949460085591e-02;\n '\
            '   fineGridProjector1d[8][2][3][1] = 8.548944317449e-02;\n '\
            '   fineGridProjector1d[8][2][3][2] = 1.292610188375e-01;\n '\
            '   fineGridProjector1d[8][2][3][3] = 7.219113702779e-02;\n '\
            '   fineGridProjector1d[8][2][3][4] = -6.934904295262e-02;\n '\
            '   fineGridProjector1d[8][2][3][5] = -9.235604337196e-02;\n '\
            '   fineGridProjector1d[8][2][3][6] = 6.919032205266e-02;\n '\
            '   fineGridProjector1d[8][2][3][7] = 9.169461454539e-02;\n '\
            '   fineGridProjector1d[8][2][3][8] = -1.637673418227e-01;\n '\
            '   fineGridProjector1d[8][2][4][0] = -6.228395877433e-02;\n '\
            '   fineGridProjector1d[8][2][4][1] = -1.705842000953e-01;\n '\
            '   fineGridProjector1d[8][2][4][2] = -2.390700449409e-01;\n '\
            '   fineGridProjector1d[8][2][4][3] = -1.240208925713e-01;\n '\
            '   fineGridProjector1d[8][2][4][4] = 1.120439925195e-01;\n '\
            '   fineGridProjector1d[8][2][4][5] = 1.424035453955e-01;\n '\
            '   fineGridProjector1d[8][2][4][6] = -1.032021296234e-01;\n '\
            '   fineGridProjector1d[8][2][4][7] = -1.338561210733e-01;\n '\
            '   fineGridProjector1d[8][2][4][8] = 2.363501286965e-01;\n '\
            '   fineGridProjector1d[8][2][5][0] = 1.000768244617e+00;\n '\
            '   fineGridProjector1d[8][2][5][1] = 9.553323694180e-01;\n '\
            '   fineGridProjector1d[8][2][5][2] = 7.368938551619e-01;\n '\
            '   fineGridProjector1d[8][2][5][3] = 2.719807184868e-01;\n '\
            '   fineGridProjector1d[8][2][5][4] = -2.006913028073e-01;\n '\
            '   fineGridProjector1d[8][2][5][5] = -2.253057212804e-01;\n '\
            '   fineGridProjector1d[8][2][5][6] = 1.512397410921e-01;\n '\
            '   fineGridProjector1d[8][2][5][7] = 1.874368719192e-01;\n '\
            '   fineGridProjector1d[8][2][5][8] = -3.234409572727e-01;\n '\
            '   fineGridProjector1d[8][2][6][0] = 5.578602610345e-02;\n '\
            '   fineGridProjector1d[8][2][6][1] = 2.060325947833e-01;\n '\
            '   fineGridProjector1d[8][2][6][2] = 5.128820931979e-01;\n '\
            '   fineGridProjector1d[8][2][6][3] = 8.871202945778e-01;\n '\
            '   fineGridProjector1d[8][2][6][4] = 9.833424775980e-01;\n '\
            '   fineGridProjector1d[8][2][6][5] = 4.796560826380e-01;\n '\
            '   fineGridProjector1d[8][2][6][6] = -2.447156358544e-01;\n '\
            '   fineGridProjector1d[8][2][6][7] = -2.674373249127e-01;\n '\
            '   fineGridProjector1d[8][2][6][8] = 4.363271678107e-01;\n '\
            '   fineGridProjector1d[8][2][7][0] = -1.766650402728e-02;\n '\
            '   fineGridProjector1d[8][2][7][1] = -5.994649707373e-02;\n '\
            '   fineGridProjector1d[8][2][7][2] = -1.199562235996e-01;\n '\
            '   fineGridProjector1d[8][2][7][3] = -1.013262605992e-01;\n '\
            '   fineGridProjector1d[8][2][7][4] = 1.789773990040e-01;\n '\
            '   fineGridProjector1d[8][2][7][5] = 7.306133731926e-01;\n '\
            '   fineGridProjector1d[8][2][7][6] = 1.039620473952e+00;\n '\
            '   fineGridProjector1d[8][2][7][7] = 4.697652355865e-01;\n '\
            '   fineGridProjector1d[8][2][7][8] = -6.187928527741e-01;\n '\
            '   fineGridProjector1d[8][2][8][0] = 4.261981043666e-03;\n '\
            '   fineGridProjector1d[8][2][8][1] = 1.416709359521e-02;\n '\
            '   fineGridProjector1d[8][2][8][2] = 2.712272997683e-02;\n '\
            '   fineGridProjector1d[8][2][8][3] = 2.100500706413e-02;\n '\
            '   fineGridProjector1d[8][2][8][4] = -3.076782684190e-02;\n '\
            '   fineGridProjector1d[8][2][8][5] = -7.084083809576e-02;\n '\
            '   fineGridProjector1d[8][2][8][6] = 1.150565385411e-01;\n '\
            '   fineGridProjector1d[8][2][8][7] = 6.887407902660e-01;\n '\
            '   fineGridProjector1d[8][2][8][8] = 1.368109827999e+00;\n '\
            '   \n '\
			'   equidistantGridProjector1d[9][0][0] = 1.588005378675e+00;\n'\
			'   equidistantGridProjector1d[9][0][1] = -6.543263366210e-02;\n'\
			'   equidistantGridProjector1d[9][0][2] = 2.661530002799e-02;\n'\
			'   equidistantGridProjector1d[9][0][3] = -1.489517532211e-02;\n'\
			'   equidistantGridProjector1d[9][0][4] = 4.621049087800e-03;\n'\
			'   equidistantGridProjector1d[9][0][5] = 3.674613048621e-03;\n'\
			'   equidistantGridProjector1d[9][0][6] = -7.298928320166e-03;\n'\
			'   equidistantGridProjector1d[9][0][7] = 7.280034450233e-03;\n'\
			'   equidistantGridProjector1d[9][0][8] = -7.326217773156e-03;\n'\
			'   equidistantGridProjector1d[9][0][9] = -2.099216577097e-02;\n'\
			'   equidistantGridProjector1d[9][1][0] = -1.016287965645e+00;\n'\
			'   equidistantGridProjector1d[9][1][1] = 4.865818879221e-01;\n'\
			'   equidistantGridProjector1d[9][1][2] = -1.190593091440e-01;\n'\
			'   equidistantGridProjector1d[9][1][3] = 5.938626508614e-02;\n'\
			'   equidistantGridProjector1d[9][1][4] = -1.750118232789e-02;\n'\
			'   equidistantGridProjector1d[9][1][5] = -1.351710804930e-02;\n'\
			'   equidistantGridProjector1d[9][1][6] = 2.634975606536e-02;\n'\
			'   equidistantGridProjector1d[9][1][7] = -2.593924774000e-02;\n'\
			'   equidistantGridProjector1d[9][1][8] = 2.585252188494e-02;\n'\
			'   equidistantGridProjector1d[9][1][9] = 7.352805218734e-02;\n'\
			'   equidistantGridProjector1d[9][2][0] = 7.575227986515e-01;\n'\
			'   equidistantGridProjector1d[9][2][1] = 7.646159098158e-01;\n'\
			'   equidistantGridProjector1d[9][2][2] = 5.268955145225e-01;\n'\
			'   equidistantGridProjector1d[9][2][3] = -1.615864863197e-01;\n'\
			'   equidistantGridProjector1d[9][2][4] = 4.111822591127e-02;\n'\
			'   equidistantGridProjector1d[9][2][5] = 2.955953583038e-02;\n'\
			'   equidistantGridProjector1d[9][2][6] = -5.521761023130e-02;\n'\
			'   equidistantGridProjector1d[9][2][7] = 5.284207829810e-02;\n'\
			'   equidistantGridProjector1d[9][2][8] = -5.161580505763e-02;\n'\
			'   equidistantGridProjector1d[9][2][9] = -1.446071081332e-01;\n'\
			'   equidistantGridProjector1d[9][3][0] = -5.836053892999e-01;\n'\
			'   equidistantGridProjector1d[9][3][1] = -2.973790015347e-01;\n'\
			'   equidistantGridProjector1d[9][3][2] = 7.273745729843e-01;\n'\
			'   equidistantGridProjector1d[9][3][3] = 7.609576748703e-01;\n'\
			'   equidistantGridProjector1d[9][3][4] = -9.872447857643e-02;\n'\
			'   equidistantGridProjector1d[9][3][5] = -5.843336578743e-02;\n'\
			'   equidistantGridProjector1d[9][3][6] = 9.930890861228e-02;\n'\
			'   equidistantGridProjector1d[9][3][7] = -8.984894071994e-02;\n'\
			'   equidistantGridProjector1d[9][3][8] = 8.455610788859e-02;\n'\
			'   equidistantGridProjector1d[9][3][9] = 2.306924543937e-01;\n'\
			'   equidistantGridProjector1d[9][4][0] = 4.466023128803e-01;\n'\
			'   equidistantGridProjector1d[9][4][1] = 1.871899444062e-01;\n'\
			'   equidistantGridProjector1d[9][4][2] = -2.511590647618e-01;\n'\
			'   equidistantGridProjector1d[9][4][3] = 4.745099964829e-01;\n'\
			'   equidistantGridProjector1d[9][4][4] = 9.685234137744e-01;\n'\
			'   equidistantGridProjector1d[9][4][5] = 1.406792970885e-01;\n'\
			'   equidistantGridProjector1d[9][4][6] = -1.815144009237e-01;\n'\
			'   equidistantGridProjector1d[9][4][7] = 1.449990620827e-01;\n'\
			'   equidistantGridProjector1d[9][4][8] = -1.270427138901e-01;\n'\
			'   equidistantGridProjector1d[9][4][9] = -3.308583679391e-01;\n'\
			'   equidistantGridProjector1d[9][5][0] = -3.308583679391e-01;\n'\
			'   equidistantGridProjector1d[9][5][1] = -1.270427138901e-01;\n'\
			'   equidistantGridProjector1d[9][5][2] = 1.449990620827e-01;\n'\
			'   equidistantGridProjector1d[9][5][3] = -1.815144009237e-01;\n'\
			'   equidistantGridProjector1d[9][5][4] = 1.406792970885e-01;\n'\
			'   equidistantGridProjector1d[9][5][5] = 9.685234137744e-01;\n'\
			'   equidistantGridProjector1d[9][5][6] = 4.745099964829e-01;\n'\
			'   equidistantGridProjector1d[9][5][7] = -2.511590647618e-01;\n'\
			'   equidistantGridProjector1d[9][5][8] = 1.871899444062e-01;\n'\
			'   equidistantGridProjector1d[9][5][9] = 4.466023128803e-01;\n'\
			'   equidistantGridProjector1d[9][6][0] = 2.306924543937e-01;\n'\
			'   equidistantGridProjector1d[9][6][1] = 8.455610788859e-02;\n'\
			'   equidistantGridProjector1d[9][6][2] = -8.984894071994e-02;\n'\
			'   equidistantGridProjector1d[9][6][3] = 9.930890861228e-02;\n'\
			'   equidistantGridProjector1d[9][6][4] = -5.843336578743e-02;\n'\
			'   equidistantGridProjector1d[9][6][5] = -9.872447857643e-02;\n'\
			'   equidistantGridProjector1d[9][6][6] = 7.609576748703e-01;\n'\
			'   equidistantGridProjector1d[9][6][7] = 7.273745729843e-01;\n'\
			'   equidistantGridProjector1d[9][6][8] = -2.973790015347e-01;\n'\
			'   equidistantGridProjector1d[9][6][9] = -5.836053892999e-01;\n'\
			'   equidistantGridProjector1d[9][7][0] = -1.446071081332e-01;\n'\
			'   equidistantGridProjector1d[9][7][1] = -5.161580505763e-02;\n'\
			'   equidistantGridProjector1d[9][7][2] = 5.284207829810e-02;\n'\
			'   equidistantGridProjector1d[9][7][3] = -5.521761023130e-02;\n'\
			'   equidistantGridProjector1d[9][7][4] = 2.955953583038e-02;\n'\
			'   equidistantGridProjector1d[9][7][5] = 4.111822591127e-02;\n'\
			'   equidistantGridProjector1d[9][7][6] = -1.615864863197e-01;\n'\
			'   equidistantGridProjector1d[9][7][7] = 5.268955145225e-01;\n'\
			'   equidistantGridProjector1d[9][7][8] = 7.646159098158e-01;\n'\
			'   equidistantGridProjector1d[9][7][9] = 7.575227986515e-01;\n'\
			'   equidistantGridProjector1d[9][8][0] = 7.352805218734e-02;\n'\
			'   equidistantGridProjector1d[9][8][1] = 2.585252188494e-02;\n'\
			'   equidistantGridProjector1d[9][8][2] = -2.593924774000e-02;\n'\
			'   equidistantGridProjector1d[9][8][3] = 2.634975606536e-02;\n'\
			'   equidistantGridProjector1d[9][8][4] = -1.351710804930e-02;\n'\
			'   equidistantGridProjector1d[9][8][5] = -1.750118232789e-02;\n'\
			'   equidistantGridProjector1d[9][8][6] = 5.938626508614e-02;\n'\
			'   equidistantGridProjector1d[9][8][7] = -1.190593091440e-01;\n'\
			'   equidistantGridProjector1d[9][8][8] = 4.865818879221e-01;\n'\
			'   equidistantGridProjector1d[9][8][9] = -1.016287965645e+00;\n'\
			'   equidistantGridProjector1d[9][9][0] = -2.099216577097e-02;\n'\
			'   equidistantGridProjector1d[9][9][1] = -7.326217773156e-03;\n'\
			'   equidistantGridProjector1d[9][9][2] = 7.280034450233e-03;\n'\
			'   equidistantGridProjector1d[9][9][3] = -7.298928320166e-03;\n'\
			'   equidistantGridProjector1d[9][9][4] = 3.674613048621e-03;\n'\
			'   equidistantGridProjector1d[9][9][5] = 4.621049087800e-03;\n'\
			'   equidistantGridProjector1d[9][9][6] = -1.489517532211e-02;\n'\
			'   equidistantGridProjector1d[9][9][7] = 2.661530002799e-02;\n'\
			'   equidistantGridProjector1d[9][9][8] = -6.543263366210e-02;\n'\
			'   equidistantGridProjector1d[9][9][9] = 1.588005378675e+00;\n'\
            '   fineGridProjector1d[9][0][0][0] = 1.369756044103e+00;\n '\
            '   fineGridProjector1d[9][0][0][1] = 6.855092197838e-01;\n '\
            '   fineGridProjector1d[9][0][0][2] = 1.095186765437e-01;\n '\
            '   fineGridProjector1d[9][0][0][3] = -7.079373286275e-02;\n '\
            '   fineGridProjector1d[9][0][0][4] = -2.443340119521e-02;\n '\
            '   fineGridProjector1d[9][0][0][5] = 2.445445929348e-02;\n '\
            '   fineGridProjector1d[9][0][0][6] = 2.136690378358e-02;\n '\
            '   fineGridProjector1d[9][0][0][7] = 1.601079014457e-03;\n '\
            '   fineGridProjector1d[9][0][0][8] = -1.053816905737e-02;\n '\
            '   fineGridProjector1d[9][0][0][9] = -1.437946599702e-02;\n '\
            '   fineGridProjector1d[9][0][1][0] = -6.246744747216e-01;\n '\
            '   fineGridProjector1d[9][0][1][1] = 4.762816438513e-01;\n '\
            '   fineGridProjector1d[9][0][1][2] = 1.042821275709e+00;\n '\
            '   fineGridProjector1d[9][0][1][3] = 7.071337534415e-01;\n '\
            '   fineGridProjector1d[9][0][1][4] = 1.400222542117e-01;\n '\
            '   fineGridProjector1d[9][0][1][5] = -1.164488255495e-01;\n '\
            '   fineGridProjector1d[9][0][1][6] = -9.316229264276e-02;\n '\
            '   fineGridProjector1d[9][0][1][7] = -6.656227420708e-03;\n '\
            '   fineGridProjector1d[9][0][1][8] = 4.267480429475e-02;\n '\
            '   fineGridProjector1d[9][0][1][9] = 5.749215321345e-02;\n '\
            '   fineGridProjector1d[9][0][2][0] = 4.477557585207e-01;\n '\
            '   fineGridProjector1d[9][0][2][1] = -2.752990080563e-01;\n '\
            '   fineGridProjector1d[9][0][2][2] = -2.425724291088e-01;\n '\
            '   fineGridProjector1d[9][0][2][3] = 5.127257837626e-01;\n '\
            '   fineGridProjector1d[9][0][2][4] = 1.000240920768e+00;\n '\
            '   fineGridProjector1d[9][0][2][5] = 8.200954832639e-01;\n '\
            '   fineGridProjector1d[9][0][2][6] = 3.598190773165e-01;\n '\
            '   fineGridProjector1d[9][0][2][7] = 2.093610904907e-02;\n '\
            '   fineGridProjector1d[9][0][2][8] = -1.221718631285e-01;\n '\
            '   fineGridProjector1d[9][0][2][9] = -1.578408210702e-01;\n '\
            '   fineGridProjector1d[9][0][3][0] = -3.408299499285e-01;\n '\
            '   fineGridProjector1d[9][0][3][1] = 1.980597343724e-01;\n '\
            '   fineGridProjector1d[9][0][3][2] = 1.535467069462e-01;\n '\
            '   fineGridProjector1d[9][0][3][3] = -2.434492941669e-01;\n '\
            '   fineGridProjector1d[9][0][3][4] = -1.775594896897e-01;\n '\
            '   fineGridProjector1d[9][0][3][5] = 3.792228598505e-01;\n '\
            '   fineGridProjector1d[9][0][3][6] = 8.672994080226e-01;\n '\
            '   fineGridProjector1d[9][0][3][7] = 1.002614711422e+00;\n '\
            '   fineGridProjector1d[9][0][3][8] = 9.093105755987e-01;\n '\
            '   fineGridProjector1d[9][0][3][9] = 7.936216710458e-01;\n '\
            '   fineGridProjector1d[9][0][4][0] = 2.594668864628e-01;\n '\
            '   fineGridProjector1d[9][0][4][1] = -1.473181745929e-01;\n '\
            '   fineGridProjector1d[9][0][4][2] = -1.090292863919e-01;\n '\
            '   fineGridProjector1d[9][0][4][3] = 1.596194319118e-01;\n '\
            '   fineGridProjector1d[9][0][4][4] = 1.017615252569e-01;\n '\
            '   fineGridProjector1d[9][0][4][5] = -1.709980271030e-01;\n '\
            '   fineGridProjector1d[9][0][4][6] = -2.371579450852e-01;\n '\
            '   fineGridProjector1d[9][0][4][7] = -2.690758442992e-02;\n '\
            '   fineGridProjector1d[9][0][4][8] = 2.509466614880e-01;\n '\
            '   fineGridProjector1d[9][0][4][9] = 4.315140742818e-01;\n '\
            '   fineGridProjector1d[9][0][5][0] = -1.917090249451e-01;\n '\
            '   fineGridProjector1d[9][0][5][1] = 1.075827166023e-01;\n '\
            '   fineGridProjector1d[9][0][5][2] = 7.787478784090e-02;\n '\
            '   fineGridProjector1d[9][0][5][3] = -1.101130037797e-01;\n '\
            '   fineGridProjector1d[9][0][5][4] = -6.674007604191e-02;\n '\
            '   fineGridProjector1d[9][0][5][5] = 1.045228242326e-01;\n '\
            '   fineGridProjector1d[9][0][5][6] = 1.319336815931e-01;\n '\
            '   fineGridProjector1d[9][0][5][7] = 1.330702412862e-02;\n '\
            '   fineGridProjector1d[9][0][5][8] = -1.092149790593e-01;\n '\
            '   fineGridProjector1d[9][0][5][9] = -1.697880472775e-01;\n '\
            '   fineGridProjector1d[9][0][6][0] = 1.334678633222e-01;\n '\
            '   fineGridProjector1d[9][0][6][1] = -7.441070596041e-02;\n '\
            '   fineGridProjector1d[9][0][6][2] = -5.321525660320e-02;\n '\
            '   fineGridProjector1d[9][0][6][3] = 7.389124278872e-02;\n '\
            '   fineGridProjector1d[9][0][6][4] = 4.369092251552e-02;\n '\
            '   fineGridProjector1d[9][0][6][5] = -6.629900659563e-02;\n '\
            '   fineGridProjector1d[9][0][6][6] = -8.060041978630e-02;\n '\
            '   fineGridProjector1d[9][0][6][7] = -7.805930664592e-03;\n '\
            '   fineGridProjector1d[9][0][6][8] = 6.170660535874e-02;\n '\
            '   fineGridProjector1d[9][0][6][9] = 9.350806765124e-02;\n '\
            '   fineGridProjector1d[9][0][7][0] = -8.358815437639e-02;\n '\
            '   fineGridProjector1d[9][0][7][1] = 4.642325174487e-02;\n '\
            '   fineGridProjector1d[9][0][7][2] = 3.296836054964e-02;\n '\
            '   fineGridProjector1d[9][0][7][3] = -4.531064015366e-02;\n '\
            '   fineGridProjector1d[9][0][7][4] = -2.643170940764e-02;\n '\
            '   fineGridProjector1d[9][0][7][5] = 3.945188500553e-02;\n '\
            '   fineGridProjector1d[9][0][7][6] = 4.707550489088e-02;\n '\
            '   fineGridProjector1d[9][0][7][7] = 4.473165280534e-03;\n '\
            '   fineGridProjector1d[9][0][7][8] = -3.477816740467e-02;\n '\
            '   fineGridProjector1d[9][0][7][9] = -5.213428845947e-02;\n '\
            '   fineGridProjector1d[9][0][8][0] = 4.247986734444e-02;\n '\
            '   fineGridProjector1d[9][0][8][1] = -2.354024109702e-02;\n '\
            '   fineGridProjector1d[9][0][8][2] = -1.665071452305e-02;\n '\
            '   fineGridProjector1d[9][0][8][3] = 2.275202396951e-02;\n '\
            '   fineGridProjector1d[9][0][8][4] = 1.317312185792e-02;\n '\
            '   fineGridProjector1d[9][0][8][5] = -1.948701682833e-02;\n '\
            '   fineGridProjector1d[9][0][8][6] = -2.302501434818e-02;\n '\
            '   fineGridProjector1d[9][0][8][7] = -2.166624302403e-03;\n '\
            '   fineGridProjector1d[9][0][8][8] = 1.670614308005e-02;\n '\
            '   fineGridProjector1d[9][0][8][9] = 2.491126198622e-02;\n '\
            '   fineGridProjector1d[9][0][9][0] = -1.212481578163e-02;\n '\
            '   fineGridProjector1d[9][0][9][1] = 6.711563352020e-03;\n '\
            '   fineGridProjector1d[9][0][9][2] = 4.737879037219e-03;\n '\
            '   fineGridProjector1d[9][0][9][3] = -6.455564911144e-03;\n '\
            '   fineGridProjector1d[9][0][9][4] = -3.724068275295e-03;\n '\
            '   fineGridProjector1d[9][0][9][5] = 5.485364430402e-03;\n '\
            '   fineGridProjector1d[9][0][9][6] = 6.451096255893e-03;\n '\
            '   fineGridProjector1d[9][0][9][7] = 6.042779234051e-04;\n '\
            '   fineGridProjector1d[9][0][9][8] = -4.641611170299e-03;\n '\
            '   fineGridProjector1d[9][0][9][9] = -6.904605374294e-03;\n '\
            '   fineGridProjector1d[9][1][0][0] = -1.525548037785e-02;\n '\
            '   fineGridProjector1d[9][1][0][1] = -1.516633714233e-02;\n '\
            '   fineGridProjector1d[9][1][0][2] = -1.020722785489e-02;\n '\
            '   fineGridProjector1d[9][1][0][3] = 5.788916826268e-04;\n '\
            '   fineGridProjector1d[9][1][0][4] = 9.570005354179e-03;\n '\
            '   fineGridProjector1d[9][1][0][5] = 8.642023174724e-03;\n '\
            '   fineGridProjector1d[9][1][0][6] = 4.293355628099e-04;\n '\
            '   fineGridProjector1d[9][1][0][7] = -6.355720189296e-03;\n '\
            '   fineGridProjector1d[9][1][0][8] = -8.237055567888e-03;\n '\
            '   fineGridProjector1d[9][1][0][9] = -7.627740188925e-03;\n '\
            '   fineGridProjector1d[9][1][1][0] = 6.065645165005e-02;\n '\
            '   fineGridProjector1d[9][1][1][1] = 5.966605371946e-02;\n '\
            '   fineGridProjector1d[9][1][1][2] = 3.953853135219e-02;\n '\
            '   fineGridProjector1d[9][1][1][3] = -2.205223668492e-03;\n '\
            '   fineGridProjector1d[9][1][1][4] = -3.589947702920e-02;\n '\
            '   fineGridProjector1d[9][1][1][5] = -3.200415238853e-02;\n '\
            '   fineGridProjector1d[9][1][1][6] = -1.574081650028e-03;\n '\
            '   fineGridProjector1d[9][1][1][7] = 2.313172793333e-02;\n '\
            '   fineGridProjector1d[9][1][1][8] = 2.983302685973e-02;\n '\
            '   fineGridProjector1d[9][1][1][9] = 2.755355732908e-02;\n '\
            '   fineGridProjector1d[9][1][2][0] = -1.636298313657e-01;\n '\
            '   fineGridProjector1d[9][1][2][1] = -1.558280441450e-01;\n '\
            '   fineGridProjector1d[9][1][2][2] = -9.871979522477e-02;\n '\
            '   fineGridProjector1d[9][1][2][3] = 5.260618101552e-03;\n '\
            '   fineGridProjector1d[9][1][2][4] = 8.231648646382e-02;\n '\
            '   fineGridProjector1d[9][1][2][5] = 7.111005287223e-02;\n '\
            '   fineGridProjector1d[9][1][2][6] = 3.415735899292e-03;\n '\
            '   fineGridProjector1d[9][1][2][7] = -4.935989761239e-02;\n '\
            '   fineGridProjector1d[9][1][2][8] = -6.296716127904e-02;\n '\
            '   fineGridProjector1d[9][1][2][9] = -5.781774214778e-02;\n '\
            '   fineGridProjector1d[9][1][3][0] = 7.267729103083e-01;\n '\
            '   fineGridProjector1d[9][1][3][1] = 5.720659732597e-01;\n '\
            '   fineGridProjector1d[9][1][3][2] = 2.942277850796e-01;\n '\
            '   fineGridProjector1d[9][1][3][3] = -1.326189386383e-02;\n '\
            '   fineGridProjector1d[9][1][3][4] = -1.839333175888e-01;\n '\
            '   fineGridProjector1d[9][1][3][5] = -1.461392425051e-01;\n '\
            '   fineGridProjector1d[9][1][3][6] = -6.630946931915e-03;\n '\
            '   fineGridProjector1d[9][1][3][7] = 9.226617511639e-02;\n '\
            '   fineGridProjector1d[9][1][3][8] = 1.149608631284e-01;\n '\
            '   fineGridProjector1d[9][1][3][9] = 1.042750862710e-01;\n '\
            '   fineGridProjector1d[9][1][4][0] = 5.169633930643e-01;\n '\
            '   fineGridProjector1d[9][1][4][1] = 6.838166551461e-01;\n '\
            '   fineGridProjector1d[9][1][4][2] = 9.019374125030e-01;\n '\
            '   fineGridProjector1d[9][1][4][3] = 9.989719425882e-01;\n '\
            '   fineGridProjector1d[9][1][4][4] = 8.175584144311e-01;\n '\
            '   fineGridProjector1d[9][1][4][5] = 4.087792072156e-01;\n '\
            '   fineGridProjector1d[9][1][4][6] = 1.501562227981e-02;\n '\
            '   fineGridProjector1d[9][1][4][7] = -1.864589911031e-01;\n '\
            '   fineGridProjector1d[9][1][4][8] = -2.181439739791e-01;\n '\
            '   fineGridProjector1d[9][1][4][9] = -1.918906045426e-01;\n '\
            '   fineGridProjector1d[9][1][5][0] = -1.918906045426e-01;\n '\
            '   fineGridProjector1d[9][1][5][1] = -2.181439739791e-01;\n '\
            '   fineGridProjector1d[9][1][5][2] = -1.864589911031e-01;\n '\
            '   fineGridProjector1d[9][1][5][3] = 1.501562227981e-02;\n '\
            '   fineGridProjector1d[9][1][5][4] = 4.087792072156e-01;\n '\
            '   fineGridProjector1d[9][1][5][5] = 8.175584144311e-01;\n '\
            '   fineGridProjector1d[9][1][5][6] = 9.989719425882e-01;\n '\
            '   fineGridProjector1d[9][1][5][7] = 9.019374125030e-01;\n '\
            '   fineGridProjector1d[9][1][5][8] = 6.838166551461e-01;\n '\
            '   fineGridProjector1d[9][1][5][9] = 5.169633930643e-01;\n '\
            '   fineGridProjector1d[9][1][6][0] = 1.042750862710e-01;\n '\
            '   fineGridProjector1d[9][1][6][1] = 1.149608631284e-01;\n '\
            '   fineGridProjector1d[9][1][6][2] = 9.226617511639e-02;\n '\
            '   fineGridProjector1d[9][1][6][3] = -6.630946931915e-03;\n '\
            '   fineGridProjector1d[9][1][6][4] = -1.461392425051e-01;\n '\
            '   fineGridProjector1d[9][1][6][5] = -1.839333175888e-01;\n '\
            '   fineGridProjector1d[9][1][6][6] = -1.326189386383e-02;\n '\
            '   fineGridProjector1d[9][1][6][7] = 2.942277850796e-01;\n '\
            '   fineGridProjector1d[9][1][6][8] = 5.720659732597e-01;\n '\
            '   fineGridProjector1d[9][1][6][9] = 7.267729103083e-01;\n '\
            '   fineGridProjector1d[9][1][7][0] = -5.781774214778e-02;\n '\
            '   fineGridProjector1d[9][1][7][1] = -6.296716127904e-02;\n '\
            '   fineGridProjector1d[9][1][7][2] = -4.935989761239e-02;\n '\
            '   fineGridProjector1d[9][1][7][3] = 3.415735899292e-03;\n '\
            '   fineGridProjector1d[9][1][7][4] = 7.111005287223e-02;\n '\
            '   fineGridProjector1d[9][1][7][5] = 8.231648646382e-02;\n '\
            '   fineGridProjector1d[9][1][7][6] = 5.260618101552e-03;\n '\
            '   fineGridProjector1d[9][1][7][7] = -9.871979522477e-02;\n '\
            '   fineGridProjector1d[9][1][7][8] = -1.558280441450e-01;\n '\
            '   fineGridProjector1d[9][1][7][9] = -1.636298313657e-01;\n '\
            '   fineGridProjector1d[9][1][8][0] = 2.755355732908e-02;\n '\
            '   fineGridProjector1d[9][1][8][1] = 2.983302685973e-02;\n '\
            '   fineGridProjector1d[9][1][8][2] = 2.313172793333e-02;\n '\
            '   fineGridProjector1d[9][1][8][3] = -1.574081650028e-03;\n '\
            '   fineGridProjector1d[9][1][8][4] = -3.200415238853e-02;\n '\
            '   fineGridProjector1d[9][1][8][5] = -3.589947702920e-02;\n '\
            '   fineGridProjector1d[9][1][8][6] = -2.205223668492e-03;\n '\
            '   fineGridProjector1d[9][1][8][7] = 3.953853135219e-02;\n '\
            '   fineGridProjector1d[9][1][8][8] = 5.966605371946e-02;\n '\
            '   fineGridProjector1d[9][1][8][9] = 6.065645165005e-02;\n '\
            '   fineGridProjector1d[9][1][9][0] = -7.627740188925e-03;\n '\
            '   fineGridProjector1d[9][1][9][1] = -8.237055567888e-03;\n '\
            '   fineGridProjector1d[9][1][9][2] = -6.355720189296e-03;\n '\
            '   fineGridProjector1d[9][1][9][3] = 4.293355628099e-04;\n '\
            '   fineGridProjector1d[9][1][9][4] = 8.642023174724e-03;\n '\
            '   fineGridProjector1d[9][1][9][5] = 9.570005354179e-03;\n '\
            '   fineGridProjector1d[9][1][9][6] = 5.788916826268e-04;\n '\
            '   fineGridProjector1d[9][1][9][7] = -1.020722785489e-02;\n '\
            '   fineGridProjector1d[9][1][9][8] = -1.516633714233e-02;\n '\
            '   fineGridProjector1d[9][1][9][9] = -1.525548037785e-02;\n '\
            '   fineGridProjector1d[9][2][0][0] = -6.904605374294e-03;\n '\
            '   fineGridProjector1d[9][2][0][1] = -4.641611170299e-03;\n '\
            '   fineGridProjector1d[9][2][0][2] = 6.042779234051e-04;\n '\
            '   fineGridProjector1d[9][2][0][3] = 6.451096255893e-03;\n '\
            '   fineGridProjector1d[9][2][0][4] = 5.485364430402e-03;\n '\
            '   fineGridProjector1d[9][2][0][5] = -3.724068275295e-03;\n '\
            '   fineGridProjector1d[9][2][0][6] = -6.455564911144e-03;\n '\
            '   fineGridProjector1d[9][2][0][7] = 4.737879037219e-03;\n '\
            '   fineGridProjector1d[9][2][0][8] = 6.711563352020e-03;\n '\
            '   fineGridProjector1d[9][2][0][9] = -1.212481578163e-02;\n '\
            '   fineGridProjector1d[9][2][1][0] = 2.491126198622e-02;\n '\
            '   fineGridProjector1d[9][2][1][1] = 1.670614308005e-02;\n '\
            '   fineGridProjector1d[9][2][1][2] = -2.166624302403e-03;\n '\
            '   fineGridProjector1d[9][2][1][3] = -2.302501434818e-02;\n '\
            '   fineGridProjector1d[9][2][1][4] = -1.948701682833e-02;\n '\
            '   fineGridProjector1d[9][2][1][5] = 1.317312185792e-02;\n '\
            '   fineGridProjector1d[9][2][1][6] = 2.275202396951e-02;\n '\
            '   fineGridProjector1d[9][2][1][7] = -1.665071452305e-02;\n '\
            '   fineGridProjector1d[9][2][1][8] = -2.354024109702e-02;\n '\
            '   fineGridProjector1d[9][2][1][9] = 4.247986734444e-02;\n '\
            '   fineGridProjector1d[9][2][2][0] = -5.213428845947e-02;\n '\
            '   fineGridProjector1d[9][2][2][1] = -3.477816740467e-02;\n '\
            '   fineGridProjector1d[9][2][2][2] = 4.473165280534e-03;\n '\
            '   fineGridProjector1d[9][2][2][3] = 4.707550489088e-02;\n '\
            '   fineGridProjector1d[9][2][2][4] = 3.945188500553e-02;\n '\
            '   fineGridProjector1d[9][2][2][5] = -2.643170940764e-02;\n '\
            '   fineGridProjector1d[9][2][2][6] = -4.531064015366e-02;\n '\
            '   fineGridProjector1d[9][2][2][7] = 3.296836054964e-02;\n '\
            '   fineGridProjector1d[9][2][2][8] = 4.642325174487e-02;\n '\
            '   fineGridProjector1d[9][2][2][9] = -8.358815437639e-02;\n '\
            '   fineGridProjector1d[9][2][3][0] = 9.350806765124e-02;\n '\
            '   fineGridProjector1d[9][2][3][1] = 6.170660535874e-02;\n '\
            '   fineGridProjector1d[9][2][3][2] = -7.805930664592e-03;\n '\
            '   fineGridProjector1d[9][2][3][3] = -8.060041978630e-02;\n '\
            '   fineGridProjector1d[9][2][3][4] = -6.629900659563e-02;\n '\
            '   fineGridProjector1d[9][2][3][5] = 4.369092251552e-02;\n '\
            '   fineGridProjector1d[9][2][3][6] = 7.389124278872e-02;\n '\
            '   fineGridProjector1d[9][2][3][7] = -5.321525660320e-02;\n '\
            '   fineGridProjector1d[9][2][3][8] = -7.441070596041e-02;\n '\
            '   fineGridProjector1d[9][2][3][9] = 1.334678633222e-01;\n '\
            '   fineGridProjector1d[9][2][4][0] = -1.697880472775e-01;\n '\
            '   fineGridProjector1d[9][2][4][1] = -1.092149790593e-01;\n '\
            '   fineGridProjector1d[9][2][4][2] = 1.330702412862e-02;\n '\
            '   fineGridProjector1d[9][2][4][3] = 1.319336815931e-01;\n '\
            '   fineGridProjector1d[9][2][4][4] = 1.045228242326e-01;\n '\
            '   fineGridProjector1d[9][2][4][5] = -6.674007604191e-02;\n '\
            '   fineGridProjector1d[9][2][4][6] = -1.101130037797e-01;\n '\
            '   fineGridProjector1d[9][2][4][7] = 7.787478784090e-02;\n '\
            '   fineGridProjector1d[9][2][4][8] = 1.075827166023e-01;\n '\
            '   fineGridProjector1d[9][2][4][9] = -1.917090249451e-01;\n '\
            '   fineGridProjector1d[9][2][5][0] = 4.315140742818e-01;\n '\
            '   fineGridProjector1d[9][2][5][1] = 2.509466614880e-01;\n '\
            '   fineGridProjector1d[9][2][5][2] = -2.690758442992e-02;\n '\
            '   fineGridProjector1d[9][2][5][3] = -2.371579450852e-01;\n '\
            '   fineGridProjector1d[9][2][5][4] = -1.709980271030e-01;\n '\
            '   fineGridProjector1d[9][2][5][5] = 1.017615252569e-01;\n '\
            '   fineGridProjector1d[9][2][5][6] = 1.596194319118e-01;\n '\
            '   fineGridProjector1d[9][2][5][7] = -1.090292863919e-01;\n '\
            '   fineGridProjector1d[9][2][5][8] = -1.473181745929e-01;\n '\
            '   fineGridProjector1d[9][2][5][9] = 2.594668864628e-01;\n '\
            '   fineGridProjector1d[9][2][6][0] = 7.936216710458e-01;\n '\
            '   fineGridProjector1d[9][2][6][1] = 9.093105755987e-01;\n '\
            '   fineGridProjector1d[9][2][6][2] = 1.002614711422e+00;\n '\
            '   fineGridProjector1d[9][2][6][3] = 8.672994080226e-01;\n '\
            '   fineGridProjector1d[9][2][6][4] = 3.792228598505e-01;\n '\
            '   fineGridProjector1d[9][2][6][5] = -1.775594896897e-01;\n '\
            '   fineGridProjector1d[9][2][6][6] = -2.434492941669e-01;\n '\
            '   fineGridProjector1d[9][2][6][7] = 1.535467069462e-01;\n '\
            '   fineGridProjector1d[9][2][6][8] = 1.980597343724e-01;\n '\
            '   fineGridProjector1d[9][2][6][9] = -3.408299499285e-01;\n '\
            '   fineGridProjector1d[9][2][7][0] = -1.578408210702e-01;\n '\
            '   fineGridProjector1d[9][2][7][1] = -1.221718631285e-01;\n '\
            '   fineGridProjector1d[9][2][7][2] = 2.093610904907e-02;\n '\
            '   fineGridProjector1d[9][2][7][3] = 3.598190773165e-01;\n '\
            '   fineGridProjector1d[9][2][7][4] = 8.200954832639e-01;\n '\
            '   fineGridProjector1d[9][2][7][5] = 1.000240920768e+00;\n '\
            '   fineGridProjector1d[9][2][7][6] = 5.127257837626e-01;\n '\
            '   fineGridProjector1d[9][2][7][7] = -2.425724291088e-01;\n '\
            '   fineGridProjector1d[9][2][7][8] = -2.752990080563e-01;\n '\
            '   fineGridProjector1d[9][2][7][9] = 4.477557585207e-01;\n '\
            '   fineGridProjector1d[9][2][8][0] = 5.749215321345e-02;\n '\
            '   fineGridProjector1d[9][2][8][1] = 4.267480429475e-02;\n '\
            '   fineGridProjector1d[9][2][8][2] = -6.656227420708e-03;\n '\
            '   fineGridProjector1d[9][2][8][3] = -9.316229264276e-02;\n '\
            '   fineGridProjector1d[9][2][8][4] = -1.164488255495e-01;\n '\
            '   fineGridProjector1d[9][2][8][5] = 1.400222542117e-01;\n '\
            '   fineGridProjector1d[9][2][8][6] = 7.071337534415e-01;\n '\
            '   fineGridProjector1d[9][2][8][7] = 1.042821275709e+00;\n '\
            '   fineGridProjector1d[9][2][8][8] = 4.762816438513e-01;\n '\
            '   fineGridProjector1d[9][2][8][9] = -6.246744747216e-01;\n '\
            '   fineGridProjector1d[9][2][9][0] = -1.437946599702e-02;\n '\
            '   fineGridProjector1d[9][2][9][1] = -1.053816905737e-02;\n '\
            '   fineGridProjector1d[9][2][9][2] = 1.601079014457e-03;\n '\
            '   fineGridProjector1d[9][2][9][3] = 2.136690378358e-02;\n '\
            '   fineGridProjector1d[9][2][9][4] = 2.445445929348e-02;\n '\
            '   fineGridProjector1d[9][2][9][5] = -2.443340119521e-02;\n '\
            '   fineGridProjector1d[9][2][9][6] = -7.079373286275e-02;\n '\
            '   fineGridProjector1d[9][2][9][7] = 1.095186765437e-01;\n '\
            '   fineGridProjector1d[9][2][9][8] = 6.855092197838e-01;\n '\
            '   fineGridProjector1d[9][2][9][9] = 1.369756044103e+00;\n ')

        l_sourceFile.write("}\n\n")
        l_sourceFile.close()
