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


import Backend
import TemplatingUtils
import Utils #matrix operation and build functions
import time #for runtime measurments


class DGMatrixGenerator:
    m_context = {}

    # name of generated output file
    m_filenameRoot = "DGMatrices"
    
    # quadrature nodes and weights mapped onto [0,1]
    m_xGPN       = []
    m_wGPN       = []
    

    def __init__(self, i_context):
        self.m_context = i_context
        self.m_wGPN, self.m_xGPN = Utils.getGaussLegendre(self.m_context['nDof'])
        

    def generateCode(self):
        l_padSize = self.m_context['nDofPad'] - self.m_context['nDof']
        self.m_context['nDofPad_seq'] = range(self.m_context['nDofPad'])
        self.m_context['nDofPadTimesnDof_seq'] = range(self.m_context['nDofPad']*self.m_context['nDof'])
        
        start = time.perf_counter()
        # [FLCoeff 0...0]; [FRCoeff 0...0];
        # right now FLCoeff, FRCoeff no pad (gives no benefit w.r.t libxsmm)
        FLCoeff, _ = Utils.BaseFunc1d(0.0, self.m_xGPN, self.m_context['nDof']) #is also F0
        FRCoeff, _ = Utils.BaseFunc1d(1.0, self.m_xGPN, self.m_context['nDof'])
        self.m_context['FLCoeff'] = Utils.vectorPad(FLCoeff, l_padSize)
        self.m_context['FRCoeff'] = Utils.vectorPad(FRCoeff, l_padSize)
        if self.m_context['runtimeDebug']:
            print("FRLCoeff: "+str(time.perf_counter()-start))
        
        # Matrices are stored in column major order (so the padding should be on the bottom rows)
        # [ Mat ]
        # [0...0]
        # [0...0]
        start = time.perf_counter()
        # Kxi
        Kxi = Utils.assembleStiffnessMatrix(self.m_xGPN, self.m_wGPN, self.m_context['nDof'])
        self.m_context['Kxi']   = Utils.matrixPadAndFlatten_ColMajor(Kxi,l_padSize)
        self.m_context['Kxi_T'] = Utils.matrixPadAndFlatten_RowMajor(Kxi,l_padSize) #transpose
        if self.m_context['runtimeDebug']:
            print("Kxi: "+str(time.perf_counter()-start))
        
        start = time.perf_counter()
        # iK1_T
        iK1 = Utils.matrixInverse(Utils.assembleK1(Kxi, self.m_xGPN, self.m_context['nDof']))
        self.m_context['iK1_T'] = Utils.matrixPadAndFlatten_RowMajor(iK1,l_padSize) #transpose
        if self.m_context['runtimeDebug']:
            print("iK1_T: "+str(time.perf_counter()-start))
        
        start = time.perf_counter()
        # dudx
        MM   = Utils.assembleMassMatrix(self.m_xGPN, self.m_wGPN, self.m_context['nDof'])
        dudx = Utils.assembleDiscreteDerivativeOperator(MM,Kxi)
        self.m_context['dudx']   = Utils.matrixPadAndFlatten_ColMajor(dudx,l_padSize)
        self.m_context['dudx_T'] = Utils.matrixPadAndFlatten_RowMajor(dudx,l_padSize) #transpose
        if self.m_context['runtimeDebug']:
            print("dudx: "+str(time.perf_counter()-start))
        
        start = time.perf_counter()
        #fineGridProjector1d
        fineGridProjector1d_0 = Utils.assembleFineGridProjector1d(self.m_xGPN, 0, self.m_context['nDof'])
        fineGridProjector1d_1 = Utils.assembleFineGridProjector1d(self.m_xGPN, 1, self.m_context['nDof'])
        fineGridProjector1d_2 = Utils.assembleFineGridProjector1d(self.m_xGPN, 2, self.m_context['nDof'])
        self.m_context['fineGridProjector1d_0']   = Utils.matrixPadAndFlatten_ColMajor(fineGridProjector1d_0,l_padSize)
        self.m_context['fineGridProjector1d_1']   = Utils.matrixPadAndFlatten_ColMajor(fineGridProjector1d_1,l_padSize)
        self.m_context['fineGridProjector1d_2']   = Utils.matrixPadAndFlatten_ColMajor(fineGridProjector1d_2,l_padSize)
        
        #fineGridProjector1d_T_weighted
        for i in range(self.m_context['nDof']):
            for j in range(self.m_context['nDof']):
                fineGridProjector1d_0[i][j] *= self.m_wGPN[i]/self.m_wGPN[j]/3.0
                fineGridProjector1d_1[i][j] *= self.m_wGPN[i]/self.m_wGPN[j]/3.0
                fineGridProjector1d_2[i][j] *= self.m_wGPN[i]/self.m_wGPN[j]/3.0
        self.m_context['fineGridProjector1d_T_weighted_0'] = Utils.matrixPadAndFlatten_RowMajor(fineGridProjector1d_0,l_padSize)
        self.m_context['fineGridProjector1d_T_weighted_1'] = Utils.matrixPadAndFlatten_RowMajor(fineGridProjector1d_1,l_padSize)
        self.m_context['fineGridProjector1d_T_weighted_2'] = Utils.matrixPadAndFlatten_RowMajor(fineGridProjector1d_2,l_padSize)
        if self.m_context['runtimeDebug']:
            print("fineGrid: "+str(time.perf_counter()-start))
        
        #generate files 
        TemplatingUtils.renderAsFile('DGMatrices_h.template',   self.m_filenameRoot+'.h',   self.m_context)
        TemplatingUtils.renderAsFile('DGMatrices_cpp.template', self.m_filenameRoot+'.cpp', self.m_context)
        
