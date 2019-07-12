#!/usr/bin/python
"""
.. module:: generateLookupTable
  :platform: Unix, Windows, Mac
  :synopsis: Generates lookup table initialisation code up to a specific order.
.. moduleauthor:: Angelika Schwarz <angelika.schwarz@tum.de>,Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Generates lookup table initialisation code up to a specific order.
"""
import mpmath as mp
import numpy as np
from os import remove
from numpy import linalg
import sys
import os
import errno

from aderdg import *
from generateBasisFunctionImplementations import *
import gaussLobattoQuadrature  as lob
import gaussLegendreQuadrature as leg

#-----------------------------------------------------------
# main
#-----------------------------------------------------------
maxOrder  = 9
printPrec = 64

outputDirectory = "generatedCode"

# create directory for output files if not existing
try:
    os.makedirs(outputDirectory)
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise


# we process one order after another
# (1) write the #ifdef construct to the output file
# (2) compute the quadrature points and weights
# (3) compute the system matrices and write them to the output file

gaulegImpl   = open("generatedCode/gaussLegendreQuadrature.csnippet", "w")
gaulobImpl   = open("generatedCode/gaussLobattoQuadrature.csnippet", "w")
gaulegHeader = open("generatedCode/gaussLegendreQuadrature.hsnippet", "w")
gaulobHeader = open("generatedCode/gaussLobattoQuadrature.hsnippet", "w")

# header
for rule, outfile in { "legendre" : gaulegHeader, "lobatto" : gaulobHeader }.items():
    def writeVector(label,suffix="\n\n"):
        outfile.write("\n".join(["extern const double %s_%d[%d+1];" % (label,i,i) for i in range(0,maxOrder+1)]))
        outfile.write("\nextern const double* const %s[%d+1];" % (label,maxOrder))
        outfile.write(suffix)
    
    def writeMatrix(label,suffix="\n\n"):
        for r in range(0,maxOrder+1):
            outfile.write("\n".join(["extern const double %s_%d_%d[%d+1];" % (label,r,i,i) for i in range(0,maxOrder+1)]))
            outfile.write("\nextern const double* const %s_%d[%d+1];" % (label,r,maxOrder))
            outfile.write("\n")
        outfile.write("extern const double** const %s[%d+1];" % (label,maxOrder))
        outfile.write(suffix)
    
    outfile.write("namespace kernels {\n")
    outfile.write("namespace gauss%s {\n" % rule)
    writeVector("weights")
    writeVector("nodes")
    writeMatrix("Kxi")
    writeMatrix("MM")
    writeMatrix("iK1")
    writeMatrix("dudx")
    writeMatrix("equidistantGridProjector")

    writeVector("FCoeff_0",suffix="\n")
    writeVector("FCoeff_1",suffix="\n")
    outfile.write("extern const double** const FCoeff[%d+1];\n\n" % (maxOrder))
    outfile.write("\n};");
    outfile.write("\n};");

    outfile.write("\n")    
    for p in range(0,maxOrder+1):
        writeMatrix("fineGridProjector_0_%d" % p,suffix="\n")
        writeMatrix("fineGridProjector_1_%d" % p,suffix="\n")
        writeMatrix("fineGridProjector_2_%d" % p,suffix="\n")
    outfile.write("extern const double*** const fineGridProjector[%d+1];\n\n" % (maxOrder))

    # basis functions    
    outfile.write("\n")
    for p in range(0,maxOrder+1):
        for i in range(0,p+1):
            outfile.write(                "double basisFunction_%d_%d(const double& s);\n" % (p,i))    
            outfile.write( "double basisFunctionFirstDerivative_%d_%d(const double& s);\n" % (p,i))    
            outfile.write("double basisFunctionSecondDerivative_%d_%d(const double& s);\n" % (p,i))    
    outfile.write("typedef double (* UnivariateFunction) (const double s);")
    for p in range(0,maxOrder+1):
        outfile.write("extern const UnivariateFunction basisFunctions_%d[%d+1];\n"                 % (p,p))
        outfile.write("extern const UnivariateFunction basisFunctionFirstDerivatives_%d[%d+1];\n"  % (p,p))
        outfile.write("extern const UnivariateFunction basisFunctionSecondDerivatives_%d[%d+1];\n" % (p,p))
    outfile.write("extern const UnivariateFunction* basisFunctions[%d+1];\n"                 % (maxOrder))
    outfile.write("extern const UnivariateFunction* basisFunctionFirstDerivatives[%d+1];\n"  % (maxOrder))
    outfile.write("extern const UnivariateFunction* basisFunctionSecondDerivatives[%d+1];\n" % (maxOrder))

# source
def quadNodes(order,rule):
    if rule is "legendre": 
        return leg.gauleg(order+1)
    elif rule is "lobatto":
        if(order+1 == 1):
            return [0.5], [1]
        elif(order+1 > 1):
            return lob.gaulob(0., 1., order+1)

for rule, outfile in { "legendre" : gaulegImpl, "lobatto" : gaulobImpl }.items():
    for order in range(0,maxOrder+1):    
        outfile.write("// order %d\n" % order)
        
        def writeVector(vector,label):
            outfile.write("const double kernels::gauss::%s::%s_%d[%d]={\n  %s\n};\n" % (rule,label,order,order+1,",\n  ".join([mp.nstr(i,printPrec) for i in vector])))
        
        def writeMatrix(matrix,label):
            for r in range(0,order+1):
                outfile.write("const double kernels::gauss::%s::%s_%d_%d[%d+1]={\n  %s\n};\n"   % (rule,label,order,r,order,
                    ",\n  ".join([mp.nstr(i,printPrec) for i in matrix[r]])))
            outfile.write("const double* const kernels::gauss::%s::%s_%d[%d+1]={\n  %s\n};\n"  % (rule,label,order,order,
                ",\n  ".join(["%s_%d_%d" % (label,order,i) for i in range(0,order+1)])))

        # quad nodes
        x,w = quadNodes(order,rule)
        writeVector(w,"weights")
        writeVector(x,"nodes")

        # reference stiffness matrix
        Kxi = assembleStiffnessMatrix(x, w, order)
        writeMatrix(Kxi,"Kxi")

        # mass matrix
        MM = assembleMassMatrix(x, w, order)
        writeMatrix(MM,"MM")

        # left-hand side matrix appearing in predictor computation 
        K1 = assembleK1(Kxi, x, order) 
        iK1 = mp.inverse(K1).tolist()
        writeMatrix(iK1,"iK1")

        # derivative Matrix
        dudx = assembleDiscreteDerivativeOperator(MM, Kxi)
        writeMatrix(dudx,"dudx")

        # equidistant grid projectors
        equidistantGridProjector1d = assembleEquidistantGridProjector1d(x, order)
        writeMatrix(equidistantGridProjector1d,"equidistantGridProjector1d")

        # lifting operators
        FLCoeff, _ = BaseFunc1d(mp.mpf(0.0), x, order) 
        FRCoeff, _ = BaseFunc1d(mp.mpf(1.0), x, order)
        writeVector(FLCoeff,"FCoeff_0")
        writeVector(FLCoeff,"FCoeff_1")
        outfile.write("const double** kernels::gauss::%s::FCoeff_%d[2]={FCoeff_0_%d,FCoeff_1_%d};\n" % (rule,order,order,order))

        # fine grid projectors
        fineGridProjector1d0 = assembleFineGridProjector1d(x, 0, order)
        fineGridProjector1d1 = assembleFineGridProjector1d(x, 1, order)
        fineGridProjector1d2 = assembleFineGridProjector1d(x, 2, order)
        writeMatrix(fineGridProjector1d0,"fineGridProjector_0")
        writeMatrix(fineGridProjector1d1,"fineGridProjector_1")
        writeMatrix(fineGridProjector1d2,"fineGridProjector_2")
        outfile.write("const double*** kernels::gauss::%s::fineGridProjector_%d[3]={fineGridProjector_0_%d,fineGridProjector_1_%d,fineGridProjector_2_%d};\n" % (rule,order,order,order,order))

        # basis functions
        outfile.write("\n")
        generateBasisFunctionDefinitions(outfile,rule,order,x,printPrec)

        def writeBasisFunctionArray(label):
            outfile.write("const UnivariateFunction* kernels::gauss::%ss::%s_%d[%d+1]={\n  %s\n};\n" % (rule,label,order,order+1,
                ",\n  ".join(["%s_%d_%d" % (label,order,i) for i in range(0,order+1)])))

        writeBasisFunctionArray("basisFunction")
        writeBasisFunctionArray("basisFunctionFirstDerivative")
        writeBasisFunctionArray("basisFunctionSecondDerivative")

        print("Done with operators for order "+ str(order))

    # collect the arrays per order
    outfile.write("\n")
    # weights,nodes
    for array in ["weights","nodes"]:
        outfile.write("const double** const kernels::gauss%s::%s[%d+1] = {\n  %s\n};\n" % (rule,array,maxOrder,
        ",\n  ".join(["%s_%d" % (array,i) for i in range(0,maxOrder+1)])))
    # Kxi,MM,iK1,dudx,equidistantGridProjector1d
    for array in ["Kxi","MM","iK1","dudx","equidistantGridProjector1d"]:
        outfile.write("const double*** const kernels::gauss%s::%s[%d+1] = {\n  %s\n};\n" % (rule,array,maxOrder,
            ",\n  ".join(["%s_%d" % (array,i) for i in range(0,maxOrder+1)])))
    
    # FCoeff
    array="FCoeff"
    outfile.write("const double*** const kernels::gauss%s::%s[%d+1] = {\n  %s\n};\n" % (rule,array,maxOrder,
        ",\n  ".join(["%s_%d" % (array,i) for i  in range(0,maxOrder+1)])))
    
    # fineGridProjector1d
    array="fineGridProjector"
    outfile.write("const double**** const kernels::gauss%s::%s[%d+1] = {\n  %s\n};\n" % (rule,array,maxOrder,
        ",\n  ".join(["%s_%d" % (array,i) for i in range(0,maxOrder+1)])))

    # basis functions
    def writeBasisFunctionArray(label):
        outfile.write("const UnivariateFunction** kernels::gauss::%ss::%s[%d+1]={\n  %s\n};\n" % (rule,label,order+1,
            ",\n  ".join(["%s_%d" % (label,order,i) for i in range(0,order+1)])))

        writeBasisFunctionArray("basisFunction")
        writeBasisFunctionArray("basisFunctionFirstDerivative")
        writeBasisFunctionArray("basisFunctionSecondDerivative")
    

gaulegHeader.close()
gaulegImpl.close()
gaulobHeader.close()
gaulobImpl.close()

#out.close()
