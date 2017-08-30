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
# Starting point of the code generator
#
# @note
# requires python3

import argparse
import os
import sys

import CodeGenArgumentParser
import Backend

# --------------------------------------------------------
# Configuration parameters
# --------------------------------------------------------

pathFromHereToExaHyPERoot = "../"                     #path to the root of ExaHyPe from this file
pathToLibxsmmGemmGenerator = "libxsmm_gemm_generator" #path to the gemm generator from this file

# --------------------------------------------------------
# Require python3
# --------------------------------------------------------
requiredVersion = (3,3)
currentVersion  = sys.version_info

if(requiredVersion > currentVersion):
    sys.exit("CodeGenerator: Requires Python 3.3 or newer. Abort.")


# --------------------------------------------------------
# Process the command line arguments
# --------------------------------------------------------
l_parser = argparse.ArgumentParser(description="This is the front end of the ExaHyPE code generator.")

l_parser.add_argument("pathToApplication",
                      help="path to the application as given by the ExaHyPE specification file (application directory as root)")
l_parser.add_argument("pathToOptKernel",
                      help="desired relative path to the generated code (application directory as root)")
l_parser.add_argument("namespace",
                      help="desired namespace for the generated code")                      
l_parser.add_argument("solverName",
                      help="name of the user-solver")
l_parser.add_argument("numberOfVariables", 
                      type=int, 
                      help="the number of quantities")
l_parser.add_argument("order", 
                      type=int, 
                      help="the order of the approximation polynomial")
l_parser.add_argument("dimension", 
                      type=int, 
                      help="number of dimensions you want to simulate")
l_parser.add_argument("numerics",
                      type=lambda numericsArg: CodeGenArgumentParser.validateNumerics(l_parser, numericsArg),
                      help="linear or nonlinear")
l_parser.add_argument("architecture",
                      type=lambda architectureArg: CodeGenArgumentParser.validateArchitecture(l_parser, architectureArg), 
                      help="the microarchitecture of the target device")
l_parser.add_argument("--deepProfiling",
                      action="store_true",
                      help="enable deep-rpofiling (use only with profiler enable)")
l_parser.add_argument("--useFlux",
                      action="store_true",
                      help="enable flux")
l_parser.add_argument("--useNCP",
                      action="store_true",
                      help="enable non conservative product")
l_parser.add_argument("--useSource",
                      action="store_true",
                      help="enable source terms")
l_parser.add_argument("--noTimeAveraging",
                      action="store_true",
                      help="disable time averaging in the spacetimepredictor (less memory usage, more computation)")
# l_parser.add_argument("--precision",
                      # type=lambda precisionArg: CodeGenArgumentParser.validatePrecision(l_parser, precisionArg),
                      # default="DP",
                      # help="SP or DP")


l_commandLineArguments = l_parser.parse_args()

config = { 
           "numerics"              : l_commandLineArguments.numerics,
           "pathToOptKernel"       : l_commandLineArguments.pathToOptKernel,
           "solverName"            : l_commandLineArguments.solverName,
           "nVar"                  : l_commandLineArguments.numberOfVariables,
           "nDof"                  : (l_commandLineArguments.order)+1,
           "nDim"                  : l_commandLineArguments.dimension,
           "nPar"                  : 0, #TODO JMG add paramters ?
           "useDeepProfiler"       : l_commandLineArguments.deepProfiling,
           "useFlux"               : l_commandLineArguments.useFlux,
           "useNCP"                : l_commandLineArguments.useNCP,
           "useSource"             : l_commandLineArguments.useSource,
           "useSourceOrNCP"        : (l_commandLineArguments.useSource or l_commandLineArguments.useNCP),
           "noTimeAveraging"       : l_commandLineArguments.noTimeAveraging,
           "codeNamespace"         : l_commandLineArguments.namespace,
           "pathToOutputDirectory" : os.path.join(os.path.dirname(__file__),pathFromHereToExaHyPERoot,l_commandLineArguments.pathToApplication,l_commandLineArguments.pathToOptKernel),
           "architecture"          : l_commandLineArguments.architecture,
           "precision"             : "DP", #l_commandLineArguments.precision,
           "pathToLibxsmmGemmGenerator"  : os.path.join(os.path.dirname(__file__),pathToLibxsmmGemmGenerator)
          }

# configure global setup of the code generator
Backend.setConfig(config)

# clean up output directory
Backend.prepareOutputDirectory(config['pathToOutputDirectory'])

# generate the compute kernels.
Backend.generateComputeKernels()



