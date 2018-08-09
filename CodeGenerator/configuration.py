import sys
import os

class Configuration:

    ######################################
    ###### Configuration parameters ######
    ######################################
    # Change them if required

    # path to the root of ExaHyPe from this file
    pathToExaHyPERoot          = os.path.join("..")

    # path to the gemm generator from this file
    pathToLibxsmmGemmGenerator = os.path.join("dependencies", "libxsmm", "bin", "libxsmm_gemm_generator")
    
    # path to jinja2
    pathToJinja2               = os.path.join("dependencies", "jinja")
    
    # path to markupsafe
    pathToMarkupsafe           = os.path.join("dependencies", "markupsafe")
    
    # simd size of the accepted architectures
    simdWidth = { "noarch" : 1,
                  "wsm"    : 2,
                  "snb"    : 4,
                  "hsw"    : 4,
                  "knc"    : 8,
                  "knl"    : 8 
                }
