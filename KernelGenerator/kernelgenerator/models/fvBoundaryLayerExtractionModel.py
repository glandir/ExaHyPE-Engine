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
# Generates the ghostLayerFillingAtBoundary Kernels
#


from .abstractModelBaseClass import AbstractModelBaseClass


class FVBoundaryLayerExtractionModel(AbstractModelBaseClass):

    def generateCode(self):
        self.render(("fv", "boundaryLayerExtraction_cpp.template"), "boundaryLayerExtraction.cpp")
