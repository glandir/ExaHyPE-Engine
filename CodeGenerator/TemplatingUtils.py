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
# @section DESCRIPTION
#
# Base class for template based generator
#


import os

# Change this boolean to change how jinja2 is imported.
#       True = jinja2 is available directly with your python installation
#       False = jinja2's source code is available locally under CodeGenerator/jinja
isJinjaAvailableAsPackage=True

if(isJinjaAvailableAsPackage): 
    # 1. if jinja2 is available with Python
    from jinja2 import Template
else:
    # 2. if jinja2 sources are locally available
    import sys
    sys.path.insert(0, 'jinja')
    from jinja.jinja2 import Template


def renderAsFile(inputFilename, outputFilename, context):
    dir = os.path.dirname(__file__)
           
    with open(os.path.join(dir,'templates',inputFilename), 'r') as input:
        template = Template(input.read(), trim_blocks=True)                
        with open(outputFilename, 'w') as output:
            output.write(template.render(context))
