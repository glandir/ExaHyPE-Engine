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
 
#ifndef _EXAHYPE_KERNELS_KERNEL_CALLS_H_
#define _EXAHYPE_KERNELS_KERNEL_CALLS_H_

#include "exahype/parser/Parser.h"
#include <ostream>

namespace kernels {
  /**
   * Is implemented within the application folder generated by the toolkit.
   */
  void registerSolvers(exahype::parser::Parser& parser);

  /**
   * This callback enables the user to deallocate memory.
   */
  void finalise();
  
  /**
   * This allows dumping information what code was generated by the toolkit. 
   */
  void toString(std::ostream& os);

  /**
   * This allows to access the specfile that was seen by the toolkit.
   **/
  const char* compiledSpecfile();

  /**
   * Allows users to implement a custom specification file parser. The
   * output string must hold a valid JSON file.
   **/
  std::string readSpecificationFileToJSON(const std::string& filename);
}

#endif
