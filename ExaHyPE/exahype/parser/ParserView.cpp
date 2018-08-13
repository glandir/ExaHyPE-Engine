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
#include "exahype/parser/ParserView.h"

#include "tarch/Assertions.h"

#include <fstream>

#include <stdio.h>
#include <string.h>
#include <string>
#include <regex>
#include <cstdlib> // getenv, exit
#include <sstream>

#include "tarch/la/ScalarOperations.h"

#include "exahype/parser/Parser.h"

tarch::logging::Log exahype::parser::ParserView::_log( "exahype::parser::ParserView" );

exahype::parser::ParserView::ParserView(const exahype::parser::Parser* parser,
                                        std::string basePath)
    : _parser(parser),
      _basePath(basePath) {}
      
std::string exahype::parser::ParserView::getPath(const std::string& key) const {
  return _basePath + "/" + key;
}

bool exahype::parser::ParserView::hasKey(const std::string& key) const {
  if (_parser!=nullptr) {
    assertion(getParser().isValid());
    return getParser().hasPath(getPath(key));
  } else {
    logError("hasKey()", "No parameters found at all!");
    std::abort();
    return 0;
  }
}

int exahype::parser::ParserView::getValueAsInt(const std::string& key) const {
  if (_parser!=nullptr) {
    return getParser().getIntFromPath(getPath(key));
  } else {
    logError("getValueAsInt()", "No parameters found at all!");
    std::abort();
    return 0;
  }
}

bool exahype::parser::ParserView::getValueAsBool(const std::string& key) const {
  if (_parser!=nullptr) {
    return getParser().getBoolFromPath(getPath(key));
  } else {
    logError("getValueAsBool()", "No parameters found at all!");
    std::abort();
    return false;
  }
}

double exahype::parser::ParserView::getValueAsDouble(
    const std::string& key) const {
  if (_parser!=nullptr) {
    return getParser().getDoubleFromPath(getPath(key));
  } else {
    logError("getValueAsDouble()", "No parameters found at all!");
    std::abort();
    return 0.0;
  }
}

std::string exahype::parser::ParserView::getValueAsString(
    const std::string& key) const {
  if (_parser!=nullptr) {
    return getParser().getStringFromPath(getPath(key));
  } else {
    logError("getValueAsString()", "No parameters found at all!");
    std::abort();
    return "";
  }
}

bool exahype::parser::ParserView::isValueValidInt(
    const std::string& key) const {
  if (_parser!=nullptr) {
    return getParser().isValueValidInt(getPath(key));
  } else {
    logError("isValueValidInt()", "No parameters found at all!");
    std::abort();
    return false;
  }
}

bool exahype::parser::ParserView::isValueValidDouble(
    const std::string& key) const {
  if (_parser!=nullptr) {
    return getParser().isValueValidDouble(getPath(key));
  } else {
    logError("isValueValidDouble()", "No parameters found at all!");
    std::abort();
    return false;
  }
}


bool exahype::parser::ParserView::isValueValidBool(
    const std::string& key) const {
  if (_parser!=nullptr) {
    return getParser().isValueValidBool(getPath(key));
  } else {
    logError("isValueValidBool()", "No parameters found at all!");
    std::abort();
    return false;
  }
}

std::vector< std::pair<std::string, std::string> > exahype::parser::ParserView::getAllAsOrderedMap() const {
  std::vector< std::pair<std::string, std::string> > retvec;

  // Well, this makes no sense in the moment and therefore I won't implement it.
  // It does not make sense anymore because constants may be any kind of nested
  // data.
  
  logError("getAllAsOrderedMap()", "Not yet implemented");
  return retvec;
}

std::map<std::string, std::string> exahype::parser::ParserView::getAllAsMap() const {
  std::map<std::string, std::string>  retmap;
	
  // Well, this makes no sense in the moment and therefore I won't implement it.
  // It does not make sense anymore because constants may be any kind of nested
  // data.
  
  logError("getAllAsMap()", "Not yet implemented");
  return retmap;
}

bool exahype::parser::ParserView::isValueValidString(
    const std::string& key) const {
  if (_parser!=nullptr) {
    return getParser().isValueValidString(getPath(key));
  } else {
    logError("isValueValidString()", "No parameters found at all!");
    std::abort();
    return false;
  }
}

const exahype::parser::Parser& exahype::parser::ParserView::getParser() const {
  if(!_parser) {
    logError("getParser()", "Trying to use a non-initialised parserView!");
    std::abort();
  }
  return *_parser;
}


std::string exahype::parser::ParserView::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::parser::ParserView::toString(std::ostream& out) const {
  if (_parser!=nullptr) {
    out << "ParserView(";
    out << "specfile:" << getParser().getSpecfileName();
    out << ",";
    out << "basePath:" <<  _basePath;
    out <<  ")";
  } {
    out << "not available";
  }
}

std::string exahype::parser::ParserView::dump(const std::string path) const {
   if (_parser!=nullptr) {
     return getParser().dumpPath( _basePath +  ( path.compare("")==0 ? "" : "/"+path ) );
   } else {
     return "not available";
   }
}

