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
 
#include "tarch/logging/Log.h"
#include "tarch/tests/TestCaseRegistry.h"
#include "tarch/logging/CommandLineLogger.h"
#include "tarch/logging/LogFilterFileReader.h"
#include "tarch/parallel/Node.h"

#include "peano/peano.h"

#include "exahype/Parser.h"
#include "exahype/Vertex.h"
#include "exahype/runners/Runner.h"
#include "buildinfo.h"

#include "kernels/KernelCalls.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

#include <vector>
#include <string>
#include <cstdlib> // getenv, exit
#include <iostream>
#include <cstdio>

tarch::logging::Log _log("");

void version(const std::string& programname); // version dumping, see below
void help(const std::string& programname);  // A help message

int main(int argc, char** argv) {
  peano::fillLookupTables();

  //
  //   Setup environment
  // =====================
  //
  int parallelSetup = peano::initParallelEnvironment(&argc, &argv);
  if (parallelSetup != 0) {
#ifdef Parallel
    // Please do not use the logging if MPI doesn't work properly.
    std::cerr << "mpi initialisation wasn't successful. Application shut down"
              << std::endl;
#else
    _log.error("main()",
               "mpi initialisation wasn't successful. Application shut down");
#endif
    return parallelSetup;
  }

  int sharedMemorySetup = peano::initSharedMemoryEnvironment();
  if (sharedMemorySetup != 0) {
    logError("main()",
             "shared memory initialisation wasn't successful. Application shut "
             "down");
    return sharedMemorySetup;
  }

  //
  //   Parse config file
  // =====================
  //

  std::string progname = argv[0];

  if (argc < 2) {
    logError("main()", "Usage: " << progname << " config-file [additional args passed to Solver...]");
    return -1;
  }

  // cmdlineargs contains all argv expect the progname.
  std::vector<std::string> cmdlineargs(argv + 1, argv + argc);
  std::string firstarg = cmdlineargs[0];

  bool showHelp    = firstarg == "-h" || firstarg == "--help";
  bool showVersion = firstarg == "-v" || firstarg == "--version";

  if(showHelp) {
    help(progname);
    return EXIT_SUCCESS;
  }

  if(showVersion) {
    version(progname);
    return EXIT_SUCCESS;
  }

  exahype::Parser parser;
  parser.readFile(firstarg);

  if (!parser.isValid()) {
    logError("main()", "invalid config file. Quit");
    return -2;
  }




  #ifdef Asserts
  logInfo( "run()", "start ping-pong test" );
  exahype::Vertex::initDatatype();
  exahype::Vertex sendVertex[5];

  if (tarch::parallel::Node::getInstance().getNumberOfNodes()>1) {
    if (tarch::parallel::Node::getInstance().getRank()==0) {
      sendVertex[0].setIsParentingRegularPersistentSubgridFlag();
      sendVertex[0].setPosition( tarch::la::Vector<DIMENSIONS,double>(2.0), 4);
      sendVertex[0].send(1,100,false,-1);
      sendVertex[1].setIsParentingRegularPersistentSubgridFlag();
      sendVertex[1].setPosition( tarch::la::Vector<DIMENSIONS,double>(3.0), 5);
      sendVertex[1].send(1,100,false,-1);
      sendVertex[2].setIsParentingRegularPersistentSubgridFlag();
      sendVertex[2].setPosition( tarch::la::Vector<DIMENSIONS,double>(4.0), 6);

      sendVertex[2].send(1,100,false,-1);
      logInfo( "run()", "vertex left system" );
      MPI_Send( sendVertex, 3, exahype::Vertex::MPIDatatypeContainer::Datatype, 1, 1, tarch::parallel::Node::getInstance().getCommunicator() );
      logInfo( "run()", "vertices left system" );
    }
    else {
      logInfo( "run()", "start to receive vertex " );
      exahype::Vertex receivedVertex;
      receivedVertex.receive(0,100,false,-1);
      logInfo( "run()", "received vertex " << receivedVertex.toString() );
      assertion1( receivedVertex.getLevel()==4, receivedVertex.toString() );
      assertion1( receivedVertex.getX()(0)==2.0, receivedVertex.toString() );
      assertion1( receivedVertex.getX()(1)==2.0, receivedVertex.toString() );
      assertion1( receivedVertex.getX()(2)==2.0, receivedVertex.toString() );


      exahype::Vertex receivedVertices[5];
      MPI_Recv( receivedVertices, 3, exahype::Vertex::MPIDatatypeContainer::Datatype, 0, 1, tarch::parallel::Node::getInstance().getCommunicator(), MPI_STATUS_IGNORE );
      logInfo( "run()", "received vertices" );
      assertion3( receivedVertices[0].getLevel()==4,  receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      assertion3( receivedVertices[0].getX()(0)==2.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      assertion3( receivedVertices[0].getX()(1)==2.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      assertion3( receivedVertices[0].getX()(2)==2.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );

      assertion3( receivedVertices[1].getLevel()==5,  receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      assertion3( receivedVertices[1].getX()(0)==3.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      assertion3( receivedVertices[1].getX()(1)==3.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      assertion3( receivedVertices[1].getX()(2)==3.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      logInfo( "run()", "ping-poing test ok" );
    }
  }
  #endif

  

  //
  //   Init solver registries
  // =====================================
  //
  kernels::initSolvers(parser, cmdlineargs);

  /*
  // We had this before to show the real solver registration, ie. like
  //    exahype::solvers::RegisteredSolvers->toString(ostream);
  // This needs kernels::initSolvers(parser...) to be called.
  if (onlyShowVersion) {
    version();
    kernels::finalise();
    return 0;
  }
  */

  //
  //   Configure the logging
  // =========================
  //
  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  #ifdef Parallel
  tarch::logging::CommandLineLogger::getInstance().setLogFormat(
      " ",    // columnSeparator
      true,   // logTimeStamp
      false,  // logTimeStampHumanReadable
      true,   // logMachineName
      true,   // logMessageType
      true,   // logTrace
      parser.getLogFileName() );
  #elif defined(Asserts) || defined(Debug)
  tarch::logging::CommandLineLogger::getInstance().setLogFormat(
      " ",    // columnSeparator
      true,   // logTimeStamp
      false,  // logTimeStampHumanReadable
      false,  // logMachineName
      true,   // logMessageType
      true,   // logTrace
      parser.getLogFileName() );
  #else
  tarch::logging::CommandLineLogger::getInstance().setLogFormat(
      " ",    // columnSeparator
      true,   // logTimeStamp
      false,  // logTimeStampHumanReadable
      false,  // logMachineName
      true,   // logMessageType
      false,   // logTrace
      parser.getLogFileName() );
  #endif

  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  if (!tarch::logging::LogFilterFileReader::parsePlainTextFile(
          "exahype.log-filter")) {
    tarch::logging::CommandLineLogger::getInstance().clearFilterList();
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry("info", false));
    #if !defined(Asserts)
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry(
            "info", -1, "peano::grid", true));
    #endif
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry("debug", true));
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry("debug", -1,
                                                             "exahype", false));
  }

//
//   Run tests
// =============
//
#if defined(Debug) || defined(Asserts)
if(! std::getenv("EXAHYPE_SKIP_TESTS")) { // cf issue #74
  tarch::tests::TestCaseRegistry::getInstance().getTestCaseCollection().run();
  int testExitCode = tarch::tests::TestCaseRegistry::getInstance()
                         .getTestCaseCollection()
                         .getNumberOfErrors();

  if (testExitCode != 0) {
    logError("main()", "unit tests failed. Quit.");
    return -2;
  }
} else {
  logInfo("main()", "Skipping tests as EXAHYPE_SKIP_TESTS is set."
     "We do so because tests are broken in the moment and nobody repairs them."); //  TODO(Sven,Dominic,JM): Fix tests.
} // end if getenv(EXAHYPE_SKIP_TESTS)
#endif

  exahype::runners::Runner runner(parser);
  int programExitCode = runner.run();

  if (programExitCode == 0) {
#ifdef Parallel
    if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      logInfo("main()", "Peano terminates successfully");
    }
#else
    logInfo("main()", "Peano terminates successfully");
#endif
  } else {
    logInfo("main()", "quit with error code " << programExitCode);
  }

  peano::shutdownParallelEnvironment();
  peano::shutdownSharedMemoryEnvironment();
  peano::releaseCachedData();

  kernels::finalise();

  return programExitCode;
}


void version(const std::string& programname) {
  std::cout << "This is " << programname << ", an ExaHyPE executable (http://exahype.eu)\n";
  std::cout << "Compiled on host " << EXAHYPE_BUILD_HOST << " at " << EXAHYPE_BUILD_DATE << "\n";
#ifdef EXAHYPE_GIT_INFO
  std::cout << "ExaHyPE git version: " << EXAHYPE_GIT_INFO << "\n";
#else
  std::cout << "ExaHyPE git version: n/a\n";
#endif
#ifdef PEANO_SVN_INFO
  std::cout << "Peano svn version:   " << PEANO_SVN_INFO << "\n";
#else
  std::cout << "Peano svn version:   n/a\n";
#endif
  std::cout << "\n";

  std::cout << "Compile time options\n";
  std::cout << "====================\n";
#ifdef DIMENSIONS
  std::cout << "Dimensions:    "<< DIMENSIONS << "\n";
#else
  std::cout << "Dimensions:    not determinable!\n";
#endif

#ifdef Debug
  std::cout << "Debug:         YES\n";
#else
  std::cout << "Debug:         no\n";
#endif
  
#ifdef Asserts
  std::cout << "Assertions:    YES\n";
#else
  std::cout << "Assertions:    no\n";
#endif

#ifdef Parallel
  std::cout << "MPI Support:   YES\n";
#else
  std::cout << "MPI Support:   no\n";
#endif
  
#ifdef EXAHYPE_CFL_FACTOR // issue #100
  std::cout << "CFL Factor:    "<< EXAHYPE_CFL_FACTOR << "\n";
#else
  std::cout << "CFL Factor:    Default\n";
#endif

  std::cout << "\n";
  std::cout << "Makesystem build options\n";
  std::cout << "========================\n";
#ifdef EXAHYPE_BUILDINFO_AVAILABLE
  std::cout << EXAHYPE_BUILD_INFO << "\n";
#else
  std::cout << "Symbols n/a" << "\n";
#endif

  std::cout << "\n";
  std::cout << "Toolkit static registry info\n";
  std::cout << "============================\n";
  kernels::toString(std::cout);
  std::cout << "\n";
}

void help(const std::string& programname) {
  std::cout << "Usage: " << programname << " <YourApplication.exahype>\n";
  std::cout << "\n";
  std::cout << "   where YourApplication.exahype is an ExaHyPE specification file.\n";
  std::cout << "   Note that you should have compiled ExaHyPE with this file as there\n";
  std::cout << "   are some compile time constants.\n";
  std::cout << "\n";
  std::cout << "   Other possible parameters:\n";
  std::cout << "\n";
  std::cout << "    --help     Show this help message\n";
  std::cout << "    --version  Show version and other hard coded information\n";
  std::cout << "\n";
}
