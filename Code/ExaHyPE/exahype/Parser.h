#ifndef EXAHYPE_PARSER
#define EXAHYPE_PARSER

namespace exahype {
  class Parser;
}


#include <vector>
#include "tarch/logging/Log.h"
#include "tarch/la/Vector.h"
#include "peano/utils/Globals.h"


/**
 * ExaHyPE command line parser
 *
 * @author Tobias Weinzierl
 */
class exahype::Parser {
  private:
    static tarch::logging::Log _log;

    std::vector< std::string > _tokenStream;

    /**
     * @return "notoken" if not found.
     */
    std::string getTokenAfter( std::string token, int additionalTokensToSkip=0 ) const;
    std::string getTokenAfter( std::string token0, std::string token1, int additionalTokensToSkip=0 ) const;
    std::string getTokenAfter( std::string token0, int occurance0, std::string token1, int occurance1, int additionalTokensToSkip=0 ) const;
  public:
    enum MulticoreOracleType {
      Dummy,
      Autotuning,
      GrainSizeSampling
      // evtl. spaeter mal InvadeSHM
    };

    void readFile( const std::string& filename );

    bool isValid() const;

    /**
     * @return How many threads is the code supposed to use?
     */
    int getNumberOfThreads() const;

    /**
     * Our domain always is cubic.
     */
    double getSize() const;

    tarch::la::Vector<DIMENSIONS,double> getOffset() const;

    std::string getMulticorePropertiesFile() const;

    MulticoreOracleType getMulticoreOracleType() const;

    double getSimulationEndTime() const;

    bool fuseAlgorithmicSteps() const;

    double getFirstSnapshotTimeForPlotter( int solverNumber, int plotterNumber ) const;
    double getRepeatTimeForPlotter( int solverNumber, int plotterNumber ) const;
    std::string getIdentifierForPlotter( int solverNumber, int plotterNumber ) const;
    std::string getFilenameForPlotter( int solverNumber, int plotterNumber ) const;
};

#endif

