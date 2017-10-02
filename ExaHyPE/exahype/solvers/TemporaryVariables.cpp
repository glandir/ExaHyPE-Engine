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

#include "exahype/solvers/TemporaryVariables.h"

#include "exahype/solvers/Solver.h"
#include "exahype/solvers/LimitingADERDGSolver.h"

#include <algorithm>
#include <mm_malloc.h> //g++
#include <cstring> //memset

#include "tarch/Assertions.h"

exahype::solvers::PredictionTemporaryVariables::PredictionTemporaryVariables() {}

void exahype::solvers::initialiseTemporaryVariables(exahype::solvers::PredictionTemporaryVariables& temporaryVariables) {
  assertion(temporaryVariables._tempSpaceTimeUnknowns    ==nullptr);
  assertion(temporaryVariables._tempSpaceTimeFluxUnknowns==nullptr);
  assertion(temporaryVariables._tempUnknowns             ==nullptr);
  assertion(temporaryVariables._tempFluxUnknowns         ==nullptr);
  assertion(temporaryVariables._tempStateSizedVectors    ==nullptr);
  assertion(temporaryVariables._tempPointForceSources    ==nullptr);

  int numberOfSolvers        = exahype::solvers::RegisteredSolvers.size();
  temporaryVariables._tempSpaceTimeUnknowns     = new double**[numberOfSolvers]; // == lQi, rhs
  temporaryVariables._tempSpaceTimeFluxUnknowns = new double**[numberOfSolvers]; // == lFi, gradQ
  temporaryVariables._tempUnknowns              = new double* [numberOfSolvers]; // == lQhi
  temporaryVariables._tempFluxUnknowns          = new double* [numberOfSolvers]; // == lFhi
  temporaryVariables._tempStateSizedVectors     = new double* [numberOfSolvers]; // == BGradQ
  temporaryVariables._tempPointForceSources     = new double* [numberOfSolvers];

  exahype::solvers::ADERDGSolver* aderdgSolver = nullptr;

  int solverNumber=0;
  for (auto solver : exahype::solvers::RegisteredSolvers) {
    switch( solver->getType() ) {
    case exahype::solvers::Solver::Type::ADERDG:
      aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
      break;
    case exahype::solvers::Solver::Type::LimitingADERDG:
      aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
      break;
    default:
      aderdgSolver = nullptr;
      break;
    }

    if (aderdgSolver!=nullptr) {
      if(aderdgSolver->alignTempArray()) {
        temporaryVariables._tempSpaceTimeUnknowns[solverNumber] = new double*[2];
        for (int i=0; i<2; ++i) { // max; see spaceTimePredictorNonlinear
          #ifdef ALIGNMENT
          temporaryVariables._tempSpaceTimeUnknowns[solverNumber][i] =
              (double *) _mm_malloc(sizeof(double)*aderdgSolver->getTempSpaceTimeUnknownsSize(), ALIGNMENT);
          #else
          temporaryVariables._tempSpaceTimeUnknowns[solverNumber][i] = new double[aderdgSolver->getTempSpaceTimeUnknownsSize()];
          #endif
          std::memset(temporaryVariables._tempSpaceTimeUnknowns[solverNumber][i], 0, sizeof(double)*aderdgSolver->getTempSpaceTimeUnknownsSize());
        }
        //
        temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber] = new double*[2];
        for (int i=0; i<2; ++i) { // max; see spaceTimePredictorNonlinear
          #ifdef ALIGNMENT
          temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber][i] =
              (double *) _mm_malloc(sizeof(double)*aderdgSolver->getTempSpaceTimeFluxUnknownsSize(), ALIGNMENT);
          #else
          temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber][i] = new double[aderdgSolver->getTempSpaceTimeFluxUnknownsSize()];
          #endif
          std::memset(temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber][i], 0, sizeof(double)*aderdgSolver->getTempSpaceTimeFluxUnknownsSize());
        }
        //
        #ifdef ALIGNMENT
        temporaryVariables._tempUnknowns    [solverNumber]      = (double *) _mm_malloc(sizeof(double)*aderdgSolver->getTempUnknownsSize(), ALIGNMENT);
        //
        temporaryVariables._tempFluxUnknowns[solverNumber]      = (double *) _mm_malloc(sizeof(double)*aderdgSolver->getTempFluxUnknownsSize(), ALIGNMENT);
         //
        temporaryVariables._tempStateSizedVectors[solverNumber] = (double *) _mm_malloc(sizeof(double)*aderdgSolver->getTempStateSizedVectorsSize(), ALIGNMENT);
        #else
        temporaryVariables._tempUnknowns    [solverNumber]      = new double[aderdgSolver->getTempUnknownsSize()];
        //
        temporaryVariables._tempFluxUnknowns[solverNumber]      = new double[aderdgSolver->getTempFluxUnknownsSize()];
         //
        temporaryVariables._tempStateSizedVectors[solverNumber] = new double[aderdgSolver->getTempStateSizedVectorsSize()];
        #endif


        #ifdef ALIGNMENT
        temporaryVariables._tempPointForceSources    [solverNumber] = (double *) _mm_malloc(sizeof(double)*aderdgSolver->getTempSpaceTimeUnknownsSize(), ALIGNMENT);
        #else
        temporaryVariables._tempPointForceSources    [solverNumber] = new double[aderdgSolver->getTempSpaceTimeUnknownsSize()];
        #endif

      } else {
        temporaryVariables._tempSpaceTimeUnknowns[solverNumber] = new double*[2];
        for (int i=0; i<2; ++i) { // max; see spaceTimePredictorNonlinear
          temporaryVariables._tempSpaceTimeUnknowns[solverNumber][i] =
              new double[aderdgSolver->getTempSpaceTimeUnknownsSize()]();
        }
        //
        temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber] = new double*[2];
        for (int i=0; i<2; ++i) { // max; see spaceTimePredictorNonlinear
          temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber][i] =
              new double[aderdgSolver->getTempSpaceTimeFluxUnknownsSize()]();
        }
        //
        temporaryVariables._tempUnknowns    [solverNumber]      = new double[aderdgSolver->getTempUnknownsSize()];
        //
        temporaryVariables._tempFluxUnknowns[solverNumber]      = new double[aderdgSolver->getTempFluxUnknownsSize()];
         //
        temporaryVariables._tempStateSizedVectors[solverNumber] = new double[aderdgSolver->getTempStateSizedVectorsSize()];

        temporaryVariables._tempPointForceSources[solverNumber] = new double[aderdgSolver->getTempSpaceTimeUnknownsSize()];
      }
    } else {
      temporaryVariables._tempSpaceTimeUnknowns    [solverNumber] = nullptr;
      temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber] = nullptr;
      temporaryVariables._tempUnknowns             [solverNumber] = nullptr;
      temporaryVariables._tempFluxUnknowns         [solverNumber] = nullptr;
      temporaryVariables._tempStateSizedVectors    [solverNumber] = nullptr;
      temporaryVariables._tempPointForceSources    [solverNumber] = nullptr;
    }

    ++solverNumber;
  }
}

void exahype::solvers::deleteTemporaryVariables(exahype::solvers::PredictionTemporaryVariables& temporaryVariables) {
  if (temporaryVariables._tempSpaceTimeUnknowns!=nullptr) {
    assertion(temporaryVariables._tempSpaceTimeFluxUnknowns!=nullptr);
    assertion(temporaryVariables._tempUnknowns             !=nullptr);
    assertion(temporaryVariables._tempFluxUnknowns         !=nullptr);
    assertion(temporaryVariables._tempStateSizedVectors    !=nullptr);
    assertion(temporaryVariables._tempPointForceSources    !=nullptr);

    int solverNumber=0;
    exahype::solvers::ADERDGSolver* aderdgSolver = nullptr;

    for (auto solver : exahype::solvers::RegisteredSolvers) {
      switch( solver->getType() ) {
      case exahype::solvers::Solver::Type::ADERDG:
        aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
        break;
      case exahype::solvers::Solver::Type::LimitingADERDG:
        aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
        break;
      default:
        aderdgSolver = nullptr;
        break;
      }

      if (aderdgSolver!=nullptr) {
        if(aderdgSolver->alignTempArray()) {
          //
          for (int i=0; i<2; ++i) {
            _mm_free(temporaryVariables._tempSpaceTimeUnknowns[solverNumber][i]);
          }
          delete[] temporaryVariables._tempSpaceTimeUnknowns[solverNumber];
          temporaryVariables._tempSpaceTimeUnknowns[solverNumber] = nullptr;
          //
          for (int i=0; i<2; ++i) {
            _mm_free(temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber][i]);
          }
          delete[] temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber];
          temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber] = nullptr;
          //
          _mm_free(temporaryVariables._tempUnknowns[solverNumber]);
          temporaryVariables._tempUnknowns[solverNumber] = nullptr;
          //
          _mm_free(temporaryVariables._tempFluxUnknowns[solverNumber]);
          temporaryVariables._tempFluxUnknowns[solverNumber] = nullptr;
          //
          _mm_free(temporaryVariables._tempStateSizedVectors[solverNumber]);
          temporaryVariables._tempStateSizedVectors[solverNumber] = nullptr;

          _mm_free(temporaryVariables._tempPointForceSources[solverNumber]);
          temporaryVariables._tempPointForceSources[solverNumber] = nullptr;
          
        } else {
          //
          for (int i=0; i<2; ++i) {
            delete[] temporaryVariables._tempSpaceTimeUnknowns[solverNumber][i];
          }
          delete[] temporaryVariables._tempSpaceTimeUnknowns[solverNumber];
          temporaryVariables._tempSpaceTimeUnknowns[solverNumber] = nullptr;
          //
          for (int i=0; i<2; ++i) {
            delete[] temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber][i];
          }
          delete[] temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber];
          temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber] = nullptr;
          //
          delete[] temporaryVariables._tempUnknowns[solverNumber];
          temporaryVariables._tempUnknowns[solverNumber] = nullptr;
          //
          delete[] temporaryVariables._tempFluxUnknowns[solverNumber];
          temporaryVariables._tempFluxUnknowns[solverNumber] = nullptr;
          //
          delete[] temporaryVariables._tempStateSizedVectors[solverNumber];
          temporaryVariables._tempStateSizedVectors[solverNumber] = nullptr;

          delete[] temporaryVariables._tempPointForceSources[solverNumber];
          temporaryVariables._tempPointForceSources[solverNumber] = nullptr;

        }

      }

      ++solverNumber;
    }

    delete[] temporaryVariables._tempSpaceTimeUnknowns;
    delete[] temporaryVariables._tempSpaceTimeFluxUnknowns;
    delete[] temporaryVariables._tempUnknowns;
    delete[] temporaryVariables._tempFluxUnknowns;
    delete[] temporaryVariables._tempStateSizedVectors;
    delete[] temporaryVariables._tempPointForceSources;
    temporaryVariables._tempSpaceTimeUnknowns     = nullptr;
    temporaryVariables._tempSpaceTimeFluxUnknowns = nullptr;
    temporaryVariables._tempUnknowns              = nullptr;
    temporaryVariables._tempFluxUnknowns          = nullptr;
    temporaryVariables._tempStateSizedVectors     = nullptr;
    temporaryVariables._tempPointForceSources     = nullptr;
  }
}

exahype::solvers::MergingTemporaryVariables::MergingTemporaryVariables() {}

void exahype::solvers::initialiseTemporaryVariables(exahype::solvers::MergingTemporaryVariables& temporaryVariables) {
  assertion(temporaryVariables._tempStateSizedVectors       ==nullptr);
  assertion(temporaryVariables._tempStateSizedSquareMatrices==nullptr);
  assertion(temporaryVariables._tempFaceUnknowns            ==nullptr);

  int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  temporaryVariables._tempStateSizedVectors        = new double**[numberOfSolvers];
  temporaryVariables._tempStateSizedSquareMatrices = new double**[numberOfSolvers];
  temporaryVariables._tempFaceUnknowns             = new double**[numberOfSolvers];
//    _tempSpaceTimeFaceUnknownsArray = new double* [numberOfSolvers]; todo

  int solverNumber=0;
  for (auto solver : exahype::solvers::RegisteredSolvers) {
    int numberOfStateSizedVectors  = 0;
    int numberOfStateSizedMatrices = 0;
    int numberOfFaceUnknowns       = 0;

    int lengthOfStateSizedVectors  = 0;
    int lengthOfFaceUnknowns       = 0;

    switch (solver->getType()) {
      case exahype::solvers::Solver::Type::ADERDG:
        numberOfStateSizedVectors  = 6; // See riemannSolverNonlinear
        numberOfStateSizedMatrices = 3; // See riemannSolverLinear
        numberOfFaceUnknowns       = 3; // See exahype::solvers::ADERDGSolver::applyBoundaryConditions
        lengthOfFaceUnknowns       =
            static_cast<exahype::solvers::ADERDGSolver*>(solver)->getBndFaceSize(); // == getDataPerFace() + eventual padding
        lengthOfStateSizedVectors       =
             static_cast<exahype::solvers::ADERDGSolver*>(solver)->getTempStateSizedVectorsSize(); // variables + parameters
        break;
      case exahype::solvers::Solver::Type::LimitingADERDG:
        // Needs the same temporary variables as the normal ADER-DG scheme.
        numberOfStateSizedVectors  = 6;
        numberOfStateSizedMatrices = 3;
        numberOfFaceUnknowns       = 3;
        lengthOfFaceUnknowns       = std::max(
            static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver()->getBndFaceSize(), // == getDataPerFace() + eventual padding
            static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiter()->getDataPerPatchFace() );
        lengthOfStateSizedVectors       =
             static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver()->getTempStateSizedVectorsSize(); // variables + parameters; same for both solvers
        break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        numberOfFaceUnknowns = 2; // See exahype::solvers::FiniteVolumesSolver::mergeWithBoundaryData
        lengthOfFaceUnknowns =
            static_cast<exahype::solvers::FiniteVolumesSolver*>(solver)->getDataPerPatchFace();
        lengthOfStateSizedVectors       =
            static_cast<exahype::solvers::FiniteVolumesSolver*>(solver)->getTempStateSizedVectorsSize(); // variables + parameters
        break;
    }

    temporaryVariables._tempStateSizedVectors[solverNumber] = nullptr;
    if (numberOfStateSizedVectors>0) {
      temporaryVariables._tempStateSizedVectors[solverNumber] = new double*[numberOfStateSizedVectors];
      temporaryVariables._tempStateSizedVectors[solverNumber][0] = new double[numberOfStateSizedVectors * lengthOfStateSizedVectors];
      for (int i=1; i<numberOfStateSizedVectors; ++i) { // see riemanSolverLinear
        temporaryVariables._tempStateSizedVectors[solverNumber][i] = temporaryVariables._tempStateSizedVectors[solverNumber][i-1] + lengthOfStateSizedVectors;
      }
    }
    //
    temporaryVariables._tempStateSizedSquareMatrices[solverNumber] = nullptr;
    if (numberOfStateSizedMatrices>0) {
      temporaryVariables._tempStateSizedSquareMatrices[solverNumber] = new double*[numberOfStateSizedMatrices];
      temporaryVariables._tempStateSizedSquareMatrices[solverNumber][0] =
          new double[numberOfStateSizedMatrices* solver->getNumberOfVariables() * solver->getNumberOfVariables()];
      for (int i=1; i<numberOfStateSizedMatrices; ++i) { // see riemanSolverLinear
        temporaryVariables._tempStateSizedSquareMatrices[solverNumber][i] =
            temporaryVariables. _tempStateSizedSquareMatrices[solverNumber][i-1] +
            solver->getNumberOfVariables() * solver->getNumberOfVariables();
      }
    }
    //
    temporaryVariables._tempFaceUnknowns[solverNumber] = nullptr;
    if (numberOfFaceUnknowns>0) {
      temporaryVariables._tempFaceUnknowns[solverNumber] = new double*[numberOfFaceUnknowns];
      temporaryVariables._tempFaceUnknowns[solverNumber][0] = new double[numberOfFaceUnknowns*lengthOfFaceUnknowns](); //initialized to 0 to ensure padding is initialized if existing
      for (int i=1; i<numberOfFaceUnknowns; ++i) { // see ADERDGSolver::applyBoundaryConditions(...)
        temporaryVariables._tempFaceUnknowns[solverNumber][i] = temporaryVariables._tempFaceUnknowns[solverNumber][i-1] + lengthOfFaceUnknowns;
      }
    }

    ++solverNumber;
  }
}

void exahype::solvers::deleteTemporaryVariables(exahype::solvers::MergingTemporaryVariables& temporaryVariables) {
  if (temporaryVariables._tempStateSizedVectors!=nullptr) {
    assertion(temporaryVariables._tempStateSizedSquareMatrices!=nullptr);

    int solverNumber=0;
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      int numberOfStateSizedVectors  = 0;
      int numberOfStateSizedMatrices = 0;
      int numberOfFaceUnknowns       = 0;
      switch (solver->getType()) {
        case exahype::solvers::Solver::Type::ADERDG:
        case exahype::solvers::Solver::Type::LimitingADERDG:
          numberOfStateSizedVectors  = 6; // See riemannSolverLinear
          numberOfStateSizedMatrices = 3; // See riemannSolverLinear
          numberOfFaceUnknowns       = 3; // See exahype::solvers::ADERDGSolver::applyBoundaryConditions
          break;
        case exahype::solvers::Solver::Type::FiniteVolumes:
          numberOfFaceUnknowns       = 2; // See exahype::solvers::FiniteVolumesSolver::mergeWithBoundaryData
          break;
      }

      if (numberOfStateSizedVectors>0) {
        delete[] temporaryVariables._tempStateSizedVectors[solverNumber][0];
        delete[] temporaryVariables._tempStateSizedVectors[solverNumber];
        temporaryVariables._tempStateSizedVectors[solverNumber] = nullptr;
      }
      //
      if (numberOfStateSizedMatrices>0) {
        delete[] temporaryVariables._tempStateSizedSquareMatrices[solverNumber][0];
        delete[] temporaryVariables._tempStateSizedSquareMatrices[solverNumber];
        temporaryVariables._tempStateSizedSquareMatrices[solverNumber] = nullptr;
      }
      //
      if (numberOfFaceUnknowns>0) {
        delete[] temporaryVariables._tempFaceUnknowns[solverNumber][0];
        delete[] temporaryVariables._tempFaceUnknowns[solverNumber];
        temporaryVariables._tempFaceUnknowns[solverNumber] = nullptr;
      }
      //
      // _tempSpaceTimeFaceUnknownsArray[solverNumber] = nullptr; // todo

      ++solverNumber;
    }

    delete[] temporaryVariables._tempStateSizedVectors;
    delete[] temporaryVariables._tempStateSizedSquareMatrices;
    delete[] temporaryVariables._tempFaceUnknowns;
//    delete[] _tempSpaceTimeFaceUnknownsArray; todo
    temporaryVariables._tempStateSizedVectors        = nullptr;
    temporaryVariables._tempStateSizedSquareMatrices = nullptr;
    temporaryVariables._tempFaceUnknowns             = nullptr;
//    _tempSpaceTimeFaceUnknownsArray  = nullptr; todo
  }
}

exahype::solvers::SolutionUpdateTemporaryVariables::SolutionUpdateTemporaryVariables() {}

void exahype::solvers::initialiseTemporaryVariables(exahype::solvers::SolutionUpdateTemporaryVariables& temporaryVariables) {
  assertion(temporaryVariables._tempStateSizedVectors==nullptr);
  assertion(temporaryVariables._tempUnknowns         ==nullptr);

  int numberOfSolvers    = exahype::solvers::RegisteredSolvers.size();
  temporaryVariables._tempStateSizedVectors = new double**[numberOfSolvers];
  temporaryVariables._tempUnknowns          = new double**[numberOfSolvers];

  int solverNumber=0;
  for (auto solver : exahype::solvers::RegisteredSolvers) {
    if  (solver->getType()==exahype::solvers::Solver::Type::FiniteVolumes ||
        solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
      const int numberOfStateSizedVectors = 1+2*DIMENSIONS; // max; see riemannSolverNonlinear(5) or kernels::finitevolumes::godunov::solutionUpdate (1+2*DIMENSIONS)
      temporaryVariables._tempStateSizedVectors[solverNumber]    = new double*[numberOfStateSizedVectors];
      temporaryVariables._tempStateSizedVectors[solverNumber][0] = new double[numberOfStateSizedVectors * solver->getNumberOfVariables()];
      for (int i=1; i<numberOfStateSizedVectors; ++i) {
        temporaryVariables._tempStateSizedVectors[solverNumber][i] = temporaryVariables._tempStateSizedVectors[solverNumber][i-1]+solver->getNumberOfVariables();
      }
      //
      // TODO(Dominic): This will change if we use a different method than a 1st order Godunov method:
      temporaryVariables._tempUnknowns[solverNumber] = nullptr;
    } else {
      temporaryVariables._tempStateSizedVectors[solverNumber] = nullptr;
      temporaryVariables._tempUnknowns         [solverNumber] = nullptr;
    }

    ++solverNumber;
  }
}

void exahype::solvers::deleteTemporaryVariables(exahype::solvers::SolutionUpdateTemporaryVariables& temporaryVariables) {
  if (temporaryVariables._tempStateSizedVectors!=nullptr) {
    assertion(temporaryVariables._tempUnknowns!=nullptr);

    int solverNumber=0;
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      if (solver->getType()==exahype::solvers::Solver::Type::FiniteVolumes ||
          solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
        //
        delete[] temporaryVariables._tempStateSizedVectors[solverNumber][0];
        delete[] temporaryVariables._tempStateSizedVectors[solverNumber];
        temporaryVariables._tempStateSizedVectors[solverNumber] = nullptr;
        // TODO(Dominic): This will change if we use a different method than a 1st order Godunov method:
        temporaryVariables._tempUnknowns[solverNumber] = nullptr;
      }

      ++solverNumber;
    }

    delete[] temporaryVariables._tempStateSizedVectors;
    delete[] temporaryVariables._tempUnknowns;
    temporaryVariables._tempStateSizedVectors = nullptr;
    temporaryVariables._tempUnknowns          = nullptr;
  }
}

exahype::solvers::TimeStepSizeComputationTemporaryVariables::TimeStepSizeComputationTemporaryVariables() {}

void exahype::solvers::initialiseTemporaryVariables(exahype::solvers::TimeStepSizeComputationTemporaryVariables& temporaryVariables) {
  int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  if (temporaryVariables._tempEigenValues==nullptr && numberOfSolvers>0) {
    temporaryVariables._tempEigenValues = new double*[numberOfSolvers];

    int solverNumber=0;
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      assertion( solver->getNumberOfVariables()>0 );
      temporaryVariables._tempEigenValues[solverNumber]  = new double[
        solver->getNumberOfVariables() + solver->getNumberOfParameters() ];
      ++solverNumber;
    }
  }
}

void exahype::solvers::deleteTemporaryVariables(exahype::solvers::TimeStepSizeComputationTemporaryVariables& temporaryVariables) {
  if(temporaryVariables._tempEigenValues!=nullptr) {
    for (unsigned int solverNumber=0;
        solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
      assertion( temporaryVariables._tempEigenValues[solverNumber]!=nullptr );

      delete[] temporaryVariables._tempEigenValues[solverNumber];
      temporaryVariables._tempEigenValues[solverNumber] = nullptr;
    }

    delete[] temporaryVariables._tempEigenValues;
    temporaryVariables._tempEigenValues = nullptr;
  }
}
