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
 
#include "exahype/solvers/Solver.h"

namespace {
constexpr const char* tags[]{
    "solutionUpdate",           "volumeIntegral",
    "surfaceIntegral",          "riemannSolver",
    "spaceTimePredictor",       "stableTimeStepSize",
    "solutionAdjustment",       "faceUnknownsProlongation",
    "faceUnknownsRestriction",  "volumeUnknownsProlongation",
    "volumeUnknownsRestriction"};
}  // namespace

// TODO: std::vector<std::unique_ptr<Solver>> ?!
std::vector<exahype::solvers::Solver*> exahype::solvers::RegisteredSolvers;

exahype::solvers::Solver::Solver(const std::string& identifier,
                                 exahype::solvers::Solver::Type type,
                                 int kernelNumber,
                                 int numberOfVariables,
                                 int numberOfParameters,
                                 int nodesPerCoordinateAxis,
                                 double maximumMeshSize,
                                 exahype::solvers::Solver::TimeStepping timeStepping,
                                 std::unique_ptr<profilers::Profiler> profiler)
    : _identifier(identifier),
      _type(type),
      _kernelNumber(kernelNumber),
      _numberOfVariables(numberOfVariables),
      _numberOfParameters(numberOfParameters),
      _nodesPerCoordinateAxis(nodesPerCoordinateAxis),
      _maximumMeshSize(maximumMeshSize),
      _unknownsPerFace(numberOfVariables * addPadding(power(nodesPerCoordinateAxis, DIMENSIONS - 1))),
      _unknownsPerCellBoundary(DIMENSIONS_TIMES_TWO * _unknownsPerFace),
      _unknownsPerCell(numberOfVariables * power(nodesPerCoordinateAxis, DIMENSIONS + 0)),
      _fluxUnknownsPerCell(addPadding(numberOfVariables) * addPadding(power(nodesPerCoordinateAxis, DIMENSIONS + 0)) * DIMENSIONS),
      _spaceTimeUnknownsPerCell(addPadding(numberOfVariables) * power(nodesPerCoordinateAxis, DIMENSIONS + 1)),
      _spaceTimeFluxUnknownsPerCell(_spaceTimeUnknownsPerCell * DIMENSIONS),
      _dataPerCell(addPadding(numberOfVariables)*power(nodesPerCoordinateAxis, DIMENSIONS + 0)),
      _timeStepping(timeStepping),
      _profiler(std::move(profiler)),
      _minCorrectorTimeStamp(std::numeric_limits<double>::max()),
      _minPredictorTimeStamp(std::numeric_limits<double>::max()),
      _minCorrectorTimeStepSize(std::numeric_limits<double>::max()),
      _minPredictorTimeStepSize(std::numeric_limits<double>::max()),
      _minNextPredictorTimeStepSize(std::numeric_limits<double>::max()) {
  assertion(numberOfParameters==0);
  // register tags with profiler
  for (const char* tag : tags) {
    _profiler->registerTag(tag);
  }
}

std::string exahype::solvers::Solver::getIdentifier() const {
  return _identifier;
}

exahype::solvers::Solver::Type exahype::solvers::Solver::getType() const {
  return _type;
}

int exahype::solvers::Solver::getNumberOfVariables() const {
  return _numberOfVariables;
}

int exahype::solvers::Solver::getNumberOfParameters() const {
  return _numberOfParameters;
}

int exahype::solvers::Solver::getNodesPerCoordinateAxis() const {
  return _nodesPerCoordinateAxis;
}

double exahype::solvers::Solver::getMaximumMeshSize() const {
  return _maximumMeshSize;
}


int exahype::solvers::Solver::getUnknownsPerFace() const {
  return _unknownsPerFace;
}

int exahype::solvers::Solver::getUnknownsPerCellBoundary() const {
  return _unknownsPerCellBoundary;
}

int exahype::solvers::Solver::getUnknownsPerCell() const {
  return _unknownsPerCell;
}

int exahype::solvers::Solver::getFluxUnknownsPerCell() const {
  return _fluxUnknownsPerCell;
}

int exahype::solvers::Solver::getSpaceTimeUnknownsPerCell() const {
  return _spaceTimeUnknownsPerCell;
}

int exahype::solvers::Solver::getSpaceTimeFluxUnknownsPerCell() const {
  return _spaceTimeFluxUnknownsPerCell;
}

int exahype::solvers::Solver::getDataPerCell() const {
  return _dataPerCell;
}

void exahype::solvers::Solver::synchroniseTimeStepping(
    exahype::records::ADERDGCellDescription& p) const {
  // todo 16/02/27:Dominic Etienne Charrier
  // in case we use optimistic time stepping:
  // if last predictor time step size is larger
  // as admissibleTimeStepSize + tolerance:
  // make sure that corrector time step size
  // will equal predictor time step size in next
  // sweep.
  // Extra attention must be paid to time stamps.
  // All this should be done by the solver.

  //  if (p.getNextPredictorTimeStepSize() < )

  if (_timeStepping == TimeStepping::Global) {
    p.setCorrectorTimeStamp(_minCorrectorTimeStamp);
    p.setCorrectorTimeStepSize(_minCorrectorTimeStepSize);
    p.setPredictorTimeStamp(_minPredictorTimeStamp);
    p.setPredictorTimeStepSize(_minPredictorTimeStepSize);

    /*
        assertionNumericalEquals1(p.getCorrectorTimeStamp()
       ,solve.getCorrectorTimeStamp(),   1e-12); // todo precision
        assertionNumericalEquals1(p.getCorrectorTimeStepSize(),solve.getCorrectorTimeStepSize(),1e-12);
        assertionNumericalEquals1(p.getPredictorTimeStamp()
       ,solve.getPredictorTimeStamp(),   1e-12);
        assertionNumericalEquals1(p.getPredictorTimeStepSize(),solve.getPredictorTimeStepSize(),1e-12);
     */
  }
/*  if (!solve.isCorrectorTimeLagging()) {*/
//    p.setCorrectorTimeStamp   (p.getPredictorTimeStamp   ());
//    p.setCorrectorTimeStepSize(p.getPredictorTimeStepSize());
//  }

#if defined(Debug) || defined(Asserts)
// @ŧodo Wieder reinnehmen
/*
if (solver.getTimeStepping()==exahype::solvers::Solver::GLOBAL &&
!solver.isCorrectorTimeLagging()) {
  // Note that the solve time stamps and time step sizes are not modified if
corrector time lagging
  // is deactivated. Thus, solve.getPredictor... and solve.getCorrector... are
not the same in general
  // for any value of solve.isCorrectorTimeLagging().
  assertionNumericalEquals1(p.getPredictorTimeStamp()
,solver.getPredictorTimeStamp(),   1e-12); // todo precision
  assertionNumericalEquals1(p.getPredictorTimeStepSize(),solver.getPredictorTimeStepSize(),1e-12);
  assertionNumericalEquals1(p.getCorrectorTimeStamp()
,solver.getPredictorTimeStamp(),   1e-12);
  assertionNumericalEquals1(p.getCorrectorTimeStepSize(),solver.getPredictorTimeStepSize(),1e-12);
}
 */
#endif
}

void exahype::solvers::Solver::startNewTimeStep() {
  _minCorrectorTimeStamp = _minPredictorTimeStamp;
  _minCorrectorTimeStepSize = _minPredictorTimeStepSize;

  _minPredictorTimeStepSize = _minNextPredictorTimeStepSize;
  _minPredictorTimeStamp =
      _minPredictorTimeStamp + _minNextPredictorTimeStepSize;

  _minNextPredictorTimeStepSize = std::numeric_limits<double>::max();
}

void exahype::solvers::Solver::updateMinNextPredictorTimeStepSize(
    const double& minNextPredictorTimeStepSize) {
  _minNextPredictorTimeStepSize =
      std::min(_minNextPredictorTimeStepSize, minNextPredictorTimeStepSize);
}

double exahype::solvers::Solver::getMinNextPredictorTimeStepSize() const {
  return _minNextPredictorTimeStepSize;
}

void exahype::solvers::Solver::setMinCorrectorTimeStamp(
    double minCorrectorTimeStamp) {
  _minCorrectorTimeStamp = minCorrectorTimeStamp;
}

double exahype::solvers::Solver::getMinCorrectorTimeStamp() const {
  return _minCorrectorTimeStamp;
}

void exahype::solvers::Solver::setMinPredictorTimeStamp(
    double minPredictorTimeStamp) {
  _minPredictorTimeStamp = minPredictorTimeStamp;
}

double exahype::solvers::Solver::getMinPredictorTimeStamp() const {
  return _minPredictorTimeStamp;
}

double exahype::solvers::Solver::getMinCorrectorTimeStepSize() const {
  return _minCorrectorTimeStepSize;
}

void exahype::solvers::Solver::setMinPredictorTimeStepSize(
    double minPredictorTimeStepSize) {
  _minPredictorTimeStepSize = minPredictorTimeStepSize;
}

double exahype::solvers::Solver::getMinPredictorTimeStepSize() const {
  return _minPredictorTimeStepSize;
}

void exahype::solvers::Solver::sendToRank(int rank, int tag) {
#ifdef Parallel
  MPI_Send(&_minCorrectorTimeStamp, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator());
  MPI_Send(&_minCorrectorTimeStepSize, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator());
  MPI_Send(&_minPredictorTimeStepSize, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator());
  MPI_Send(&_minPredictorTimeStamp, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator());
  MPI_Send(&_minNextPredictorTimeStepSize, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator());
#endif
}

void exahype::solvers::Solver::receiveFromRank(int rank, int tag) {
#ifdef Parallel
  MPI_Recv(&_minCorrectorTimeStamp, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator(),
           MPI_STATUS_IGNORE);
  MPI_Recv(&_minCorrectorTimeStepSize, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator(),
           MPI_STATUS_IGNORE);
  MPI_Recv(&_minPredictorTimeStepSize, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator(),
           MPI_STATUS_IGNORE);
  MPI_Recv(&_minPredictorTimeStamp, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator(),
           MPI_STATUS_IGNORE);
  MPI_Recv(&_minNextPredictorTimeStepSize, 1, MPI_DOUBLE, rank, tag,
           tarch::parallel::Node::getInstance().getCommunicator(),
           MPI_STATUS_IGNORE);
#endif
}

double exahype::solvers::Solver::getMinSolverTimeStamp() {
  double currentMinTimeStamp = std::numeric_limits<double>::max();

  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    currentMinTimeStamp =
        std::min(currentMinTimeStamp, p->getMinCorrectorTimeStamp());
  }

  return currentMinTimeStamp;
}

double exahype::solvers::Solver::getMinSolverTimeStepSize() {
  double currentMinTimeStepSize = std::numeric_limits<double>::max();

  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    currentMinTimeStepSize =
        std::min(currentMinTimeStepSize, p->getMinCorrectorTimeStepSize());
  }

  return currentMinTimeStepSize;
}
