#include "exahype/records/ADERDGCellDescription.h"

exahype::records::ADERDGCellDescription::PersistentRecords::
    PersistentRecords() {}

exahype::records::ADERDGCellDescription::PersistentRecords::PersistentRecords(
    const int& solverNumber, const bool& isParent, const bool& isGhost,
    const std::bitset<DIMENSIONS_TIMES_TWO>& riemannSolvePerformed,
    const double& correctorTimeStepSize, const double& correctorTimeStamp,
    const double& predictorTimeStepSize, const double& predictorTimeStamp,
    const double& nextPredictorTimeStepSize, const int& spaceTimePredictor,
    const int& spaceTimeVolumeFlux, const int& solution, const int& update,
    const int& predictor, const int& volumeFlux,
    const int& extrapolatedPredictor, const int& fluctuation, const int& level,
    const tarch::la::Vector<DIMENSIONS, double>& offset,
    const tarch::la::Vector<DIMENSIONS, double>& size)
    : _solverNumber(solverNumber),
      _isParent(isParent),
      _isGhost(isGhost),
      _riemannSolvePerformed(riemannSolvePerformed),
      _correctorTimeStepSize(correctorTimeStepSize),
      _correctorTimeStamp(correctorTimeStamp),
      _predictorTimeStepSize(predictorTimeStepSize),
      _predictorTimeStamp(predictorTimeStamp),
      _nextPredictorTimeStepSize(nextPredictorTimeStepSize),
      _spaceTimePredictor(spaceTimePredictor),
      _spaceTimeVolumeFlux(spaceTimeVolumeFlux),
      _solution(solution),
      _update(update),
      _predictor(predictor),
      _volumeFlux(volumeFlux),
      _extrapolatedPredictor(extrapolatedPredictor),
      _fluctuation(fluctuation),
      _level(level),
      _offset(offset),
      _size(size) {}

exahype::records::ADERDGCellDescription::ADERDGCellDescription() {}

exahype::records::ADERDGCellDescription::ADERDGCellDescription(
    const PersistentRecords& persistentRecords)
    : _persistentRecords(
          persistentRecords._solverNumber, persistentRecords._isParent,
          persistentRecords._isGhost, persistentRecords._riemannSolvePerformed,
          persistentRecords._correctorTimeStepSize,
          persistentRecords._correctorTimeStamp,
          persistentRecords._predictorTimeStepSize,
          persistentRecords._predictorTimeStamp,
          persistentRecords._nextPredictorTimeStepSize,
          persistentRecords._spaceTimePredictor,
          persistentRecords._spaceTimeVolumeFlux, persistentRecords._solution,
          persistentRecords._update, persistentRecords._predictor,
          persistentRecords._volumeFlux,
          persistentRecords._extrapolatedPredictor,
          persistentRecords._fluctuation, persistentRecords._level,
          persistentRecords._offset, persistentRecords._size) {}

exahype::records::ADERDGCellDescription::ADERDGCellDescription(
    const int& solverNumber, const bool& isParent, const bool& isGhost,
    const std::bitset<DIMENSIONS_TIMES_TWO>& riemannSolvePerformed,
    const double& correctorTimeStepSize, const double& correctorTimeStamp,
    const double& predictorTimeStepSize, const double& predictorTimeStamp,
    const double& nextPredictorTimeStepSize, const int& spaceTimePredictor,
    const int& spaceTimeVolumeFlux, const int& solution, const int& update,
    const int& predictor, const int& volumeFlux,
    const int& extrapolatedPredictor, const int& fluctuation, const int& level,
    const tarch::la::Vector<DIMENSIONS, double>& offset,
    const tarch::la::Vector<DIMENSIONS, double>& size)
    : _persistentRecords(
          solverNumber, isParent, isGhost, riemannSolvePerformed,
          correctorTimeStepSize, correctorTimeStamp, predictorTimeStepSize,
          predictorTimeStamp, nextPredictorTimeStepSize, spaceTimePredictor,
          spaceTimeVolumeFlux, solution, update, predictor, volumeFlux,
          extrapolatedPredictor, fluctuation, level, offset, size) {}

exahype::records::ADERDGCellDescription::~ADERDGCellDescription() {}

std::string exahype::records::ADERDGCellDescription::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::records::ADERDGCellDescription::toString(
    std::ostream& out) const {
  out << "(";
  out << "solverNumber:" << getSolverNumber();
  out << ",";
  out << "isParent:" << getIsParent();
  out << ",";
  out << "isGhost:" << getIsGhost();
  out << ",";
  out << "riemannSolvePerformed:[";
  for (int i = 0; i < DIMENSIONS_TIMES_TWO - 1; i++) {
    out << getRiemannSolvePerformed(i) << ",";
  }
  out << getRiemannSolvePerformed(DIMENSIONS_TIMES_TWO - 1) << "]";
  out << ",";
  out << "correctorTimeStepSize:" << getCorrectorTimeStepSize();
  out << ",";
  out << "correctorTimeStamp:" << getCorrectorTimeStamp();
  out << ",";
  out << "predictorTimeStepSize:" << getPredictorTimeStepSize();
  out << ",";
  out << "predictorTimeStamp:" << getPredictorTimeStamp();
  out << ",";
  out << "nextPredictorTimeStepSize:" << getNextPredictorTimeStepSize();
  out << ",";
  out << "spaceTimePredictor:" << getSpaceTimePredictor();
  out << ",";
  out << "spaceTimeVolumeFlux:" << getSpaceTimeVolumeFlux();
  out << ",";
  out << "solution:" << getSolution();
  out << ",";
  out << "update:" << getUpdate();
  out << ",";
  out << "predictor:" << getPredictor();
  out << ",";
  out << "volumeFlux:" << getVolumeFlux();
  out << ",";
  out << "extrapolatedPredictor:" << getExtrapolatedPredictor();
  out << ",";
  out << "fluctuation:" << getFluctuation();
  out << ",";
  out << "level:" << getLevel();
  out << ",";
  out << "offset:[";
  for (int i = 0; i < DIMENSIONS - 1; i++) {
    out << getOffset(i) << ",";
  }
  out << getOffset(DIMENSIONS - 1) << "]";
  out << ",";
  out << "size:[";
  for (int i = 0; i < DIMENSIONS - 1; i++) {
    out << getSize(i) << ",";
  }
  out << getSize(DIMENSIONS - 1) << "]";
  out << ")";
}

exahype::records::ADERDGCellDescription::PersistentRecords
exahype::records::ADERDGCellDescription::getPersistentRecords() const {
  return _persistentRecords;
}

exahype::records::ADERDGCellDescriptionPacked
exahype::records::ADERDGCellDescription::convert() const {
  return ADERDGCellDescriptionPacked(
      getSolverNumber(), getIsParent(), getIsGhost(),
      getRiemannSolvePerformed(), getCorrectorTimeStepSize(),
      getCorrectorTimeStamp(), getPredictorTimeStepSize(),
      getPredictorTimeStamp(), getNextPredictorTimeStepSize(),
      getSpaceTimePredictor(), getSpaceTimeVolumeFlux(), getSolution(),
      getUpdate(), getPredictor(), getVolumeFlux(), getExtrapolatedPredictor(),
      getFluctuation(), getLevel(), getOffset(), getSize());
}

#ifdef Parallel
tarch::logging::Log exahype::records::ADERDGCellDescription::_log(
    "exahype::records::ADERDGCellDescription");

MPI_Datatype exahype::records::ADERDGCellDescription::Datatype = 0;
MPI_Datatype exahype::records::ADERDGCellDescription::FullDatatype = 0;

void exahype::records::ADERDGCellDescription::initDatatype() {
  {
    ADERDGCellDescription dummyADERDGCellDescription[2];

    const int Attributes = 21;
    MPI_Datatype subtypes[Attributes] = {
        MPI_INT,  // solverNumber
        MPI_CHAR,  // isParent
        MPI_CHAR,  // isGhost
        MPI_INT,  // riemannSolvePerformed
        MPI_DOUBLE,  // correctorTimeStepSize
        MPI_DOUBLE,  // correctorTimeStamp
        MPI_DOUBLE,  // predictorTimeStepSize
        MPI_DOUBLE,  // predictorTimeStamp
        MPI_DOUBLE,  // nextPredictorTimeStepSize
        MPI_INT,  // spaceTimePredictor
        MPI_INT,  // spaceTimeVolumeFlux
        MPI_INT,  // solution
        MPI_INT,  // update
        MPI_INT,  // predictor
        MPI_INT,  // volumeFlux
        MPI_INT,  // extrapolatedPredictor
        MPI_INT,  // fluctuation
        MPI_INT,  // level
        MPI_DOUBLE,  // offset
        MPI_DOUBLE,  // size
        MPI_UB  // end/displacement flag
    };

    int blocklen[Attributes] = {
        1,  // solverNumber
        1,  // isParent
        1,  // isGhost
        DIMENSIONS_TIMES_TWO,  // riemannSolvePerformed
        1,  // correctorTimeStepSize
        1,  // correctorTimeStamp
        1,  // predictorTimeStepSize
        1,  // predictorTimeStamp
        1,  // nextPredictorTimeStepSize
        1,  // spaceTimePredictor
        1,  // spaceTimeVolumeFlux
        1,  // solution
        1,  // update
        1,  // predictor
        1,  // volumeFlux
        1,  // extrapolatedPredictor
        1,  // fluctuation
        1,  // level
        DIMENSIONS,  // offset
        DIMENSIONS,  // size
        1  // end/displacement flag
    };

    MPI_Aint disp[Attributes];

    MPI_Aint base;
    MPI_Address(const_cast<void*>(
                    static_cast<const void*>(&(dummyADERDGCellDescription[0]))),
                &base);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._solverNumber))),
        &disp[0]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._isParent))),
        &disp[1]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._isGhost))),
        &disp[2]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._riemannSolvePerformed))),
                &disp[3]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._correctorTimeStepSize))),
                &disp[4]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._correctorTimeStamp))),
                &disp[5]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._predictorTimeStepSize))),
                &disp[6]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._predictorTimeStamp))),
                &disp[7]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._nextPredictorTimeStepSize))),
                &disp[8]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._spaceTimePredictor))),
                &disp[9]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._spaceTimeVolumeFlux))),
                &disp[10]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._solution))),
        &disp[11]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._update))),
        &disp[12]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._predictor))),
        &disp[13]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._volumeFlux))),
        &disp[14]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._extrapolatedPredictor))),
                &disp[15]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._fluctuation))),
        &disp[16]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._level))),
        &disp[17]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))),
        &disp[18]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._size[0]))),
        &disp[19]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[1]._persistentRecords._solverNumber))),
        &disp[20]);

    for (int i = 1; i < Attributes; i++) {
      assertion1(disp[i] > disp[i - 1], i);
    }
    for (int i = 0; i < Attributes; i++) {
      disp[i] -= base;
    }
    MPI_Type_struct(Attributes, blocklen, disp, subtypes,
                    &ADERDGCellDescription::Datatype);
    MPI_Type_commit(&ADERDGCellDescription::Datatype);
  }
  {
    ADERDGCellDescription dummyADERDGCellDescription[2];

    const int Attributes = 21;
    MPI_Datatype subtypes[Attributes] = {
        MPI_INT,  // solverNumber
        MPI_CHAR,  // isParent
        MPI_CHAR,  // isGhost
        MPI_INT,  // riemannSolvePerformed
        MPI_DOUBLE,  // correctorTimeStepSize
        MPI_DOUBLE,  // correctorTimeStamp
        MPI_DOUBLE,  // predictorTimeStepSize
        MPI_DOUBLE,  // predictorTimeStamp
        MPI_DOUBLE,  // nextPredictorTimeStepSize
        MPI_INT,  // spaceTimePredictor
        MPI_INT,  // spaceTimeVolumeFlux
        MPI_INT,  // solution
        MPI_INT,  // update
        MPI_INT,  // predictor
        MPI_INT,  // volumeFlux
        MPI_INT,  // extrapolatedPredictor
        MPI_INT,  // fluctuation
        MPI_INT,  // level
        MPI_DOUBLE,  // offset
        MPI_DOUBLE,  // size
        MPI_UB  // end/displacement flag
    };

    int blocklen[Attributes] = {
        1,  // solverNumber
        1,  // isParent
        1,  // isGhost
        DIMENSIONS_TIMES_TWO,  // riemannSolvePerformed
        1,  // correctorTimeStepSize
        1,  // correctorTimeStamp
        1,  // predictorTimeStepSize
        1,  // predictorTimeStamp
        1,  // nextPredictorTimeStepSize
        1,  // spaceTimePredictor
        1,  // spaceTimeVolumeFlux
        1,  // solution
        1,  // update
        1,  // predictor
        1,  // volumeFlux
        1,  // extrapolatedPredictor
        1,  // fluctuation
        1,  // level
        DIMENSIONS,  // offset
        DIMENSIONS,  // size
        1  // end/displacement flag
    };

    MPI_Aint disp[Attributes];

    MPI_Aint base;
    MPI_Address(const_cast<void*>(
                    static_cast<const void*>(&(dummyADERDGCellDescription[0]))),
                &base);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._solverNumber))),
        &disp[0]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._isParent))),
        &disp[1]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._isGhost))),
        &disp[2]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._riemannSolvePerformed))),
                &disp[3]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._correctorTimeStepSize))),
                &disp[4]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._correctorTimeStamp))),
                &disp[5]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._predictorTimeStepSize))),
                &disp[6]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._predictorTimeStamp))),
                &disp[7]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._nextPredictorTimeStepSize))),
                &disp[8]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._spaceTimePredictor))),
                &disp[9]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._spaceTimeVolumeFlux))),
                &disp[10]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._solution))),
        &disp[11]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._update))),
        &disp[12]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._predictor))),
        &disp[13]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._volumeFlux))),
        &disp[14]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescription[0]
                          ._persistentRecords._extrapolatedPredictor))),
                &disp[15]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._fluctuation))),
        &disp[16]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._level))),
        &disp[17]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))),
        &disp[18]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[0]._persistentRecords._size[0]))),
        &disp[19]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescription[1]._persistentRecords._solverNumber))),
        &disp[20]);

    for (int i = 1; i < Attributes; i++) {
      assertion1(disp[i] > disp[i - 1], i);
    }
    for (int i = 0; i < Attributes; i++) {
      disp[i] -= base;
    }
    MPI_Type_struct(Attributes, blocklen, disp, subtypes,
                    &ADERDGCellDescription::FullDatatype);
    MPI_Type_commit(&ADERDGCellDescription::FullDatatype);
  }
}

void exahype::records::ADERDGCellDescription::shutdownDatatype() {
  MPI_Type_free(&ADERDGCellDescription::Datatype);
  MPI_Type_free(&ADERDGCellDescription::FullDatatype);
}

void exahype::records::ADERDGCellDescription::send(
    int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise,
    int communicateSleep) {
  if (communicateSleep < 0) {
    const int result = MPI_Send(
        this, 1,
        exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype,
        destination, tag,
        tarch::parallel::Node::getInstance().getCommunicator());
    if (result != MPI_SUCCESS) {
      std::ostringstream msg;
      msg << "was not able to send message "
             "exahype::records::ADERDGCellDescription " << toString()
          << " to node " << destination << ": "
          << tarch::parallel::MPIReturnValueToString(result);
      _log.error("send(int)", msg.str());
    }

  } else {
    MPI_Request* sendRequestHandle = new MPI_Request();
    MPI_Status status;
    int flag = 0;
    int result;

    clock_t timeOutWarning = -1;
    clock_t timeOutShutdown = -1;
    bool triggeredTimeoutWarning = false;

    if (exchangeOnlyAttributesMarkedWithParallelise) {
      result = MPI_Isend(this, 1, Datatype, destination, tag,
                         tarch::parallel::Node::getInstance().getCommunicator(),
                         sendRequestHandle);

    } else {
      result = MPI_Isend(this, 1, FullDatatype, destination, tag,
                         tarch::parallel::Node::getInstance().getCommunicator(),
                         sendRequestHandle);
    }
    if (result != MPI_SUCCESS) {
      std::ostringstream msg;
      msg << "was not able to send message "
             "exahype::records::ADERDGCellDescription " << toString()
          << " to node " << destination << ": "
          << tarch::parallel::MPIReturnValueToString(result);
      _log.error("send(int)", msg.str());
    }
    result = MPI_Test(sendRequestHandle, &flag, &status);
    while (!flag) {
      if (timeOutWarning == -1)
        timeOutWarning =
            tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
      if (timeOutShutdown == -1)
        timeOutShutdown =
            tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
      result = MPI_Test(sendRequestHandle, &flag, &status);
      if (result != MPI_SUCCESS) {
        std::ostringstream msg;
        msg << "testing for finished send task for "
               "exahype::records::ADERDGCellDescription " << toString()
            << " sent to node " << destination
            << " failed: " << tarch::parallel::MPIReturnValueToString(result);
        _log.error("send(int)", msg.str());
      }

      // deadlock aspect
      if (tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
          (clock() > timeOutWarning) && (!triggeredTimeoutWarning)) {
        tarch::parallel::Node::getInstance().writeTimeOutWarning(
            "exahype::records::ADERDGCellDescription", "send(int)", destination,
            tag, 1);
        triggeredTimeoutWarning = true;
      }
      if (tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
          (clock() > timeOutShutdown)) {
        tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "exahype::records::ADERDGCellDescription", "send(int)", destination,
            tag, 1);
      }
      tarch::parallel::Node::getInstance().receiveDanglingMessages();
      usleep(communicateSleep);
    }

    delete sendRequestHandle;
#ifdef Debug
    _log.debug("send(int,int)", "sent " + toString());
#endif
  }
}

void exahype::records::ADERDGCellDescription::receive(
    int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise,
    int communicateSleep) {
  if (communicateSleep < 0) {
    MPI_Status status;
    const int result = MPI_Recv(
        this, 1,
        exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype,
        source, tag, tarch::parallel::Node::getInstance().getCommunicator(),
        &status);
    if (result != MPI_SUCCESS) {
      std::ostringstream msg;
      msg << "failed to start to receive "
             "exahype::records::ADERDGCellDescription from node " << source
          << ": " << tarch::parallel::MPIReturnValueToString(result);
      _log.error("receive(int)", msg.str());
    }

  } else {
    MPI_Request* sendRequestHandle = new MPI_Request();
    MPI_Status status;
    int flag = 0;
    int result;

    clock_t timeOutWarning = -1;
    clock_t timeOutShutdown = -1;
    bool triggeredTimeoutWarning = false;

    if (exchangeOnlyAttributesMarkedWithParallelise) {
      result = MPI_Irecv(this, 1, Datatype, source, tag,
                         tarch::parallel::Node::getInstance().getCommunicator(),
                         sendRequestHandle);

    } else {
      result = MPI_Irecv(this, 1, FullDatatype, source, tag,
                         tarch::parallel::Node::getInstance().getCommunicator(),
                         sendRequestHandle);
    }
    if (result != MPI_SUCCESS) {
      std::ostringstream msg;
      msg << "failed to start to receive "
             "exahype::records::ADERDGCellDescription from node " << source
          << ": " << tarch::parallel::MPIReturnValueToString(result);
      _log.error("receive(int)", msg.str());
    }

    result = MPI_Test(sendRequestHandle, &flag, &status);
    while (!flag) {
      if (timeOutWarning == -1)
        timeOutWarning =
            tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
      if (timeOutShutdown == -1)
        timeOutShutdown =
            tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
      result = MPI_Test(sendRequestHandle, &flag, &status);
      if (result != MPI_SUCCESS) {
        std::ostringstream msg;
        msg << "testing for finished receive task for "
               "exahype::records::ADERDGCellDescription failed: "
            << tarch::parallel::MPIReturnValueToString(result);
        _log.error("receive(int)", msg.str());
      }

      // deadlock aspect
      if (tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
          (clock() > timeOutWarning) && (!triggeredTimeoutWarning)) {
        tarch::parallel::Node::getInstance().writeTimeOutWarning(
            "exahype::records::ADERDGCellDescription", "receive(int)", source,
            tag, 1);
        triggeredTimeoutWarning = true;
      }
      if (tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
          (clock() > timeOutShutdown)) {
        tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "exahype::records::ADERDGCellDescription", "receive(int)", source,
            tag, 1);
      }
      tarch::parallel::Node::getInstance().receiveDanglingMessages();
      usleep(communicateSleep);
    }

    delete sendRequestHandle;

#ifdef Debug
    _log.debug("receive(int,int)", "received " + toString());
#endif
  }
}

bool exahype::records::ADERDGCellDescription::isMessageInQueue(
    int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
  MPI_Status status;
  int flag = 0;
  MPI_Iprobe(MPI_ANY_SOURCE, tag,
             tarch::parallel::Node::getInstance().getCommunicator(), &flag,
             &status);
  if (flag) {
    int messageCounter;
    if (exchangeOnlyAttributesMarkedWithParallelise) {
      MPI_Get_count(&status, Datatype, &messageCounter);
    } else {
      MPI_Get_count(&status, FullDatatype, &messageCounter);
    }
    return messageCounter > 0;
  } else
    return false;
}

#endif

exahype::records::ADERDGCellDescriptionPacked::PersistentRecords::
    PersistentRecords() {}

exahype::records::ADERDGCellDescriptionPacked::PersistentRecords::
    PersistentRecords(
        const int& solverNumber, const bool& isParent, const bool& isGhost,
        const std::bitset<DIMENSIONS_TIMES_TWO>& riemannSolvePerformed,
        const double& correctorTimeStepSize, const double& correctorTimeStamp,
        const double& predictorTimeStepSize, const double& predictorTimeStamp,
        const double& nextPredictorTimeStepSize, const int& spaceTimePredictor,
        const int& spaceTimeVolumeFlux, const int& solution, const int& update,
        const int& predictor, const int& volumeFlux,
        const int& extrapolatedPredictor, const int& fluctuation,
        const int& level, const tarch::la::Vector<DIMENSIONS, double>& offset,
        const tarch::la::Vector<DIMENSIONS, double>& size)
    : _solverNumber(solverNumber),
      _isParent(isParent),
      _isGhost(isGhost),
      _riemannSolvePerformed(riemannSolvePerformed),
      _correctorTimeStepSize(correctorTimeStepSize),
      _correctorTimeStamp(correctorTimeStamp),
      _predictorTimeStepSize(predictorTimeStepSize),
      _predictorTimeStamp(predictorTimeStamp),
      _nextPredictorTimeStepSize(nextPredictorTimeStepSize),
      _spaceTimePredictor(spaceTimePredictor),
      _spaceTimeVolumeFlux(spaceTimeVolumeFlux),
      _solution(solution),
      _update(update),
      _predictor(predictor),
      _volumeFlux(volumeFlux),
      _extrapolatedPredictor(extrapolatedPredictor),
      _fluctuation(fluctuation),
      _level(level),
      _offset(offset),
      _size(size) {}

exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked() {}

exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked(
    const PersistentRecords& persistentRecords)
    : _persistentRecords(
          persistentRecords._solverNumber, persistentRecords._isParent,
          persistentRecords._isGhost, persistentRecords._riemannSolvePerformed,
          persistentRecords._correctorTimeStepSize,
          persistentRecords._correctorTimeStamp,
          persistentRecords._predictorTimeStepSize,
          persistentRecords._predictorTimeStamp,
          persistentRecords._nextPredictorTimeStepSize,
          persistentRecords._spaceTimePredictor,
          persistentRecords._spaceTimeVolumeFlux, persistentRecords._solution,
          persistentRecords._update, persistentRecords._predictor,
          persistentRecords._volumeFlux,
          persistentRecords._extrapolatedPredictor,
          persistentRecords._fluctuation, persistentRecords._level,
          persistentRecords._offset, persistentRecords._size) {}

exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked(
    const int& solverNumber, const bool& isParent, const bool& isGhost,
    const std::bitset<DIMENSIONS_TIMES_TWO>& riemannSolvePerformed,
    const double& correctorTimeStepSize, const double& correctorTimeStamp,
    const double& predictorTimeStepSize, const double& predictorTimeStamp,
    const double& nextPredictorTimeStepSize, const int& spaceTimePredictor,
    const int& spaceTimeVolumeFlux, const int& solution, const int& update,
    const int& predictor, const int& volumeFlux,
    const int& extrapolatedPredictor, const int& fluctuation, const int& level,
    const tarch::la::Vector<DIMENSIONS, double>& offset,
    const tarch::la::Vector<DIMENSIONS, double>& size)
    : _persistentRecords(
          solverNumber, isParent, isGhost, riemannSolvePerformed,
          correctorTimeStepSize, correctorTimeStamp, predictorTimeStepSize,
          predictorTimeStamp, nextPredictorTimeStepSize, spaceTimePredictor,
          spaceTimeVolumeFlux, solution, update, predictor, volumeFlux,
          extrapolatedPredictor, fluctuation, level, offset, size) {}

exahype::records::ADERDGCellDescriptionPacked::~ADERDGCellDescriptionPacked() {}

std::string exahype::records::ADERDGCellDescriptionPacked::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::records::ADERDGCellDescriptionPacked::toString(
    std::ostream& out) const {
  out << "(";
  out << "solverNumber:" << getSolverNumber();
  out << ",";
  out << "isParent:" << getIsParent();
  out << ",";
  out << "isGhost:" << getIsGhost();
  out << ",";
  out << "riemannSolvePerformed:[";
  for (int i = 0; i < DIMENSIONS_TIMES_TWO - 1; i++) {
    out << getRiemannSolvePerformed(i) << ",";
  }
  out << getRiemannSolvePerformed(DIMENSIONS_TIMES_TWO - 1) << "]";
  out << ",";
  out << "correctorTimeStepSize:" << getCorrectorTimeStepSize();
  out << ",";
  out << "correctorTimeStamp:" << getCorrectorTimeStamp();
  out << ",";
  out << "predictorTimeStepSize:" << getPredictorTimeStepSize();
  out << ",";
  out << "predictorTimeStamp:" << getPredictorTimeStamp();
  out << ",";
  out << "nextPredictorTimeStepSize:" << getNextPredictorTimeStepSize();
  out << ",";
  out << "spaceTimePredictor:" << getSpaceTimePredictor();
  out << ",";
  out << "spaceTimeVolumeFlux:" << getSpaceTimeVolumeFlux();
  out << ",";
  out << "solution:" << getSolution();
  out << ",";
  out << "update:" << getUpdate();
  out << ",";
  out << "predictor:" << getPredictor();
  out << ",";
  out << "volumeFlux:" << getVolumeFlux();
  out << ",";
  out << "extrapolatedPredictor:" << getExtrapolatedPredictor();
  out << ",";
  out << "fluctuation:" << getFluctuation();
  out << ",";
  out << "level:" << getLevel();
  out << ",";
  out << "offset:[";
  for (int i = 0; i < DIMENSIONS - 1; i++) {
    out << getOffset(i) << ",";
  }
  out << getOffset(DIMENSIONS - 1) << "]";
  out << ",";
  out << "size:[";
  for (int i = 0; i < DIMENSIONS - 1; i++) {
    out << getSize(i) << ",";
  }
  out << getSize(DIMENSIONS - 1) << "]";
  out << ")";
}

exahype::records::ADERDGCellDescriptionPacked::PersistentRecords
exahype::records::ADERDGCellDescriptionPacked::getPersistentRecords() const {
  return _persistentRecords;
}

exahype::records::ADERDGCellDescription
exahype::records::ADERDGCellDescriptionPacked::convert() const {
  return ADERDGCellDescription(
      getSolverNumber(), getIsParent(), getIsGhost(),
      getRiemannSolvePerformed(), getCorrectorTimeStepSize(),
      getCorrectorTimeStamp(), getPredictorTimeStepSize(),
      getPredictorTimeStamp(), getNextPredictorTimeStepSize(),
      getSpaceTimePredictor(), getSpaceTimeVolumeFlux(), getSolution(),
      getUpdate(), getPredictor(), getVolumeFlux(), getExtrapolatedPredictor(),
      getFluctuation(), getLevel(), getOffset(), getSize());
}

#ifdef Parallel
tarch::logging::Log exahype::records::ADERDGCellDescriptionPacked::_log(
    "exahype::records::ADERDGCellDescriptionPacked");

MPI_Datatype exahype::records::ADERDGCellDescriptionPacked::Datatype = 0;
MPI_Datatype exahype::records::ADERDGCellDescriptionPacked::FullDatatype = 0;

void exahype::records::ADERDGCellDescriptionPacked::initDatatype() {
  {
    ADERDGCellDescriptionPacked dummyADERDGCellDescriptionPacked[2];

    const int Attributes = 21;
    MPI_Datatype subtypes[Attributes] = {
        MPI_INT,  // solverNumber
        MPI_CHAR,  // isParent
        MPI_CHAR,  // isGhost
        MPI_INT,  // riemannSolvePerformed
        MPI_DOUBLE,  // correctorTimeStepSize
        MPI_DOUBLE,  // correctorTimeStamp
        MPI_DOUBLE,  // predictorTimeStepSize
        MPI_DOUBLE,  // predictorTimeStamp
        MPI_DOUBLE,  // nextPredictorTimeStepSize
        MPI_INT,  // spaceTimePredictor
        MPI_INT,  // spaceTimeVolumeFlux
        MPI_INT,  // solution
        MPI_INT,  // update
        MPI_INT,  // predictor
        MPI_INT,  // volumeFlux
        MPI_INT,  // extrapolatedPredictor
        MPI_INT,  // fluctuation
        MPI_INT,  // level
        MPI_DOUBLE,  // offset
        MPI_DOUBLE,  // size
        MPI_UB  // end/displacement flag
    };

    int blocklen[Attributes] = {
        1,  // solverNumber
        1,  // isParent
        1,  // isGhost
        DIMENSIONS_TIMES_TWO,  // riemannSolvePerformed
        1,  // correctorTimeStepSize
        1,  // correctorTimeStamp
        1,  // predictorTimeStepSize
        1,  // predictorTimeStamp
        1,  // nextPredictorTimeStepSize
        1,  // spaceTimePredictor
        1,  // spaceTimeVolumeFlux
        1,  // solution
        1,  // update
        1,  // predictor
        1,  // volumeFlux
        1,  // extrapolatedPredictor
        1,  // fluctuation
        1,  // level
        DIMENSIONS,  // offset
        DIMENSIONS,  // size
        1  // end/displacement flag
    };

    MPI_Aint disp[Attributes];

    MPI_Aint base;
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]))),
                &base);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._solverNumber))),
                &disp[0]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._isParent))),
                &disp[1]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._isGhost))),
                &disp[2]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._riemannSolvePerformed))),
                &disp[3]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._correctorTimeStepSize))),
                &disp[4]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._correctorTimeStamp))),
                &disp[5]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._predictorTimeStepSize))),
                &disp[6]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._predictorTimeStamp))),
                &disp[7]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._nextPredictorTimeStepSize))),
                &disp[8]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._spaceTimePredictor))),
                &disp[9]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._spaceTimeVolumeFlux))),
                &disp[10]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._solution))),
                &disp[11]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescriptionPacked[0]._persistentRecords._update))),
        &disp[12]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._predictor))),
                &disp[13]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._volumeFlux))),
                &disp[14]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._extrapolatedPredictor))),
                &disp[15]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._fluctuation))),
                &disp[16]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))),
        &disp[17]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._offset[0]))),
                &disp[18]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._size[0]))),
                &disp[19]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[1]
                          ._persistentRecords._solverNumber))),
                &disp[20]);

    for (int i = 1; i < Attributes; i++) {
      assertion1(disp[i] > disp[i - 1], i);
    }
    for (int i = 0; i < Attributes; i++) {
      disp[i] -= base;
    }
    MPI_Type_struct(Attributes, blocklen, disp, subtypes,
                    &ADERDGCellDescriptionPacked::Datatype);
    MPI_Type_commit(&ADERDGCellDescriptionPacked::Datatype);
  }
  {
    ADERDGCellDescriptionPacked dummyADERDGCellDescriptionPacked[2];

    const int Attributes = 21;
    MPI_Datatype subtypes[Attributes] = {
        MPI_INT,  // solverNumber
        MPI_CHAR,  // isParent
        MPI_CHAR,  // isGhost
        MPI_INT,  // riemannSolvePerformed
        MPI_DOUBLE,  // correctorTimeStepSize
        MPI_DOUBLE,  // correctorTimeStamp
        MPI_DOUBLE,  // predictorTimeStepSize
        MPI_DOUBLE,  // predictorTimeStamp
        MPI_DOUBLE,  // nextPredictorTimeStepSize
        MPI_INT,  // spaceTimePredictor
        MPI_INT,  // spaceTimeVolumeFlux
        MPI_INT,  // solution
        MPI_INT,  // update
        MPI_INT,  // predictor
        MPI_INT,  // volumeFlux
        MPI_INT,  // extrapolatedPredictor
        MPI_INT,  // fluctuation
        MPI_INT,  // level
        MPI_DOUBLE,  // offset
        MPI_DOUBLE,  // size
        MPI_UB  // end/displacement flag
    };

    int blocklen[Attributes] = {
        1,  // solverNumber
        1,  // isParent
        1,  // isGhost
        DIMENSIONS_TIMES_TWO,  // riemannSolvePerformed
        1,  // correctorTimeStepSize
        1,  // correctorTimeStamp
        1,  // predictorTimeStepSize
        1,  // predictorTimeStamp
        1,  // nextPredictorTimeStepSize
        1,  // spaceTimePredictor
        1,  // spaceTimeVolumeFlux
        1,  // solution
        1,  // update
        1,  // predictor
        1,  // volumeFlux
        1,  // extrapolatedPredictor
        1,  // fluctuation
        1,  // level
        DIMENSIONS,  // offset
        DIMENSIONS,  // size
        1  // end/displacement flag
    };

    MPI_Aint disp[Attributes];

    MPI_Aint base;
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]))),
                &base);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._solverNumber))),
                &disp[0]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._isParent))),
                &disp[1]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._isGhost))),
                &disp[2]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._riemannSolvePerformed))),
                &disp[3]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._correctorTimeStepSize))),
                &disp[4]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._correctorTimeStamp))),
                &disp[5]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._predictorTimeStepSize))),
                &disp[6]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._predictorTimeStamp))),
                &disp[7]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._nextPredictorTimeStepSize))),
                &disp[8]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._spaceTimePredictor))),
                &disp[9]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._spaceTimeVolumeFlux))),
                &disp[10]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._solution))),
                &disp[11]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescriptionPacked[0]._persistentRecords._update))),
        &disp[12]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._predictor))),
                &disp[13]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._volumeFlux))),
                &disp[14]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._extrapolatedPredictor))),
                &disp[15]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._fluctuation))),
                &disp[16]);
    MPI_Address(
        const_cast<void*>(static_cast<const void*>(
            &(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))),
        &disp[17]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._offset[0]))),
                &disp[18]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[0]
                          ._persistentRecords._size[0]))),
                &disp[19]);
    MPI_Address(const_cast<void*>(static_cast<const void*>(
                    &(dummyADERDGCellDescriptionPacked[1]
                          ._persistentRecords._solverNumber))),
                &disp[20]);

    for (int i = 1; i < Attributes; i++) {
      assertion1(disp[i] > disp[i - 1], i);
    }
    for (int i = 0; i < Attributes; i++) {
      disp[i] -= base;
    }
    MPI_Type_struct(Attributes, blocklen, disp, subtypes,
                    &ADERDGCellDescriptionPacked::FullDatatype);
    MPI_Type_commit(&ADERDGCellDescriptionPacked::FullDatatype);
  }
}

void exahype::records::ADERDGCellDescriptionPacked::shutdownDatatype() {
  MPI_Type_free(&ADERDGCellDescriptionPacked::Datatype);
  MPI_Type_free(&ADERDGCellDescriptionPacked::FullDatatype);
}

void exahype::records::ADERDGCellDescriptionPacked::send(
    int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise,
    int communicateSleep) {
  if (communicateSleep < 0) {
    const int result = MPI_Send(
        this, 1,
        exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype,
        destination, tag,
        tarch::parallel::Node::getInstance().getCommunicator());
    if (result != MPI_SUCCESS) {
      std::ostringstream msg;
      msg << "was not able to send message "
             "exahype::records::ADERDGCellDescriptionPacked " << toString()
          << " to node " << destination << ": "
          << tarch::parallel::MPIReturnValueToString(result);
      _log.error("send(int)", msg.str());
    }

  } else {
    MPI_Request* sendRequestHandle = new MPI_Request();
    MPI_Status status;
    int flag = 0;
    int result;

    clock_t timeOutWarning = -1;
    clock_t timeOutShutdown = -1;
    bool triggeredTimeoutWarning = false;

    if (exchangeOnlyAttributesMarkedWithParallelise) {
      result = MPI_Isend(this, 1, Datatype, destination, tag,
                         tarch::parallel::Node::getInstance().getCommunicator(),
                         sendRequestHandle);

    } else {
      result = MPI_Isend(this, 1, FullDatatype, destination, tag,
                         tarch::parallel::Node::getInstance().getCommunicator(),
                         sendRequestHandle);
    }
    if (result != MPI_SUCCESS) {
      std::ostringstream msg;
      msg << "was not able to send message "
             "exahype::records::ADERDGCellDescriptionPacked " << toString()
          << " to node " << destination << ": "
          << tarch::parallel::MPIReturnValueToString(result);
      _log.error("send(int)", msg.str());
    }
    result = MPI_Test(sendRequestHandle, &flag, &status);
    while (!flag) {
      if (timeOutWarning == -1)
        timeOutWarning =
            tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
      if (timeOutShutdown == -1)
        timeOutShutdown =
            tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
      result = MPI_Test(sendRequestHandle, &flag, &status);
      if (result != MPI_SUCCESS) {
        std::ostringstream msg;
        msg << "testing for finished send task for "
               "exahype::records::ADERDGCellDescriptionPacked " << toString()
            << " sent to node " << destination
            << " failed: " << tarch::parallel::MPIReturnValueToString(result);
        _log.error("send(int)", msg.str());
      }

      // deadlock aspect
      if (tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
          (clock() > timeOutWarning) && (!triggeredTimeoutWarning)) {
        tarch::parallel::Node::getInstance().writeTimeOutWarning(
            "exahype::records::ADERDGCellDescriptionPacked", "send(int)",
            destination, tag, 1);
        triggeredTimeoutWarning = true;
      }
      if (tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
          (clock() > timeOutShutdown)) {
        tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "exahype::records::ADERDGCellDescriptionPacked", "send(int)",
            destination, tag, 1);
      }
      tarch::parallel::Node::getInstance().receiveDanglingMessages();
      usleep(communicateSleep);
    }

    delete sendRequestHandle;
#ifdef Debug
    _log.debug("send(int,int)", "sent " + toString());
#endif
  }
}

void exahype::records::ADERDGCellDescriptionPacked::receive(
    int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise,
    int communicateSleep) {
  if (communicateSleep < 0) {
    MPI_Status status;
    const int result = MPI_Recv(
        this, 1,
        exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype,
        source, tag, tarch::parallel::Node::getInstance().getCommunicator(),
        &status);
    if (result != MPI_SUCCESS) {
      std::ostringstream msg;
      msg << "failed to start to receive "
             "exahype::records::ADERDGCellDescriptionPacked from node "
          << source << ": " << tarch::parallel::MPIReturnValueToString(result);
      _log.error("receive(int)", msg.str());
    }

  } else {
    MPI_Request* sendRequestHandle = new MPI_Request();
    MPI_Status status;
    int flag = 0;
    int result;

    clock_t timeOutWarning = -1;
    clock_t timeOutShutdown = -1;
    bool triggeredTimeoutWarning = false;

    if (exchangeOnlyAttributesMarkedWithParallelise) {
      result = MPI_Irecv(this, 1, Datatype, source, tag,
                         tarch::parallel::Node::getInstance().getCommunicator(),
                         sendRequestHandle);

    } else {
      result = MPI_Irecv(this, 1, FullDatatype, source, tag,
                         tarch::parallel::Node::getInstance().getCommunicator(),
                         sendRequestHandle);
    }
    if (result != MPI_SUCCESS) {
      std::ostringstream msg;
      msg << "failed to start to receive "
             "exahype::records::ADERDGCellDescriptionPacked from node "
          << source << ": " << tarch::parallel::MPIReturnValueToString(result);
      _log.error("receive(int)", msg.str());
    }

    result = MPI_Test(sendRequestHandle, &flag, &status);
    while (!flag) {
      if (timeOutWarning == -1)
        timeOutWarning =
            tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
      if (timeOutShutdown == -1)
        timeOutShutdown =
            tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
      result = MPI_Test(sendRequestHandle, &flag, &status);
      if (result != MPI_SUCCESS) {
        std::ostringstream msg;
        msg << "testing for finished receive task for "
               "exahype::records::ADERDGCellDescriptionPacked failed: "
            << tarch::parallel::MPIReturnValueToString(result);
        _log.error("receive(int)", msg.str());
      }

      // deadlock aspect
      if (tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
          (clock() > timeOutWarning) && (!triggeredTimeoutWarning)) {
        tarch::parallel::Node::getInstance().writeTimeOutWarning(
            "exahype::records::ADERDGCellDescriptionPacked", "receive(int)",
            source, tag, 1);
        triggeredTimeoutWarning = true;
      }
      if (tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
          (clock() > timeOutShutdown)) {
        tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "exahype::records::ADERDGCellDescriptionPacked", "receive(int)",
            source, tag, 1);
      }
      tarch::parallel::Node::getInstance().receiveDanglingMessages();
      usleep(communicateSleep);
    }

    delete sendRequestHandle;

#ifdef Debug
    _log.debug("receive(int,int)", "received " + toString());
#endif
  }
}

bool exahype::records::ADERDGCellDescriptionPacked::isMessageInQueue(
    int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
  MPI_Status status;
  int flag = 0;
  MPI_Iprobe(MPI_ANY_SOURCE, tag,
             tarch::parallel::Node::getInstance().getCommunicator(), &flag,
             &status);
  if (flag) {
    int messageCounter;
    if (exchangeOnlyAttributesMarkedWithParallelise) {
      MPI_Get_count(&status, Datatype, &messageCounter);
    } else {
      MPI_Get_count(&status, FullDatatype, &messageCounter);
    }
    return messageCounter > 0;
  } else
    return false;
}

#endif
