#include "SHMMultipleRanksPerNodeStrategy.h"
#include "SHMController.h"
#include "SHMSharedMemoryBetweenTasks.h"
#include "SHMMacros.h"


#include <iostream>


shminvade::SHMMultipleRanksPerNodeStrategy::SHMMultipleRanksPerNodeStrategy() {
}


shminvade::SHMMultipleRanksPerNodeStrategy::~SHMMultipleRanksPerNodeStrategy() {
}


std::set<pid_t> shminvade::SHMMultipleRanksPerNodeStrategy::invadeThreads(int wantedNumberOfThreads) {
  std::set<pid_t> bookedCores;

  if (SHMController::getInstance()._switchedOn) {
    typedef tbb::spin_mutex CoreTraversalMutex;
    static CoreTraversalMutex invadeMutex;
    CoreTraversalMutex::scoped_lock lock(invadeMutex);

    for (auto p: SHMController::getInstance()._threads) {
      SHMController::ThreadTable::accessor a;
      SHMController::getInstance()._threads.find(a,p.first);
      if (
        wantedNumberOfThreads>0
        &&
        ( a->second.type==SHMController::ThreadType::NotOwned )
      ) {
        bool success = SHMSharedMemoryBetweenTasks::getInstance().tryToBookThreadForProcess(a->first);
        if (success) {
          wantedNumberOfThreads--;
          a->second.type = SHMController::ThreadType::ExclusivelyOwned;
          bookedCores.insert(a->first);
          #if SHM_INVADE_DEBUG>=4
          std::cout << SHM_DEBUG_PREFIX <<  "invade thread " << a->first << " (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
          #endif
        }
      }
    }

//    Geht er wirklich hier in die Knie?
// -> ja!


    lock.release();
    #if SHM_INVADE_DEBUG>=4
    std::cout << SHM_DEBUG_PREFIX <<  "invaded " << bookedCores.size() << " thread(s) in total with " << wantedNumberOfThreads << " open requests (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
    std::cout << SHM_DEBUG_PREFIX <<  "known thread-process association: " << SHMSharedMemoryBetweenTasks::getInstance().getThreadProcessAssociation() << std::endl;
    #endif
  }

  return bookedCores;
}


void shminvade::SHMMultipleRanksPerNodeStrategy::cleanUp() {
  SHMSharedMemoryBetweenTasks::getInstance().cleanUp();
}


void shminvade::SHMMultipleRanksPerNodeStrategy::retreat(const std::set<pid_t>& threadIds) {
  for (auto p: threadIds) {
    SHMSharedMemoryBetweenTasks::getInstance().freeThread(p);
    SHMController::getInstance().retreatFromCore(p);
  }
  #if SHM_INVADE_DEBUG>=4
  std::cout << SHM_DEBUG_PREFIX <<  "retreated from " << threadIds.size() << " thread(s) (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
  std::cout << SHM_DEBUG_PREFIX <<  "known thread-process association: " << SHMSharedMemoryBetweenTasks::getInstance().getThreadProcessAssociation() << std::endl;
  #endif
}
