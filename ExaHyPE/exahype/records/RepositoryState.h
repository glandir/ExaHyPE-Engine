#ifndef _EXAHYPE_RECORDS_REPOSITORYSTATE_H
#define _EXAHYPE_RECORDS_REPOSITORYSTATE_H

#include "tarch/compiler/CompilerSpecificSettings.h"
#include "peano/utils/PeanoOptimisations.h"
#ifdef Parallel
	#include "tarch/parallel/Node.h"
#endif
#ifdef Parallel
	#include <mpi.h>
#endif
#include "tarch/logging/Log.h"
#include "tarch/la/Vector.h"
#include <bitset>
#include <complex>
#include <string>
#include <iostream>

namespace exahype {
   namespace records {
      class RepositoryState;
      class RepositoryStatePacked;
   }
}

/**
 * @author This class is generated by DaStGen
 * 		   DataStructureGenerator (DaStGen)
 * 		   2007-2009 Wolfgang Eckhardt
 * 		   2012      Tobias Weinzierl
 *
 * 		   build date: 09-02-2014 14:40
 *
 * @date   05/03/2018 14:11
 */
class exahype::records::RepositoryState { 
   
   public:
      
      typedef exahype::records::RepositoryStatePacked Packed;
      
      enum Action {
         WriteCheckpoint = 0, ReadCheckpoint = 1, Terminate = 2, RunOnAllNodes = 3, UseAdapterMeshRefinement = 4, UseAdapterMeshRefinementAndPlotTree = 5, UseAdapterFinaliseMeshRefinement = 6, UseAdapterFinaliseMeshRefinementOrLocalRollback = 7, UseAdapterFusedTimeStep = 8, UseAdapterPredictionRerun = 9, UseAdapterBroadcastAndDropNeighbourMessages = 10, UseAdapterLimiterStatusSpreading = 11, UseAdapterPredictionOrLocalRecomputation = 12, UseAdapterGlobalRollback = 13, UseAdapterBroadcastAndMergeNeighbours = 14, UseAdapterUpdateAndReduce = 15, UseAdapterPrediction = 16, NumberOfAdapters = 17
      };
      
      struct PersistentRecords {
         Action _action;
         int _numberOfIterations;
         bool _exchangeBoundaryVertices;
         /**
          * Generated
          */
         PersistentRecords();
         
         /**
          * Generated
          */
         PersistentRecords(const Action& action, const int& numberOfIterations, const bool& exchangeBoundaryVertices);
         
         
         inline Action getAction() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _action;
         }
         
         
         
         inline void setAction(const Action& action) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _action = action;
         }
         
         
         
         inline int getNumberOfIterations() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _numberOfIterations;
         }
         
         
         
         inline void setNumberOfIterations(const int& numberOfIterations) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _numberOfIterations = numberOfIterations;
         }
         
         
         
         inline bool getExchangeBoundaryVertices() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _exchangeBoundaryVertices;
         }
         
         
         
         inline void setExchangeBoundaryVertices(const bool& exchangeBoundaryVertices) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _exchangeBoundaryVertices = exchangeBoundaryVertices;
         }
         
         
         
      };
      private: 
         PersistentRecords _persistentRecords;
         
      public:
         /**
          * Generated
          */
         RepositoryState();
         
         /**
          * Generated
          */
         RepositoryState(const PersistentRecords& persistentRecords);
         
         /**
          * Generated
          */
         RepositoryState(const Action& action, const int& numberOfIterations, const bool& exchangeBoundaryVertices);
         
         /**
          * Generated
          */
         virtual ~RepositoryState();
         
         
         inline Action getAction() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _persistentRecords._action;
         }
         
         
         
         inline void setAction(const Action& action) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _persistentRecords._action = action;
         }
         
         
         
         inline int getNumberOfIterations() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _persistentRecords._numberOfIterations;
         }
         
         
         
         inline void setNumberOfIterations(const int& numberOfIterations) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _persistentRecords._numberOfIterations = numberOfIterations;
         }
         
         
         
         inline bool getExchangeBoundaryVertices() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _persistentRecords._exchangeBoundaryVertices;
         }
         
         
         
         inline void setExchangeBoundaryVertices(const bool& exchangeBoundaryVertices) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _persistentRecords._exchangeBoundaryVertices = exchangeBoundaryVertices;
         }
         
         
         /**
          * Generated
          */
         static std::string toString(const Action& param);
         
         /**
          * Generated
          */
         static std::string getActionMapping();
         
         /**
          * Generated
          */
         std::string toString() const;
         
         /**
          * Generated
          */
         void toString(std::ostream& out) const;
         
         
         PersistentRecords getPersistentRecords() const;
         /**
          * Generated
          */
         RepositoryStatePacked convert() const;
         
         
      #ifdef Parallel
         protected:
            static tarch::logging::Log _log;
            
            int _senderDestinationRank;
            
         public:
            
            /**
             * Global that represents the mpi datatype.
             * There are two variants: Datatype identifies only those attributes marked with
             * parallelise. FullDatatype instead identifies the whole record with all fields.
             */
            static MPI_Datatype Datatype;
            static MPI_Datatype FullDatatype;
            
            /**
             * Initializes the data type for the mpi operations. Has to be called
             * before the very first send or receive operation is called.
             */
            static void initDatatype();
            
            static void shutdownDatatype();
            
            /**
             * @param communicateSleep -1 Data exchange through blocking mpi
             * @param communicateSleep  0 Data exchange through non-blocking mpi, i.e. pending messages are received via polling until MPI_Test succeeds
             * @param communicateSleep >0 Same as 0 but in addition, each unsuccessful MPI_Test is follows by an usleep
             */
            void send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep);
            
            void receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep);
            
            static bool isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise);
            
            int getSenderRank() const;
            #endif
   
};

#ifndef DaStGenPackedPadding
  #define DaStGenPackedPadding 1      // 32 bit version
  // #define DaStGenPackedPadding 2   // 64 bit version
#endif


#ifdef PackedRecords
   #pragma pack (push, DaStGenPackedPadding)
#endif

/**
 * @author This class is generated by DaStGen
 * 		   DataStructureGenerator (DaStGen)
 * 		   2007-2009 Wolfgang Eckhardt
 * 		   2012      Tobias Weinzierl
 *
 * 		   build date: 09-02-2014 14:40
 *
 * @date   05/03/2018 14:11
 */
class exahype::records::RepositoryStatePacked { 
   
   public:
      
      typedef exahype::records::RepositoryState::Action Action;
      
      struct PersistentRecords {
         Action _action;
         int _numberOfIterations;
         bool _exchangeBoundaryVertices;
         /**
          * Generated
          */
         PersistentRecords();
         
         /**
          * Generated
          */
         PersistentRecords(const Action& action, const int& numberOfIterations, const bool& exchangeBoundaryVertices);
         
         
         inline Action getAction() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _action;
         }
         
         
         
         inline void setAction(const Action& action) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _action = action;
         }
         
         
         
         inline int getNumberOfIterations() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _numberOfIterations;
         }
         
         
         
         inline void setNumberOfIterations(const int& numberOfIterations) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _numberOfIterations = numberOfIterations;
         }
         
         
         
         inline bool getExchangeBoundaryVertices() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _exchangeBoundaryVertices;
         }
         
         
         
         inline void setExchangeBoundaryVertices(const bool& exchangeBoundaryVertices) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _exchangeBoundaryVertices = exchangeBoundaryVertices;
         }
         
         
         
      };
      private: 
         PersistentRecords _persistentRecords;
         
      public:
         /**
          * Generated
          */
         RepositoryStatePacked();
         
         /**
          * Generated
          */
         RepositoryStatePacked(const PersistentRecords& persistentRecords);
         
         /**
          * Generated
          */
         RepositoryStatePacked(const Action& action, const int& numberOfIterations, const bool& exchangeBoundaryVertices);
         
         /**
          * Generated
          */
         virtual ~RepositoryStatePacked();
         
         
         inline Action getAction() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _persistentRecords._action;
         }
         
         
         
         inline void setAction(const Action& action) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _persistentRecords._action = action;
         }
         
         
         
         inline int getNumberOfIterations() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _persistentRecords._numberOfIterations;
         }
         
         
         
         inline void setNumberOfIterations(const int& numberOfIterations) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _persistentRecords._numberOfIterations = numberOfIterations;
         }
         
         
         
         inline bool getExchangeBoundaryVertices() const 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            return _persistentRecords._exchangeBoundaryVertices;
         }
         
         
         
         inline void setExchangeBoundaryVertices(const bool& exchangeBoundaryVertices) 
 #ifdef UseManualInlining
 __attribute__((always_inline))
 #endif 
 {
            _persistentRecords._exchangeBoundaryVertices = exchangeBoundaryVertices;
         }
         
         
         /**
          * Generated
          */
         static std::string toString(const Action& param);
         
         /**
          * Generated
          */
         static std::string getActionMapping();
         
         /**
          * Generated
          */
         std::string toString() const;
         
         /**
          * Generated
          */
         void toString(std::ostream& out) const;
         
         
         PersistentRecords getPersistentRecords() const;
         /**
          * Generated
          */
         RepositoryState convert() const;
         
         
      #ifdef Parallel
         protected:
            static tarch::logging::Log _log;
            
            int _senderDestinationRank;
            
         public:
            
            /**
             * Global that represents the mpi datatype.
             * There are two variants: Datatype identifies only those attributes marked with
             * parallelise. FullDatatype instead identifies the whole record with all fields.
             */
            static MPI_Datatype Datatype;
            static MPI_Datatype FullDatatype;
            
            /**
             * Initializes the data type for the mpi operations. Has to be called
             * before the very first send or receive operation is called.
             */
            static void initDatatype();
            
            static void shutdownDatatype();
            
            /**
             * @param communicateSleep -1 Data exchange through blocking mpi
             * @param communicateSleep  0 Data exchange through non-blocking mpi, i.e. pending messages are received via polling until MPI_Test succeeds
             * @param communicateSleep >0 Same as 0 but in addition, each unsuccessful MPI_Test is follows by an usleep
             */
            void send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep);
            
            void receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep);
            
            static bool isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise);
            
            int getSenderRank() const;
            #endif
   
};

#ifdef PackedRecords
#pragma pack (pop)
#endif


#endif

