// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STACK_H_ 
#define _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STACK_H_ 


#include "EulerFlow/repositories/Repository.h"
#include "EulerFlow/records/RepositoryState.h"

#include "EulerFlow/State.h"
#include "EulerFlow/Vertex.h"
#include "EulerFlow/Cell.h"

#include "peano/grid/Grid.h"
#include "peano/stacks/CellArrayStack.h"
#include "peano/stacks/VertexArrayStack.h"


 #include "EulerFlow/adapters/InitialGrid.h" 
 #include "EulerFlow/adapters/GridExport.h" 
 #include "EulerFlow/adapters/PatchInitialisation.h" 
 #include "EulerFlow/adapters/PatchInitialisationAndExport.h" 
 #include "EulerFlow/adapters/FaceDataExchange.h" 
 #include "EulerFlow/adapters/InitialConditionAndGlobalTimeStepComputation.h" 
 #include "EulerFlow/adapters/InitialConditionAndExportAndGlobalTimeStepComputation.h" 
 #include "EulerFlow/adapters/PredictorAndGlobalTimeStepComputation.h" 
 #include "EulerFlow/adapters/CorrectorAndPredictorAndGlobalTimeStepComputation.h" 
 #include "EulerFlow/adapters/CorrectorAndPredictorAndGlobalTimeStepComputationAndExport.h" 



namespace exahype {
      namespace repositories {
        class RepositoryArrayStack;  
      }
}


class exahype::repositories::RepositoryArrayStack: public exahype::repositories::Repository {
  private:
    static tarch::logging::Log _log;
  
    peano::geometry::Geometry& _geometry;
    
    typedef peano::stacks::CellArrayStack<exahype::Cell>       CellStack;
    typedef peano::stacks::VertexArrayStack<exahype::Vertex>   VertexStack;

    CellStack    _cellStack;
    VertexStack  _vertexStack;
    exahype::State          _solverState;
    peano::grid::RegularGridContainer<exahype::Vertex,exahype::Cell>  _regularGridContainer;
    peano::grid::TraversalOrderOnTopLevel                                         _traversalOrderOnTopLevel;

    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::InitialGrid> _gridWithInitialGrid;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::GridExport> _gridWithGridExport;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PatchInitialisation> _gridWithPatchInitialisation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PatchInitialisationAndExport> _gridWithPatchInitialisationAndExport;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::FaceDataExchange> _gridWithFaceDataExchange;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::InitialConditionAndGlobalTimeStepComputation> _gridWithInitialConditionAndGlobalTimeStepComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::InitialConditionAndExportAndGlobalTimeStepComputation> _gridWithInitialConditionAndExportAndGlobalTimeStepComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictorAndGlobalTimeStepComputation> _gridWithPredictorAndGlobalTimeStepComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::CorrectorAndPredictorAndGlobalTimeStepComputation> _gridWithCorrectorAndPredictorAndGlobalTimeStepComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::CorrectorAndPredictorAndGlobalTimeStepComputationAndExport> _gridWithCorrectorAndPredictorAndGlobalTimeStepComputationAndExport;

  
   exahype::records::RepositoryState               _repositoryState;
   
    tarch::timing::Measurement _measureInitialGridCPUTime;
    tarch::timing::Measurement _measureGridExportCPUTime;
    tarch::timing::Measurement _measurePatchInitialisationCPUTime;
    tarch::timing::Measurement _measurePatchInitialisationAndExportCPUTime;
    tarch::timing::Measurement _measureFaceDataExchangeCPUTime;
    tarch::timing::Measurement _measureInitialConditionAndGlobalTimeStepComputationCPUTime;
    tarch::timing::Measurement _measureInitialConditionAndExportAndGlobalTimeStepComputationCPUTime;
    tarch::timing::Measurement _measurePredictorAndGlobalTimeStepComputationCPUTime;
    tarch::timing::Measurement _measureCorrectorAndPredictorAndGlobalTimeStepComputationCPUTime;
    tarch::timing::Measurement _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCPUTime;

    tarch::timing::Measurement _measureInitialGridCalendarTime;
    tarch::timing::Measurement _measureGridExportCalendarTime;
    tarch::timing::Measurement _measurePatchInitialisationCalendarTime;
    tarch::timing::Measurement _measurePatchInitialisationAndExportCalendarTime;
    tarch::timing::Measurement _measureFaceDataExchangeCalendarTime;
    tarch::timing::Measurement _measureInitialConditionAndGlobalTimeStepComputationCalendarTime;
    tarch::timing::Measurement _measureInitialConditionAndExportAndGlobalTimeStepComputationCalendarTime;
    tarch::timing::Measurement _measurePredictorAndGlobalTimeStepComputationCalendarTime;
    tarch::timing::Measurement _measureCorrectorAndPredictorAndGlobalTimeStepComputationCalendarTime;
    tarch::timing::Measurement _measureCorrectorAndPredictorAndGlobalTimeStepComputationAndExportCalendarTime;


  public:
    RepositoryArrayStack(
      peano::geometry::Geometry&                   geometry,
      const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
      const tarch::la::Vector<DIMENSIONS,double>&  computationalDomainOffset,
      int                                          maximumSizeOfCellInOutStack,
      int                                          maximumSizeOfVertexInOutStack,
      int                                          maximumSizeOfVertexTemporaryStack
    );
    
    /**
     * Parallel Constructor
     *
     * Used in parallel mode only where the size of the domain is not known 
     * when the type of repository is determined.  
     */
    RepositoryArrayStack(
      peano::geometry::Geometry&                   geometry,
      int                                          maximumSizeOfCellInOutStack,
      int                                          maximumSizeOfVertexInOutStack,
      int                                          maximumSizeOfVertexTemporaryStack
    );
    
    virtual ~RepositoryArrayStack();

    virtual void restart(
      const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
      const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
      int                                          domainLevel,
      const tarch::la::Vector<DIMENSIONS,int>&     positionOfCentralElementWithRespectToCoarserRemoteLevel
    );
         
    virtual void terminate();
        
    virtual exahype::State& getState();
    virtual const exahype::State& getState() const;

    virtual void iterate(int numberOfIterations=1, bool exchangeBoundaryVertices=true);
    
    virtual void writeCheckpoint(peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> * const checkpoint); 
    virtual void readCheckpoint( peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> const * const checkpoint );
    virtual peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>* createEmptyCheckpoint(); 

    virtual void switchToInitialGrid();    
    virtual void switchToGridExport();    
    virtual void switchToPatchInitialisation();    
    virtual void switchToPatchInitialisationAndExport();    
    virtual void switchToFaceDataExchange();    
    virtual void switchToInitialConditionAndGlobalTimeStepComputation();    
    virtual void switchToInitialConditionAndExportAndGlobalTimeStepComputation();    
    virtual void switchToPredictorAndGlobalTimeStepComputation();    
    virtual void switchToCorrectorAndPredictorAndGlobalTimeStepComputation();    
    virtual void switchToCorrectorAndPredictorAndGlobalTimeStepComputationAndExport();    

    virtual bool isActiveAdapterInitialGrid() const;
    virtual bool isActiveAdapterGridExport() const;
    virtual bool isActiveAdapterPatchInitialisation() const;
    virtual bool isActiveAdapterPatchInitialisationAndExport() const;
    virtual bool isActiveAdapterFaceDataExchange() const;
    virtual bool isActiveAdapterInitialConditionAndGlobalTimeStepComputation() const;
    virtual bool isActiveAdapterInitialConditionAndExportAndGlobalTimeStepComputation() const;
    virtual bool isActiveAdapterPredictorAndGlobalTimeStepComputation() const;
    virtual bool isActiveAdapterCorrectorAndPredictorAndGlobalTimeStepComputation() const;
    virtual bool isActiveAdapterCorrectorAndPredictorAndGlobalTimeStepComputationAndExport() const;

     
    #ifdef Parallel
    virtual ContinueCommand continueToIterate();
    virtual void runGlobalStep();
    #endif

    virtual void setMaximumMemoryFootprintForTemporaryRegularGrids(double value);
    virtual void logIterationStatistics() const;
    virtual void clearIterationStatistics();
};


#endif
