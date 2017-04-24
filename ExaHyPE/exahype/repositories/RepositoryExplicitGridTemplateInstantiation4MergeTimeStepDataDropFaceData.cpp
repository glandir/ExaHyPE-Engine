#include "exahype/repositories/Repository.h"
#include "exahype/records/RepositoryState.h"

#include "exahype/State.h"
#include "exahype/Vertex.h"
#include "exahype/Cell.h"

#include "peano/grid/Grid.h"

#include "peano/stacks/CellArrayStack.h"
#include "peano/stacks/CellSTDStack.h"

#include "peano/stacks/VertexArrayStack.h"
#include "peano/stacks/VertexSTDStack.h"

 #include "exahype/adapters/MeshRefinement.h" 
 #include "exahype/adapters/MeshRefinementAndPlotGrid.h" 
 #include "exahype/adapters/PlotAugmentedAMRGrid.h" 
 #include "exahype/adapters/InitialConditionAndTimeStepSizeComputation.h" 
 #include "exahype/adapters/PredictionAndFusedTimeSteppingInitialisation.h" 
 #include "exahype/adapters/PredictionAndFusedTimeSteppingInitialisationAndPlot.h" 
 #include "exahype/adapters/PredictionAndFusedTimeSteppingInitialisationAndPlot2d.h" 
 #include "exahype/adapters/GridErasing.h" 
 #include "exahype/adapters/ADERDGTimeStep.h" 
 #include "exahype/adapters/PlotAndADERDGTimeStep.h" 
 #include "exahype/adapters/Reinitialisation.h" 
 #include "exahype/adapters/SolutionRecomputationAndTimeStepSizeComputation.h" 
 #include "exahype/adapters/NeighbourDataMerging.h" 
 #include "exahype/adapters/SolutionUpdate.h" 
 #include "exahype/adapters/TimeStepSizeComputation.h" 
 #include "exahype/adapters/Prediction.h" 
 #include "exahype/adapters/PredictionAndPlot.h" 
 #include "exahype/adapters/PredictionAndPlot2d.h" 
 #include "exahype/adapters/FinaliseMeshRefinementAndTimeStepSizeComputation.h" 
 #include "exahype/adapters/MergeTimeStepData.h" 
 #include "exahype/adapters/MergeTimeStepDataDropFaceData.h" 
 #include "exahype/adapters/FinaliseMeshRefinementAndSubcellSending.h" 


namespace peano {
  namespace grid {
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexArrayStack<exahype::Vertex> ,peano::stacks::CellArrayStack<exahype::Cell> ,exahype::adapters::MergeTimeStepDataDropFaceData>;
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexSTDStack<  exahype::Vertex> ,peano::stacks::CellSTDStack<  exahype::Cell> ,exahype::adapters::MergeTimeStepDataDropFaceData>;
  }
}

#include "peano/grid/Grid.cpph"
