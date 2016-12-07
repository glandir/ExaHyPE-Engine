// Do not modify any part of this file. It will be overwritten throughout the 
// next pdt run.


#include "exahype/VertexOperations.h"
#include "peano/utils/Loop.h"
#include "peano/grid/Checkpoint.h"


exahype::VertexOperations::VertexOperations() { 
}



























































 tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D,int>  exahype::VertexOperations::readCellDescriptionsIndex(const peano::grid::VertexEnumerator& enumerator, const Vertex* const vertices)  { tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D,int> result; dfor2(x) tarch::la::slice(result,vertices[ enumerator(x) ]._vertexData.getCellDescriptionsIndex(),xScalar*TWO_POWER_D); enddforx return result; }











 tarch::la::Vector<TWO_POWER_D,int>  exahype::VertexOperations::readCellDescriptionsIndex(const Vertex& vertex)  { return vertex._vertexData.getCellDescriptionsIndex(); }














 void exahype::VertexOperations::writeCellDescriptionsIndex(const peano::grid::VertexEnumerator& enumerator, Vertex* const vertices, const tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D,int>& values) { dfor2(x) tarch::la::Vector<TWO_POWER_D,int> temp = tarch::la::slice<TWO_POWER_D>(values,xScalar*TWO_POWER_D); vertices[ enumerator(x) ]._vertexData.setCellDescriptionsIndex( temp ); enddforx }











 void exahype::VertexOperations::writeCellDescriptionsIndex(Vertex& vertex, const tarch::la::Vector<TWO_POWER_D,int>& values) { vertex._vertexData.setCellDescriptionsIndex(values ); }






















































