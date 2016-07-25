// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef _EXAHYPE_VERTEX_OPERATIONS_H_ 
#define _EXAHYPE_VERTEX_OPERATIONS_H_


#include "exahype/Vertex.h"
#include "peano/grid/Vertex.h"
#include "peano/grid/VertexEnumerator.h"
#include "peano/utils/Globals.h"


namespace exahype { 
      class VertexOperations; 
}


/**
 * This class comprises a collection of static operations on sets of vertices. 
 * It is generated by the Peano Development Toolkit (PDT) and will be 
 * overwritten by every rerun. Thus, do edit neither the class header nor its 
 * implementation. 
 */
class exahype::VertexOperations { 
  private: 
    /**
     * One should not instantiate this class as it is a collection of static operations.
     */
    VertexOperations();
  public:





























    




    static tarch::la::Vector<TWO_POWER_D,int> readCellDescriptionsIndex(const Vertex& vertex);






    static tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D,int> readCellDescriptionsIndex(const peano::grid::VertexEnumerator& enumerator, const Vertex* const vertices);









    static void writeCellDescriptionsIndex(const peano::grid::VertexEnumerator& enumerator, Vertex* const vertices, const tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D,int>& values);






    static void writeCellDescriptionsIndex(Vertex&  vertex, const tarch::la::Vector<TWO_POWER_D,int>& values);

































};


#endif
