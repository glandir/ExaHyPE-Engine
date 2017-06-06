// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _MPIBALANCING_SFC_DIFFUSION_NODE_POOL_STRATEGY_H_
#define _MPIBALANCING_SFC_DIFFUSION_NODE_POOL_STRATEGY_H_


#ifdef Parallel
#include <mpi.h>
#endif

#include "tarch/parallel/NodePoolStrategy.h"
#include "tarch/logging/Log.h"

#include <vector>


namespace mpibalancing {
  class SFCDiffusionNodePoolStrategy;
}


/**
 * SFC Diffusion Node Pool Strategy
 *
 * <h2> User pitfalls </h2>
 * The strategy tries to book solely the primary MPI nodes. Therefore,
 * statements alike

  assertion( tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()==0 );

 * are doomed to fail.
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.1 $
 */
class mpibalancing::SFCDiffusionNodePoolStrategy: public tarch::parallel::NodePoolStrategy {
  private:
    /**
     * Copy from FCFS but enriched by counter how many rank have already
     * requested an update.
     *
     * @author Tobias Weinzierl
     * @version $Revision: 1.1 $
     */
    class NodePoolListEntry {
      public:
        /**
         * Represents the state of the worker, i.e. whether it is idle or busy.
         */
        enum class State {
          Undef,
          IdlePrimaryNode,
          IdleSecondaryNode,
          WorkingPrimaryNode,
          WorkingSecondaryNode
        };

      private:
        /**
         * Holds the state of the process.
         */
        State             _state;

        /**
         * Machine name
         *
         * Should be const but then using the operator= becomes a mess.
         */
        std::string       _name;

      public:
        /**
         * Construct one entry. By default this entry corresponds to an idle worker.
         */
        NodePoolListEntry( const std::string& name, bool isPrimaryNode );

        /**
         * I need a default constructor for some resorting, but it is not
         * available from outside.
         */
        NodePoolListEntry();

        virtual ~NodePoolListEntry();

        /**
         * Activates the node. Precondition: Node is idle. Thus, the local min level
         * is overwritten by the argument level and the state is set to working.
         */
        void activate();

        /**
         * The state is switched to idle.
         */
        void deActivate();

        bool isRegistered() const;

        /**
         * @return Rank of process.
         */
        int getRank() const;

        /**
         * @return Name of the node the process is running on.
         */
        std::string getNodeName() const;

        /**
         * Create string representation.
         */
        void toString(std::ostream& out) const;

        /**
         * Return string representation.
         */
        std::string toString() const;

        bool isIdlePrimaryNode() const;
        bool isIdleSecondaryNode() const;
    };

    enum class NodePoolState {
      /**
       * Standard mode. As soon as we run out of primary nodes, we switch it to
       * DeployingAlsoSecondaryNodes.
       */
      DeployingIdlePrimaryNodes,
      /**
       * If this flag is set, we also deploy secondary nodes. However, as soon
       * as our request queue is empty, i.e. as soon as we have answered one
       * batch of requests, we switch into NoNodesLeft.
       */
      DeployingAlsoSecondaryNodes,
      /**
       * We do not hand out any nodes anymore.
       */
      NoNodesLeft
    };


    /**
     * Is ordered along ranks.
     */
    typedef std::vector<NodePoolListEntry>   NodeContainer;

    /**
     * Logging Device
     */
    static tarch::logging::Log _log;

    /**
     * Tag on which the node pool works
     */
    int _tag;

    /**
     * The list the list of active nodes. Every entry corresponds to one node.
     * If the entry is set, the node is working already and the server is not
     * allowed to deploy another job on this node. If the entry isn't set, there
     * is a job request message in the queue and the server is allowed to send
     * a job. Therefore, in the beginning, all the entries are set. For the very
     * first entry, corresponding to the server node, the invariant holds, that
     * this entry is set always.
     */
    NodeContainer _nodes;

    double        _waitTimeOut;

    /**
     * is not used at the moment, but should be used in first mode.
     */
    const int     _mpiRanksPerNode;

    const int     _primaryMPIRanksPerNode;

    NodePoolState _nodePoolState;

    bool continueToFillRequestQueue(int queueSize) const;

    int getNumberOfIdlePrimaryNodes() const;

    /**
     * Take queue and return re-sorted queue. Input and output queue contain
     * the same elements in different ordering though.
     */
    RequestQueue sortRequestQueue( const RequestQueue& queue );

    bool isPrimaryMPIRank(int rank) const;
    bool isFirstOrLastPrimaryMPIRankOnANode(int rank) const;
  public:
  /**
   * Constructor
   *
   * @param mpiRanksPerNode       Number of ranks per node.
   * @param waitTimeOutSec        How long shall the node wait for more
   *   messages dropping in before it starts to answer them.
   */
    SFCDiffusionNodePoolStrategy(int mpiRanksPerNode, int primaryMPIRanksPerNode, double waitTimeOutSec = 1e-5);
    virtual ~SFCDiffusionNodePoolStrategy();

    void setNodePoolTag(int tag) override;

    /**
     * Here's the magic to get the balancing right.
     */
    tarch::parallel::messages::WorkerRequestMessage extractElementFromRequestQueue(RequestQueue& queue) override;
    void fillWorkerRequestQueue(RequestQueue& queue) override;
    void addNode(const tarch::parallel::messages::RegisterAtNodePoolMessage& node ) override;
    void removeNode( int rank ) override;

    /**
     * This routine is used by the NodePool to compare to the total number of
     * ranks. The calling operation therefrom derives whether all nodes have
     * successfully registered. So we may not only return the number of primary
     * ranks, but we have to return the total number of ranks.
     */
    int getNumberOfIdleNodes() const override;
    void setNodeIdle( int rank ) override;
    int reserveNode(int forMaster) override;
    void reserveParticularNode(int rank) override;
    bool isRegisteredNode(int rank) const override;
    bool isIdleNode(int rank) const override;
    int getNumberOfRegisteredNodes() const override;
    std::string toString() const override;
    bool hasIdleNode(int forMaster) const override;
};

#endif
