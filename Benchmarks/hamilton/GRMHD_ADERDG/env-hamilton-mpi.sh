source /etc/profile.d/modules.sh

module purge
module load java
module load slurm
module load intel/xe_2017.2
module load intelmpi/intel/2017.2
module load gcc

export TBB_SHLIB="-L/ddn/apps/Cluster-Apps/intel/xe_2017.2/tbb/lib/intel64/gcc4.7 -ltbb"

export I_MPI_FABRICS="shm:dapl"

export EXAHYPE_CC="mpicc -DnoParallelExchangePackedRecordsAtBoundary -DnoParallelExchangePackedRecordsBetweenMasterAndWorker -DnoParallelExchangePackedRecordsInHeaps -DnoParallelExchangePackedRecordsThroughoutJoinsAndForks"
export COMPILER_LFLAGS="-lgsl -lgslcblas"

export MODE=Release
export COMPILER=Intel
export DISTRIBUTEDMEM=MPI
export ARCHITECTURE=hsw
export GPROF=off
