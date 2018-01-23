module purge
module load python
module load java
module load slurm
module load intel/xe_2017.2
module load intelmpi/intel/2017.2
module load gcc/4.9.1

export TBB_SHLIB="-L/ddn/apps/Cluster-Apps/intel/xe_2017.2/tbb/lib/intel64/gcc4.7 -ltbb"

export EXAHYPE_CC="mpicc"
export COMPILER_CFLAGS=" DnoParallelExchangePackedRecordsAtBoundary -DnoParallelExchangePackedRecordsBetweenMasterAndWorker -DnoParallelExchangePackedRecordsInHeaps -DnoParallelExchangePackedRecordsThroughoutJoinsAndForks "

export MODE=Release
export COMPILER=Intel
export DISTRIBUTEDMEM=MPI
export GPROF=off

# optimised kernels
export USE_IPO=on
