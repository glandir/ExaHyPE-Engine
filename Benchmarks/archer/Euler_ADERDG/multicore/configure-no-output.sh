echo "Configure project for multicore scaling test (no output)."
mkdir multicore/results
rm -rf *.o cfiles.mk ffiles.mk kernels cipofiles.mk
( cd ../../../ && java -XX:MaxHeapSize=512m -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/archer/Euler_ADERDG/multicore/Euler_ADERDG-no-output.exahype )
