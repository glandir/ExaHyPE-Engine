echo "Configure project for multicore scaling test (output)."
mkdir multicore/results
rm -rf *.o cfiles.mk ffiles.mk kernels cipofiles.mk
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler_FV/multicore/Euler_FV-output.exahype )
