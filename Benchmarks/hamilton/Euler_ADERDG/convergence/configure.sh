echo "Configure project for convergence test."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler_ADERDG/convergence/Euler_ADERDG.exahype )
mkdir convergence/results
rm *.o
