# EulerFV with Delta #



## Prepare code ##

Type in

Toolkit/toolkit.sh Demonstrators/EulerFV-with-Delta/EulerFV-with-Delta.exahype

If you haven't updated Peano, you have to use

cd Submodules
./updateSubmodules.sh

before



## Translate ## 

Change into the directory Demonstrators/EulerFV-with-Delta and configure your
build. I enlist the requirements below and give examples in brackets if stuff
has to be done manually, i.e. if your system is not configured properly yet. 

- Ensure that TBB_INC is set properly 
  (export TBB_INC=-I/opt/intel/tbb/include)
- Ensure that TBB_SHLIB is set properly
  (export TBB_SHLIB="-L/opt/intel/tbb/lib/intel64/gcc4.7 -ltbb")
  Some systems also require -lpthread to be in the variable.
- Ensure the Delta sources path is added to the project's build flags
  (export PROJECT_CFLAGS=-I/home/.../git/Delta/src)
- Ensure the Delta library path is added to the project's linker flags
  (export PROJECT_LFLAGS=-L/home/.../git/Delta/src/lib -ldelta)
 
