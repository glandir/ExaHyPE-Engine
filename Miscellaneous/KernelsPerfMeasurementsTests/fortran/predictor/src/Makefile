# *********************************************************************
#
#  Makefile for ADERDG3D 
#  
# *********************************************************************
#
.SUFFIXES : .o .f90 .f .c .F90
.KEEP_STATE:
#
EXE=ADERDG  
NVAR=9
ORDER=3
# *********************************************************************

F90 = ifort

#------------------------------------------------------------------------------
# LIBRARIES
#------------------------------------------------------------------------------

#-LIBDIR      = ./lib
#-LTEC   = -ltecio

# Settings for SuperMUC
# Run "module load tecplot" before compiling!

#LTEC   = /lrz/sys/applications/tecplot/2012r3/tecplot360/lib/libtecio.a

#LFLAGS = $(LTEC) -cxxlib

# COMPILATION OPTIMIZATION FLAGS for ifort
F90FLAGS = -c -O2 -r8 -align -fpp -xCORE-AVX2  
#F90FLAGS = -c -fast -r8 -align -fpp -xCORE-AVX2

# PREPROCESSOR DEFINITIONS 
FCFLAGS = -DORDER__=$(ORDER) -DNVAR__=$(NVAR) 

#
.SUFFIXES: .f90 .o .f77 .F90

.f90.o:
	@echo "   "
	@echo " Compiling Source File: " $<
	@echo " ---------------------- ===================="
	$(F90) $(F90FLAGS) $(FCFLAGS) $(INCLUDES) $<

.F90.o:
	@echo "   "
	@echo " Compiling Source File: " $<
	@echo " ---------------------- ===================="
	$(F90) $(F90FLAGS) $(FCFLAGS) $(INCLUDES) $<

.f.o:
	@echo "   "
	@echo " Compiling Source File: " $<
	@echo " ---------------------- ===================="
	$(F90) $(F90FLAGS) $(FCFLAGS) $(INCLUDES) $<

# Fortran files with no optimization
NO_OPT_OBJ = 

# Compiled object files

F90_OBJ = gauleg.o                  \
          GAUSS.o                   \
          typesDef.o                \
          DGBasis.o                 \
          PDE.o                     \
          Init.o                    \
          SpaceTimePredictor.o      \
		  predictorTest.o           \
          ADERDG3D.o                

OBJS       = $(F90_OBJ)
FOR_FILES  = $(OBJS:.o=.f90)

Executable: $(EXE)

$(EXE): Makefile $(OBJS)
	@echo "   "
	@echo " Now Linking..."
	@echo " --------------"
	$(F90) $(OBJS) $(LFLAGS)  -o $(EXE)
	@echo "   "
	@echo " Compiled ! "
	@echo " ---------- "

#
clean:
	rm -f core* *.M *.o $(EXE) *.kmo *.mod 


