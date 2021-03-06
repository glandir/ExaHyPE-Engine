! DIM Parameters
  
  MODULE Parameters 
    IMPLICIT NONE  
    PUBLIC  

    ! This typesDef.f90 is a relict still from the old Fortran interface in Exahype.
    ! However, the following two variables are needed in the Fortran code. They
    ! should provided by the glue code generator in an appropriate way.

    ! If you modify the SRHD.exahype, please do the following mapping by hand:
    !
    ! solver ADER-DG SRHDSolver / unknowns   ->  goes to ->  nVar
    ! computational-domain / dimension       ->  goes to ->  nDim
    !

    ! Here, we obtain DIMENSIONS from the C++ preprocessor
#if defined(Dim3)
    INTEGER, PARAMETER             :: nDim = 3                   ! The number of space dimensions
	!CHARACTER(LEN=10), PARAMETER	:: ICFlag='CGeom'
	CHARACTER(LEN=10), PARAMETER	:: ICFlag='LOH1'
	!CHARACTER(LEN=10), PARAMETER	:: ICFlag='PlaneWave'
#elif defined(Dim2)
    INTEGER, PARAMETER             :: nDim = 2                   ! The number of space dimensions
	CHARACTER(LEN=10), PARAMETER	:: ICFlag='PlaneWave'
#endif
	!CHARACTER(LEN=10), PARAMETER				:: ICFlag='CGeom'
    INTEGER, PARAMETER             :: nAux = 0
	INTEGER, PARAMETER             :: nVar = 14                           ! The number of variables of the PDE system  
    INTEGER, PARAMETER             :: nLin = 12
	INTEGER, PARAMETER 				:: nPar=0
	!REAL, PARAMETER             :: epsilon1 =1.e-3
	REAL, PARAMETER             :: epsilon1 =1.e-20
  END MODULE Parameters  
