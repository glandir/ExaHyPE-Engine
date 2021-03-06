SUBROUTINE AnalyseDG 
    USE typesDef
    IMPLICIT NONE
    INTEGER :: iElem, iVar, ii, jj, kk, iErr, iInfElem(nVar)  
    REAL    :: ldx(d), lx0(d), xGP(d), par0(nParam), aux(d) 
    REAL    :: u0(nVar), state(nVar), pstate(nVar), pexact(nVar)   
    REAL    :: L1norm(nVar), L2norm(nVar), Linfnorm(nVar), locError(nVar) 
    CHARACTER(LEN=10)  :: varName(nVar) 
    !
    IF(AnalyseType==0) THEN
        RETURN
    ENDIF  
    !
    L1norm = 0. 
    L2norm = 0. 
    Linfnorm = 0. 
    ! 
    ldx = dx 
    DO iElem = 1, nElem
      lx0 = x(:,tri(1,iElem))
      !  
      DO kk = 1, nDOF(3) 
       DO jj = 1, nDOF(2)
        DO ii = 1, nDOF(1) 
            xGP = lx0 + (/ xiGPN(ii), xiGPN(jj), xiGPN(kk) /)*ldx
            SELECT CASE(AnalyseType)
            CASE(1) ! Initialcondition
                CALL InitialField(u0,par0,xGP,0.0)
                CALL PDECons2Prim(pexact,u0,iErr)
            CASE(2) ! Time-shifted Initialcondition
                CALL InitialField(u0,par0,xGP,time)
                CALL PDECons2Prim(pexact,u0,iErr)
            END SELECT 
            state(:) = uh(:,ii,jj,kk,iElem)
            CALL PDECons2Prim(pstate,state,iErr)

            locError(:) = ABS(pexact(:)-pstate(:))
            aux = (/ wGPN(ii), wGPN(jj), wGPN(kk) /)
            ! 
            L1norm(:)   = L1norm(:) + locError(:)    * PRODUCT(aux(1:nDim))*PRODUCT(ldx(1:nDim)) 
            L2norm(:)   = L2norm(:) + locError(:)**2 * PRODUCT(aux(1:nDim))*PRODUCT(ldx(1:nDim))   
            !
            DO iVar = 1, nVar
                IF(locError(iVar).GT.Linfnorm(iVar)) THEN
                    Linfnorm(iVar) = locError(iVar)
                    iInfElem(iVar) = iElem
                ENDIF
            ENDDO

       ENDDO
      ENDDO
     ENDDO
     ! 
    ENDDO 
    !
    L2norm = SQRT(L2norm) 
    ! 
    CALL PDEVarName(VarName) 
    DO iVar = 1, nVar                                                         
        WRITE(*,*) '| Error analysis of variable ', TRIM(VarName(iVar)) 
        WRITE(*,*) '| ======================================'     !
        WRITE(*,*) '|   L1_norm   : ', L1norm(iVar)               !
        WRITE(*,*) '|   L2_norm   : ', L2norm(iVar)               !
        WRITE(*,*) '|   Linf_norm : ', Linfnorm(iVar), iInfElem(iVar)
        WRITE(*,*) '| ======================================'     !
    END DO
    !
END SUBROUTINE AnalyseDG 