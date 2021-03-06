SUBROUTINE BoundaryConditions 
    USE typesDef
    ! Local variables 
    REAL :: j,k,iFace
    REAL :: Qbc(nVar),Fbc(nVar,d),Vbc(nVar) 
    !
    ! Fix boundary data  
    Vbc = (/ 1., 0., 0., 0., 1. /)   ! primitive variables 
    CALL PDEPrim2Cons(qBC,Vbc)        ! convert into conservative variables    
    !
    DO iFace = 1, nFace
        ! Here, we need to take care of the boundary conditions 
        ! For the moment, we use simple extrapolation (copy from inside the domain) 
        IF(Face(iFace)%Left.EQ.0) THEN
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2) 
                    Face(iFace)%qL(:,j,k) = qbc  
                    CALL PDEFlux(Fbc,qbc) 
                    Face(iFace)%FL(:,j,k) = MATMUL(Fbc(:,:), Face(iFace)%nv) 
                ENDDO
            ENDDO 
        ENDIF
        IF(Face(iFace)%Right.EQ.0) THEN 
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2) 
                    Face(iFace)%qR(:,j,k) = qbc  
                    CALL PDEFlux(Fbc,qbc) 
                    Face(iFace)%FR(:,j,k) = MATMUL(Fbc(:,:), Face(iFace)%nv) 
                ENDDO
            ENDDO 
        ENDIF            
    ENDDO    
    ! 
END SUBROUTINE BoundaryConditions 
    
    