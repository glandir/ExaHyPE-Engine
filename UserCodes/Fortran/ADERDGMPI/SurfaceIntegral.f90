SUBROUTINE ADERSurfaceIntegral(lduh,lFbnd)
    USE typesDef
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)    :: lFbnd(nVar,nDOF(2),nDOF(3),6)            ! nonlinear flux tensor in each space-time DOF 
    REAL, INTENT(INOUT) :: lduh(nVar,nDOF(1),nDOF(2),nDOF(3))       ! spatial degrees of freedom 
    ! Local variables 
    INTEGER           :: i,j,k,l,iVar,sig  
    REAL              :: aux(d) 
    !
#ifdef LINEAR  
    sig = +1
#else
    sig = -1 
#endif
    !
    ! Now multiply the numerical fluxes on the surfaces with the test functions and compute the surface integrals 
    ! 
    ! x faces
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2) 
            aux = (/ 1., wGPN(j), wGPN(k) /)
            DO iVar = 1, nVar 
                lduh(iVar,:,j,k) = lduh(iVar,:,j,k) - PRODUCT(aux(1:nDim))/dx(1)*( lFbnd(iVar,j,k,2)*FRCoeff + sig*lFbnd(iVar,j,k,1)*FLCoeff )      ! left flux minus right flux 
            ENDDO                                                             
        ENDDO
    ENDDO 
    IF(nDim>=2) THEN
        ! y faces
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1) 
                aux = (/ 1., wGPN(i), wGPN(k) /) 
                DO iVar = 1, nVar 
                    lduh(iVar,i,:,k) = lduh(iVar,i,:,k) - PRODUCT(aux(1:nDim))/dx(2)*( lFbnd(iVar,i,k,4)*FRCoeff + sig*lFbnd(iVar,i,k,3)*FLCoeff )  ! left flux minus right flux  
                ENDDO                                                             
            ENDDO
        ENDDO 
    ENDIF
    IF(nDim>=3) THEN
        ! z faces
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1) 
                aux = (/ 1., wGPN(i), wGPN(j) /) 
                DO iVar = 1, nVar 
                    lduh(iVar,i,j,:) = lduh(iVar,i,j,:) - PRODUCT(aux(1:nDim))/dx(3)*( lFbnd(iVar,i,j,6)*FRCoeff + sig*lFbnd(iVar,i,j,5)*FLCoeff )  ! left flux minus right flux  
                ENDDO                                                             
            ENDDO
        ENDDO 
    ENDIF
    !
END SUBROUTINE ADERSurfaceIntegral 
    
    