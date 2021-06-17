MODULE GMRESUTIL
  ! subroutines for solving system of linear equations by GMRES iterative method
  ! Based on Saad, 2007 "Iterative Methods" book, Ch. 6, algorithm 6.9 and 6.11
  IMPLICIT NONE

  LOGICAL :: GMRESVERBOSE = .TRUE.
CONTAINS


  SUBROUTINE RUNZGMRES(AMAT,BVEC,X0,N,GMM,NTRY,TOL,XVEC,RINFO)
    ! see Fraysse, 1997 for implementation details
    ! run iterative GMRES method from CERFACS code
    ! to solve linear equation AMAT*X = BVEC
    ! in complex variables
    ! X0 = initial guess
    ! N = problem dimension
    ! GMM = restart parameter, controls memory required for matrix H
    ! ITERM = maximal number of iterations
    IMPLICIT NONE
    COMPLEX*16, INTENT(IN) :: AMAT(N,N),BVEC(N),X0(N)
    INTEGER, INTENT(IN) :: N, GMM, NTRY
    DOUBLE PRECISION, INTENT(IN) :: TOL
    COMPLEX*16, INTENT(OUT) :: XVEC(N)
    DOUBLE PRECISION, INTENT(OUT) :: RINFO(2)
    INTEGER :: LWORK, M
    COMPLEX*16, ALLOCATABLE :: WORK(:)
    COMPLEX*16 :: WORK0
    INTEGER :: IRC(5), ICNTL(7), INFO(3), J, REVCOM, COLX, COLY, COLZ, NBSCAL
    DOUBLE PRECISION :: CNTL(5)

    INFO = -10

    ! default GMM = N
    IF (GMM.LE.0) THEN
       M = N
    ELSE
       M = GMM
    ENDIF

    IF (M.GT.N) THEN
       PRINT*, 'ERROR IN RUNZGMRES: memory parameter M must be at most N', N, GMM
       STOP 1
    ENDIF

    LWORK = M*M+M*(N+5)+6*N+M+10
    !LWORK = 0
    
    OPEN(UNIT=33,FILE='gmresconv.out')

    IRC = 0
    ICNTL = 0
    CNTL = 0D0

    !print*, 'testx0:', ntry, LWORK

    ! control parameters for DRIVE_ZGMRES
    ! Initialize control parameters to default values
    CALL INIT_ZGMRES(ICNTL,CNTL)
    ! ignore warning messages
    ICNTL(2) = 0
    ! save convergence history to file
    ICNTL(3) = 33
    ! No preconditioning
    ICNTL(4) = 0
    ICNTL(5) = 1 ! user-supplied initial guess
    ICNTL(7) = NTRY ! max number of iterations

    CNTL(1) = TOL ! tolerance for convergence

    IRC(1) = -1

    ! ! find optimal size of work array
    ! CALL DRIVE_ZGMRES(N,N,M,LWORK,WORK0,IRC,ICNTL,CNTL,INFO,RINFO)
    ! IF (INFO(1).EQ.-3) THEN
    !    PRINT*, 'TESTX1:', INFO
    !    LWORK = INFO(2)      
    !    ALLOCATE(WORK(LWORK))
    ! ELSE
    !    PRINT*, 'ERROR IN RUNZGMRES: failed to get optimal array size.', INFO, LWORK
    !    STOP 1
    ! ENDIF

    ! set up work array
    ALLOCATE(WORK(LWORK))
    WORK = 0D0
    ! initial guess
    WORK(1:N) = X0    
    ! right-hand side
    WORK(N+1:2*N) = BVEC

    DO WHILE (IRC(1).NE.0)
       CALL DRIVE_ZGMRES(N,N,M,LWORK,WORK,IRC,ICNTL,CNTL,INFO,RINFO)          
      ! PRINT*, 'TESTX1:', IRC(1), info(1)
       ! REVCOM: what to do next (0 means done)
       ! 0: done
       ! 1: do matrix-vector product z<- Ax
       ! 2: do left preconditioning z<- M1^-1 *x
       ! 3: do right preconditioning z<-M2^-1*x
       ! 4: do NBSCAL scalar products z<-X*y
       REVCOM = IRC(1) ! what to do next (0 means done)
       ! index to read off x
       COLX = IRC(2) !
       ! index to read off y
       COLY = IRC(3)
       ! index to write Z
       COLZ = IRC(4)
       ! number of scalar products
       NBSCAL = IRC(5)

       IF (REVCOM.EQ.1) THEN
          ! do matrix vector product WORK(COLZ) <-- A*WORK(COLX)
          CALL ZGEMV('N',N,N,DCMPLX(1D0),AMAT,N,WORK(COLX),1,DCMPLX(0D0),WORK(COLZ),1)
       ELSEIF (REVCOM.EQ.4) THEN
          CALL ZGEMV('C',N,NBSCAL,DCMPLX(1D0),WORK(COLX),N,WORK(COLY),1,DCMPLX(0D0),WORK(COLZ),1)
       ELSE IF (REVCOM.NE.0) THEN
          PRINT*, 'ERROR IN RUNZGMRES: REVCOM=',revcom
          STOP 1
       ENDIF
     
    END DO

    ! Check exit conditions   
    IF (INFO(1).EQ.0) THEN
       ! PRINT*, 'normal exit'
       ! print*, 'convereged after ', info(2), ' iterations'
       ! print*, 'optimal and actual workspace size', info(3), lwork
    ELSEIF (INFO(1).EQ.-4) THEN
       PRINT*, 'ERROR IN RUNZGMRES: No convergence after ', icntl(7), ' iteractions'
       STOP 1
    ELSE
       PRINT*, 'ERROR IN RUNZGMRES: ', INFO
       print*, 'n, m:', N,M
       STOP 1
    ENDIF

    XVEC = WORK(1:N)

    ! PRINT*, 'TESTX6:'
    ! DO J = 1,N
    !    PRINT*, J, XVEC(J)
    ! ENDDO

    DEALLOCATE(WORK)
    CLOSE(33)
  END SUBROUTINE RUNZGMRES


  ! SUBROUTINE RUNDGMRES(AMAT,BVEC,X0,N,GMM,NTRY,TOL,XVEC,RINFO)
  !   ! see Fraysse, 1997 for implementation details
  !   ! run iterative GMRES method from CERFACS code
  !   ! to solve linear equation AMAT*X = BVEC
  !   ! in complex variables
  !   ! X0 = initial guess
  !   ! N = problem dimension
  !   ! GMM = restart parameter, controls memory required for matrix H
  !   ! ITERM = maximal number of iterations
  !   IMPLICIT NONE
  !   DOUBLE PRECISION, INTENT(IN) :: AMAT(N,N),BVEC(N),X0(N)
  !   INTEGER, INTENT(IN) :: N, GMM, NTRY
  !   DOUBLE PRECISION, INTENT(IN) :: TOL
  !   DOUBLE PRECISION, INTENT(OUT) :: XVEC(N)
  !   DOUBLE PRECISION, INTENT(OUT) :: RINFO(2)
  !   INTEGER :: LWORK, M
  !   DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
  !   DOUBLE PRECISION :: WORK0
  !   INTEGER :: IRC(5), ICNTL(7), INFO(3), J, REVCOM, COLX, COLY, COLZ, NBSCAL
  !   DOUBLE PRECISION :: CNTL(5)

  !   INFO = -10

  !   ! default GMM = N
  !   IF (GMM.LE.0) THEN
  !      M = N
  !   ELSE
  !      M = GMM
  !   ENDIF

  !   IF (M.GT.N) THEN
  !      PRINT*, 'ERROR IN RUNZGMRES: memory parameter M must be at most N', N, GMM
  !      STOP 1
  !   ENDIF

  !   LWORK = M*M+M*(N+5)+6*N+M+10
  !   !LWORK = 0
    
  !   OPEN(UNIT=33,FILE='gmresconv.out')

  !   IRC = 0
  !   ICNTL = 0
  !   CNTL = 0D0

  !   !print*, 'testx0:', ntry, LWORK

  !   ! control parameters for DRIVE_ZGMRES
  !   ! Initialize control parameters to default values
  !   CALL INIT_DGMRES(ICNTL,CNTL)
  !   ! ignore warning messages
  !   ICNTL(2) = 0
  !   ! save convergence history to file
  !   ICNTL(3) = 33
  !   ! No preconditioning
  !   ICNTL(4) = 0
  !   ICNTL(5) = 1 ! user-supplied initial guess
  !   ICNTL(7) = NTRY ! max number of iterations

  !   CNTL(1) = TOL ! tolerance for convergence

  !   IRC(1) = -1

  !   ! ! find optimal size of work array
  !   ! CALL DRIVE_ZGMRES(N,N,M,LWORK,WORK0,IRC,ICNTL,CNTL,INFO,RINFO)
  !   ! IF (INFO(1).EQ.-3) THEN
  !   !    PRINT*, 'TESTX1:', INFO
  !   !    LWORK = INFO(2)      
  !   !    ALLOCATE(WORK(LWORK))
  !   ! ELSE
  !   !    PRINT*, 'ERROR IN RUNZGMRES: failed to get optimal array size.', INFO, LWORK
  !   !    STOP 1
  !   ! ENDIF

  !   ! set up work array
  !   ALLOCATE(WORK(LWORK))
  !   WORK = 0D0
  !   ! initial guess
  !   WORK(1:N) = X0    
  !   ! right-hand side
  !   WORK(N+1:2*N) = BVEC

  !   DO WHILE (IRC(1).NE.0)
  !      CALL DRIVE_DGMRES(N,N,M,LWORK,WORK,IRC,ICNTL,CNTL,INFO,RINFO)          
  !     ! PRINT*, 'TESTX1:', IRC(1), info(1)
  !      ! REVCOM: what to do next (0 means done)
  !      ! 0: done
  !      ! 1: do matrix-vector product z<- Ax
  !      ! 2: do left preconditioning z<- M1^-1 *x
  !      ! 3: do right preconditioning z<-M2^-1*x
  !      ! 4: do NBSCAL scalar products z<-X*y
  !      REVCOM = IRC(1) ! what to do next (0 means done)
  !      ! index to read off x
  !      COLX = IRC(2) !
  !      ! index to read off y
  !      COLY = IRC(3)
  !      ! index to write Z
  !      COLZ = IRC(4)
  !      ! number of scalar products
  !      NBSCAL = IRC(5)

  !      IF (REVCOM.EQ.1) THEN
  !         ! do matrix vector product WORK(COLZ) <-- A*WORK(COLX)
  !         CALL DGEMV('N',N,N,DCMPLX(1D0),AMAT,N,WORK(COLX),1,DCMPLX(0D0),WORK(COLZ),1)
  !      ELSEIF (REVCOM.EQ.4) THEN
  !         CALL DGEMV('C',N,NBSCAL,DCMPLX(1D0),WORK(COLX),N,WORK(COLY),1,DCMPLX(0D0),WORK(COLZ),1)
  !      ELSE IF (REVCOM.NE.0) THEN
  !         PRINT*, 'ERROR IN RUNZGMRES: REVCOM=',revcom
  !         STOP 1
  !      ENDIF
     
  !   END DO

  !   ! Check exit conditions   
  !   IF (INFO(1).EQ.0) THEN
  !      ! PRINT*, 'normal exit'
  !      ! print*, 'convereged after ', info(2), ' iterations'
  !      ! print*, 'optimal and actual workspace size', info(3), lwork
  !   ELSEIF (INFO(1).EQ.-4) THEN
  !      PRINT*, 'ERROR IN RUNZGMRES: No convergence after ', icntl(7), ' iteractions'
  !      STOP 1
  !   ELSE
  !      PRINT*, 'ERROR IN RUNZGMRES: ', INFO
  !      print*, 'n, m:', N,M
  !      STOP 1
  !   ENDIF

  !   XVEC = WORK(1:N)

  !   ! PRINT*, 'TESTX6:'
  !   ! DO J = 1,N
  !   !    PRINT*, J, XVEC(J)
  !   ! ENDDO

  !   DEALLOCATE(WORK)
  !   CLOSE(33)
  ! END SUBROUTINE RUNDGMRES

  SUBROUTINE GMRES_RESTART(AMAT,BVEC,X0,M,N,ITERM,MAXTRY,TOL,XVEC,TRY,RES)
    ! solve AMAT*x = BVEC using iterative GMRES method
    ! run through ITERM cycles of GMRES before testing residual
    ! repeat  MAXTRY times before abandoning
    ! TOL is the tolerance on the residual to declare success
    ! WARNING: this is an inefficient implementation, should really test
    ! residual within inner loop as we go
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, N,ITERM,MAXTRY
    DOUBLE PRECISION, INTENT(IN) :: AMAT(M,N), BVEC(M), X0(N), TOL
    DOUBLE PRECISION, INTENT(OUT) :: XVEC(N), RES
    INTEGER, INTENT(OUT) :: TRY
    INTEGER :: MOUT
    DOUBLE PRECISION ::  XSTART(N)

    XSTART = X0
    DO TRY = 1,MAXTRY
       ! run inner gmres loop
       CALL GMRES(AMAT,BVEC,XSTART,M,N,ITERM,TOL,XVEC,MOUT,RES)

       IF (MOUT.LE.ITERM) THEN ! successfull completion
       ! compute residual
!       RES = BVEC
!       CALL DGEMV('N',M,N,1D0,AMAT,M,XVEC,1,-1D0,RES,1)

!       NRES = SQRT(DOT_PRODUCT(RES,RES))/sqrt(DOT_PRODUCT(BVEC,BVEC))
!       PRINT*, 'TESTX2:', NRES
!       IF (NRES.LE.TOL) THEN
          RETURN
       ELSE
          XSTART = XVEC
       ENDIF
    ENDDO
    
  END SUBROUTINE GMRES_RESTART

  SUBROUTINE ZGMRES_RESTART(AMAT,BVEC,X0,M,N,ITERM,MAXTRY,TOL,XVEC,TRY,RES)
    ! solve AMAT*x = BVEC using iterative GMRES method
    ! run through ITERM cycles of GMRES before testing residual
    ! repeat  MAXTRY times before abandoning
    ! TOL is the tolerance on the residual to declare success
    ! WARNING: this is an inefficient implementation, should really test
    ! residual within inner loop as we go
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, N,ITERM,MAXTRY
    COMPLEX*16, INTENT(IN) :: AMAT(M,N), BVEC(M), X0(N)
    DOUBLE PRECISION, INTENT(IN) :: TOL
    COMPLEX*16, INTENT(OUT) :: XVEC(N), RES
    INTEGER, INTENT(OUT) :: TRY
    INTEGER :: MOUT
    COMPLEX*16 ::  XSTART(M)

    XSTART = X0
    DO TRY = 1,MAXTRY
       ! run inner gmres loop
       CALL MYZGMRES(AMAT,BVEC,XSTART,M,N,ITERM,TOL,XVEC,MOUT,RES)

       IF (MOUT.LE.ITERM) THEN ! successfull completion
       ! compute residual
!       RES = BVEC
!       CALL DGEMV('N',M,N,1D0,AMAT,M,XVEC,1,-1D0,RES,1)

!       NRES = SQRT(DOT_PRODUCT(RES,RES))/sqrt(DOT_PRODUCT(BVEC,BVEC))
!       PRINT*, 'TESTX2:', NRES
!       IF (NRES.LE.TOL) THEN
          RETURN
       ELSE
          XSTART = XVEC
       ENDIF
    ENDDO
    
  END SUBROUTINE ZGMRES_RESTART

  SUBROUTINE GMRES(AMAT,BVEC,X0,M,N,ITERM,TOL,XVEC,MOUT,RES)
    ! algorithm 6.9 in Saad, 2007 book 
    ! solve AMAT*x = BVEC using iterative GMRES method
    ! M,N are matrix dimensions
    ! ITERM is the number of iterations of the method
    ! X = approximate solution
    ! MOUT = final inner iteration reached (will be ITERM unless exact solution found)

    INTEGER, INTENT(IN) :: M, N,ITERM
    DOUBLE PRECISION, INTENT(IN) :: AMAT(M,N), BVEC(M), X0(N), TOL
    DOUBLE PRECISION, INTENT(OUT) :: XVEC(N), RES
    INTEGER, INTENT(OUT) :: MOUT
    
    INTEGER :: LWORK
    DOUBLE PRECISION :: WORK(3*N)
    DOUBLE PRECISION :: VMAT(M,ITERM+1),BE1(ITERM+1),YVEC(ITERM), GVEC(ITERM+1)
    DOUBLE PRECISION :: HMAT(ITERM+1,ITERM), WJ(M),R0(M),BETA
    DOUBLE PRECISION :: ROTCOS(ITERM),ROTSIN(ITERM), GAMJ, NORMB, HJJ, HJ1J
    INTEGER :: I,J, INFO

    LWORK = N*3

    XVEC = X0

    ! R0 = -A*X0+B
    R0 = BVEC
    CALL DGEMV('N',M,N,-1D0,AMAT,M,X0,1,1D0,R0,1)
    
    BETA = SQRT(SUM(R0*R0))

    IF (ABS(BETA).LT.EPSILON(1D0)) THEN
       XVEC = 0D0
       MOUT = 0
       RES = 0
       RETURN
    ENDIF

    NORMB = SQRT(DOT_PRODUCT(BVEC,BVEC))

    HMAT = 0D0
    VMAT = 0D0
    GVEC = 0d0
    GVEC(1) = BETA

    VMAT(:,1) = R0/BETA    
    MOUT = ITERM

    
    DO J = 1,ITERM
       WJ = 0D0
       CALL DGEMV('N',M,N,1D0,AMAT,M,VMAT(:,J),1,0D0,WJ,1)
       
       DO I = 1,J
          HMAT(I,J) = DOT_PRODUCT(WJ,VMAT(:,I))
          WJ = WJ - HMAT(I,J)*VMAT(:,I)
       ENDDO

       HMAT(J+1,J) = SQRT(DOT_PRODUCT(WJ,WJ))
       IF (HMAT(J+1,J).NE.0D0) THEN
          VMAT(:,J+1) = WJ/HMAT(J+1,J)
       ENDIF

       ! Apply previous rotation matrices to Jth column (in place)
       DO I = 1,J-1
          CALL DROT(1,HMAT(I,J),1,HMAT(I+1,J),1,ROTCOS(I),ROTSIN(I))
       ENDDO

       ! get current rotation
       HJJ = HMAT(J,J); HJ1J = HMAT(J+1,J)

       CALL DROTG(HJJ,HJ1J,ROTCOS(J),ROTSIN(J))       

       ! apply current rotation to final column 
       CALL DROT(1,HMAT(J,J),1,HMAT(J+1,J),1,ROTCOS(J),ROTSIN(J))

       ! apply current rotation to RHS
       GAMJ = GVEC(J)
       GVEC(J) = ROTCOS(J)*GAMJ
       RES = -ROTSIN(J)*GAMJ
       GVEC(J+1) = RES       

       IF (GMRESVERBOSE) THEN
          PRINT*, "GMRES residual:", J, RES
       ENDIF       

       IF (ABS(RES).LT.TOL/NORMB) THEN
          ! reached termination 
          MOUT = J
          EXIT       
       ENDIF
              
    ENDDO

    ! minimize |beta*e1-H*y|
    !BE1=0D0; BE1(1) = BETA
    !CALL DGELS('N',MOUT+1,MOUT,1,HMAT,ITERM+1,BE1,ITERM+1,WORK,LWORK,INFO)

    ! Solve triangular system to get estimate of y
    YVEC = GVEC(1:ITERM)

    !print*, 'TESTX5:', MOUT, ITERM
    CALL DTRSV('U','N','N',MOUT,HMAT,ITERM+1,YVEC,1)
        

    ! calculate the overall approximate solution 
    XVEC = X0
    CALL DGEMV('N',N,MOUT,1D0,VMAT(:,1:MOUT),N,YVEC(1:MOUT),1,1D0,XVEC,1)
    

  END SUBROUTINE GMRES

    SUBROUTINE MYZGMRES(AMAT,BVEC,X0,M,N,ITERM,TOL,XVEC,MOUT,RES)
      ! iterative solution of linear system, for complex variables
    ! algorithm 6.9 in Saad, 2007 book 
    ! solve AMAT*x = BVEC using iterative GMRES method
    ! M,N are matrix dimensions
    ! ITERM is the number of iterations of the method
    ! X = approximate solution
    ! MOUT = final inner iteration reached (will be ITERM unless exact solution found)

    INTEGER, INTENT(IN) :: M, N,ITERM
    COMPLEX*16, INTENT(IN) :: AMAT(M,N),BVEC(M),X0(N)
    DOUBLE PRECISION, INTENT(IN) :: TOL
    COMPLEX*16, INTENT(OUT) :: XVEC(N), RES
    INTEGER, INTENT(OUT) :: MOUT
    
    INTEGER :: LWORK
    COMPLEX*16 :: WORK(3*N)
    COMPLEX*16 :: VMAT(N,ITERM+1),BE1(ITERM+1),YVEC(ITERM), GVEC(ITERM+1)
    COMPLEX*16 :: HMAT(ITERM+1,ITERM), WJ(M),R0(M)
    DOUBLE PRECISION :: BETA, NORMB
    COMPLEX*16 :: ROTCOS(ITERM),ROTSIN(ITERM), GAMJ, HJJ, HJ1J
    INTEGER :: I,J, INFO

    LWORK = N*3

    XVEC = X0

    ! R0 = -A*X0+B
    R0 = BVEC
    CALL ZGEMV('N',M,N,-1D0,AMAT,M,X0,1,1D0,R0,1)
    
    BETA = SQRT(DOT_PRODUCT(R0,R0))

    IF (ABS(BETA).LT.EPSILON(1D0)) THEN
       XVEC = 0D0
       MOUT = 0
       RES = 0
       RETURN
    ENDIF

    NORMB = SQRT(DOT_PRODUCT(BVEC,BVEC))

    HMAT = 0D0
    VMAT = 0D0
    GVEC = 0d0
    GVEC(1) = BETA

    VMAT(:,1) = R0/BETA    
    MOUT = ITERM

    
    DO J = 1,ITERM
       WJ = 0D0
       CALL ZGEMV('N',M,N,1D0,AMAT,M,VMAT(:,J),1,0D0,WJ,1)
       
       DO I = 1,J
          HMAT(I,J) = DOT_PRODUCT(WJ,VMAT(:,I))
          WJ = WJ - HMAT(I,J)*VMAT(:,I)
       ENDDO

       HMAT(J+1,J) = SQRT(DOT_PRODUCT(WJ,WJ))
       IF (HMAT(J+1,J).NE.0D0) THEN
          VMAT(:,J+1) = WJ/HMAT(J+1,J)
       ENDIF

       ! Apply previous rotation matrices to Jth column (in place)
       DO I = 1,J-1
          CALL ZDROT(1,HMAT(I,J),1,HMAT(I+1,J),1,ROTCOS(I),DBLE(ROTSIN(I)))
       ENDDO

       ! get current rotation
       HJJ = HMAT(J,J); HJ1J = HMAT(J+1,J)

       CALL ZROTG(HJJ,HJ1J,ROTCOS(J),ROTSIN(J))       

       ! apply current rotation to final column 
       CALL ZDROT(1,HMAT(J,J),1,HMAT(J+1,J),1,ROTCOS(J),DBLE(ROTSIN(J)))

       ! apply current rotation to RHS
       GAMJ = GVEC(J)
       GVEC(J) = ROTCOS(J)*GAMJ
       RES = -ROTSIN(J)*GAMJ
       GVEC(J+1) = RES       

       IF (GMRESVERBOSE) THEN
          PRINT*, "GMRES residual:", J, RES
       ENDIF       

       IF (ABS(RES).LT.TOL/NORMB) THEN
          ! reached termination 
          MOUT = J
          EXIT       
       ENDIF
              
    ENDDO

    ! minimize |beta*e1-H*y|
    !BE1=0D0; BE1(1) = BETA
    !CALL DGELS('N',MOUT+1,MOUT,1,HMAT,ITERM+1,BE1,ITERM+1,WORK,LWORK,INFO)

    ! Solve triangular system to get estimate of y
    YVEC = GVEC(1:ITERM)

    !print*, 'TESTX5:', MOUT, ITERM
    CALL ZTRSV('U','N','N',MOUT,HMAT,ITERM+1,YVEC,1)
        

    ! calculate the overall approximate solution 
    XVEC = X0
    CALL ZGEMV('N',N,MOUT,1D0,VMAT(:,1:MOUT),N,YVEC(1:MOUT),1,1D0,XVEC,1)
    

  END SUBROUTINE MYZGMRES
END MODULE GMRESUTIL
