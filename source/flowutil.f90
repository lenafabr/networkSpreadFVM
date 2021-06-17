MODULE FLOWUTIL
  ! utilities for dealing with fluid flow within the network
  USE NETWORKUTIL, ONLY : NETWORK
  IMPLICIT NONE
  
CONTAINS

  SUBROUTINE GETEDGEFLOWVELS(NETP,FLOWNODES,FLOWVELS)   
    ! given the total flow into (positive) or out of (negative) each node   
    ! calculate the flow velocities along the edges and `pressures' at nodes
    ! by solving Kirchoff's equations for the network
    ! FLOWNODES: flux into each node divided by cross-sectional area
    ! (equal to sum of flow velocities in edges conneted to that node)
    ! FLOWVELS: output flow velocity on each edge
    USE GENUTIL, ONLY : PRINTMAT
    USE NETWORKUTIL, ONLY : OUTPUTNETWORK
    IMPLICIT NONE

    TYPE(NETWORK), POINTER :: NETP
    DOUBLE PRECISION, INTENT(IN) :: FLOWNODES(NETP%NNODE)
    DOUBLE PRECISION, INTENT(OUT) :: FLOWVELS(NETP%NEDGE)    
    
    ! ALPHA = resistivity connecting pressure to flow; assumed to be 1 throughout
    ! this just sets the pressure units; in principle it should depend on the thickness of each tube
    DOUBLE PRECISION, PARAMETER :: ALPHA = 1D0
    DOUBLE PRECISION :: MAT(NETP%NEDGE,NETP%NEDGE)
    DOUBLE PRECISION :: MAT0(NETP%NEDGE,NETP%NEDGE)
    DOUBLE PRECISION :: VEC(NETP%NEDGE), VEC0(NETP%NEDGE), VECGUESS(NETP%NEDGE)
    
    INTEGER :: EC, START, START2, LC, EDGE
    INTEGER :: MATSIZE
    INTEGER :: INFO, TRY, MAXTRY, ITERM
    integer:: LWORK, IPIV(NETP%NEDGE)
    DOUBLE PRECISION :: WORK(3*(NETP%NEDGE)**2), RES
    DOUBLE PRECISION :: TOL, RINFO(2)

!    DOUBLE PRECISION :: SVALS(NETP%NNODE+NETP%NEDGE)
!    DOUBLE PRECISION :: UMAT(NETP%NNODE+NETP%NEDGE,NETP%NNODE+NETP%NEDGE)
!    DOUBLE PRECISION :: VTMAT(NETP%NNODE+NETP%NEDGE,NETP%NNODE+NETP%NEDGE)

    IF (NETP%NLOOP.NE.NETP%NEDGE-NETP%NNODE+1) THEN
       PRINT*, 'ERROR IN GETEDGEFLOW: &
            & currently only set up for networks with one connected component',&
            & NETP%NLOOP, NETP%NEDGE, NETP%NNODE
       STOP 1
    ENDIF
    
    LWORK =     3*(NETP%NEDGE)**2
    
    IF (SUM(FLOWNODES).GT.100D-16) THEN
       PRINT*, 'ERROR IN GETEDGEFLOWVELS: total flow into system does not sum to 0', SUM(FLOWNODES)
       
       STOP 1
    ENDIF

    MATSIZE = NETP%NEDGE
    
    ! Set up the matrix of coefficients for Kirchoff's equations
    MAT = 0D0
       
    ! (equations for total flow into each node)
    ! Mij = -1 if edge j starts at node i (flowing out)
    ! Mij = 1 if edge j ends at node i (flowing in)
    ! Mij = 0 otherwise
    ! skip the last node since that constraint is redundant
    DO EC = 1,NETP%NEDGE
      ! print*, 'testx1:', NETP%NNODE, EC, NETP%EDGENODE(EC,:)
       IF (NETP%EDGENODE(EC,1).LT.NETP%NNODE) MAT(NETP%EDGENODE(EC,1),EC) = -1D0
       IF (NETP%EDGENODE(EC,2).LT.NETP%NNODE) MAT(NETP%EDGENODE(EC,2),EC) = 1D0
    ENDDO

    ! equations for zero pressure change around each loop
    ! leave off the last flow conservation equation
    ! as it gives an overdetermined system
    ! (assume total volume change in whole network is zero)
    START = NETP%NNODE-1
    DO LC = 1,NETP%NLOOP
       DO EC= 1,NETP%LOOPLENS(LC)
          EDGE = ABS(NETP%LOOPEDGES(LC,EC))
          MAT(START+LC,EDGE) = ALPHA*SIGN(1,NETP%LOOPEDGES(LC,EC))*NETP%EDGELEN(EDGE)
       ENDDO
    END DO
    
    ! set up the result vector for Kirchoffs equations
    ! (total flow OUT each node in the top N spots, then zeros)
    VEC = 0d0
    VEC(1:NETP%NNODE-1) = -FLOWNODES(1:NETP%NNODE-1)

    ! add constraint on total pressures (sum to 0)
    !MAT(MATSIZE+1,NETP%NEDGE+1:MATSIZE) = 1D0
    !VEC(MATSIZE+1) = 0D0

    MAT0 = MAT;
    VEC0 = VEC;

    
    
    ! solve the equations
     CALL DGESV(MATSIZE,1,MAT,MATSIZE,IPIV,VEC,MATSIZE,INFO)

    
  !print*, 'done with dgesv', info
    
    
     CALL DGEMV('N',MATSIZE,MATSIZE,1D0,MAT0,MATSIZE,VEC,1,-1D0,VEC0,1)
    ! PRINT*, 'TESTXY:', maxval(VEC0)
    ! PRINT*, 'TESTXZ:', shape(mat)
    
    
    IF (INFO.NE.0) THEN
       PRINT*, 'ERROR IN GETEDGEFLOWVELS: failed to solve linear system', INFO

       PRINT*, 'VEC:', VEC0
       PRINT*, 'MAT:'
       CALL PRINTMAT(MAT0)
       STOP 1
    ENDIF

    IF (MAXVAL(VEC).GT.1D8) THEN
       PRINT*, 'ERROR IN GETEDGEFLOWVELS: bad solution for flow', INFO
       PRINT*, VEC
       STOP 1
    ENDIF
    
    FLOWVELS = VEC(1:NETP%NEDGE)
  END SUBROUTINE GETEDGEFLOWVELS
  
END MODULE FLOWUTIL
