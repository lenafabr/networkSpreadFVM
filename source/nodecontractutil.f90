MODULE NODECONTRACTUTIL
  ! utilities for dealing with network contractions
  USE NETWORKUTIL, ONLY : NETWORK
  USE DYNSYSUTIL, ONLY : DYNSYSTEM
  IMPLICIT NONE

  
CONTAINS

  SUBROUTINE GETVELS_NODECONTRACT(DSP,NETP,PSTART,DV,FLOWSPEED)
    ! update velocities on mesh boundaries, assuming contracting nodes
    ! PSTART: probability of a stationary node starting contraction on a time interval
    ! DV: change in volume (V/(pi a^2)) in each time step. Equal to flowspeed*delt
    ! FLOWSPEED: total flow velocity out of node during contraction
    ! updates NETP%contnodes, NETP%NODESTATE, NETP%NODEVOLS, DSP%VEL
    USE FLOWUTIL, ONLY : GETEDGEFLOWVELS
    USE MESHUTIL, ONLY : EDGETOBOUND
    
    IMPLICIT NONE
    TYPE(DYNSYSTEM), POINTER :: DSP
    TYPE(NETWORK), POINTER :: NETP
    DOUBLE PRECISION, INTENT(IN) :: PSTART, DV, FLOWSPEED
    DOUBLE PRECISION :: NODEFLOW(NETP%NNODE), EDGEFLOWS(NETP%NEDGE)
    LOGICAL :: FLOWCHANGE
    
    ! update contraction state, volume, etc.
    ! And get the total flow into each node
    CALL NODECONTRACTUPDATE(NETP,PSTART,DV,FLOWSPEED,NODEFLOW,FLOWCHANGE)

    IF (FLOWCHANGE) THEN
       ! some change to node state has occured, calculate change in flow pattern
       ! Calculate flow on each edge
       CALL GETEDGEFLOWVELS(NETP,NODEFLOW,EDGEFLOWS)   
       
       ! Set flow on each mesh boundary
       CALL EDGETOBOUND(NETP,DSP%MESHP,EDGEFLOWS,DSP%VEL)
    ENDIF
  END SUBROUTINE GETVELS_NODECONTRACT
  
  SUBROUTINE NODECONTRACTUPDATE(NETP,PSTART,DV,FLOWSPEED,NODEFLOW,FLOWCHANGE)
    ! Update flows and volumes of contracting and expanding network nodes

    ! PSTART: probability of a stationary node starting contraction on a time interval
    ! DV: change in volume (V/(pi a^2)) in each time step. Equal to flowspeed*delt
    ! FLOWSPEED: total flow velocity out of node during contraction
    ! NODEVOLS: current node volumes (in units of tube length: V/(pi a^2))
    ! nodes stop contracting when they reach 0 volume
    ! each contracting node has a partner node that expands while it contracts
    ! partner nodes remain the same while a closed node is expanding
    ! then when the partner is fully closed, it finds a new partner upon starting its expansion
    ! NCONT: number of contracted nodes
    ! CONTNODES(:,1): nodes either fully contracted or in the process of expanding
    ! CONTNODES(:,2): partner nodes for those in the process of contracting
    ! NODESTATE = 1 if contracting / expanding, 0 otherwise
    ! FLOWCHANGE = has there been any change in node state during this step?
    
    USE KEYS, ONLY : VERBOSE
    USE mt19937, ONLY : GRND
    USE GENUTIL, ONLY : RANDSELECT1_INT

    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    DOUBLE PRECISION, INTENT(IN) :: PSTART,FLOWSPEED, DV
   
    DOUBLE PRECISION, INTENT(OUT) :: NODEFLOW(NETP%NNODE) ! total flow out of each node
    LOGICAL, INTENT(OUT) :: FLOWCHANGE
    INTEGER :: CC, NC, NC2, NODELIST(NETP%NNODE), PICKNODE, IND
    DOUBLE PRECISION :: U, VOVERSHOOT, SUMFLOW
    LOGICAL :: JUSTSTOPPED(NETP%NNODE)

    NODELIST = (/ (NC, NC = 1,NETP%NNODE) /)
    NODEFLOW = 0D0

    ! has there been any change to the flows during this step?
    ! ie: anything start or stop contracting?
    FLOWCHANGE = .FALSE.

    JUSTSTOPPED = .FALSE.
    DO CC = 1,NETP%NCONT
       !PRINT*, 'TESTX2:',CC, NETP%CONTNODES(CC,:)
       !PRINT*, 'TESTX3:', cc, NETP%NODESTATE
       
       ! go through each current node that is fully contracted or in the process of expanding
       NC = NETP%CONTNODES(CC,1)
       NC2 = NETP%CONTNODES(CC,2) ! partner node
       IF (NETP%NODESTATE(NC).EQ.-2) THEN
          ! node is fully contracted. Decide whether to start expansion
          U = GRND()          
          IF (U.LT.PSTART) THEN             
             NETP%NODESTATE(NC) = 1
             NODEFLOW(NC) = -FLOWSPEED
             NETP%NODEVOLS(NC) = NETP%NODEVOLS(NC) + DV

             ! reset partner node, which will be contracting
             ! pick a new partner node from the list of fully expanded nodes
             CALL RANDSELECT1_INT(PACK(NODELIST,NETP%NODESTATE.EQ.2.and..not.JUSTSTOPPED),PICKNODE,IND)
             NETP%CONTNODES(CC,2) = PICKNODE

             IF (VERBOSE) PRINT*, 'Starting expansion: cont index, node, partner', CC, NC, PICKNODE
             
             ! start contracting the partner, flow out
             NETP%NODEVOLS(PICKNODE) = NETP%NODEVOLS(PICKNODE)-DV
             NODEFLOW(PICKNODE) = FLOWSPEED

             NETP%NODESTATE(NC) = 1 ! expanding node
             NETP%NODESTATE(PICKNODE) = -1 ! contracting node

             FLOWCHANGE = .TRUE.
          ELSE
             NODEFLOW(NC) = 0D0
             IF (NC2.GT.0) NODEFLOW(NC2) = 0D0
          ENDIF
       ELSEIF (NETP%NODESTATE(NC).EQ.1) THEN
          ! node is in the process of expanding
          NODEFLOW(NC) = -FLOWSPEED
          NETP%NODEVOLS(NC) = NETP%NODEVOLS(NC)+DV
          ! partner continues contracting
          NODEFLOW(NC2) = FLOWSPEED
          NETP%NODEVOLS(NC2) = NETP%NODEVOLS(NC2)-DV

          ! Decide if you're done expanding
          IF (NETP%NODEVOLS(NC2).LE.0) THEN             
             ! partner node fully contracted
             NETP%NODESTATE(NC2) = -2
             ! fully expanded node
             NETP%NODESTATE(NC) = 2

             ! mark these nodes as just stopped
             JUSTSTOPPED(NC) = .TRUE.
             JUSTSTOPPED(NC2)=.TRUE.
             
             
             ! Cut off contraction to give zero volume
             VOVERSHOOT = -NETP%NODEVOLS(NC2)
             NETP%NODEVOLS(NC2) = 0D0
             NODEFLOW(NC) = (1-VOVERSHOOT/DV)*NODEFLOW(NC)
             NODEFLOW(NC2) = (1-VOVERSHOOT/DV)*NODEFLOW(NC2)
             
             ! Cut off expanding node similarly
             NETP%NODEVOLS(NC) = NETP%NODEVOLS(NC)-VOVERSHOOT

             IF(VERBOSE) print*, 'Finished expanding. Index, node, partner ', CC, NC, NC2
             
             ! reset partner as the main contraction
             NETP%CONTNODES(CC,1) = NC2
             NETP%CONTNODES(CC,2) = 0 ! partner not set yet
             ! pick a new partner node from the list of fully expanded nodes
            ! CALL RANDSELECT1_INT(PACK(NODELIST,NETP%NODESTATE.EQ.2),PICKNODE,IND)
             !NETP%CONTNODES(CC,2) = PICKNODE

            FLOWCHANGE = .TRUE.
          ENDIF
       ELSE
          PRINT*, 'ERROR IN NODECONTRACTUPDATE: bad state for contracting node', CC, NC, NETP%NODESTATE(NC)
          STOP 1
       ENDIF
    END DO

    SUMFLOW = SUM(NODEFLOW)
    IF (SUMFLOW.GT.100E-16) THEN
       PRINT*, 'ERROR IN NODECONTRACTUPDATE: total flow does not sum to 0', SUMFLOW
       PRINT*, 'NODE STATE, VOL, FLOW:'
       DO NC = 1,NETP%NNODE
          PRINT*, NC, NETP%NODESTATE(NC), NETP%NODEVOLS(NC), NODEFLOW(NC)
       ENDDO
       STOP 1
    ENDIF
    
  
  END SUBROUTINE NODECONTRACTUPDATE
  
  SUBROUTINE OUTPUT_NODECONTRACT(NETP,OUTFILE,INFO,APPEND)
    ! output information on node volumes
    ! info is an additional list of floats to include (eg: for time information)
    ! APPEND = whether to append vs rewriting the file
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    CHARACTER(LEN=*) :: OUTFILE
    DOUBLE PRECISION, INTENT(IN) :: INFO(:)
    LOGICAL, INTENT(IN) :: APPEND
    INTEGER, PARAMETER :: FU=88
    
    IF (APPEND) THEN
       OPEN(UNIT=FU, FILE=OUTFILE,STATUS='UNKNOWN',ACCESS='APPEND')
    ELSE
       OPEN(UNIT=FU, FILE=OUTFILE,STATUS='UNKNOWN')
    ENDIF

    ! write number of nodes, followed by extra information
    WRITE(FU,*) NETP%NNODE, INFO
    ! write states of each node
    WRITE(FU,*) NETP%NODESTATE
    ! write volume of each node
    WRITE(FU,*) NETP%NODEVOLS

    CLOSE(FU)
    
  END SUBROUTINE OUTPUT_NODECONTRACT

  SUBROUTINE INITIALIZE_NODECONTRACT(NETP,NCONT)
    ! initialize network node contractions
    USE GENUTIL, ONLY : RANDSELECT_INT
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: NCONT
    INTEGER :: INDS(2*NCONT), VALS(2*NCONT), NODELIST(NETP%NNODE)      
    INTEGER :: I


    NODELIST = (/(I, I = 1,NETP%NNODE)/)

    IF (.NOT.NETP%CONTSET) THEN
       ! allocate arrays
       ALLOCATE(NETP%CONTNODES(NCONT,2))         
       NETP%CONTSET = .TRUE.
       NETP%NCONT = NCONT
    ENDIF

    ! pick contracting nodes and partners at random, without replacement
    CALL RANDSELECT_INT(NODELIST,2*NCONT,.FALSE.,VALS,INDS)
    NETP%CONTNODES(:,1) = INDS(1:NCONT)
    NETP%CONTNODES(:,2) = INDS(NCONT+1:2*NCONT)

    ! start with contracted nodes fully closed, everything else fully open
    NETP%NODESTATE = 2
    NETP%NODESTATE(INDS(1:NCONT)) = -2

    NETP%NODEVOLS(INDS(1:NCONT)) = 0D0
  END SUBROUTINE INITIALIZE_NODECONTRACT

END MODULE NODECONTRACTUTIL
