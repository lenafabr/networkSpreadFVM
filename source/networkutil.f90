MODULE NETWORKUTIL
  ! Utilities for dealing with connectivity and geometry of a network
  ! including defining, manipulating, input/output

  IMPLICIT NONE
  
  TYPE NETWORK
     ! define a network structure (connectivity and geometry)

     ! dimension of the space the network is embedded in
     INTEGER :: DIM 

     INTEGER :: NFIELD ! number of fields that will be evolved on this network     
     ! ----------------------
     ! information on network nodes
     ! ----------------------
     INTEGER :: NNODE ! number of nodes
     ! list of other node indices each node connects to
     INTEGER, POINTER :: NODENODE(:,:)
     ! degree (number of branches) of each node
     INTEGER, POINTER :: NODEDEG(:)
     ! list of branch indices each node connects to
     INTEGER, POINTER :: NODEEDGE(:,:)          
     ! lengths of all edges leaving the node
     DOUBLE PRECISION, POINTER :: NODELEN(:,:)
     ! spatial location of node
     DOUBLE PRECISION, POINTER :: NODEPOS(:,:)
     ! which nodes are absorbing (conc set to 0 and flux tracked)
     LOGICAL, POINTER :: NODEABS(:,:)
     ! which nodes have fixed conc, and value of that fixed concentration
     !LOGICAL, POINTER :: NODEFIX(:,:)
     !DOUBLE PRECISION, POINTER :: NODEFIXVAL(:,:)
     ! width of nodes (for use in keeping track of contractions)
     DOUBLE PRECISION, POINTER :: NODEWIDTH(:)
     
     ! ------------------
     ! information on network branches
     ! ------------------
     INTEGER :: NEDGE ! Number of branches
     ! nodes at the start and end of each branch
     INTEGER, POINTER :: EDGENODE(:,:)
     ! spatial position of branch starting points
     ! branch direction and length
     DOUBLE PRECISION, POINTER :: EDGESTART(:,:), EDGEDIR(:,:), EDGELEN(:)
     ! number of independent cycles
     ! for each independent cycle, list edge indices
     ! negative values for going around the edge backward
     INTEGER :: NLOOP, MAXLOOPLEN
     INTEGER, POINTER :: LOOPEDGES(:,:), LOOPLENS(:)
     
     ! -------------
     ! Map nodes and edges to a mesh object     
     ! -------------
     ! max number of cells on an edge
     INTEGER :: MAXCELLEDGE
     ! EDGENCELL: number of cells along each edge
     ! EDGECELLS: list of cell indices for each edge (1st dim is edge number)
     ! NODECELLS: list of cell index for each node
     !    lists zero when the node is terminal and has no cell
     INTEGER, POINTER :: EDGENCELL(:), EDGECELLS(:,:), NODECELLS(:)
     
     ! ----------------------
     ! Define contraction points along the edges
     ! ----------------------
     INTEGER :: NCONT ! number of contraction points
     ! For node contractions
     INTEGER, POINTER :: NODESTATE(:), CONTNODES(:,:)
     DOUBLE PRECISION, POINTER :: NODEVOLS(:)
     
          
     ! ARRAYSET: Have arrays been allocated?
     ! STRUCTURESET: has node and edge information been set up?
     ! MESHSET: has mesh array been set up yet?
     LOGICAL :: ARRAYSET = .FALSE., STRUCTURESET = .FALSE., MESHSET=.FALSE.
     ! have contraction arrays been set up yet?
     LOGICAL :: CONTSET = .FALSE.
     
     ! --------
     ! Reservoirs
     ! ---------
     INTEGER :: NRESV ! number of reservoirs     
     ! reservoir associated with each node (could be 0 if node is not part of reservoir)
     INTEGER, POINTER :: NODERESV(:)
     ! which nodes are associated with each reservoir
     ! Also number of nodes for each reservoir
     INTEGER, POINTER :: RESVNODES(:,:), RESVNNODE(:)
     ! volume and surface area and effective length of each reservoir
     DOUBLE PRECISION, POINTER :: RESVVOL(:), RESVSA(:), RESVLEN(:)
     INTEGER, POINTER :: RESVCELLS(:) ! which cell each reservoir corresponds to
     ! is each reservoir well-mixed?
     LOGICAL, POINTER :: RESVMIX(:)

     ! permeability of nodes
     LOGICAL, POINTER :: ISPERM(:)

     ! Evolve buffer protein concentrations
     ! buffertype = 0: no buffer prots
     ! buffertype = 1: explicit buffer prots
     ! buffertype = 2: rapidly equilibrated buffer prots
     INTEGER :: BUFFERTYPE
     DOUBLE PRECISION :: KON, KOFF, KDEQUIL

     ! Fixing edge velocity
     INTEGER :: NFIXEDGEVEL
     LOGICAL, POINTER :: FIXEDGEVEL(:)
     DOUBLE PRECISION, POINTER :: FIXEDGEVELVAL(:)
     
  END type NETWORK

CONTAINS
  

  SUBROUTINE INSERTINTERNODE(NETP,EDGE,FRAC,NID,EIDIN,WIPEPREV)
    ! insert intermediate node along an edge
    ! reconnecting everything appropriately
    ! place new node in index NID, new edge in index EIDIN
    ! new node will be along edge EDGE, a fraction FRAC of the way along it
    ! WIPEPREV: if true, then the node NID was already set up
    ! as a degree 2 node and needs to be removed from the network
    ! in this case, use the 2nd edge attached to it instead of EIDIN
    ! So when WIPEPREV is true, the value of EIDIN does not matter
    
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: EDGE, NID, EIDIN
    DOUBLE PRECISION, INTENT(IN) :: FRAC
    LOGICAL, INTENT(IN) :: WIPEPREV

    INTEGER :: EID, NODE1, NODE2, NC
    DOUBLE PRECISION :: LEN1, LEN2, POS(NETP%DIM)

    IF (FRAC.LT.-TINY(0D0).OR.FRAC.GT.1D0+TINY(1D0)) THEN
       PRINT*, 'ERROR IN INSERTINTERNODE: frac out of range', FRAC
       STOP 1
    ENDIF
    
    IF (WIPEPREV) THEN
       IF (NETP%NODEDEG(NID).NE.2) THEN
          PRINT*, 'ERROR IN INSERTINTERNODE: previous node does not have degree 2'
          STOP 1
       ENDIF
       PRINT*, 'ERROR IN INSERTINTERNODE: wipeprev not set up yet'
       STOP 1
    ELSE
       EID = EIDIN
    ENDIF
    
    ! edge where contraction point is       
    NODE1 = NETP%EDGENODE(EDGE,1) ! first node for this edge
    NODE2 = NETP%EDGENODE(EDGE,2) ! second node for this edge
       
    ! length of edge up to contraction and after it
    LEN1 = FRAC*NETP%EDGELEN(EDGE)
    LEN2 = NETP%EDGELEN(EDGE)-LEN1
       
    ! new node position in space       
    POS = NETP%EDGESTART(EDGE,:) &
            & + LEN1*NETP%EDGEDIR(EDGE,:)
       
    NETP%NODEPOS(NID,:) = POS
    NETP%NODEDEG(NID) = 2; ! extra node is inserted into edge, deg=2
    ! connect this node to the nodes belonging to this edge
    NETP%NODENODE(NID,1:2) = NETP%EDGENODE(EDGE,:)
    NETP%NODELEN(NID,:) = (/LEN1,LEN2/) ! total length connected to node

    ! break up edge into 2. Current EDGE index keeps the first half
    ! EID index gets the second half
    NETP%EDGELEN(EDGE) = LEN1   
    NETP%EDGENODE(EDGE,2) = NID
       
    NETP%EDGELEN(EID) = LEN2
    NETP%EDGENODE(EID,:) = (/NID,NODE2/)
    NETP%EDGESTART(EID,:) = POS
    NETP%EDGEDIR(EID,:) = NETP%EDGEDIR(EDGE,:)

    NETP%NODEEDGE(NID,1:2) = (/EDGE,EID/)

    ! map from new edge indices to original edges
   ! NETP%EDGEORIG(EID) = NETP%EDGEORIG(EDGE)
    
    ! for first node on contraction edge, update connectivity
    DO NC = 1,NETP%NODEDEG(NODE1)
       IF (NETP%NODENODE(NODE1,NC).EQ.NODE2) THEN
          NETP%NODENODE(NODE1,NC) = NID
       ENDIF
    ENDDO
       
    ! for last node on contraction edge, update connectivity
    DO NC = 1,NETP%NODEDEG(NODE2)
       IF (NETP%NODENODE(NODE2,NC).EQ.NODE1) THEN
          NETP%NODENODE(NODE2,NC) = NID
          NETP%NODEEDGE(NODE2,NC) = EID
       ENDIF
    ENDDO

    ! update loops to deal with the split edge
    CALL UPDATELOOP_EDGESPLIT(NETP,EDGE,EID)

  END SUBROUTINE INSERTINTERNODE
  
  SUBROUTINE UPDATELOOP_EDGESPLIT(NETP,TARGETEDGE,EDGE2)
    ! update loop structures to account for a new degree 2 node being inserted into an edge
    ! TARGETEDGE = index of edge being split
    ! EDGE2 = index for second edge after the split
    ! The first edge in the split keeps the original index
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: TARGETEDGE, EDGE2
    INTEGER :: LC, EC, LOOPEDGES(NETP%MAXLOOPLEN)

    ! Update any loops containing the newly split edge
    DO LC = 1,NETP%NLOOP ! for each loop
       LOOPEDGES = NETP%LOOPEDGES(LC,:)

       DO EC = 1,NETP%LOOPLENS(LC) ! for each edge in loop
          IF (ABS(LOOPEDGES(EC)).EQ.TARGETEDGE) THEN ! loop includes edge
             NETP%LOOPLENS(LC) = NETP%LOOPLENS(LC)+1
             IF (NETP%LOOPLENS(LC).GT.NETP%MAXLOOPLEN) THEN
                PRINT*, 'ERROR IN JUMPNODE: loop too long', LC, EC
                STOP 1
             ENDIF
             IF (LOOPEDGES(EC).GT.0) THEN !  forward along edges                
                NETP%LOOPEDGES(LC,EC+2:NETP%LOOPLENS(LC)) = &
                     & NETP%LOOPEDGES(LC,EC+1:NETP%LOOPLENS(LC)-1)
                NETP%LOOPEDGES(LC,EC+1) = EDGE2
             ELSE !  backward along edge
                NETP%LOOPEDGES(LC,EC+1:NETP%LOOPLENS(LC)) = &
                     & NETP%LOOPEDGES(LC,EC:NETP%LOOPLENS(LC)-1)
                NETP%LOOPEDGES(LC,EC) = -EDGE2
             ENDIF
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE UPDATELOOP_EDGESPLIT

  SUBROUTINE UPDATELOOP_EDGEMERGE(NETP,EDGE1)
    ! update loop structures to account for a new degree 2 node being inserted into an edge
    ! EDGE1 = index of first piece of edge that will be merged
    ! newly merged edge keeps this index
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: EDGE1
    INTEGER :: LC, EC, LOOPEDGES(NETP%MAXLOOPLEN)


    ! update any loops containing the now-merged edge
    DO LC = 1,NETP%NLOOP ! for each loop
       LOOPEDGES = NETP%LOOPEDGES(LC,:)

       DO EC = 1,NETP%LOOPLENS(LC) ! for each edge in loop
          IF (ABS(LOOPEDGES(EC)).EQ.EDGE1) THEN ! loop includes edge
             NETP%LOOPLENS(LC) = NETP%LOOPLENS(LC)-1
             IF (NETP%LOOPLENS(LC).LE.0) THEN
                PRINT*, 'ERROR IN JUMPNODE: loop too short', LC, EC
                STOP 1
             ENDIF
             IF (LOOPEDGES(EC).GT.0) THEN !  forward along edges                
                NETP%LOOPEDGES(LC,EC+1:NETP%LOOPLENS(LC)) = &
                     & NETP%LOOPEDGES(LC,EC+2:NETP%LOOPLENS(LC)+1)
             ELSE !  backward along edge
                NETP%LOOPEDGES(LC,EC-1:NETP%LOOPLENS(LC)) = &
                     & NETP%LOOPEDGES(LC,EC:NETP%LOOPLENS(LC)+1)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE UPDATELOOP_EDGEMERGE
  
  ! SUBROUTINE SETUPFIXNODES(NETP,NFIX,FIXVALS,FIXNODES)
  !   ! set up fixed nodes on a network
  !   ! either from a given list, or set randomly
  !   USE mt19937, ONLY : GRND
  !   IMPLICIT NONE
  !   TYPE(NETWORK), POINTER :: NETP
  !   INTEGER, INTENT(IN) :: NFIX(:)
  !   DOUBLE PRECISION, INTENT(IN) :: FIXVALS(:,:)
  !   INTEGER, INTENT(IN), OPTIONAL :: FIXNODES(:,:)
  !   INTEGER :: FC, A, CC, TRY, NID
  !   INTEGER, PARAMETER :: MAXNTRY= 100
    
  !   IF (.NOT.NETP%ARRAYSET) THEN
  !      PRINT*, 'ERROR IN SETFIXNODES: network not set up yet'
  !      STOP 1    
  !   ENDIF
    
  !   IF (PRESENT(FIXNODES)) THEN ! set the desired fixed nodes
  !      NETP%NODEFIX = .FALSE.
  !      DO FC = 1,NETP%NFIELD
  !         DO A = 1,NFIX(FC)
  !            IF (FIXNODES(A,FC).GT.NETP%NNODE.OR.FIXNODES(A,FC).LT.1) THEN
  !               PRINT*, 'ERROR IN SETUPFIXNODES: bad fixed node', A, FIXNODES(A,FC)
  !               STOP 1
  !            ENDIF
  !            NETP%NODEFIX(FIXNODES(A,FC),FC) = .TRUE.             
  !            NETP%NODEFIXVAL(FIXNODES(A,FC),FC) = FIXVALS(A,FC)              
  !         ENDDO
  !      ENDDO
  !   ELSE
  !      ! set random fixed nodes
  !      DO FC = 1,NETP%NFIELD
  !         NETP%NODEFIX(:,FC) = .FALSE.
  !         ! Sample fix nodes for this field
  !         DO CC = 1,NFIX(FC)
  !            DO TRY = 1,MAXNTRY ! try picking a new node
  !               NID = FLOOR(GRND()*NETP%NNODE)+1
  !               IF (NID.GT.NETP%NNODE) CYCLE
  !               ! cannot fix reservoir nodes
  !               IF (.not.NETP%NODEFIX(NID,FC).and.NETP%NODERESV(NID).EQ.0) EXIT                
  !            ENDDO
  !            IF (TRY.GE.MAXNTRY) THEN
  !               PRINT*, 'ERROR IN SETTING FIXED NODES:   failed to find a new node', FC, CC, NID
  !               STOP 1
  !            ENDIF

  !            PRINT*, 'FIXED NODE:', NID, NETP%NODERESV(NID)
  !            NETP%NODEFIX(NID,FC) = .TRUE.
  !            NETP%NODEFIXVAL(NID,FC) = FIXVALS(CC,FC)             
  !         ENDDO
  !      ENDDO
  !   ENDIF

    
  ! END SUBROUTINE SETUPFIXNODES

  
  SUBROUTINE NETWORKFROMFILE(NETP,NETFILE)
    USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI
    USE KEYS, ONLY : MAXBRANCH,  NABS, ABSORBERS, &
         & MAXLOOPLEN, DORESERVOIRS, NFIELD, RANDFIXNODES,&
         & NODEVOL,RESVVOL, RESVSA,RESVLEN,FIXEDGEVEL,FIXEDGEVELVAL,&
         & NFIXEDGEVEL, FIXNODES, FIXVALS, NFIX, &
         & FIXNODEFROMNETFILE, NETWORKDIM
    
    USE GENUTIL, ONLY : NORMALIZE
    ! Set up (allocate) a network structure, reading in connectivity and
    ! geometry from an input file    
    
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    CHARACTER(LEN=*), INTENT(IN) :: NETFILE
    LOGICAL :: LDUM, CASESET
    LOGICAL :: FILEEND=.FALSE., DIMSET = .FALSE.
    CHARACTER(LEN=100) :: WORD, LBL, STR
    INTEGER :: NITEMS, NNODE, NEDGE, NLOOP, NODE1, NODE2, DIM, EID, NID, NE
    INTEGER :: LC, WHICHLOOP,I, CC
    DOUBLE PRECISION :: LEN, LENS(MAXBRANCH), LENSAVE(MAXBRANCH)
    INTEGER :: TMPARRAY(MAXBRANCH)
    LOGICAL, ALLOCATABLE :: EDGELENSET(:)
    INTEGER :: NODERESV, MAXRESV
    INTEGER :: NC, RC, FC
    LOGICAL, ALLOCATABLE :: ISFIXED(:)
    DOUBLE PRECISION :: SAVEFIXVAL
    
    INTEGER, PARAMETER :: NF = 55 ! input file unit number

    
    MAXRESV = 0
    
    ! deallocate any previously set arrays
    IF (NETP%ARRAYSET) CALL CLEANUPNETWORK(NETP)

    
    ! go through file and count nodes and branches (total and from each node)
    PRINT*, 'Reading network structure file: ', NETFILE
    INQUIRE(FILE=NETFILE,EXIST=LDUM)
    IF (.NOT.LDUM) THEN
       PRINT*, 'ERROR in NETWORKFROMFILE: network file ', TRIM(ADJUSTL(NETFILE)), ' does not exist.'
        STOP 1
     ENDIF
     OPEN(UNIT=NF, FILE=NETFILE, STATUS='OLD')

     ! Count number of nodes and edges, set spatial dimension
     DIMSET = .FALSE.
     NNODE = 0; NEDGE = 0; NLOOP = 0

     ! predefine network dimension
     IF (NETWORKDIM.GT.0) THEN
        DIMSET = .TRUE.
        DIM = NETWORKDIM
     END IF    
     
     DO 
        CALL READLINE(NF,FILEEND,NITEMS)
        
        IF (FILEEND.and.nitems.eq.0) EXIT
        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE
        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)

        IF (WORD.EQ.'NODE') THEN
           NNODE = NNODE+1
           ! set spatial dimension of network
           ! IF (DIMSET.AND.DIM.NE.NITEMS-2) THEN
           !    PRINT*, 'ERROR IN NETWORKFROMFILE: dimension of node positions is inconsistent', DIM, NITEMS
           !    STOP 1
           ! ENDIF     
           
           IF (NETWORKDIM.LE.0) DIM = NITEMS-2

           ! loop through remaining items looking for reservoir label
           DO I = 1,NITEMS-1
              CALL READA(LBL)
              IF ((LBL(1:1).EQ.'R').OR.LBL(1:1).EQ.'P'.OR.LBL(1:1).EQ.'F') THEN
                 IF (NETWORKDIM.LE.0) DIM = NITEMS - 3
              ENDIF
              IF (LBL(1:1).EQ.'R'.AND.DORESERVOIRS) THEN ! attach reservoir label                 
                 READ(LBL(2:100),*) NODERESV
                 IF (NODERESV.GT.MAXRESV) MAXRESV = NODERESV
              ENDIF
           ENDDO        
           
        ELSEIF (WORD.EQ.'EDGE') THEN
           NEDGE = NEDGE + 1
        ELSEIF (WORD.EQ.'LOOP') THEN
           NLOOP = NLOOP+1
        ENDIF
     ENDDO
     CLOSE(NF)

     PRINT*, 'Number of nodes, edges, loops, dimension: ', NNODE, NEDGE, NLOOP,DIM
     ! allocate arrays         
     CALL SETUPNETWORK(NETP,NNODE,NEDGE,MAXRESV,NLOOP,DIM,NFIELD,MAXBRANCH,MAXLOOPLEN,&
          & ABSORBERS,NABS)
          
     ALLOCATE(EDGELENSET(NEDGE), ISFIXED(NNODE))
     EDGELENSET=.FALSE. ! which edges are set directly from file

     PRINT*, 'Reading in node positions and edge connectivity...'
     OPEN(UNIT=NF, FILE=NETFILE, STATUS='OLD')
     ! Get node and edge information
     WHICHLOOP = 0
     IF (FIXNODEFROMNETFILE) ISFIXED = .FALSE.
     DO 
        CALL READLINE(NF,FILEEND,NITEMS)
        IF (FILEEND.and.nitems.eq.0) EXIT
        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE
        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)
        ! Skip any comment lines
        IF (WORD(1:1).EQ.'#') CYCLE

        IF (WORD.EQ.'NODE') THEN
           CALL READI(NID)
           IF (NID.LT.1.OR.NID.GT.NNODE) THEN
              PRINT*, 'ERROR IN NETWORKFROMFILE: invalid  node index', NNODE, NID
              STOP 1
           ENDIF
           CALL READF(NETP%NODEPOS(NID,1))
           CALL READF(NETP%NODEPOS(NID,2))
           IF (NETP%DIM.GT.2) THEN
              CALL READF(NETP%NODEPOS(NID,3))!this line had to be added for 3d
           ENDIF
           
           ! loop through remaining items looking for reservoir label or permeability label          
           DO I = 1,NITEMS-2-NETP%DIM
              CALL READA(LBL)
              IF (LBL(1:1).EQ.'R'.AND.DORESERVOIRS) THEN ! attach reservoir label
                 READ(LBL(2:100),*) NETP%NODERESV(NID)
                 
              ELSE IF(LBL(1:1).EQ.'P') THEN ! set node as permeable for a given field
                 READ(LBL(2:100),*) FC
                 NETP%ISPERM(NID) = .TRUE.
              ELSEIF(LBL(1:1).EQ.'F'.AND.FIXNODEFROMNETFILE) THEN ! set this as a fixed node
                 ISFIXED(NID) = .TRUE.
              ENDIF              
           ENDDO

        ELSEIF (WORD.EQ.'EDGE') THEN
           CALL READI(EID) ! edge id
           CALL READI(NODE1)
           CALL READI(NODE2)
           IF (NITEMS.GT.4) THEN ! read in edge length directly
              CALL READF (NETP%EDGELEN(EID))
              EDGELENSET(EID) = .TRUE.
           ENDIF
           
           IF (NODE1.LT.1.OR.NODE2.LT.1.OR.EID.LT.1&
                & .OR.NODE1.GT.NNODE.OR.NODE2.GT.NNODE.OR.EID.GT.NEDGE) THEN
              PRINT*, 'ERROR IN NETWORK FROM FILE: &
                   & bad edge or node indices while reading edge.', &
                   ' EID, NODE1, NODE2, NNODE, NEDGE:', &
                   & EID, NODE1, NODE2, NNODE, NEDGE
              STOP 1
           ENDIF

           ! nodes connected to this edge           
           NETP%EDGENODE(EID,:) = (/NODE1,NODE2/)
        ELSEIF (WORD.EQ.'LOOP') THEN
           WHICHLOOP = WHICHLOOP+1           
           IF (WHICHLOOP.LE.NETP%NLOOP) THEN
 !             IF (NITEMS.GT.MAXLOOPLEN+1) THEN
                 !PRINT*, 'ERROR IN NETWORK FROM FILE: loop is too long', NITEMS, MAXLOOPLEN
                 !STOP 1
!              ENDIF
                 NETP%LOOPLENS(WHICHLOOP) = NITEMS-1
                 DO LC = 1,NITEMS-1
                    IF (LC.GT.MAXLOOPLEN) THEN
                       PRINT*, 'WARNING: LOOP IS TOO LONG. IGNORING THE REST'
                       EXIT
                    ENDIF
                    CALL READI(NETP%LOOPEDGES(WHICHLOOP,LC) )
                 ENDDO
           ELSE
              PRINT*, 'WARNING: too many loops, ignoring the rest', WHICHLOOP, NETP%NLOOP
           ENDIF
        ENDIF        
     END DO
     CLOSE(NF)

     IF (FIXNODEFROMNETFILE) THEN
        ! update fixed nodes from network file (assume field 1 fixed only)
        IF (NFIX(1).LT.1) THEN
           PRINT*, 'ERROR: cannot set up fix nodes from network if no FIXVAL set in param file.'
           STOP 1
        ENDIF
        IF (NFIELD.GT.1.AND.NFIX(2).GT.0) THEN
           PRINT*, 'ERROR: cannot get fixed nodes beyond field 1 from network file.'
           STOP 1
        ENDIF
        
        SAVEFIXVAL = FIXVALS(1,1)
        
        NFIX = 0
        FIXNODES = 0
        DO NC = 1,NETP%NNODE
           IF (ISFIXED(NC)) THEN
              NFIX(1) = NFIX(1) + 1
              FIXNODES(NFIX(1),1) = NC
              FIXVALS(NFIX(1),1) = SAVEFIXVAL
           ENDIF
        ENDDO
        
        PRINT*, 'Fixed nodes for field 1 only:', NFIX(1), FIXNODES(1:NFIX(1),1)
        
     END IF
     
     PRINT*, 'Setting up connectivity and edge lengths...'

     DO EID = 1,NEDGE
        NODE1 = NETP%EDGENODE(EID,1)
        NODE2 = NETP%EDGENODE(EID,2)

        ! edge start, length, and (normalized) direction
        NETP%EDGESTART(EID,:) = NETP%NODEPOS(NODE1,:)
        NETP%EDGEDIR(EID,:) = NETP%NODEPOS(NODE2,:)-NETP%NODEPOS(NODE1,:)   
        CALL NORMALIZE(NETP%EDGEDIR(EID,:), LEN)
        IF (.NOT.EDGELENSET(EID))  THEN ! edge length was not set in network file
           NETP%EDGELEN(EID) = LEN
        ENDIF

        ! increment degrees of the nodes
        NETP%NODEDEG(NODE1) = NETP%NODEDEG(NODE1)+1
        NETP%NODEDEG(NODE2) = NETP%NODEDEG(NODE2)+1

        IF ( MAX(NETP%NODEDEG(NODE1),NETP%NODEDEG(NODE2)).GT.MAXBRANCH) THEN
           PRINT*, 'ERROR IN NETWORKFROMFILE: node degree exceeds maximum.',&
                & NODE1, NODE2, MAXBRANCH, NETP%NODEDEG(NODE1), NETP%NODEDEG(NODE2)
           STOP 1
        ENDIF



        ! edges connected to each node
        NETP%NODEEDGE(NODE1,NETP%NODEDEG(NODE1)) = EID
        NETP%NODEEDGE(NODE2,NETP%NODEDEG(NODE2)) = EID

        ! nodes connected to each node
        NETP%NODENODE(NODE1,NETP%NODEDEG(NODE1)) = NODE2
        NETP%NODENODE(NODE2,NETP%NODEDEG(NODE2)) = NODE1              
     END DO

     ! for each node, sort edges by length
     DO NID = 1,NETP%NNODE
        NE = NETP%NODEDEG(NID)
        LENS(1:NE) = NETP%EDGELEN(NETP%NODEEDGE(NID,1:NE))
        LENSAVE = LENS
        TMPARRAY = NETP%NODEEDGE(NID,:)        
        CALL SORT2INT(NE,LENS,TMPARRAY)
        NETP%NODEEDGE(NID,1:NE) = TMPARRAY(1:NE)
        
        TMPARRAY = NETP%NODENODE(NID,:)
        CALL SORT2INT(NE,LENSAVE,TMPARRAY)
        NETP%NODENODE(NID,1:NE) = TMPARRAY(1:NE)

        ! length of all edges leaving the node
        NETP%NODELEN(NID,1:NE)= LENS(1:NE)
     ENDDO

     ! set up reservoir info
     NETP%RESVNNODE = 0
     DO NC = 1,NETP%NNODE
        RC = NETP%NODERESV(NC)
        IF (RC.GT.0) THEN ! how many and which nodes attached to this reservoir
           NETP%RESVNNODE(RC) = NETP%RESVNNODE(RC)+1
           NETP%RESVNODES(RC,NETP%RESVNNODE(RC)) = NC
        ENDIF
     ENDDO

     PRINT*, 'Reservoir nodes:', NETP%NRESV
     DO RC = 1,NETP%NRESV
        PRINT*, RC, NETP%RESVNNODE(RC),NETP%RESVNODES(RC,1:NETP%RESVNNODE(RC))
        IF (NETP%RESVNNODE(RC).LT.1) THEN
           PRINT*, 'ERROR: reservoir with no nodes', RC
           STOP 1
        ENDIF
     ENDDO

     
     DEALLOCATE(EDGELENSET)
     NETP%STRUCTURESET = .TRUE.

     ! set all nodes to default volume
     ! Volume is actually in terms of tubule length (ie NODEVOL = V/(pi*a^2))
     ! where a is the tube radius
     NETP%NODEVOLS = NODEVOL

     ! set up edges with fixed velocities      
     NETP%FIXEDGEVEL = .FALSE.
     NETP%FIXEDGEVELVAL = 0D0
     DO CC = 1,NFIXEDGEVEL
        NETP%FIXEDGEVEL(FIXEDGEVEL(CC)) = .TRUE.
        NETP%FIXEDGEVELVAL(FIXEDGEVEL(CC)) = FIXEDGEVELVAL(CC)
    ENDDO
    
  END SUBROUTINE NETWORKFROMFILE

  SUBROUTINE OUTPUTNETWORK(NETP,NETOUTFILE,APPEND,NODEVAL,EDGEVAL)
    ! output a network file describing network structure
    ! if NODEVAL is present: write an additional value for each node
    ! if EDGEVAL is present: write an additional value for each edge
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    CHARACTER (LEN=*), INTENT(IN) :: NETOUTFILE
    LOGICAL, INTENT(IN) :: APPEND
    DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: NODEVAL(NETP%NNODE), EDGEVAL(NETP%NEDGE)
    INTEGER, PARAMETER :: FU=51
    INTEGER :: NC, EC, LC
    CHARACTER(LEN=100) :: NODEFMTSTR, NODEFMTSTRRESV
    CHARACTER(LEN=6) :: NUMSTRING

    IF (APPEND) THEN
       OPEN(UNIT=FU,FILE=NETOUTFILE,STATUS='UNKNOWN',ACCESS='APPEND')
    ELSE
       OPEN(UNIT=FU,FILE=NETOUTFILE,STATUS='UNKNOWN')
    ENDIF
    WRITE(FU,'(A)')'# network structure output by fortran code'
!    WRITE(FU,'(A)') '#list of node indices and positions'

    
    IF (PRESENT(NODEVAL)) THEN
       WRITE(NODEFMTSTR,'(A,I1,A)') '(A,1X,I5,1X,',NETP%DIM,'F20.10,F20.10)'
       WRITE(NODEFMTSTRRESV,'(A,I1,A)') '(A,1X,I5,1X,',NETP%DIM,'F20.10,F20.10,A,I1)'
    ELSE
       WRITE(NODEFMTSTR,'(A,I1,A)') '(A,1X,I5,1X,',NETP%DIM,'F20.10)'
       WRITE(NODEFMTSTRRESV,'(A,I1,A)') '(A,1X,I5,1X,',NETP%DIM,'F20.10,A,I1)'
    ENDIF
    
    DO NC = 1,NETP%NNODE
       IF (PRESENT(NODEVAL)) THEN
          IF (NETP%NODERESV(NC).GT.0) THEN
             WRITE(FU,NODEFMTSTRRESV) 'NODE ',  NC, NETP%NODEPOS(NC,:), NODEVAL(NC),' R',NETP%NODERESV(NC)
          ELSE
             WRITE(FU,NODEFMTSTR) 'NODE ',  NC, NETP%NODEPOS(NC,:), NODEVAL(NC)
          ENDIF
       ELSE
          IF (NETP%NODERESV(NC).GT.0) THEN
             WRITE(FU,NODEFMTSTRRESV) 'NODE ',  NC, NETP%NODEPOS(NC,:), ' R',NETP%NODERESV(NC)
          ELSE
             WRITE(FU,NODEFMTSTR) 'NODE ',  NC, NETP%NODEPOS(NC,:)
          ENDIF
       ENDIF
    ENDDO


    !    WRITE(FU,'(A,/)') '#list of edges connecting nodes'
    
    DO EC = 1,NETP%NEDGE
       IF (PRESENT(EDGEVAL)) THEN
          WRITE(FU,'(A,1X,I6,1X,I5,1X,I5,2F20.10)') 'EDGE ', EC, NETP%EDGENODE(EC,:), NETP%EDGELEN(EC),EDGEVAL(EC)
       ELSE
          WRITE(FU,'(A,1X,I6,1X,I5,1X,I5,F20.10)') 'EDGE ', EC, NETP%EDGENODE(EC,:), NETP%EDGELEN(EC)
       END IF

       IF (NETP%EDGENODE(EC,1).EQ.NETP%EDGENODE(EC,2)) THEN
          PRINT*, 'ERROR IN OUTPUTNETWORK: loop edge found.'
          STOP 1
       ENDIF
    ENDDO  

    DO LC = 1,NETP%NLOOP
       WRITE(NUMSTRING,'(I6)') NETP%LOOPLENS(LC)
       WRITE(FU,'(A,1X,'// TRIM(ADJUSTL(NUMSTRING))//'I6)') 'LOOP ',NETP%LOOPEDGES(LC,1:NETP%LOOPLENS(LC))
    ENDDO
    
    WRITE(FU,*) ''
    CLOSE(FU)
    
  END SUBROUTINE OUTPUTNETWORK

 
  
  SUBROUTINE SETUPNETWORK(NETP,NNODE,NEDGE,NRESV,NLOOP,DIM,NFIELD,&
       & MAXBRANCH,MAXLOOPLEN,ABSORBERS,NABS)
    ! set up a network by allocating arrays
    ! NETP: Pointer to a network
    ! NNODE: number of nodes
    ! NEDGE: number of branches
    ! NLOOP: Number of independent loops listed (should be < nedge-nnode+1 for fully connected network)
    ! DIM: spatial dimensionality where network resides
    ! MAXBRANCH: maximum branches attached to each node
    
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: NNODE, NEDGE, NRESV, NLOOP, DIM, MAXBRANCH, MAXLOOPLEN, NFIELD
    INTEGER, INTENT(IN), OPTIONAL :: ABSORBERS(:,:), NABS(:)
    INTEGER :: A, FC

    NETP%NNODE = NNODE
    NETP%NEDGE = NEDGE
    NETP%DIM = DIM
    NETP%NFIELD = NFIELD
    ! WARNING: network is assumed to have a single connected component !
    NETP%NLOOP = NLOOP ! number of independent cycles
    NETP%MAXLOOPLEN = MAXLOOPLEN    
    
    ! allocate node data
    ALLOCATE(NETP%NODENODE(NNODE,MAXBRANCH), NETP%NODEEDGE(NNODE,MAXBRANCH))
    ALLOCATE(NETP%NODEPOS(NNODE,DIM), NETP%NODEDEG(NNODE),&
         & NETP%NODELEN(NNODE,MAXBRANCH), NETP%NODEWIDTH(NNODE))    
    
    ! allocate branch data
    ALLOCATE(NETP%EDGENODE(NEDGE,2), NETP%EDGESTART(NEDGE,DIM),&
         & NETP%EDGEDIR(NEDGE,DIM), NETP%EDGELEN(NEDGE))
    ! allocate loop data
    ALLOCATE(NETP%LOOPEDGES(NETP%NLOOP,MAXLOOPLEN), NETP%LOOPLENS(NETP%NLOOP))
    
    NETP%LOOPLENS = 0
    NETP%LOOPEDGES = 0

    ! allocate reservoir data
    NETP%NRESV = NRESV
    ALLOCATE(NETP%NODERESV(NNODE),NETP%RESVNODES(NRESV,NNODE),&
         & NETP%RESVNNODE(NRESV),NETP%RESVVOL(NRESV),NETP%RESVSA(NRESV),&
         & NETP%RESVLEN(NRESV),NETP%RESVMIX(NRESV))
    NETP%NODERESV = 0
    IF (NRESV.GT.0) THEN      
       NETP%RESVNNODE = 0
       NETP%RESVMIX = .TRUE.
    ENDIF

    ALLOCATE(NETP%ISPERM(NETP%NNODE))
    NETP%ISPERM = .FALSE.
    
    ! for tracking contracting nodes
    ALLOCATE(NETP%NODEVOLS(NNODE),NETP%NODESTATE(NNODE))
    NETP%NODEVOLS = 0D0
    NETP%NODESTATE = 2 ! all nodes start fully open
    NETP%CONTSET = .FALSE. ! contractions not yet set up
 
    NETP%ARRAYSET = .TRUE.
    NETP%NODELEN = HUGE(1D0)
    NETP%NODEDEG = 0
    NETP%NODEEDGE = 0; NETP%EDGENODE = 0; NETP%NODENODE = 0
    NETP%NODEWIDTH=1D0
    
    ! absorber and fixed nodes

    !PRINT*, 'TESTX1:', NETP%nNODE
    ALLOCATE(NETP%NODEABS(NETP%NNODE,NFIELD))
    ! allocate(NETP%NODEFIX(NETP%NNODE,NFIELD),NETP%NODEFIXVAL(NETP%NNODE,NFIELD)
    !PRINT*, 'TESTX2:', NETP%nNODE
    
    NETP%NODEABS = .FALSE.
    !NETP%NODEFIX = .FALSE.
    !NETP%NODEFIXVAL = 0D0
    
    IF (PRESENT(ABSORBERS)) THEN
       DO FC = 1,NFIELD
          DO A = 1,NABS(FC)
             IF (ABSORBERS(A,FC).GT.NETP%NNODE.OR.ABSORBERS(A,FC).LT.1) THEN
                PRINT*, 'ERROR IN SETUPNETWORK: bad absorber', ABSORBERS(A,FC)
                STOP 1
             ENDIF
             NETP%NODEABS(ABSORBERS(A,FC),FC) = .TRUE.
          ENDDO
       ENDDO
    END IF
    ! IF (PRESENT(FIXNODES)) THEN
    !    DO FC = 1,NFIELD
    !       DO A = 1,NFIX(FC)
    !          IF (FIXNODES(A,FC).GT.NETP%NNODE.OR.FIXNODES(A,FC).LT.1) THEN
    !             PRINT*, 'ERROR IN SETUPNETWORK: bad fixed node', A, FIXNODES(A,FC)
    !             STOP 1
    !          ENDIF
    !          NETP%NODEFIX(FIXNODES(A,FC),FC) = .TRUE.
    !          IF (PRESENT(FIXVALS)) THEN
    !             NETP%NODEFIXVAL(FIXNODES(A,FC),FC) = FIXVALS(A,FC)            
    !          ENDIF
    !       ENDDO
    !    ENDDO
    ! ENDIF

    ! Array for mapping to mesh objects later
    ALLOCATE(NETP%EDGENCELL(NETP%NEDGE), NETP%NODECELLS(NETP%NNODE), NETP%RESVCELLS(NETP%NRESV))

    ! For fixing edge velocities
    ALLOCATE(NETP%FIXEDGEVEL(NETP%NEDGE), NETP%FIXEDGEVELVAL(NETP%NEDGE))
  END SUBROUTINE SETUPNETWORK

  SUBROUTINE CLEANUPNETWORK(NETP)
    ! deallocate arrays for the network structure
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP

    
    IF (NETP%ARRAYSET) THEN
       DEALLOCATE(NETP%NODENODE, NETP%NODEEDGE, NETP%NODEPOS, NETP%NODELEN, &
            & NETP%NODEDEG, NETP%NODEABS,NETP%NODEWIDTH)
!       DEALLOCATE(NETP%NODEFIX, NETP%NODEFIXVAL)
       DEALLOCATE(NETP%EDGENODE, NETP%EDGESTART, NETP%EDGEDIR, NETP%EDGELEN)
       DEALLOCATE(NETP%LOOPEDGES,NETP%LOOPLENS)
       DEALLOCATE(NETP%NODESTATE, NETP%NODEVOLS)
       IF (NETP%CONTSET) THEN
          DEALLOCATE(NETP%CONTNODES)
       ENDIF

       DEALLOCATE(NETP%NODERESV,NETP%RESVNODES,NETP%RESVVOL, NETP%RESVSA,&
            & NETP%RESVLEN,NETP%RESVMIX, NETP%ISPERM, NETP%RESVCELLS)

       DEALLOCATE(NETP%NODECELLS)
       DEALLOCATE(NETP%EDGENCELL)
       DEALLOCATE(NETP%FIXEDGEVEL,NETP%FIXEDGEVELVAL)
    ENDIF

    ! arrays mapping to mesh object    
    IF (NETP%MESHSET) THEN       
       DEALLOCATE(NETP%EDGECELLS)
       NETP%MESHSET = .FALSE.
    ENDIF
    
    NETP%ARRAYSET = .FALSE.
    NETP%STRUCTURESET = .FALSE.
  END SUBROUTINE CLEANUPNETWORK
END MODULE NETWORKUTIL
