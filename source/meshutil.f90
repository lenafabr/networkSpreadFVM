MODULE MESHUTIL
  ! utilities for dealing with a cell-centered mesh on a network
  USE NETWORKUTIL, ONLY : NETWORK
  USE RESVUTIL, ONLY : RESERVOIRS
  USE GENUTIL, ONLY : PI
  
  IMPLICIT NONE
  
  TYPE MESH
     ! object defining the mesh

     INTEGER :: NCELL ! number of cells
     INTEGER :: DIM ! dimensionality of space where mesh is embedded
     INTEGER :: MAXDEG ! maximum cell connectivity degree, sets size of BOUNDS array
     
     ! POS: absolute positions of cell centers
     ! LEN: total length of tubule within cell, used for velocity weighted average
     ! VOL: volume of cell, used to convert flux to conc
     ! for non-reservoir cell, LEN and VOL should be same
     ! for a reservoir cell, this is the volume of the cell divided by pi a^2
     ! SA: metric for surface area of cell, used for computing permeability and
     ! for weighting selection probability in RANDFIXCELLS
     ! SA is in units of length. For a tube cell, SA=LEN
     ! for reservoir is specified directly in terms of SA = (surface area)/(2*pi*a)
     ! LENPM: length to the center of each adjacent cell (h_plus, h_minus in the math notes motation)
     ! RAD: radius for each mesh cell
     ! NOTE: radius is only used if USEVARRAD keyword is on
     DOUBLE PRECISION, POINTER :: POS(:,:), LEN(:), LENPM(:,:), VOL(:)
     DOUBLE PRECISION, POINTER :: SA(:), RAD(:)
     ! number of neighbors for each cell
     INTEGER, POINTER :: DEG(:)
     ! for each edge (boundary) of the cell, list the type of boundary
     ! -1 = dirichlet (fixed value) boundary
     ! 0 = reflecting (no-flux) boundary
     ! >0 = connected to another cell 
     INTEGER, POINTER :: BOUNDS(:,:) ! list of neighbors for each cell     
     ! Bounddir is +1 if the edge direction across this boundary points
     ! outward from the cell; -1 if it points inward
     INTEGER, POINTER :: BOUNDDIR(:,:)

     ! keep track of which boundaries are closed off
     LOGICAL, POINTER :: BOUNDCLOSED(:,:)

     ! cross-sectional area of each boundary between cells
     DOUBLE PRECISION, POINTER :: BOUNDAREA(:,:)

     ! -----------
     ! mapping to a network structure
     ! -----------
     ! celltype = 0 for a cell cell corresponding to a network node
     ! celltype = 1 for a cell lying along a network edge
     ! celltype = 2 for a cell corresponding to a reservoir
     INTEGER, POINTER :: CELLTYPE(:)
     ! index of node or edge on which the cell is located
     ! EDGEIND(:,1) gives which edge, EDGEIND(:,2) gives which cell along the edge
     ! RESVIND lists which reservoir each cell belongs to; lists 0 otherwise
     INTEGER, POINTER :: NODEIND(:), EDGEIND(:,:), RESVIND(:)
     ! Boundedge gives the edge index for each boundary from this cell
     INTEGER, POINTER :: BOUNDEDGE(:,:)
     ! for cells with a reflective boundary, what terminal node does it correspond to?
     INTEGER, POINTER :: TERMNODE(:,:)
     
     ! have arrays been allocated?
     LOGICAL :: ARRAYSET = .FALSE.
     ! have the cell positions and connectivities been set?
     LOGICAL :: CELLSET = .FALSE.
     
  END type MESH
  
CONTAINS
  SUBROUTINE EDGETOBOUND(NETP,MESHP,EDGEVALS,BOUNDVALS)
    ! for a mesh on a network, copy value from network edges to boundaries of mesh cells
    IMPLICIT NONE

    TYPE(NETWORK), POINTER :: NETP
    TYPE(MESH), POINTER :: MESHP
    DOUBLE PRECISION, INTENT(IN) :: EDGEVALS(NETP%NEDGE)
    DOUBLE PRECISION, INTENT(OUT) :: BOUNDVALS(MESHP%NCELL,MESHP%MAXDEG)
    INTEGER :: CC, BC, EC
    
    BOUNDVALS = 0D0
    DO CC= 1,MESHP%NCELL       
       DO BC = 1,MESHP%DEG(CC)
          EC = MESHP%BOUNDEDGE(CC,BC)
          IF (EC.GT.0) THEN
             BOUNDVALS(CC,BC) = EDGEVALS(EC)
          ENDIF
       ENDDO
    ENDDO
    
  END SUBROUTINE EDGETOBOUND

  SUBROUTINE CLOSEBOUNDS(MESHP,CLOSERATEPERLEN,OPENRATE,EDGECLOSED,DELT)
    ! for each edge, randomly close of a mesh cell boundary along that edge    
    ! PCLOSEPERLEN = closing probability per length (use LENPM for the length surrounding each boundary)
    ! each mesh boundary is able to close or open independently
    ! POPEN = probability of opening an already closed edge
    ! DELT: if provided, then treat CLOSERATEPERLEN, OPENRATE as rates and DELT as the relevant timestep
    ! if DELT not provided, then CLOSERATEPERLEN, OPENRATE are treated directly as probabilities
    USE MT19937, ONLY : GRND
    IMPLICIT NONE
    TYPE(MESH), POINTER :: MESHP
    DOUBLE PRECISION, INTENT(IN) :: CLOSERATEPERLEN, OPENRATE
    LOGICAL, INTENT(OUT) :: EDGECLOSED(:)
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: DELT
    LOGICAL :: ISDONE(MESHP%NCELL,MESHP%MAXDEG)
    INTEGER :: CC, BC, BCC, BCC2
    DOUBLE PRECISION :: PCLOSE, POPEN
    
    ! go through and decide which boundaries to close
    
    ISDONE = .FALSE. ! Which boundaries have already been checked
    EDGECLOSED = .FALSE. ! which edges have a closure along them
    
    DO CC = 1,MESHP%NCELL
       ! only close boundaries between edge-edge or edge-node cells
       IF (MESHP%CELLTYPE(CC).EQ.2) CYCLE
       
       DO BCC = 1,MESHP%DEG(CC)          
          ! this boundary has already been checked for closing
          IF (ISDONE(CC,BCC)) CYCLE

          ! cell connected to through this boundary
          BC = MESHP%BOUNDS(CC,BCC)
          IF (BC.GT.0) THEN ! there is a neighbor cell across this boundary
             IF (MESHP%CELLTYPE(CC).EQ.2) CYCLE ! do not close boundary to reservoir
             DO BCC2 = 1,MESHP%DEG(BC) ! find back-connecting index from the neighbor cell
                IF (MESHP%BOUNDS(BC,BCC2).EQ.CC) EXIT
             ENDDO
             ISDONE(BC,BCC2) = .TRUE. ! already checked boundary for connected cell
          ENDIF

          !IF (CC.GT.1600) PRINT*, 'TESTX2:', CC, BC, BCC, BCC2, ISDONE(16347,1)
          
          IF (MESHP%BOUNDCLOSED(CC,BCC)) THEN
             ! already closed, decide whether to open
             IF (PRESENT(DELT)) THEN
                ! input rates treated as rates
                POPEN = 1D0 - EXP(-DELT*OPENRATE)
             ELSE
                ! close/open rates are actually probabilities
                POPEN = OPENRATE
             ENDIF
             
             IF (GRND().LT.POPEN) THEN
                MESHP%BOUNDCLOSED(CC,BCC) = .FALSE.
                IF (BC.GT.0) MESHP%BOUNDCLOSED(BC,BCC2) = .FALSE.
             ENDIF
          ELSE ! Currently open, decide whether to close
             ! probability of closing for this boundary
             IF (PRESENT(DELT)) THEN
                ! input rates treated as rates
                PCLOSE = 1D0 - EXP(-DELT*CLOSERATEPERLEN*MESHP%LENPM(CC,BCC))
             ELSE
                ! close/open rates are actually probabilities
                PCLOSE = CLOSERATEPERLEN*MESHP%LENPM(CC,BCC)
                !PRINT*, 'TESTX1:', CC, BCC,  BC, meshp%edgeind(cc,:), PCLOSE, MESHP%LENPM(CC,BCC)
             ENDIF

             IF (GRND().LT.PCLOSE) THEN
                MESHP%BOUNDCLOSED(CC,BCC) = .TRUE.
                ! close the boundary for the connected cell as well
                IF (BC.GT.0) MESHP%BOUNDCLOSED(BC,BCC2) = .TRUE.
             ENDIF
          ENDIF

          ! count up total boundaries and number that are closed
          !IF (MESHP%BOUNDCLOSED(CC,BCC)) NCLOSED = NCLOSED+1
          !NTOT = NTOT+1

          ! keep track of which edges have a closed boundary          
          IF (MESHP%BOUNDCLOSED(CC,BCC)) THEN
             IF (MESHP%CELLTYPE(CC).EQ.1) THEN
                EDGECLOSED(MESHP%EDGEIND(CC,1)) = .TRUE.
             ELSEIF (MESHP%CELLTYPE(BC).EQ.1) THEN
                EDGECLOSED(MESHP%EDGEIND(BC,1)) = .TRUE.
             ENDIF
          ENDIF             
       ENDDO
    ENDDO

!    FRACCLOSED = DBLE(NCLOSED)/NTOT
  END SUBROUTINE CLOSEBOUNDS
  
  SUBROUTINE CLOSEBOUNDSOLD(MESHP, PCLOSE)    
    ! randomly close off some of the mesh-cell boundaries
    ! pclose = probability each cell boundary is closed
    USE MT19937, ONLY : GRND
    IMPLICIT NONE
    TYPE(MESH), POINTER :: MESHP
    DOUBLE PRECISION, INTENT(IN) :: PCLOSE
    INTEGER :: CC, BC, BCC, BCC2, BC2, NCLOSE, NTOT, MINCC
    DOUBLE PRECISION :: U, PTRY

    ! each boundary will be counted twice for the 2 adjacent mesh cells,
    ! this is the right probability to sample so that overall prob
    ! of the boundary closing is PCLOSE
    PTRY = 1 - SQRT(1D0-PCLOSE)

    NCLOSE = 0
    NTOT = 0
    DO CC = 1,MESHP%NCELL       
       DO BCC = 1,MESHP%DEG(CC)
          BC = MESHP%BOUNDS(CC,BCC)
          IF (BC.GT.0) THEN ! this is a real boundary between cells
             NTOT = NTOT+1
             DO BCC2 = 1,MESHP%DEG(BC) ! boundary must be closed for both adjacent cells
                BC2 = MESHP%BOUNDS(BC,BCC2)
                IF (BC2.EQ.CC) EXIT
             ENDDO
             IF (.NOT.MESHP%BOUNDCLOSED(BC,BCC2)) THEN ! not already closed
                U = GRND()
                IF (U.LT.PTRY) THEN
                   MESHP%BOUNDCLOSED(CC,BCC) = .TRUE.
                   MESHP%BOUNDCLOSED(BC,BCC2) = .TRUE.

                   NCLOSE = NCLOSE + 1
                ENDIF               
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    PRINT*, 'Closed N boundaries out of total:', NCLOSE, NTOT/2, PTRY, PCLOSE
  END SUBROUTINE CLOSEBOUNDSOLD
  
  SUBROUTINE SETUPNETWORKMESH(MESHP,NETP,MAXCELLLEN,MINNCELL,CONC3D,RESVP)
    IMPLICIT NONE
    ! using a network object, set up a mesh on it for FVM simulations
    ! * CELL-CENTERED mesh. Same distance from cell position to each boundary *
    ! all nodes of deg > 1 are treated as cell centers
    ! all nodes of deg=1 are boundaries of a cell
    ! MAXCELLLEN: maximum tube length in a cell
    ! MINNCELL: minimum number of internal cells on an edge; does not include nodal cells
    ! this subroutine allocates all arrays, including the mapping from network to mesh indices
    ! Does NOT set diffusivity or velocities on the mesh
    ! Network object must be fully set up already
    ! CONC3D: work with 3D concentrations? default is 1D
    ! OPTIONAL: RESVP input is a pointer to a reservoirs object desribing interconnected
    ! reservoir elements to be included in the mesh
    
    TYPE(MESH), POINTER :: MESHP
    TYPE(NETWORK), POINTER :: NETP
    TYPE(RESERVOIRS), POINTER, OPTIONAL :: RESVP
    DOUBLE PRECISION, INTENT(IN) :: MAXCELLLEN
    INTEGER, INTENT(IN) :: MINNCELL
    logical, INTENT(IN) :: CONC3D
    
    INTEGER :: EC, NC, CC, CT, N1, N2, D1, D2, BC, RC, bct, CCT
    INTEGER :: NCELL
    DOUBLE PRECISION :: CX1, CX2, ELEN, MINELEN, CELLLEN
    DOUBLE PRECISION :: EDGECELLLEN(NETP%NEDGE), NODECELLLEN(NETP%NNODE)

    INTEGER :: NINTNODES, NCELLTOT, DIM, MAXDEG, INTNODES(NETP%NEDGE,2), ESHORT, MINCC

    DOUBLE PRECISION :: EXTLEN
    
    IF (.NOT.NETP%ARRAYSET.OR..NOT.NETP%STRUCTURESET) THEN
       PRINT*, 'ERROR IN SETUPMESH: Network not fully set up', NETP%ARRAYSET, NETP%STRUCTURESET
       STOP 1
    ENDIF

    IF (NETP%MESHSET) THEN
       PRINT*, 'ERROR IN SETUPMESH: mesh mapping for this network has already been set up. Need to erase.'
       STOP 1
    ENDIF

    NODECELLLEN = 0D0
    EDGECELLLEN = 0D0

    ! label each edge according to whether its nodes are terminal or not
    INTNODES = 0
    DO EC = 1,NETP%NEDGE
       N1 = NETP%EDGENODE(EC,1)
       N2 = NETP%EDGENODE(EC,2)
    
       IF (NETP%NODEDEG(N1).GT.1) INTNODES(EC,1) = 1
       IF (NETP%NODEDEG(N2).GT.1) INTNODES(EC,2) = 1
    ENDDO
       
    ! Go through each node, and decide on the length the corresponding cell
    ! will extend along each edge
    DO NC = 1,NETP%NNODE
       D1 = NETP%NODEDEG(NC)
       ! ignore terminal nodes; they do not get their own cell
       IF (D1.EQ.1) CYCLE

       ! find shortest edge adjacent to this node
       MINELEN =HUGE(1D0)
       ESHORT = 0
       DO EC = 1,D1
          ELEN = NETP%EDGELEN(NETP%NODEEDGE(NC,EC))
          IF (ELEN.LT.MINELEN) THEN
             MINELEN = ELEN
             ESHORT = EC
          ENDIF
       ENDDO

       ! decide on cell size based on the shortest attached edge
       ! CELLLEN is the approximated cell size on this edge
       ! Use it to set the total cell size at the node, equally split over attached edges
       IF (SUM(INTNODES(ESHORT,:)).EQ.0.AND.MINNCELL.EQ.0) THEN
          ! two terminal nodes; must have at least one internal cell
          CELLLEN = MINELEN
       ELSE
          CELLLEN = MINELEN/(MINNCELL + 0.5*SUM(INTNODES(ESHORT,:)))
       ENDIF
       IF (CELLLEN.GT.MAXCELLLEN) CELLLEN=MAXCELLLEN
       
       NODECELLLEN(NC) = CELLLEN       
    ENDDO
    
    ! go through each edge, decide how many cells and of what type and size
    DO EC = 1,NETP%NEDGE
       N1 = NETP%EDGENODE(EC,1)
       N2 = NETP%EDGENODE(EC,2)
       D1 = NETP%NODEDEG(N1); D2 = NETP%NODEDEG(N2)

       ! how much of the edge length goes into nodal cells?
       CX1 = 0D0; CX2 = 0D0       
       IF (D1.GT.1) CX1 = NODECELLLEN(N1)/D1
       IF (D2.GT.1) CX2 = NODECELLLEN(N2)/D2          
       
       ! how many internal cells could fit in here (assuming largest allowed)
       NCELL = CEILING((NETP%EDGELEN(EC)-CX1-CX2)/MAXCELLLEN)      
       NCELL = MAX(NCELL,MINNCELL)
       
       
       ! number of internal cells for each edge
       NETP%EDGENCELL(EC) = NCELL       

       IF (NCELL.GT.0) THEN
          ! determine internal cell size on this edge       
          EDGECELLLEN(EC) = (NETP%EDGELEN(EC)-CX1-CX2)/NCELL
       ELSE
          EDGECELLLEN(EC) = 0D0 ! No internal cells
       ENDIF

    ENDDO

    
    ! ----------
    ! allocate arrays for mesh object
    ! ----------
    ! number of nodes of degree > 1 (interior network nodes)
    NINTNODES = 0    
    DO NC = 1,NETP%NNODE
       IF (NETP%NODEDEG(NC).GT.1) NINTNODES = NINTNODES + 1
    ENDDO

    DIM = NETP%DIM
    IF (PRESENT(RESVP)) THEN
       IF (DIM.NE.RESVP%DIM) THEN
          PRINT*, 'ERROR: network and reservoirs are in different dimensions', DIM, RESVP%DIM
       ENDIF
    ENDIF
    
    ! each reservoir gets its own cell
    IF (PRESENT(RESVP)) THEN
       ! reservoirs are defined in separate data structures
       NCELLTOT = SUM(NETP%EDGENCELL) + NINTNODES + RESVP%NRESV
       MAXDEG = MAX(MAXVAL(NETP%NODEDEG),MAXVAL(NETP%RESVNNODE),MAXVAL(RESVP%RESVDEG))
    ELSE  
       NCELLTOT = SUM(NETP%EDGENCELL) + NINTNODES + NETP%NRESV
       MAXDEG = MAX(MAXVAL(NETP%NODEDEG),MAXVAL(NETP%RESVNNODE))
    ENDIF
    
    CALL ALLOCATEMESH(MESHP,NCELLTOT,DIM,MAXDEG)

    ! Allocate network arrays mapping to cells
    NETP%MAXCELLEDGE = MAXVAL(NETP%EDGENCELL)
    ALLOCATE(NETP%EDGECELLS(NETP%NEDGE,NETP%MAXCELLEDGE))

    NETP%NODECELLS = 0; NETP%EDGECELLS= 0
    NETP%RESVCELLS = 0
    MESHP%NODEIND = 0; MESHP%EDGEIND = 0; MESHP%RESVIND = 0

    ! ----------
    ! set up reservoir cells defined in network file that are NOT defined in reservoir file
    ! ----------
    CT= 0 ! total cell counter

    IF (PRESENT(RESVP)) THEN
       ! leave space for explicitly defined reservoirs
       ! WARNING: cannot have connection defined in network file to a reservoir
       ! defined in the reservoir file
       CT = RESVP%NRESV
       MESHP%CELLTYPE(1:CT) = 2
       MESHP%DEG(1:CT) = 0
    ENDIF

    ! these are *separate* implicitly defined reservoirs in network file
    ! set up based on network nodes marked as connecting to reservoir
    ! This does *not* include connections defined in the explicit reservoir file
    DO RC = 1,NETP%NRESV

       IF (RC.LE.CT) THEN ! this reservoir will be explicitly defined in reservoir file
          CYCLE
       ENDIF
       
       CT = CT+1;

       NETP%RESVCELLs(RC) = CT ! map from reservoir to cell


       MESHP%CELLTYPE(CT) = 2 ! reservoir cell
       MESHP%RESVIND(CT) = RC ! which reservoir it belongs to

       ! volume of reservoir (in terms of edge length; this is really V/(pi a^2)       
       MESHP%VOL(CT) = NETP%RESVVOL(RC)

       ! surface area of reservoir (in terms of edge length; this is really V/(2*pi*a)
       MESHP%SA(CT) = NETP%RESVSA(RC)

       ! number of connected cells
       D1 = NETP%RESVNNODE(RC)
       MESHP%DEG(CT) = D1
!       print*, 'TESTX1:', CT, D1
       ! place position of cell at average of connected nodes
       MESHP%POS(CT,:)= SUM(NETP%NODEPOS(NETP%RESVNODES(RC,1:D1),:),1)/D1
    ENDDO    
    
    ! -------------
    ! Set up cells on nodes
    ! --------------
    DO NC = 1,NETP%NNODE
       D1 = NETP%NODEDEG(NC)
       IF (D1.EQ.1) CYCLE ! no cell on terminal nodes

       CT = CT+1
       ! set node cell positions and lengths
       MESHP%POS(CT,:) = NETP%NODEPOS(NC,:)
       MESHP%LEN(CT) = NODECELLLEN(NC)
       MESHP%VOL(CT) = MESHP%LEN(CT)
       MESHP%SA(CT) = MESHP%LEN(CT)
       
       ! Mapping from mesh to network indices
       MESHP%CELLTYPE(CT) = 0 ! this is a cell on a node
       MESHP%DEG(CT) = D1 ! how many other cells it connects to
       ! increase degree if connected to reservoir
       IF (NETP%NODERESV(NC).GT.0) MESHP%DEG(CT) = MESHP%DEG(CT)+1
       MESHP%NODEIND(CT) = NC ! corresponding network node index

       ! Mapping from network to mesh indices
       NETP%NODECELLS(NC) = CT
    ENDDO

    ! ----------
    ! set up cells on edges
    ! -----------
    DO EC = 1,NETP%NEDGE       
       N1 = NETP%EDGENODE(EC,1); N2 = NETP%EDGENODE(EC,2)
       
       DO CC = 1,NETP%EDGENCELL(EC)

          CT = CT+1      
          ! cell position and length
          CX1 = NODECELLLEN(N1)/NETP%NODEDEG(N1)
          MESHP%POS(CT,:) =  NETP%NODEPOS(N1,:) + &
               & (CX1+((CC-1) + 0.5D0)*EDGECELLLEN(EC))*NETP%EDGEDIR(EC,:)
          MESHP%LEN(CT) = EDGECELLLEN(EC)
          MESHP%VOL(CT) = MESHP%LEN(CT)
          MESHP%SA(CT) = MESHP%LEN(CT)
          
          ! distance to center of neighbor cell
         ! MESHP%LENPM(CT,1:2) = EDGECELLLEN(EC)
          
          ! Mapping from mesh to network indices
          MESHP%CELLTYPE(CT) = 1 ! this is a cell on an edge
          MESHP%DEG(CT) = 2
          MESHP%EDGEIND(CT,1) = EC
          MESHP%EDGEIND(CT,2) = CC

          ! Mapping from network to mesh indices
          NETP%EDGECELLS(EC,CC) = CT

          ! boundaries for the cell
          IF (CC.EQ.1) THEN
             IF (NETP%NODEDEG(N1).EQ.1) THEN
                RC = NETP%NODERESV(N1)
                IF (RC.GT.0) THEN
                   ! reservoir node                   
                   MESHP%BOUNDS(CT,1) = NETP%RESVCELLS(RC)                   
                ELSE
                   ! terminal node, no flux boundary
                   MESHP%BOUNDS(CT,1) = 0
                  ! MESHP%DEG(CT) = 1
                ENDIF
                MESHP%TERMNODE(CT,1) = N1
             ELSE
                ! connect to a node cell
                MESHP%BOUNDS(CT,1) = NETP%NODECELLS(N1)
!                MESHP%LENPM(CT,1) = EDGECELLLEN(EC)/2 + NODECELLLEN(NETP%NODECELLS(N1))/D1
             ENDIF
          ELSE ! connect interior cell to the previous one
             MESHP%BOUNDS(CT,1) = CT-1
!             MESHP%LENPM(CT,1) = EDGECELLLEN(EC)
          ENDIF          

          IF (CC.EQ.NETP%EDGENCELL(EC)) THEN
             !IF (Ct.EQ.199) PRINT*, 'TESTX1:', N2, NETP%NODECELLS(N2)
             IF (NETP%NODEDEG(N2).EQ.1) THEN
                RC = NETP%NODERESV(N2)
                IF (RC.GT.0) THEN
                   ! reservoir node
                   MESHP%BOUNDS(CT,2) = NETP%RESVCELLS(RC)
                ELSE
                   ! no flux boundary
                   MESHP%BOUNDS(CT,2) = 0
                  ! MESHP%DEG(CT) = 1
                END IF
                MESHP%TERMNODE(CT,2) = N2
             ELSE
                ! connect to a node cell
                MESHP%BOUNDS(CT,2) = NETP%NODECELLS(N2)
 !               MESHP%LENPM(CT,2) = EDGECELLLEN(EC)/2 + NODECELLLEN(NETP%NODECELLS(N2))/NETP%NODEDEG(N2)
             ENDIF
          ELSE ! connect interior cell to next one
             MESHP%BOUNDS(CT,2) = CT+1
 !            MESHP%LENPM(CT,2) = EDGECELLLEN(EC)
          ENDIF          
          ! edge direction points inward along first boundary, outward along second
          MESHP%BOUNDDIR(CT,1:2) = (/-1,1/)
          MESHP%BOUNDEDGE(CT,1:2) = EC
       ENDDO
    ENDDO

   
    IF (CT.NE.MESHP%NCELL) THEN
       PRINT*, 'ERROR IN SETUPMESH: wrong number of cells', MESHP%NCELL, CT
       STOP 1
    ENDIF

    IF (.NOT.PRESENT(RESVP)) THEN
    ! set up connectivity for reservoirs defined implicitly in the network
    DO RC = 1,NETP%NRESV
       ! mesh cell belonging to this reservoir
       CT = NETP%RESVCELLS(RC)

       ! Cells bordering the reservoir
       DO CC = 1,NETP%RESVNNODE(RC)
          ! node belonging to this reservoir
          NC = NETP%RESVNODES(RC,CC)

          ! which network node connects to this reservoir mesh cell
          MESHP%TERMNODE(CT,CC) = NC
          
          IF (NETP%NODECELLS(NC).EQ.0) THEN
             ! terminal node bordering reservoir

             IF (NETP%NODEDEG(NC).NE.1) THEN
                PRINT*, 'ERROR IN SETUPMESH: bad terminal node at reservoir.'
                PRINT*, RC, CC, NC, NETP%NODEDEG(NC)
                STOP 1
             ENDIF
             
             ! attach first/last cell on the associated edge
             EC = NETP%NODEEDGE(NC,1)
             MESHP%BOUNDEDGE(CT,CC) = EC
             IF (NETP%EDGENODE(EC,1).EQ.NC) THEN
                MESHP%BOUNDS(CT,CC) = NETP%EDGECELLS(EC,1)
                MESHP%BOUNDDIR(CT,CC) = 1 ! edge points out of reservoir
             ELSEIF (NETP%EDGENODE(EC,2).EQ.NC) THEN
                MESHP%BOUNDS(CT,CC) = NETP%EDGECELLS(EC,NETP%EDGENCELL(EC))
                MESHP%BOUNDDIR(CT,CC) = -1 ! edge points into reservoir
             ELSE
                PRINT*, 'ERROR IN SETUPMESH: bad node/edge connectivity'
                PRINT*, RC, CC, NC, EC, NETP%NODEEDGE(NC,1), NETP%EDGENODE(EC,:)
                STOP 1                
             ENDIF
             
          ELSE
             print*, 'ERROR IN SETUP MESH: currently cannot handle reservoir nodes &
               & that have degree > 1. Not sure how to define flow velocities &
               & on the connection to the reservoir in that case.'
             print*, RC, CT, NC, CC
             STOP 1
          ENDIF
       ENDDO       
    ENDDO
    ENDIF
    
    ! set up connectivity for node cells
    DO NC = 1,NETP%NNODE
       IF (NETP%NODEDEG(NC).EQ.1) CYCLE ! no separate cell for terminal nodes
       CT = NETP%NODECELLS(NC)
       DO CC = 1,NETP%NODEDEG(NC)
          EC = NETP%NODEEDGE(NC,CC)
          IF (NETP%EDGENCELL(EC).GE.1) THEN
             ! there is at least one internal cell; connect the nodal cell to it
             IF (NETP%EDGENODE(EC,1).EQ.NC) THEN
                ! this is the first node of an edge
                MESHP%BOUNDS(CT,CC) = NETP%EDGECELLS(EC,1)
             ELSEIF (NETP%EDGENODE(EC,2).EQ.NC) THEN
                ! this is the second node of an edge
                MESHP%BOUNDS(CT,CC) = NETP%EDGECELLS(EC,NETP%EDGENCELL(EC))
             ELSE
                PRINT*, 'ERROR IN SETUPMESH: nodes and edges do not match'
                PRINT*, 'NODEEDGE:', NETP%NODEEDGE(NC,:)
                PRINT*, 'EDGENODE: ', NETP%EDGENODE(EC,:)
                STOP 1
             ENDIF
 !            MESHP%LENPM(CT,CC) = MESHP%LEN(CT)/NETP%NODEDEG(NC) + EDGECELLLEN(EC)/2
          ELSE
             ! there are no internal cells
             N2 = NETP%NODENODE(NC,CC)
             IF (NETP%NODEDEG(N2).EQ.1) THEN
                ! no-flux boundary
                MESHP%BOUNDS(CT,CC) = 0
                MESHP%TERMNODE(CT,CC) = N2
             ELSE
                ! connect directly to the next node
                MESHP%BOUNDS(CT,CC) = NETP%NODECELLS(N2)
             ENDIF

             !            MESHP%LENPM(CT,CC) = NETP%EDGELEN(EC)
          ENDIF

          ! Set up boundary directions
          MESHP%BOUNDEDGE(CT,CC) = EC
          IF (NETP%EDGENODE(EC,1).EQ.NC) THEN
             MESHP%BOUNDDIR(CT,CC) = 1 ! edge points out from node
          ELSEIF (NETP%EDGENODE(eC,2).EQ.NC) THEN
             MESHP%BOUNDDIR(CT,CC) = -1 ! edge points into node
          ELSE
             PRINT*, 'ERROR IN SETUPMESH: EDGES AND NODES DONT MATCH', NETP%EDGENODE(EC,:), NC
             STOP 1
          ENDIF
       ENDDO
       IF (MESHP%RESVIND(CT).GT.0) THEN
          ! this node is associated with a reservoir
          print*, 'ERROR IN SETUP MESH: currently cannot handle reservoir nodes &
               & that have degree > 1. Not sure how to define flow velocities &
               & on the connection to the reservoir in that case.'
          PRINT*, CT, MESHP%RESVIND(CT), NC, NETP%NODEDEG(NC)
          STOP 1
          !CC = NETP%NODEDEG(NC)+1
          !RC = MESHP%RESVIND(CT)
          !MESHP%BOUNDS(CT,CC) = NETP%RESVCELLS(RC)
       ENDIF
    ENDDO

    ! get appropriate length (beyond end nodes) for reservoir
    ! for reservoir cell, length = avg cell lngth on all attached edges
    ! NOTE: this length is not the cell volume, but rather a length
    ! used in weighted avg for velocity and flux calculations

    ! from diffusive flux calculations to narrow exits (1/19/2023 notes)
    ! effective RESVLEN should be:
    ! for spheres: L = pi*r/2 (r = tube radius)
    ! for sheets: L = r^2/h*ln(R/r) (h=sheet thickness, pi*R^2 = A = sheet area)
    IF (PRESENT(RESVP)) THEN
       MINCC = RESVP%NRESV+1
    ELSE
       MINCC = 1
    ENDIF
     
    DO CC = MINCC,MESHP%NCELL
       IF (MESHP%CELLTYPE(CC).EQ.2) THEN ! implicit reservoir cells only
          RC = MESHP%RESVIND(CC) ! which reservoir is this?
          MESHP%LEN(CC) = NETP%RESVLEN(RC)/2*MESHP%DEG(CC)
       ENDIF
    ENDDO

    ! check that length is defined for all cells
    DO CC = 1,MESHP%NCELL
       IF (MESHP%LEN(CC).LE.0D0) THEN
          PRINT*, 'ERROR: negative length found.', CC, MESHP%CELLTYPE(CC)
          STOP 1
       ENDIF
    ENDDO   
    
    ! For each cell, get the length to each neighboring cell
    DO CT = 1,MESHP%NCELL
       DO CC = 1,MESHP%DEG(CT)
          BC = MESHP%BOUNDS(CT,CC)
          !IF (MESHP%CELLTYPE(CT).EQ.2) CYCLE


          IF (BC.LE.0) THEN
             MESHP%LENPM(CT,CC) = MESHP%LEN(CT)/MESHP%DEG(CT) ! length to reflecting bound
          !ELSEIF (MESHP%CELLTYPE(BC).EQ.2) THEN ! boundary to reservoir     
          !   MESHP%LENPM(CT,CC) = 2*MESHP%LEN(CT)/MESHP%DEG(CT)
          ELSE ! boundary to other cell OR to reservoir
             MESHP%LENPM(CT,CC) = MESHP%LEN(CT)/MESHP%DEG(CT) + MESHP%LEN(BC)/MESHP%DEG(BC)
          ENDIF
       ENDDO
    ENDDO

    IF (PRESENT(RESVP)) THEN
       ! Update explicit reservoir elements
       ! this will assume the first NRESV mesh cells are reservoirs
       CALL RESERVOIRSTOMESH(RESVP,MESHP,NETP)
    ENDIF

    ! default mesh cell radius
    !MESHP%RAD = SQRT(1D0/PI)

    ! set mesh radii
    CALL SETMESHRADII(MESHP,NETP)        
    CALL SETMESHBOUNDAREAS(MESHP,NETP)
    
    ! set up to work with 3D concentrations
    IF (CONC3D) THEN
       CALL UPDATEMESH3D(MESHP)
    ENDIF
    
    MESHP%CELLSET = .TRUE.
    NETP%MESHSET = .TRUE.    
    
    ! DEBUG: double check the LENPM matches up to LENs for bounding cells
    PRINT*, 'Testing that LENPM are self-consistent...'
    DO CC = 1,MESHP%NCELL
       DO BCT = 1,MESHP%DEG(CC)
          BC = MESHP%BOUNDS(CC,BCT)
          IF (BC.EQ.0) CYCLE

          CT =0
          DO CCT = 1,MESHP%DEG(BC)
             IF (MESHP%BOUNDS(BC,CCT).EQ.CC) THEN
                CT = CCT
                EXIT
             ENDIF
          ENDDO

          IF (CT.EQ.0) THEN
             PRINT*, 'BOUNDARY MISMATCH'
             print*, CC, MESHP%DEG(CC), MESHP%CELLTYPE(CC), BCT
             PRINT*, BC, MESHP%DEG(BC), MESHP%CELLTYPE(BC), MESHP%BOUNDS(BC,:)
             STOP 1
          ENDIF
          IF (ABS(MESHP%LENPM(CC,BCT)-MESHP%LENPM(BC,CT)).GT.1D-10) THEN
             PRINT*, 'ERROR IN CHECKING LENPM:', CC, BC
             PRINT*, MESHP%BOUNDS(CC,BCT),MESHP%LENPM(CC,BCT)
             PRINT*, MESHP%BOUNDS(BC,CT), MESHP%LENPM(BC,CT)
             STOP 1
          ENDIF
       ENDDO
    ENDDO

    ! DEBUG: check that mesh boundary areas are self-consistent (same from both sides of boundary)
    PRINT*, 'Testing that mesh boundary areas are self-consistent'
    DO CC = 1,MESHP%NCELL
       DO BCT = 1,MESHP%DEG(CC)
          BC = MESHP%BOUNDS(CC,BCT)
          IF (BC.EQ.0) CYCLE

          CT =0
          DO CCT = 1,MESHP%DEG(BC)
             IF (MESHP%BOUNDS(BC,CCT).EQ.CC) THEN
                CT = CCT
                EXIT
             ENDIF
          ENDDO
          
          IF (ABS(MESHP%BOUNDAREA(CC,BCT) - MESHP%BOUNDAREA(BC,CCT)).GT.1D-10) THEN
             PRINT*, 'ERROR: boundary area mismatch', CC, BCT, BC, CCT, MESHP%BOUNDAREA(CC,BCT),  MESHP%BOUNDAREA(BC,CCT)
             STOP 1
          ENDIF
       ENDDO
    ENDDO

    
  END SUBROUTINE SETUPNETWORKMESH  

  SUBROUTINE RESERVOIRSTOMESH(RESVP,MESHP,NETP)
    ! Update mesh object info for reservoir cells, based on explicit reservoirs
    IMPLICIT NONE

    TYPE(RESERVOIRS), POINTER :: RESVP
    TYPE(MESH), POINTER :: MESHP
    TYPE(NETWORK), POINTER :: NETP
    INTEGER :: CC, BC, IND, RC, RC2, NETDEG, IC, CC2, TERM, ct
    INTEGER :: EC, V1, V2, MEC, NC
    DOUBLE PRECISION :: ECENT(MESHP%DIM), INTERLEN
    DOUBLE PRECISION :: DIFF(RESVP%DIM)

    ! update info for the reservoir cells
    DO RC = 1,RESVP%NRESV
       CC = RC ! index among mesh cells
       MESHP%RESVIND(CC) = RC ! index for this reservoir
       RESVP%RESVCELL(RC) = CC

       ! volume of reservoir (in terms of edge length; this is really V/(pi a^2)
       MESHP%VOL(CC) = RESVP%RESVVOL(RC)

       ! surface area of reservoir (in terms of edge length; this is really V/(2*pi*a)
       MESHP%SA(CC) = RESVP%RESVSA(RC)

       ! effective length used to calculate narrow escape flux into network tube
       ! for spheres: L = pi*r/2 (r = tube radius)
       ! for sheets: L = r^2/h*ln(R/r) (h=sheet thickness, pi*R^2 = A = sheet area)
       ! NOTE: this should never actually be used, since connecting nodes explicitly through edge
       MESHP%LEN(CC) = RESVP%RESVLENEFF(RC)

       ! number of connected cells
       ! Add on resv-resv connections on top of node-resv connections
       MESHP%DEG(CC) = RESVP%RESVDEG(RC)

       ! place position of cell at centroid
       MESHP%POS(CC,:)= RESVP%RESVCENT(RC,1:MESHP%DIM)

       ! Set up boundaries to connected nodes
       DO IC = 1,RESVP%RESVNNODE(RC)

          ! this node is connected to our reservoir
          NC = RESVP%RESVCONNODE(RC,IC,1)
          IF (NETP%NODEDEG(NC).GT.1) THEN
             PRINT*, 'ERROR: can only connect terminal nodes to explicit reservoir'
             STOP 1
          ENDIF
          EC = NETP%NODEEDGE(NC, 1) ! edge going out of this node

          ! find the mesh cell leading to this node
          DO CC2 = 1,MESHP%NCELL
             IF (MESHP%CELLTYPE(CC2).EQ.1.and.MESHP%TERMNODE(CC2,1).EQ.NC) THEN
                TERM = 1 ! which side of mesh cell is the terminal node
                EXIT
             ELSEIF (MESHP%CELLTYPE(CC2).EQ.1.and.MESHP%TERMNODE(CC2,2).EQ.NC) THEN
                TERM = 2
                EXIT
             ENDIF
          ENDDO
          IF (CC2.GT.MESHP%NCELL) THEN
             PRINT*, 'Failed to find mesh cell terminating at connected node'
             STOP 1
          ENDIF

          ! find the reservoir mesh edge for the connection
          DO BC = 1,RESVP%RESVDEG(RC)
             IF (RESVP%RESVEDGE(RC,BC).EQ.RESVP%RESVCONNODE(RC,IC,2)) EXIT
          ENDDO
          IF (BC.GT.RESVP%RESVDEG(RC)) THEN
             PRINT*, 'ERROR: reservoir does not actually abut this edge', RC, IC,&
                  & RESVP%RESVCONNODE(RC,IC,2), RESVP%RESVEDGE(RC,:)
             STOP 1
          ENDIF

          ! link the network mesh cell to the reservoir mesh cell
          MESHP%BOUNDS(CC,BC) = CC2
          MESHP%BOUNDS(CC2,TERM) = CC

          ! boundary direction always outward from reservoir, into network
          MESHP%BOUNDDIR(CC,BC) = 1
          MESHP%BOUNDDIR(CC2,TERM) = -1
          MESHP%BOUNDEDGE(CC,BC) = EC
          MESHP%BOUNDEDGE(CC2,TERM) = EC

          ! length between cell centers. Actually, take length to the point halfway
          ! along the connecting mesh-edge + length along the network tube
          MEC = RESVP%RESVEDGE(RC,BC)
          V1 = RESVP%EDGEVERT(MEC,1); V2 = RESVP%EDGEVERT(MEC,2)
          ECENT = (RESVP%VERTPOS(V1,:) + RESVP%VERTPOS(V2,:))/2

          INTERLEN = SQRT(SUM((RESVP%RESVCENT(RC,:)-ECENT))**2)
          INTERLEN = INTERLEN + MESHP%LEN(CC2)/2                    
          MESHP%LENPM(CC,BC) = INTERLEN
          MESHP%LENPM(CC2,TERM) = INTERLEN
       ENDDO
       
       ! set up boundaries between reservoir cells
       DO BC = 1,RESVP%RESVDEG(RC)
          RC2 = RESVP%RESVRESV(RC,BC)

          IF (RC2.EQ.0) THEN
             CYCLE ! reflecting boundary or connected to node
          ENDIF
          MESHP%BOUNDS(CC,BC) = RC2
          ! area of the boundary to the reservoir edge
          MESHP%BOUNDAREA(CC,BC) = RESVP%EDGEAREA(RESVP%RESVEDGE(CC,BC))

          ! boundary directions defined to always point towards higher index
          IF (RC2.GT.RC) THEN
             MESHP%BOUNDDIR(CC,BC) = 1
          ELSE
             MESHP%BOUNDDIR(CC,BC) = -1
          ENDIF

          ! no edge associated with resv-resv boundary
          MESHP%BOUNDEDGE(CC,BC) = 0

          ! length between cell centers          
          DIFF = RESVP%RESVCENT(RC,:) - RESVP%RESVCENT(RC2,:)
          MESHP%LENPM(CC,BC) = SQRT(SUM(DIFF**2))
       ENDDO
       
    ENDDO

    RESVP%MESHSET = .TRUE.
    
  END SUBROUTINE RESERVOIRSTOMESH


  
  SUBROUTINE SETMESHRADII(MESHP,NETP)
    ! Set up radii for each mesh cell.
    ! Also change volume of each non-reservoir mesh cell to be in terms of actual volume units.
    ! node and reservoir radii are set to average of the connected edge radii
    ! for now, sets radii based on radius of network edge
    ! TODO: UPDATE TO ALLOW SINUSOIDALLY VARYING RADII ALONG AN EDGE
    
    IMPLICIT NONE
    TYPE(MESH), POINTER :: MESHP
    TYPE(NETWORK), POINTER :: NETP
    INTEGER :: CC, NC, EC, RC, ECC, NCC, CT
    DOUBLE PRECISION :: TOTRAD
    LOGICAL :: NETCON
    INTEGER :: BC, BCT

    DO CC = 1,MESHP%NCELL
       SELECT CASE (MESHP%CELLTYPE(CC))
       CASE (0) ! node cell radius (avg of surrounding edges)
          NC = MESHP%NODEIND(CC)

          ! get average radii of nearby edges
          TOTRAD = 0D0
          DO ECC = 1,NETP%NODEDEG(NC)
             TOTRAD = TOTRAD + NETP%EDGERAD(NETP%NODEEDGE(NC,ECC))
          ENDDO
          MESHP%RAD(CC) = TOTRAD/NETP%NODEDEG(NC)
       CASE (1) ! edge cell
          EC = MESHP%EDGEIND(CC,1)
          MESHP%RAD(CC) = NETP%EDGERAD(EC)
       CASE (2) ! reservoir cell
          RC = MESHP%RESVIND(CC)

          ! check if this reservoir is connected to any nodes or edges of the network
          NETCON = .FALSE.
          DO BCT = 1,MESHP%DEG(CC)
             BC = MESHP%BOUNDS(CC,BCT)
             IF (BC.GT.0) THEN
                IF (MESHP%CELLTYPE(BC).EQ.1.OR.MESHP%CELLTYPE(BC).EQ.0) THEN
                   NETCON = .TRUE.
                   EXIT
                ENDIF
             ENDIF
          ENDDO

          IF (NETCON) THEN
             ! get average radii of nearby edges
             TOTRAD = 0D0
             CT = 0
             DO NCC = 1,NETP%RESVNNODE(RC)
                NC = NETP%RESVNODES(RC,NCC) ! node belonging to this reservoir
                DO ECC = 1,NETP%NODEDEG(NC) ! go over all edges from this node
                   TOTRAD = TOTRAD + NETP%EDGERAD(NETP%NODEEDGE(NC,ECC))
                   CT = CT+1
                ENDDO
             ENDDO
             MESHP%RAD(CC) = TOTRAD/CT
          ELSE
             ! no network nodes are attached to this reservoir
             ! set radius to arbitrary value, should never be used
             MESHP%RAD(CC) = 0D0
          ENDIF

       END SELECT
       
    ENDDO
       
  END SUBROUTINE SETMESHRADII

  SUBROUTINE SETMESHBOUNDAREAS(MESHP,NETP)
    ! set areas of boundaries between mesh elements
    ! assumes radii of edge cells have been predefined
    ! ignores boundaries between two reservoir elements
    ! (these should be set in RESERVOIRSTOMESH)

    IMPLICIT NONE
    
    TYPE(MESH), POINTER :: MESHP
    TYPE(NETWORK), POINTER :: NETP
    INTEGER :: CC, CC2, NC, NC2, BC, DC, EC
    
    DO CC = 1,MESHP%NCELL
       SELECT CASE (MESHP%CELLTYPE(CC))
       CASE (0) ! node cell
          NC = MESHP%NODEIND(CC)
          
          ! set each boundary to cross-sectional area of corresponding edge
          DO BC = 1,MESHP%DEG(CC)
             CC2 = MESHP%BOUNDS(CC,BC) ! boundary cell

             IF (CC2.EQ.0) THEN ! reflecting boundary
                ! find edge leading to the right terminal node
                NC= MESHP%TERMNODE(CC,BC)
                EC = NETP%NODEEDGE(NC,1)
                MESHP%BOUNDAREA(CC,BC) = PI*NETP%EDGERAD(EC)**2
                CYCLE
             ENDIF
                
             IF (MESHP%CELLTYPE(CC2).EQ.1) THEN ! boundary leads to edge                
                MESHP%BOUNDAREA(CC,BC) = PI*MESHP%RAD(CC2)**2
             ELSEIF (MESHP%CELLTYPE(CC2).EQ.0) THEN ! boundary leads to another node
                NC2 = MESHP%NODEIND(CC2)
                
                DO DC = 1,NETP%NODEDEG(NC)
                   EC = NETP%NODEEDGE(NC,DC) ! edge adjacent to this node
                   IF (NETP%EDGENODE(EC,1).EQ.NC2.OR.NETP%EDGENODE(EC,2).EQ.NC2) THEN
                      ! this edge leads to the correct bordering node
                      MESHP%BOUNDAREA(CC,BC) = PI*NETP%EDGERAD(EC)**2
                      EXIT
                   ENDIF
                ENDDO
             ELSEIF (MESHP%CELLTYPE(CC2).EQ.2) THEN ! boundary leads to a reservoir
                PRINT*, 'ERROR IN SETMESHBOUNDAREAS: should not have network node directly connected to reservoir', CC, BC, CC2
                STOP 1
             ELSE
                PRINT*, 'ERROR IN SETMESHBOUNDAREAS: this should never happen', CC, BC, CC2
                STOP 1
             ENDIF
          ENDDO
       CASE (1) ! edge cell
          ! set both boundaries to cross-sectional area of the edge
          MESHP%BOUNDAREA(CC,1) = PI*MESHP%RAD(CC)**2
          MESHP%BOUNDAREA(CC,2) = MESHP%BOUNDAREA(CC,1)
       CASE (2) ! reservoir cell
          DO BC = 1,MESHP%DEG(CC)
             CC2 = MESHP%BOUNDS(CC,BC) ! boundary cell
             IF (CC2.EQ.0) CYCLE
             
             IF (MESHP%CELLTYPE(CC2).EQ.1) THEN
                ! edge cell connected to this reservoir
                MESHP%BOUNDAREA(CC,BC) = PI*MESHP%RAD(CC2)**2
             ELSEIF (MESHP%CELLTYPE(CC2).EQ.0) THEN
                PRINT*, 'ERROR IN SETMESHBOUNDAREAS 2: should not have network node directly connected to reservoir', CC, BC, CC2
                STOP 1
             ENDIF
          ENDDO
       END SELECT
    ENDDO
    
    
  END SUBROUTINE SETMESHBOUNDAREAS

  SUBROUTINE UPDATEMESH3D(MESHP)
    ! update node and edge cells on mesh to deal with 3D concentrations
    ! just changes the values of VOL and SA for each mesh cell on network node or edge
    ! does not do anything to reservoir cells
    IMPLICIT NONE
    TYPE(MESH), POINTER :: MESHP
    INTEGER :: CC

    DO CC = 1,MESHP%NCELL
       IF (MESHP%CELLTYPE(CC).EQ.0.OR.MESHP%CELLTYPE(CC).EQ.1) THEN
          MESHP%VOL(CC) = MESHP%LEN(CC)*PI*MESHP%RAD(CC)**2
          MESHP%SA(CC) = MESHP%LEN(CC)*2*PI*MESHP%RAD(CC)
       ENDIF
    END DO
    
  END SUBROUTINE UPDATEMESH3D
  
  SUBROUTINE OUTPUTMESH(MESHP,OUTFILE)
    ! output mesh structure to a file (for loading in matlab)
    IMPLICIT NONE

    TYPE(MESH), POINTER :: MESHP
    CHARACTER(LEN=*), INTENT(IN) :: OUTFILE
    INTEGER, PARAMETER :: OU=99
    INTEGER :: CT, DEG
    
    OPEN(UNIT=OU, FILE=OUTFILE)

    ! output general info
    WRITE(OU, *) MESHP%NCELL, MESHP%DIM
    
    DO CT = 1,MESHP%NCELL
       DEG = MESHP%DEG(CT)
       ! for each cell output: index, position, degree, total cell length, bounding cells, reservoir
       WRITE(OU,*) CT, MESHP%POS(CT,:), MESHP%VOL(CT), DEG, &
            & MESHP%BOUNDS(CT,1:DEG), MESHP%NODEIND(CT), MESHP%EDGEIND(CT,1:2), MESHP%RESVind(CT), MESHP%RAD(CT)
    ENDDO
    
    CLOSE(OU)
    
  END SUBROUTINE OUTPUTMESH
  
  SUBROUTINE ALLOCATEMESH(MESHP,NCELLTOT,DIM,MAXDEG)
    ! Allocate arrays for a mesh object
    ! NCELLTOT = total cells
    ! DIM = spatial dimension in which its embedded
    ! MAXDEG = maximum connection degree of the mesh

    IMPLICIT NONE
    TYPE(MESH), POINTER :: MESHP
    INTEGER, INTENT(IN) :: NCELLTOT, DIM, MAXDEG

    MESHP%NCELL = NCELLTOT
    MESHP%DIM = DIM
    MESHP%MAXDEG = MAXDEG
    
    ALLOCATE(MESHP%POS(NCELLTOT,DIM), MESHP%LEN(NCELLTOT), MESHP%VOL(NCELLTOT),&
         & MESHP%LENPM(NCELLTOT,MAXDEG), MESHP%SA(NCELLTOT))
    ALLOCATE(MESHP%DEG(NCELLTOT),MESHP%BOUNDS(NCELLTOT,MAXDEG))
    ALLOCATE(MESHP%CELLTYPE(NCELLTOT),MESHP%NODEIND(NCELLTOT), MESHP%EDGEIND(NCELLTOT,2), MESHP%TERMNODE(NCELLTOT,MAXDEG))
    ALLOCATE(MESHP%BOUNDDIR(NCELLTOT,MAXDEG), MESHP%BOUNDEDGE(NCELLTOT,MAXDEG), MESHP%BOUNDCLOSED(NCELLTOT,MAXDEG))
    ALLOCATE(MESHP%RESVIND(NCELLTOT), MESHP%RAD(NCELLTOT), MESHP%BOUNDAREA(NCELLTOT,MAXDEG))
    
    MESHP%ARRAYSET = .TRUE.
    
    MESHP%BOUNDS = 0
    MESHP%BOUNDEDGE = 0
    MESHP%BOUNDDIR = 0
    MESHP%DEG = 0
    MESHP%LEN=-1D0

    MESHP%BOUNDCLOSED = .FALSE.
  END SUBROUTINE ALLOCATEMESH

   SUBROUTINE CLEANUPMESH(MESHP)
    ! Deallocate arrays for a mesh object

    IMPLICIT NONE
    TYPE(MESH), POINTER :: MESHP
    
    DEALLOCATE(MESHP%POS, MESHP%LEN, MESHP%LENPM, MESHP%VOL, MESHP%SA)
    DEALLOCATE(MESHP%DEG,MESHP%BOUNDS)
    DEALLOCATE(MESHP%CELLTYPE,MESHP%NODEIND, MESHP%EDGEIND)
    DEALLOCATE(MESHP%TERMNODE, MESHP%BOUNDDIR, MESHP%BOUNDEDGE, MESHP%BOUNDCLOSED)
    DEALLOCATE(MESHP%RESVIND, MESHP%RAD, MESHP%BOUNDAREA)
    
    MESHP%ARRAYSET = .FALSE.
    MESHP%CELLSET = .FALSE.
  END SUBROUTINE CLEANUPMESH
END MODULE MESHUTIL
