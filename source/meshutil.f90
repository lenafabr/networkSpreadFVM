MODULE MESHUTIL
  ! utilities for dealing with a cell-centered mesh on a network
  USE NETWORKUTIL, ONLY : NETWORK
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

  SUBROUTINE CLOSEBOUNDS(MESHP, PCLOSE)    
    ! randomly close off some of the mesh-cell boundaries
    ! pclose = probability each cell boundary is closed
    USE MT19937, ONLY : GRND
    IMPLICIT NONE
    TYPE(MESH), POINTER :: MESHP
    DOUBLE PRECISION, INTENT(IN) :: PCLOSE
    INTEGER :: CC, BC, BCC, BCC2, BC2, NCLOSE, NTOT
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
  END SUBROUTINE CLOSEBOUNDS
  
  SUBROUTINE SETUPNETWORKMESH(MESHP,NETP,MAXCELLLEN,MINNCELL)
    IMPLICIT NONE
    ! using a network object, set up a mesh on it for FVM simulations
    ! * CELL-CENTERED mesh. Same distance from cell position to each boundary *
    ! all nodes of deg > 1 are treated as cell centers
    ! all nodes of deg=1 are boundaries of a cell
    ! MAXCELLLEN: maximum tube length in a cell
    ! MINNCELL: minimum number of internal cells on an edge; does not include nodal cells
    ! this subroutine allocates all arrays, including the mapping from network to mesh indices
    ! Does NOT set diffusivity of velocities on the mesh
    ! Network object must be fully set up already
    
    TYPE(MESH), POINTER :: MESHP
    TYPE(NETWORK), POINTER :: NETP
    DOUBLE PRECISION, INTENT(IN) :: MAXCELLLEN
    INTEGER, INTENT(IN) :: MINNCELL
    
    INTEGER :: EC, NC, CC, CT, N1, N2, D1, D2, BC, RC, bct, CCT
    INTEGER :: NCELL
    DOUBLE PRECISION :: CX1, CX2, ELEN, MINELEN, CELLLEN
    DOUBLE PRECISION :: EDGECELLLEN(NETP%NEDGE), NODECELLLEN(NETP%NNODE)

    INTEGER :: NINTNODES, NCELLTOT, DIM, MAXDEG, INTNODES(NETP%NEDGE,2), ESHORT

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

    ! each reservoir gets its own cell
    NCELLTOT = SUM(NETP%EDGENCELL) + NINTNODES + NETP%NRESV
    DIM = NETP%DIM
    MAXDEG = MAX(MAXVAL(NETP%NODEDEG),MAXVAL(NETP%RESVNNODE))
    
    CALL ALLOCATEMESH(MESHP,NCELLTOT,DIM,MAXDEG)

    ! Allocate network arrays mapping to cells
    NETP%MAXCELLEDGE = MAXVAL(NETP%EDGENCELL)
    ALLOCATE(NETP%EDGECELLS(NETP%NEDGE,NETP%MAXCELLEDGE))

    NETP%NODECELLS = 0; NETP%EDGECELLS= 0
    NETP%RESVCELLS = 0
    MESHP%NODEIND = 0; MESHP%EDGEIND = 0; MESHP%RESVIND = 0

    ! ----------
    ! set up reservoir cells
    ! ----------
    CT= 0 ! total cell counter
    
    DO RC = 1,NETP%NRESV
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
    DO CC = 1,MESHP%NCELL
       IF (MESHP%CELLTYPE(CC).EQ.2) THEN ! reservoir cells only
          RC = MESHP%RESVIND(CC) ! which reservoir is this?
          MESHP%LEN(CT) = NETP%RESVLEN(RC)/2*MESHP%DEG(CC)
       ENDIF
    ENDDO

    ! --------------------
    ! OLD inaccurate calculations
   !  DO CC = 1,MESHP%NCELL
!        IF (MESHP%CELLTYPE(CC).EQ.2) THEN
! !          MESHP%LEN(CC) = MESHP%VOL(CC)
          
!           EXTLEN = 0D0; CT = 0
!           DO BCT = 1,MESHP%DEG(CC) ! look over all boundary cells
!              BC = MESHP%BOUNDS(CC,BCT)             
!              IF (BC.GT.0) THEN
!                 EXTLEN = EXTLEN + MESHP%LEN(BC)/MESHP%DEG(BC)
!                 CT = CT+1 
!              ENDIF
!           ENDDO
!           ! multiply by degree so that in the LENPM calculation that cancels out
!           MESHP%LEN(CC) = MESHP%DEG(CC)*EXTLEN/CT                  
!        ENDIF       
!     ENDDO
    ! -----------------------
    
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

    ! default mesh cell radius
    MESHP%RAD = SQRT(1D0/PI)
    
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
  END SUBROUTINE SETUPNETWORKMESH  
  
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
       END SELECT
       
    ENDDO
       
  END SUBROUTINE SETMESHRADII
  
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
            & MESHP%BOUNDS(CT,1:DEG), MESHP%NODEIND(CT), MESHP%EDGEIND(CT,1:2), MESHP%RESVind(CT)
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
    ALLOCATE(MESHP%RESVIND(NCELLTOT), MESHP%RAD(NCELLTOT))
    
    MESHP%ARRAYSET = .TRUE.
    
    MESHP%BOUNDS = 0
    MESHP%BOUNDEDGE = 0
    MESHP%BOUNDDIR = 0

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
    DEALLOCATE(MESHP%RESVIND, MESHP%RAD)
    
    MESHP%ARRAYSET = .FALSE.
    MESHP%CELLSET = .FALSE.
  END SUBROUTINE CLEANUPMESH
END MODULE MESHUTIL
