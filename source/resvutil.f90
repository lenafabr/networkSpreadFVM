MODULE RESVUTIL
  ! utilities for dealing with many reservoirs, mesh elements in 2D
  ! currently, only 2-dim case set up. 
  
  TYPE RESERVOIRS
     ! object defining a set of interconnected reservoirs (for meshing higher dimensional domains)
     
     ! number of reservoir elements, vertices edges
     INTEGER :: NRESV, NVERT, NEDGE
     INTEGER :: DIM ! dimension in space
     ! vertex position in space
     DOUBLE PRECISION, POINTER :: VERTPOS(:,:)
     ! vertices and reservoirs bordering each edge
     INTEGER, POINTER :: EDGEVERT(:,:)
     INTEGER, POINTER :: EDGERESV(:,:)
     ! cross-sectional area  of the edges
     DOUBLE PRECISION, POINTER :: EDGEAREA(:)
     ! verticies and edges bordering the elements
     INTEGER, POINTER :: RESVVERT(:,:)     
     INTEGER, POINTER :: RESVEDGE(:,:)
     ! volume and surface area (for release) of reservoir elements
     DOUBLE PRECISION, POINTER :: RESVVOL(:), RESVSA(:)
     ! degree of connection for each reservoir
     INTEGER, POINTER :: RESVDEG(:)
     ! mesh elements corresponding to each reservoir
     ! neighboring reservoirs connected to each reservoir
     INTEGER, POINTER :: RESVCELL(:), RESVRESV(:,:)

     ! for each reservoir: number of nodes connected to it
     ! which node, and through which edge (numbered as 1,2,3 around reservoir)
     INTEGER, POINTER :: RESVNNODE(:)
     INTEGER, POINTER :: RESVCONNODE(:,:,:)
     
     ! centroid position for each reservoir element
     DOUBLE PRECISION, POINTER :: RESVCENT(:,:)

     ! effective length used to calculate narrow escape flux into network tube
     ! for spheres: L = pi*r/2 (r = tube radius)
     ! for sheets: L = r^2/h*ln(R/r) (h=sheet thickness, pi*R^2 = A = sheet area)
     DOUBLE PRECISION, POINTER :: RESVLENEFF(:)
     
     ! arrays are set up
     LOGICAL :: ARRAYSET = .FALSE.

     ! mesh elements set up?
     LOGICAL :: MESHSET = .FALSE.
  END TYPE RESERVOIRS

CONTAINS
 
  SUBROUTINE RESERVOIRSFROMFILE(RESVP,FILENAME)
    ! read a text file containing info on reservoirs corresponding to mesh elements (for meshing 2D domains)
    USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI
    
    IMPLICIT NONE
    TYPE(RESERVOIRS), POINTER :: RESVP
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    INTEGER :: NF, I, NID, VID, EID, RID
    LOGICAL :: LDUM
    LOGICAL :: FILEEND=.FALSE.
    CHARACTER(LEN=100) :: WORD, LBL, STR
    INTEGER :: NITEMS, NVERT, NEDGE,NRESV, DIM, CT, IC
    ! max number of reservoirs defined for this system. Can set huge, since only used once
    INTEGER, PARAMETER :: TMPMAXNRESV = 1D8
    INTEGER :: NNODE(TMPMAXNRESV)

    NNODE = 0
    
    ! go through file and count number of nodes, edges, reservoirs
    PRINT*, 'Reading reservoir structure file: ', FILENAME
    INQUIRE(FILE=FILENAME,EXIST=LDUM)
    IF (.NOT.LDUM) THEN
       PRINT*, 'ERROR in RESERVOIRFROMFILE: file ', TRIM(ADJUSTL(FILENAME)), ' does not exist.'
        STOP 1
     ENDIF
     OPEN(UNIT=NF, FILE=FILENAME, STATUS='OLD')

     NVERT = 0; NEDGE = 0; NRESV = 0; DIM = 0
     DO 
        CALL READLINE(NF,FILEEND,NITEMS)        
        
        IF (FILEEND.and.nitems.eq.0) EXIT
        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE
        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)

        ! Skip any comment lines
        IF (WORD(1:1).EQ.'#') CYCLE
        
        IF (WORD.EQ.'VERT') THEN
           CALL READI(VID)
           IF (VID.GT.NVERT) NVERT = VID

           DIM = NITEMS-2
        ELSEIF (WORD.EQ.'EDGE') THEN
           CALL READI(EID)
           IF (EID.GT.NEDGE) NEDGE = EID
        ELSEIF (WORD.EQ.'RESV') THEN
           CALL READI(RID)
           IF (RID.GT.NRESV) NRESV = RID
        ELSEIF (WORD.EQ.'NODECON') THEN
           ! count how many nodes connected to each reservoir
           CALL READI(NID)
           CALL READI(RID)
           IF (RID.GT.TMPMAXNRESV) THEN
              PRINT*, 'Too many reservoirs, need to increase TMPMAXNRESV'
              STOP 1
           ENDIF
           NNODE(RID) = NNODE(RID)+1
        ELSE
           PRINT*, 'Unknown keyword in reservoir file:', WORD
        ENDIF
     ENDDO
     CLOSE(NF)

     PRINT*, 'Number of vertices, edges, reservoirs, dim: ', NVERT, NEDGE, NRESV, DIM
     ! allocate arrays         
     CALL SETUPRESERVOIRS(RESVP,NVERT,NEDGE,NRESV,DIM,MAXVAL(NNODE))

     ! Now go through and read in all the reservoir info
     OPEN(UNIT=NF, FILE=FILENAME, STATUS='OLD')     
     DO 
        CALL READLINE(NF,FILEEND,NITEMS)
        IF (FILEEND.and.nitems.eq.0) EXIT
        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE
        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)
        ! Skip any comment lines
        IF (WORD(1:1).EQ.'#') CYCLE

        IF (WORD.EQ.'VERT') THEN
           CALL READI(VID)
           DO I = 1,2              
              CALL READF(RESVP%VERTPOS(VID,I))
           ENDDO
        ELSEIF (WORD.EQ.'EDGE') THEN
           CALL READI(EID)
           DO I = 1,2
              CALL READI(RESVP%EDGEVERT(EID,I))
           ENDDO
           DO I = 1,2
              CALL READI(RESVP%EDGERESV(EID,I))
           ENDDO
           CALL READF(RESVP%EDGEAREA(EID))
        ELSEIF (WORD.EQ.'RESV') THEN
           CALL READI(RID)
           DO I = 1,3
              CALL READI(RESVP%RESVVERT(RID,I))
           ENDDO
           DO I = 1,3
              CALL READI(RESVP%RESVEDGE(RID,I))
           ENDDO

           CALL READF(RESVP%RESVVOL(RID))
           CALL READF(RESVP%RESVSA(RID))
           CALL READF(RESVP%RESVLENEFF(RID))
           
            ! define degree for each reservoir (number of other reservoirs connected to it)
           ! For now this is always set to 3 (2D triangles).          
           ! should generalize at some point
           ! boundaries may lead to nodes, other reservoirs, or be reflective
           RESVP%RESVDEG(RID) = 3

        ELSEIF (WORD.EQ.'NODECON') THEN
           ! network node to be connected to reservoir through a full mesh edge
           CALL READI(NID)
           CALL READI(RID)
           CALL READI(EID)

           RESVP%RESVNNODE(RID) =  RESVP%RESVNNODE(RID) + 1
           ! save the node for this connection
           RESVP%RESVCONNODE(RID,RESVP%RESVNNODE(RID),1) = NID
           ! save the mesh edge for this connection
           RESVP%RESVCONNODE(RID,RESVP%RESVNNODE(RID),2) = EID
        ELSE
           PRINT*, 'Unknown keyword in reservoir file:', WORD
        ENDIF
     END DO

     CALL SETCENTROIDS(RESVP)
     CALL SETRESVNEIGHB(RESVP)
   END SUBROUTINE RESERVOIRSFROMFILE

   SUBROUTINE SETRESVNEIGHB(RESVP)
     ! Set up the neighbor reservoir through each edge
     IMPLICIT NONE
     TYPE(RESERVOIRS), POINTER :: RESVP
     INTEGER :: RC, ECC, EC, I, RC2

     ! 0 by default if there is no neighbor
     RESVP%RESVRESV = 0
         
     DO RC = 1,RESVP%NRESV
        DO ECC = 1,3
           EC = RESVP%RESVEDGE(RC,ECC)
           DO I = 1,2
              RC2 = RESVP%EDGERESV(EC,I)
              IF (RC2.GT.0.AND.RC2.NE.RC) THEN
                 RESVP%RESVRESV(RC,ECC) = RC2
              ENDIF
           ENDDO          
        ENDDO        
     ENDDO
   END SUBROUTINE SETRESVNEIGHB
   
   SUBROUTINE SETCENTROIDS(RESVP)
     ! set centroids for reservoir elements
     IMPLICIT NONE
     TYPE(RESERVOIRS), POINTER :: RESVP
     INTEGER ::  RC, IC
     ! number of vertices for each element
     INTEGER, PARAMETER :: NVE = 3
     DOUBLE PRECISION :: POS(NVE,RESVP%DIM)
     
     DO RC = 1,RESVP%NRESV
        DO IC = 1,NVE
           POS(IC,:) = RESVP%VERTPOS(RESVP%RESVVERT(RC,IC),:)
        ENDDO
        
        RESVP%RESVCENT(RC,:) = SUM(POS,1)/NVE        
     ENDDO
   END SUBROUTINE SETCENTROIDS
  
   SUBROUTINE SETUPRESERVOIRS(RESVP, NVERT,NEDGE,NRESV,DIM,MAXNNODE)
     ! set up arrays for reservoir object
     ! MAXNNODE: max number of nodes connected to a single reservoir

    IMPLICIT NONE
    TYPE(RESERVOIRS), POINTER :: RESVP
    INTEGER, INTENT(IN) :: NVERT,NEDGE,NRESV, DIM, MAXNNODE

    ALLOCATE(RESVP%VERTPOS(NVERT,DIM), RESVP%EDGEVERT(NEDGE,2), RESVP%EDGERESV(NEDGE,2),RESVP%EDGEAREA(NEDGE))
    ALLOCATE(RESVP%RESVVERT(NRESV,3),RESVP%RESVEDGE(NRESV,3),&
         & RESVP%RESVVOL(NRESV),RESVP%RESVSA(NRESV), RESVP%RESVCENT(NRESV,DIM))
    ALLOCATE(RESVP%RESVDEG(NRESV), RESVP%RESVCELL(NRESV), RESVP%RESVLENEFF(NRESV))
    ALLOCATE(RESVP%RESVRESV(NRESV,3))
    ALLOCATE(RESVP%RESVNNODE(NRESV), RESVP%RESVCONNODE(NRESV,MAXNNODE,2))
    
    RESVP%NVERT = NVERT
    RESVP%NEDGE = NEDGE
    RESVP%NRESV = NRESV
    RESVP%DIM = DIM

    RESVP%ARRAYSET = .TRUE.
    
  END SUBROUTINE SETUPRESERVOIRS

  SUBROUTINE CLEANUPRESERVOIRS(RESVP)
    IMPLICIT NONE
    TYPE(RESERVOIRS), POINTER :: RESVP

    DEALLOCATE(RESVP%VERTPOS, RESVP%EDGEVERT, RESVP%EDGERESV, RESVP%EDGEAREA)
    DEALLOCATE(RESVP%RESVVERT, RESVP%RESVEDGE, RESVP%RESVVOL, RESVP%RESVSA, RESVP%RESVCENT)
    DEALLOCATE(RESVP%RESVDEG, RESVP%RESVCELL, RESVP%RESVLENEFF, RESVP%RESVRESV)
    DEALLOCATE(RESVP%RESVNNODE, RESVP%RESVCONNODE)
    
    RESVP%ARRAYSET = .FALSE.
    RESVP%MESHSET = .FALSE.
  END SUBROUTINE CLEANUPRESERVOIRS
END MODULE RESVUTIL
