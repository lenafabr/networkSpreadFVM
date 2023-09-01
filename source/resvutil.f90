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
     ! area (length in 2D of the edges)
     DOUBLE PRECISION, POINTER :: EDGEAREA(:)
     ! verticies and edges bordering the elements
     INTEGER, POINTER :: RESVVERT(:,:)     
     INTEGER, POINTER :: RESVEDGE(:,:)
     ! volume and surface area (for release) of reservoir elements
     DOUBLE PRECISION, POINTER :: RESVVOL(:), RESVSA(:)

     ! centroid position for each reservoir element
     DOUBLE PRECISION, POINTER :: RESVCENT(:,:)
     
     ! arrays are set up
     LOGICAL :: ARRAYSET = .FALSE.
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
    INTEGER :: NITEMS, NVERT, NEDGE,NRESV, DIM
    
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
        ENDIF
     ENDDO
     CLOSE(NF)

     PRINT*, 'Number of vertices, edges, reservoirs, dim: ', NVERT, NEDGE, NRESV, DIM
     ! allocate arrays         
     CALL SETUPRESERVOIRS(RESVP,NVERT,NEDGE,NRESV,DIM)

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
        ENDIF
     END DO

     CALL SETCENTROIDS(RESVP)
     
   END SUBROUTINE RESERVOIRSFROMFILE

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
           POS(IC,:) = RESVP%VERTPOS(RESVP%RESVVERT(RC,1),:)
        ENDDO
        
        RESVP%RESVCENT(RC,:) = SUM(POS,1)/NVE        
     ENDDO
   END SUBROUTINE SETCENTROIDS
  
   SUBROUTINE SETUPRESERVOIRS(RESVP, NVERT,NEDGE,NRESV,DIM)
    ! set up arrays for reservoir object

    IMPLICIT NONE
    TYPE(RESERVOIRS), POINTER :: RESVP
    INTEGER, INTENT(IN) :: NVERT,NEDGE,NRESV, DIM

    ALLOCATE(RESVP%VERTPOS(NVERT,DIM), RESVP%EDGEVERT(NEDGE,2), RESVP%EDGERESV(NEDGE,2),RESVP%EDGEAREA(NEDGE))
    ALLOCATE(RESVP%RESVVERT(NRESV,3),RESVP%RESVEDGE(NRESV,3),&
         & RESVP%RESVVOL(NRESV),RESVP%RESVSA(NRESV), RESVP%RESVCENT(NRESV,DIM))

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
    
    RESVP%ARRAYSET = .FALSE.
  END SUBROUTINE CLEANUPRESERVOIRS
END MODULE RESVUTIL
