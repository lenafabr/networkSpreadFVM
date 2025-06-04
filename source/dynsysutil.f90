MODULE DYNSYSUTIL
  ! utilities for defining the dynamical system we are trying to evolve

  USE NETWORKUTIL, ONLY : NETWORK
  USE MESHUTIL, ONLY : MESH
  
  IMPLICIT NONE

  TYPE DYNSYSTEM
     INTEGER :: NFIELD=0 ! number of fields
     ! spatially uniform diffusion coefficient for each field
     DOUBLE PRECISION, POINTER :: DCOEFF(:)
               
     ! pointer to a mesh object
     TYPE(MESH), POINTER :: MESHP
     ! flow velocities across each cell boundary
     ! dimensions: which cell, which boundary
     DOUBLE PRECISION, POINTER :: VEL(:,:)

     ! FIELDS(CC,FC) = value of field FC on cell CC
     DOUBLE PRECISION, POINTER :: FIELDS(:,:)

     ! Fixed values for specific cells
     LOGICAL, POINTER :: ISFIXED(:,:)
     DOUBLE PRECISION, POINTER :: FIXVALS(:,:)

     ! field is allowed to move with flow and diffusion
     LOGICAL, POINTER :: MOBILEFIELD(:)
     
     ! ------------
     ! For dealing with buffers
     ! ------------

     ! buffer type: 0 = no buffer, 1 = explicit on/off rates,
     ! 2 = equilibrated with total buffer sites tracked as a field
     INTEGER :: BUFFERTYPE=0
     ! Assume spatially uniform buffer?
     LOGICAL :: UNIFORMBUFFER
     LOGICAL :: DOACTIVATION 
     ! explicit on/off rates or equilibrium dissociation constant
     DOUBLE PRECISION :: KON, KOFF, KDEQUIL

     ! arrays have been allocated
     LOGICAL :: ARRAYSET = .FALSE.

     ! if >0, use varying radii for mesh elements
     ! integer value controls the correction factor for diffusivities
     ! 0: no correction
     ! 1: Eq 1.18 in Berezhkovskii 2007; from Reguerra, 2001
     INTEGER :: VARRAD
     
     ! parameters have been set
     LOGICAL :: PARAMSET = .FALSE.

     ! ---------
     ! local activation and depletion of tracked molecules
     ! ---------
     DOUBLE PRECISION, POINTER :: ACTRATE(:), DEPRATE(:)
     
     ! --------
     ! For dealing with permeable nodes
     ! ---------
     LOGICAL, POINTER :: ISPERM(:)
     DOUBLE PRECISION, POINTER :: PERM(:,:), CEXT(:)
     INTEGER :: NPERM

     ! Include flows along edges?
     LOGICAL :: USEEDGEFLOW
     ! treat concentrations as 3D?
     ! if true, will use the RAD field of the mesh to track
     ! tubule radius
     LOGICAL :: CONC3D

     ! periodic global permeability
     DOUBLE PRECISION :: PERIODGLOBALPERM, DURGLOBALPERM

     ! dealing with global reservoir with pumping kinetics
     ! rate and saturation conc for recovery from global reservoir
     DOUBLE PRECISION :: GLOBALRESVKR, GLOBALRESVKMR
     ! rate and saturation conc for pumping out of global reservoir
     DOUBLE PRECISION :: GLOBALRESVKOUT, GLOBALRESVKMOUT
     LOGICAL :: PERMTOGLOBALRESV

     ! track total ion concentration directly and back-solve for free ions
     LOGICAL :: TRACKDCDT
  END type DYNSYSTEM
  
CONTAINS

  SUBROUTINE GETINITCONC(NFIELD,STARTCONC,STARTEQUIL,KD,INITCONC)
    INTEGER, INTENT(IN) :: NFIELD, STARTEQUIL
    DOUBLE PRECISION, INTENT(IN) :: STARTCONC(:), KD
    DOUBLE PRECISION, INTENT(OUT) :: INITCONC(NFIELD)
    DOUBLE PRECISION :: FREELIG, TOTPROT
    
    ! get initial concentrations:
    ! if STARTEQUIL=0, get directly from startconc
    ! if STARTEQUIL=1, start equilibrated, with free ligand and total protein fixed

    IF (STARTEQUIL.EQ.0) THEN
       INITCONC = STARTCONC(1:NFIELD)
    ELSEIF (STARTEQUIL.EQ.1) THEN
       FREELIG = STARTCONC(1)
       TOTPROT = STARTCONC(2)+STARTCONC(3)

       INITCONC(1) = FREELIG
       INITCONC(2) = FREELIG*TOTPROT/(KD + FREELIG)
       INITCONC(3) = TOTPROT - INITCONC(2)
    ELSE
       PRINT*, 'ERROR IN GETINITCONC: STARTEQUIL must be 0 or 1', STARTEQUIL
       STOP 1              
    ENDIF

    PRINT*, 'Initial concentrations:', initconc
    
  END SUBROUTINE GETINITCONC
  
  SUBROUTINE INTEGRATEFIELD(DSP, INTFIELD,TOTLEN,INCLUDERESV)
    ! get integral of field over all cells in network
    ! TOTLEN = total length
    ! INCLUDERESV = include reservoirs in calculation?
    IMPLICIT NONE
    TYPE(DYNSYSTEM), POINTER :: DSP
    DOUBLE PRECISION, INTENT(OUT) :: INTFIELD(DSP%NFIELD)
    DOUBLE PRECISION, INTENT(OUT) :: TOTLEN
    LOGICAL, INTENT(IN) :: INCLUDERESV
    INTEGER :: FC, CC
    
    DO FC = 1,DSP%NFIELD
       INTFIELD(FC) = SUM(DSP%FIELDS(:,FC)*DSP%MESHP%VOL)       
    END DO
    TOTLEN = SUM(DSP%MESHP%VOL)

    
    IF (.NOT.INCLUDERESV.AND.ANY(DSP%MESHP%CELLTYPE.EQ.2)) THEN
       ! do not integrate over reservoirs
       DO CC= 1,DSP%MESHP%NCELL
          IF (DSP%MESHP%CELLTYPE(CC).EQ.2) THEN
             INTFIELD = INTFIELD - DSP%FIELDS(CC,:)*DSP%MESHP%VOL(CC)
             TOTLEN = TOTLEN - DSP%MESHP%VOL(CC)
          ENDIF
       ENDDO
    ENDIF
    
    
  END SUBROUTINE INTEGRATEFIELD
  
  SUBROUTINE INITIALIZEFIELDNODES(NETP,DSP,STARTNODES,STARTCONC)
    ! initialize field to a constant nonzero value
    ! over specific nodes
    ! if startconc < 0, normalize to make field integrate to 1
    
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    TYPE(DYNSYSTEM), POINTER :: DSP
    INTEGER, INTENT(IN) :: STARTNODES(:)
    DOUBLE PRECISION, INTENT(IN) :: STARTCONC(DSP%NFIELD)
    INTEGER :: NSTARTNODE, FC, NC, CC
    DOUBLE PRECISION :: TOTLEN
    
    NSTARTNODE = SIZE(STARTNODES)
    !DSP%FIELDS = 0D0

   
    IF (ANY(NETP%NODEDEG(STARTNODES).EQ.1)) THEN
       PRINT*, 'ERROR IN INITIALIZEFIELDNODES: cannot initialize on terminal nodes'
       PRINT*, 'STARTNODES:', STARTNODES
       STOP 1
    ENDIF
    
    DO FC = 1,DSP%NFIELD ! for each field
       IF (STARTCONC(FC).GE.0) THEN         
          DSP%FIELDS(NETP%NODECELLS(STARTNODES),FC) = STARTCONC(FC)
       ELSE
          ! total volume of starting cells
          TOTLEN = SUM(DSP%MESHP%VOL(NETP%NODECELLS(STARTNODES)))
          IF (DSP%CONC3D) THEN
             DSP%FIELDS(NETP%NODECELLS(STARTNODES),FC) = 1D0/TOTLEN
          ELSE
             DSP%FIELDS(NETP%NODECELLS(STARTNODES),FC) = &
                  & DSP%MESHP%VOL(NETP%NODECELLS(STARTNODES))/DSP%MESHP%LEN(NETP%NODECELLS(STARTNODES))/TOTLEN
          ENDIF
       ENDIF
    ENDDO

    ! Set fixed cell values
    DO FC = 1,DSP%NFIELD
       DO CC = 1,DSP%MESHP%NCELL
          IF (DSP%ISFIXED(CC,FC)) DSP%FIELDS(CC,FC) = DSP%FIXVALS(CC,FC)
       ENDDO
    ENDDO
    
  END SUBROUTINE INITIALIZEFIELDNODES

  SUBROUTINE INITIALIZEFIELDNEARPos(NETP,DSP,STARTPOS,STARTCONC, RAD,BACKGROUNDCONC)
    ! initialize field to a constant nonzero value
    ! on all cells within a particular radius of a set of starting positions
    ! if startconc < 0, normalize to make field integrate to 1
    
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    TYPE(DYNSYSTEM), POINTER :: DSP
    DOUBLE PRECISION, INTENT(IN) :: STARTPOS(:,:)    
    DOUBLE PRECISION, INTENT(IN) :: STARTCONC(DSP%NFIELD), RAD
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: BACKGROUNDCONC(DSP%NFIELD)
    INTEGER :: NSTARTPOS, FC, NC, CC
    DOUBLE PRECISION :: TOTLEN, DIST2
    LOGICAL :: STARTCELLS(DSP%MESHP%NCELL)
    
    
    NSTARTPOS = SIZE(STARTPOS,1)
    DSP%FIELDS = 0D0

    IF (PRESENT(BACKGROUNDCONC)) THEN
       DO FC = 1,DSP%NFIELD
          DSP%FIELDS(:,FC) = BACKGROUNDCONC(FC)
       ENDDO
    ENDIF
    
    ! IF (ANY(NETP%NODEDEG(STARTNODES).EQ.1)) THEN
    !    PRINT*, 'ERROR IN INITIALIZEFIELDNODES: cannot initialize on terminal nodes'
    !    STOP 1
    ! ENDIF

    ! mark which cells to include at the start
    TOTLEN = 0
    STARTCELLS = .FALSE.
    DO CC = 1,DSP%MESHP%NCELL
       DO NC = 1,NSTARTPOS
          DIST2 = SUM((DSP%MESHP%POS(CC,:) - STARTPOS(NC,:))**2)         
          IF (DIST2.LT.RAD**2) THEN
             STARTCELLS(CC) = .TRUE.
             PRINT*, 'Cell included in starting concentration:', CC
             TOTLEN = TOTLEN + DSP%MESHP%VOL(CC)
          ENDIF
       ENDDO
    ENDDO

        
    DO FC = 1,DSP%NFIELD ! for each field
       IF (STARTCONC(FC).GE.0) THEN
          DO CC = 1,DSP%MESHP%NCELL
             IF (STARTCELLS(CC)) DSP%FIELDS(CC,FC) = STARTCONC(FC)
          ENDDO
       ELSE
          DO CC = 1,DSP%MESHP%NCELL
             IF (STARTCELLS(CC)) THEN
                IF (DSP%CONC3D) THEN
                   DSP%FIELDS(CC,FC) = 1D0/TOTLEN
                ELSE
                   DSP%FIELDS(CC,FC) = DSP%MESHP%VOL(CC)/DSP%MESHP%LEN(CC)/TOTLEN
                ENDIF
             ENDIF
          ENDDO        
       ENDIF
    ENDDO

    ! Set fixed cell values
    DO FC = 1,DSP%NFIELD
       DO CC = 1,DSP%MESHP%NCELL
          IF (DSP%ISFIXED(CC,FC)) DSP%FIELDS(CC,FC) = DSP%FIXVALS(CC,FC)
       ENDDO
    ENDDO
    
  END SUBROUTINE INITIALIZEFIELDNEARPOS
  
  SUBROUTINE INITIALIZEFIELDEDGES(NETP,DSP,STARTEDGES,STARTCONC)
    ! initialize field to a constant nonzero value
    ! on all cells along particular edges
    ! if startconc < 0, normalize to make field integrate to 1
    
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    TYPE(DYNSYSTEM), POINTER :: DSP
    INTEGER, INTENT(IN) :: STARTEDGES(:)    
    DOUBLE PRECISION, INTENT(IN) :: STARTCONC(DSP%NFIELD)
    INTEGER :: NSTARTEDGE, FC, NC, CC, EC, NCELL, BC, CC2
    DOUBLE PRECISION :: TOTLEN, DIST2, TMP(DSP%NFIELD), FRACCELLS(DSP%MESHP%NCELL)
    LOGICAL :: STARTCELLS(DSP%MESHP%NCELL), WHICHEDGES(NETP%NEDGE), INCLUDERESV
    TYPE(MESH), POINTER :: MESHP

    MESHP=>DSP%MESHP
    
    NSTARTEDGE = SIZE(STARTEDGES)
    !DSP%FIELDS = 0D0


    ! mark which cells to include at the start
    TOTLEN = 0
    STARTCELLS = .FALSE.
    FRACCELLS = 0D0 ! fraction of cell on a starting edge

    ! which edges to actually start on
    WHICHEDGES = .FALSE.
    INCLUDERESV = .FALSE.
    IF (STARTEDGES(1).LT.-1) THEN
       WHICHEDGES = .TRUE.
       INCLUDERESV = .TRUE.
    ELSEIF (STARTEDGES(1).EQ.-1) THEN
       WHICHEDGES = .TRUE.
    ELSE
       WHICHEDGES(STARTEDGES(1:NSTARTEDGE)) = .TRUE.
    ENDIF

    TOTLEN = 0D0
    DO CC = 1,MESHP%NCELL
       IF (MESHP%CELLTYPE(CC).EQ.0) THEN
          ! nodal cell
          ! partially fill cell for every edge it is attached to
          DO CC2 = 1,MESHP%DEG(CC)
             BC = MESHP%BOUNDS(CC,CC2)
             IF (BC.EQ.0) CYCLE
             IF (MESHP%CELLTYPE(BC).EQ.1.AND.WHICHEDGES(MESHP%EDGEIND(BC,1))) THEN
                ! cell contains piece of a starting edge
                STARTCELLS(CC) = .TRUE.
                FRACCELLS(CC) = FRACCELLS(CC) + 1D0/MESHP%DEG(CC)
                TOTLEN = TOTLEN + MESHP%LEN(CC)/MESHP%DEG(CC)
             ENDIF
          ENDDO
       ELSEIF (DSP%MESHP%CELLTYPE(CC).EQ.1) THEN
          IF (WHICHEDGES(DSP%MESHP%EDGEIND(CC,1))) THEN
             ! Cell is entirely on a starting edge
             STARTCELLS(CC) = .TRUE.
             FRACCELLS(CC) = 1
             TOTLEN = TOTLEN + MESHP%LEN(CC)
          ENDIF
       ELSEIF (DSP%MESHP%CELLTYPE(CC).EQ.2.AND.INCLUDERESV) THEN
          ! fill reservoir as well
          STARTCELLS(CC) = .TRUE.
          FRACCELLS(CC) = 1
          TOTLEN = TOTLEN + MESHP%LEN(CC)
       ENDIF
    ENDDO

    DO FC = 1,DSP%NFIELD ! for each field
       IF (STARTCONC(FC).GE.0) THEN
          ! start with starting concentration
          DSP%FIELDS(:,FC) = STARTCONC(FC)*FRACCELLS 
       ELSE
          ! start with concentration integrating to 1
          DSP%FIELDS(:,FC) = FRACCELLS*MESHP%LEN/TOTLEN
       ENDIF
    ENDDO

    ! Set fixed cell values
    DO FC = 1,DSP%NFIELD
       DO CC = 1,DSP%MESHP%NCELL
          IF (DSP%ISFIXED(CC,FC)) DSP%FIELDS(CC,FC) = DSP%FIXVALS(CC,FC)
       ENDDO
    ENDDO
    
    CALL INTEGRATEFIELD(DSP,TMP,TOTLEN,.FALSE.)

  END SUBROUTINE INITIALIZEFIELDEDGES
  
  SUBROUTINE OUTPUTFIELDS(DSP,OUTFILE,INFO,APPEND,OUTPUT1FIELD)
    ! output field values
    ! INFO is a list of floats containing additional information that gets dumped before the fields
    ! append: whether or not to append to previous file
    ! Snapshotvel=1 means include boundary velocities in the snapshot as well
    ! optionally: OUTPUT1FIELD = integer setting single field to output.
    ! if not provided (or if negative), outputs all
    USE KEYS, ONLY : SNAPSHOTVEL
    
    
    IMPLICIT NONE

    TYPE(DYNSYSTEM), POINTER :: DSP
    CHARACTER(LEN=*), INTENT(IN) :: OUTFILE
    DOUBLE PRECISION, INTENT(IN) :: INFO(:)
    LOGICAL, INTENT(IN) :: APPEND
    INTEGER, INTENT(IN), OPTIONAL :: OUTPUT1FIELD
    INTEGER, PARAMETER :: OU = 99
    INTEGER :: FC, CELLTYPE(DSP%MESHP%NCELL), CC, NFIELD
    LOGICAL :: DOOUTPUT1
    
    IF (APPEND) THEN
       OPEN(UNIT=OU, FILE=OUTFILE, ACCESS='APPEND')
    ELSE
       OPEN(UNIT=OU, FILE=OUTFILE)
    ENDIF

    ! write negative index for reservoir cells
    ! cell type for other cell types
    
    ! CELLTYPE = DSP%MESHP%CELLTYPE
    ! DO CC = 1, DSP%MESHP%NCELL
    !    IF (CELLTYPE(CC).EQ.2) CELLTYPE(CC) = -DSP%MESHP%RESVIND(CC)
    ! ENDDO

    

    DOOUTPUT1 = .FALSE.; NFIELD = DSP%NFIELD
    
    IF (PRESENT(OUTPUT1FIELD)) THEN       
       IF (OUTPUT1FIELD.GE.0) THEN
          DOOUTPUT1 = .TRUE.
          NFIELD = 1
       ENDIF
    ENDIF
          
    ! write info line: number of fields, number of cells, max degree,
    ! followed by extra info
    WRITE(OU,*) NFIELD, DSP%MESHP%NCELL, DSP%MESHP%MAXDEG, INFO, SNAPSHOTVEL
    IF (SNAPSHOTVEL.GT.0) THEN
       ! write velocities across each boundary
       WRITE(OU,*) DSP%VEL
    END IF
    IF (DOOUTPUT1) THEN
       ! write values for 1 field
       WRITE(OU,*) DSP%FIELDS(:,OUTPUT1FIELD)
    ELSE
       ! for each field, write field values for all cells
       DO FC = 1,DSP%NFIELD
          WRITE(OU,*) DSP%FIELDS(:,FC)
       ENDDO
    ENDIF
    
    CLOSE(OU)
    
  END SUBROUTINE OUTPUTFIELDS
  
  SUBROUTINE SETPARAMDYNSYS(DSP,NETP)
    ! set parameters of dynamical system using global keyword arguments
    ! NNODE = number of network nodes, to allow for fixing random nodes
    USE KEYS, ONLY : KON, KOFF, KDEQUIL, DOBUFFER, DCOEFF, FASTEQUIL, &
         & FIXNODES,FIXVALS, NFIX,NFIXRESV,MOBILEFIELD, RANDFIXNODES, &
         & FIXEDGEVEL,FIXEDGEVELVAL,NFIXEDGEVEL,&
         & CEXT, NPERM, PERMEABILITY, TRACKFLUXPERM, PERMNODES, PERMFROMFILE, &
         & RANDPERMNODES, USEPERMPREFACTOR, FIXNEARNODEDIST, NACT, ACTRATE, ACTNODES,&
         & ACTNEARNODEDIST, DEPRATE, PERMNEARNODEDIST, FIXRECTANGLE, &
         & ALLOWFIXEDRESV, DORESERVOIRS, ALLOWRESVFIX, &
         & NFIXCELL, RANDFIXCELLS,FIXCELLS, RANDFIXPTS,FIXPTS,NFIXPT,FIXPTCENT,&
         & FIXPTRAD, MAXNABSORBER, MAXNFIELD,NPERMPOS,PERMPOS,POSPERMEABILITY,VARRAD, &
         & RANDFIXRESV, FIXPTMAXDIST, FIXPTEXCENT, FIXPTEXRAD, DOFLOW, USERESVELEMENTS, &
         & UNIFORMBUFFER, CONCENTRATIONS3D, PERIODGLOBALPERM, DURGLOBALPERM, &
         & USEGLOBALRESV,GLOBALRESVKR,GLOBALRESVKMR, GLOBALRESVKOUT, GLOBALRESVKMOUT, PERMTOGLOBALRESV, &
         & TRACKDCDT
    USE NETWORKUTIL, ONLY : NETWORK
    USE GENUTIL, ONLY : RANDSELECT_INT
    USE MT19937, ONLY : RANDUNIFCIRCLE
    IMPLICIT NONE
    TYPE(DYNSYSTEM), POINTER :: DSP
    TYPE(NETWORK), POINTER :: NETP
    INTEGER :: NFIELD
    TYPE(MESH), POINTER :: MESHP
    INTEGER :: NC, FC, CC, CT, CT2, CC2, EC, RC, I
    INTEGER :: NODELIST(netp%NNODE), TMP(netp%NNODE), ALLNODELIST(NETP%NNODE)
    INTEGER :: NODEAVAIL, CELLLIST(DSP%MESHP%NCELL)
    DOUBLE PRECISION :: DIST, MINX, MINY, MAXX, MAXY, POS(NETP%DIM)
    DOUBLE PRECISION :: WEIGHTS(DSP%MESHP%NCELL), COORDS(MAXNABSORBER,2)
    DOUBLE PRECISION :: DIFFS(DSP%MESHP%NCELL,2), DISTS(DSP%MESHP%NCELL)
    INTEGER :: TMPARR(1), FIXRESV(MAXNABSORBER,MAXNFIELD)
    LOGICAL :: SUCCESS   
    
    FIXRESV = 0
    
    IF (NACT.GT.0) THEN
       DSP%DOACTIVATION = .TRUE.
    ELSE
       DSP%DOACTIVATION = .FALSE.
    ENDIF
    
    IF (ALLOWFIXEDRESV) THEN
       ! allow reservoirs to be fixed except the ones explicitly excluded
       CT = 0
       DO NC = 1,NETP%NNODE
          IF (NETP%NODERESV(NC).GT.0) THEN ! this node is attached to a reservoir
             IF (ALLOWRESVFIX(NETP%NODERESV(NC))) THEN ! that reservoir can be fixed
                CT = CT+1
                NODELIST(CT) = NC
             ENDIF
          ENDIF
       ENDDO
       NODEAVAIL = CT
    ELSE ! fixed nodes cannot be attached to reservoirs
       ! list of nodes not attached to reservoirs
       CT = 0
       DO NC = 1,NETP%NNODE
          ALLNODELIST(NC) = NC; ! all nodes
          IF (NETP%NODERESV(NC).EQ.0) THEN
             CT = CT+1
             NODELIST(CT) = NC ! only those nodes not attached to reservoirs
          ENDIF
       ENDDO
       NODEAVAIL = CT
    ENDIF    
 
    IF (.NOT.DSP%ARRAYSET) THEN
       PRINT*, 'ERROR IN SETPARAMSDYNSYS: dynamic system not yet allocated'
       STOP 1
    ENDIF

    MESHP=>DSP%MESHP
    
    DSP%KON = KON
    DSP%KOFF = KOFF
    IF (KDEQUIL.LT.0) THEN
       DSP%KDEQUIL = DSP%KOFF/DSP%KON
    ELSE
       DSP%KDEQUIL = KDEQUIL
    ENDIF
    DSP%DCOEFF = DCOEFF(1:DSP%NFIELD)
    DSP%VEL = 0D0
    DSP%MOBILEFIELD = MOBILEFIELD(1:DSP%NFIELD)       

    DSP%UNIFORMBUFFER = UNIFORMBUFFER
    
    ! ------------------
    ! Check buffer type matches number of fields
    ! ------------------    
    IF (DOBUFFER) THEN
       IF (FASTEQUIL) THEN          
          DSP%BUFFERTYPE = 2
          NFIELD = 2
       ELSE
          DSP%BUFFERTYPE = 1
          NFIELD = 3         
       ENDIF
    ELSE
       IF (DSP%DOACTIVATION) THEN
          NFIELD = 2
       ELSE
          NFIELD = 1
       ENDIF
       DSP%BUFFERTYPE = 0
    ENDIF

    IF (NFIELD.NE.DSP%NFIELD) THEN
       PRINT*, 'ERROR IN SETPARAMDYNSYS: system has been allocated with wrong number of fields', &
            & NFIELD, DSP%NFIELD, dsp%doactivation
       STOP 1
    ENDIF

    ! set up fixed cell values
    DSP%FIXVALS = 0D0
    DSP%ISFIXED = .FALSE.

    ! track total concentrations and solve for free ligand on each step
    DSP%TRACKDCDT = TRACKDCDT
    
    IF (RANDFIXPTS) THEN
       ! pick randomly selected points in a disc-shaped region
       ! fix the mesh cell nearest each of them

       DO FC = 1,DSP%NFIELD
          ! select points in the circle
          COORDS(1:NFIXPT(FC),:) = RANDUNIFCIRCLE(NFIXPT(FC),FIXPTCENT,FIXPTRAD)
          
          DO CC = 1,NFIXPT(FC)
             SUCCESS = .FALSE.

             DO WHILE (.NOT.SUCCESS) ! keep trying until we get a good sample
                ! find the nearest mesh cell             
                DO I = 1,2
                   DIFFS(:,I) = COORDS(CC,I) -  MESHP%POS(:,I)
                ENDDO
                DISTS = SUM(DIFFS**2,2)
                TMPARR = MINLOC(DISTS)             

                ! CHECK: did we successfully sample?
                IF (MESHP%CELLTYPE(TMPARR(1)).EQ.2.AND..NOT.ALLOWFIXEDRESV) THEN
                   ! nearest point is a reservoir mesh cell, not allowed to be fixed
                   SUCCESS = .FALSE.
                ELSEIF (DISTS(TMPARR(1)).GT.FIXPTMAXDIST**2) THEN
                   ! nearest mesh cell is too far
                   SUCCESS = .FALSE.
                ELSE
                   SUCCESS = .TRUE.
                ENDIF

                ! Check if point is inside excluded region
                IF (SUM((MESHP%POS(TMPARR(1),:) - FIXPTEXCENT(1:NETP%DIM))**2)&
                     & .LT.FIXPTEXRAD**2) THEN
                   PRINT*, 'point is in excluded region. Reselecting', CC, &
                        & FIXPTEXRAD, TMPARR(1), MESHP%POS(TMPARR(1),:)
                   SUCCESS = .FALSE.
                ENDIF
                
                IF (.NOT.SUCCESS) THEN ! try sampling again
                   COORDS(CC:CC,:) = RANDUNIFCIRCLE(1,FIXPTCENT,FIXPTRAD)          
                ENDIF
             ENDDO

             FIXCELLS(CC,FC) = TMPARR(1)

             ! set it as fixed
             DSP%ISFIXED(FIXCELLS(CC,FC),FC) = .TRUE.
             DSP%FIXVALS(FIXCELLS(CC,FC),FC) = FIXVALS(CC,FC)

             PRINT*, 'Fixed cell in position: ', FC, CC, FIXCELLS(CC,FC), MESHP%POS(TMPARR(1),:)
          ENDDO
       END DO
    END IF
    
    IF (RANDFIXCELLS) THEN
       ! fix randomly selected mesh cells, weighted by SA
       CELLLIST = (/(CC,cc=1,MESHP%NCELL)/)
       
       WEIGHTS = MESHP%SA
       IF (.NOT.ALLOWFIXEDRESV) THEN
          ! do not allow reservoirs to be fixed
          DO CC = 1,MESHP%NCELL
             IF (MESHP%CELLTYPE(CC).EQ.2) WEIGHTS(CC) = 0D0
          ENDDO
       ENDIF       
       ! select random cells to fix, weighted by cell SA
       DO FC = 1,DSP%NFIELD
          CALL RANDSELECT_INT(CELLLIST,NFIXCELL(FC),.FALSE.,FIXCELLS(1:NFIXCELL(FC),FC),TMP,WEIGHTS)

          DO CC = 1,NFIXCELL(FC)
             DSP%ISFIXED(FIXCELLS(CC,FC),FC) = .TRUE.
             DSP%FIXVALS(FIXCELLS(CC,FC),FC) = FIXVALS(CC,FC)
          ENDDO
       ENDDO
       
    ENDIF
    
    IF (RANDFIXNODES) THEN
       ! Pick network nodes to fix without replacement
       DO FC = 1,DSP%NFIELD
          CALL RANDSELECT_INT(NODELIST(1:NODEAVAIL),NFIX(FC),.FALSE.,FIXNODES(1:NFIX(FC),FC),TMP)
          PRINT*, 'Field ', FC, NFIX(FC), ' fixed nodes:', FIXNODES(1:NFIX(FC),FC)
       ENDDO
    ENDIF
    
    
    IF (RANDFIXRESV) THEN
       ! randomly select (without replacement) reservoirs to be fixed (excluding those that are not allowed to be fixed
       DO FC = 1,DSP%NFIELD
          CALL RANDSELECT_INT( PACK((/(RC, RC=1,NETP%NRESV)/),ALLOWRESVFIX(1:NETP%NRESV)), &
               & NFIXRESV(FC),.FALSE.,FIXRESV(1:NFIXRESV(FC),FC),TMP)
          PRINT*, 'Field ', FC, NFIXRESV(FC), ' fixed reservoirs:', FIXRESV(1:NFIX(FC),FC)

          ! update the total number of fixed mesh cells to include the fixed reservoirs
          NFIX(FC) = NFIX(FC) + NFIXRESV(FC)
          DO CT = 1,NFIXRESV(FC)
             DSP%FIXVALS(FIXRESV(CT,FC),FC) = FIXVALS(CT,FC)
          ENDDO
       ENDDO       
       
    END IF
    
    ! fix the cells corresponding to the fixed network nodes
    DO CC = 1,MESHP%NCELL
       IF (MESHP%CELLTYPE(CC).EQ.0) THEN
          ! this is a nodal cell
          NC = MESHP%NODEIND(CC)          
          DO FC = 1,DSP%NFIELD
             DO CT = 1,NFIX(FC)
                IF (FIXNODES(CT,FC).EQ.NC) THEN
                   ! cell is fixed
                   DSP%ISFIXED(CC,FC) = .TRUE.
                   DSP%FIXVALS(CC,FC) = FIXVALS(CT,FC)                   
                ENDIF
             ENDDO
          ENDDO
       ELSEIF (MESHP%CELLTYPE(CC).EQ.1) THEN
          ! this is an edge cell
          ! if it hits a fixed terminal node, fix its value instead
          ! value fixed to that of LAST fixed copy of this terminal node
          DO FC = 1,DSP%NFIELD
             DO CT= 1,MESHP%DEG(CC)
                IF (MESHP%BOUNDS(CC,CT).EQ.0) THEN
                   NC = MESHP%TERMNODE(CC,CT)
                   DO CT2 = 1,NFIX(FC)
                      IF (FIXNODES(CT2,FC).EQ.NC) THEN                         
                         DSP%ISFIXED(CC,FC)=.TRUE.
                         DSP%FIXVALS(CC,FC) = FIXVALS(CT2,FC)
                         EXIT
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ELSEIF (MESHP%CELLTYPE(CC).EQ.2) THEN
         
          
          ! reservoir cell
          Rc = MESHP%RESVIND(CC)          
          DO FC = 1,DSP%NFIELD

             ! check if reservoir is within range of those that have nodes attached to them
             IF (RC.LE.NETP%NRESV) THEN
                ! fix if any attached node is in the fix list
                DO CT= 1,NETP%RESVNNODE(RC)
                   NC = NETP%RESVNODES(RC,CT)

                   DO CT2 = 1,NFIX(FC)
                      IF (FIXNODES(CT2,FC).EQ.NC) THEN
                         IF (.NOT.ALLOWFIXEDRESV.and.DORESERVOIRS) THEN
                            PRINT*, 'ERROR IN SETTING UP FIXED NODES: reservoir fixed nodes not allowed'
                            STOP 1
                         ENDIF

                         ! set to first fixed value among attached nodes
                         DSP%ISFIXED(CC,FC)=.TRUE.
                         DSP%FIXVALS(CC,FC) = FIXVALS(CT2,FC)
                         EXIT
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

             ! also fix if the reservoir itself is in the fix list
             IF (ALLOWFIXEDRESV) THEN
                DO CT = 1,NFIXRESV(FC)
                   IF (FIXRESV(CT,FC).EQ.RC) THEN                     
                      DSP%ISFIXED(CC,FC)=.TRUE.
                      DSP%FIXVALS(CC,FC) = FIXVALS(CT,FC)
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ENDDO

    IF (FIXNEARNODEDIST.GT.0) THEN
       ! fix all cells whose centers are near enough to fixed nodes
       DO FC = 1,DSP%NFIELD
          DO CT = 1,NFIX(FC)          
             NC = FIXNODES(CT,FC) ! which node
             ! distances to the node
             DO CC = 1,MESHP%NCELL
                DIST = SQRT(SUM((MESHP%POS(CC,:)-NETP%NODEPOS(NC,:))**2))
                IF (DIST.LT.FIXNEARNODEDIST) THEN
                   IF (.NOT.ALLOWFIXEDRESV.AND.DORESERVOIRS.AND.MESHP%CELLTYPE(CC).EQ.2) THEN
                     PRINT*, 'ERROR IN SETTING UP FIXED NODES: reservoir fixed nodes not allowed'
                     STOP 1
                   ENDIF
                   
                   DSP%ISFIXED(CC,FC) = .TRUE.
                   DSP%FIXVALS(CC,FC) = FIXVALS(CT,FC)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    
    IF (ANY(FIXRECTANGLE(:,4).GT.0)) THEN
       ! fix all cells whose center lies within a certain rectangle
       ! rectangle format: minx, miny, width, height

       DO FC = 1,DSP%NFIELD
          IF (FIXRECTANGLE(FC,4).GT.0) THEN ! fix this field in the rectangle
             MINX = FIXRECTANGLE(FC,1)
             MINY= FIXRECTANGLE(FC,2)
             MAXX = MINX+FIXRECTANGLE(FC,3)
             MAXY = MINY + FIXRECTANGLE(FC,4)
             DO CC = 1,MESHP%NCELL
                POS = MESHP%POS(CC,:)
                IF (POS(1).GT.MINX.AND.POS(1).LT.MAXX &
                     & .AND.POS(2).GT.MINY.AND.POS(2).LT.MAXY) THEN
                   DSP%ISFIXED(CC,FC) = .TRUE.
                   DSP%FIXVALS(CC,FC) = FIXRECTANGLE(FC,5)
                ENDIF
             END DO
          ENDIF
       ENDDO
       
    END IF
    
    ! set up local activation in a region around some center
    DSP%DEPRATE = DEPRATE ! depletion everywhere
    
   
    DO CT = 1,NACT
       NC = ACTNODES(CT) ! which node
       ! Distances to the node
       DO CC = 1,MESHP%NCELL
          DIST = SQRT(SUM((MESHP%POS(CC,:)-NETP%NODEPOS(NC,:))**2))
          IF (DIST.LT.ACTNEARNODEDIST(ct)) THEN             
             DSP%ACTRATE(CC) = ACTRATE(CT)
          ENDIF
       ENDDO
    ENDDO
                
    ! Set up permeable cells (on permeable network nodes)
    DSP%NPERM = COUNT(NETP%ISPERM)
    DSP%CEXT = CEXT(1:DSP%NFIELD)

    IF (PERMFROMFILE) THEN       
       ! get permeability nodes from the network file (data stored in network object)
       NPERM = COUNT(NETP%ISPERM)
       PERMNODES(1:NPERM) = PACK(ALLNODELIST,NETP%ISPERM)
       DO FC = 1,DSP%NFIELD
          PERMEABILITY(1:NPERM,FC) = PERMEABILITY(1,FC)
       ENDDO
    ELSEIF (RANDPERMNODES) THEN
       ! Pick random network nodes to be permeable       
        CALL RANDSELECT_INT(NODELIST(1:NODEAVAIL),NPERM,.FALSE.,PERMNODES(1:NPERM),TMP)
          PRINT*, NPERM, ' permeable nodes:', PERMNODES(1:NPERM)
    ENDIF
        
    
    DO CT = 1,NPERM
       NC = PERMNODES(CT)
       ! Update network node permeability
       NETP%ISPERM(NC) = .TRUE.
       
       CC = NETP%NODECELLS(NC)

       IF (CC.EQ.0) THEN
          IF (NETP%NODEDEG(NC).NE.1) THEN
             PRINT*, 'ERROR IN SETPARAMDYNSYS. Is this a terminal node?'
             STOP 1
          ENDIF
          
          ! terminal node. Find a cell that borders this node
          EC = NETP%NODEEDGE(NC,1)
          IF (NETP%EDGENODE(EC,1).EQ.NC) THEN
             CC2 = NETP%EDGECELLS(EC,1) ! first cell on this edge
          ELSEIF( NETP%EDGENODE(EC,2).EQ.NC) THEN
             CC2 = NETP%EDGECELLS(EC,NETP%EDGENCELL(EC)) ! last cell on this edge
          ELSE
             PRINT*, 'ERROR IN SETPARAMDYNSYS. Bad edge at terminal node', NC, EC, NETP%EDGENODE(EC,:)
             STOP 1
          ENDIF

          ! make the cell touching this terminal node permeable          
          DSP%ISPERM(CC2) = .TRUE.
          DSP%PERM(CC2,:) = PERMEABILITY(CT,1:DSP%NFIELD)                       
       ELSE
          ! make the nodal cell itself permeable
          DSP%ISPERM(CC) = .TRUE.
          DSP%PERM(CC,:) = PERMEABILITY(CT,1:DSP%NFIELD)
       ENDIF
    ENDDO

    ! Only track flux from permeable nodes if there are some
    TRACKFLUXPERM = TRACKFLUXPERM.AND.(NPERM.GT.0.OR.NPERMPOS.GT.0)   
    
     IF (PERMNEARNODEDIST.GT.0.AND.NPERM.GT.0) THEN
        ! make permeable all cells whose centers are near enough to fixed nodes 
        DO CT = 1,NPERM         
           NC = PERMNODES(CT) ! which node
           ! distances to the node
           DO CC = 1,MESHP%NCELL
              DIST = SQRT(SUM((MESHP%POS(CC,:)-NETP%NODEPOS(NC,:))**2))
              IF (DIST.LT.PERMNEARNODEDIST) THEN
                 DSP%ISPERM(CC) = .TRUE.
                 DSP%PERM(CC,:) = PERMEABILITY(CT,1:DSP%NFIELD)                 
              ENDIF
           ENDDO
        ENDDO
     ENDIF

    ! Make permeable all cells within a given distance of a particular position
    IF (NPERMPOS.GT.0) THEN
       DO FC = 1,DSP%NFIELD
          DO NC = 1,NPERMPOS ! go over all permeable centers
             ! distances to the position
             DO CC = 1,MESHP%NCELL
                DIST = SQRT(SUM((MESHP%POS(CC,:)-PERMPOS(NC, 1:NETP%DIM))**2))
                IF (DIST.LT.PERMNEARNODEDIST) THEN
                   DSP%ISPERM(CC) = .TRUE.
                   DSP%PERM(CC,:) = POSPERMEABILITY(CT,1:DSP%NFIELD)
                ENDIF
             ENDDO
          ENDDO
       ENDDO    
    END IF

    IF (PERIODGLOBALPERM.GT.0) THEN
       ! set everything to permeable if using periodic global permeability
       DO FC = 1,DSP%NFIELD
          DSP%PERM(:,FC) = PERMEABILITY(1,FC)
       ENDDO
       DSP%PERIODGLOBALPERM = PERIODGLOBALPERM
       DSP%DURGLOBALPERM = DURGLOBALPERM
    ELSE
       DSP%PERIODGLOBALPERM = -1D0
       DSP%DURGLOBALPERM = 0D0
    ENDIF
    
    PRINT*, 'TRACKFLUXPERM:', TRACKFLUXPERM
    
    IF (USEPERMPREFACTOR) THEN
       ! the provided permeability is actually a prefactor that should
       ! be multiplied by mesh cell length (actually [surface area/(2*pi*a)])
       ! and by species diffusivity
       DO FC = 1, DSP%NFIELD          
          DSP%PERM(:,FC) = DSP%PERM(:,FC)*MESHP%SA*DCOEFF(FC)
       ENDDO
    END IF

    print*, 'Permeable cells:'
    DO CC = 1,MESHP%NCELL
       IF (DSP%ISPERM(CC)) THEN
          SELECT CASE (MESHP%CELLTYPE(CC))
          CASE(0)
             PRINT*, 'Node ', MESHP%NODEIND(CC), ', cell ', CC
          CASE(1)
             PRINT*, 'Edge ', MESHP%EDGEIND(CC,1), ', cell ', CC
          CASE(2)
             PRINT*, 'Reservoir ', MESHP%RESVIND(CC), ', cell ', CC
          CASE DEFAULT
             PRINT*, 'this mesh cell is an unknown case:', CC, MESHP%CELLTYPE(CC)
          END SELECT
       ENDIF
    ENDDO

    print*, 'Fixed cells for field 1:', NFIX(1)
    DO CC = 1,MESHP%NCELL
       IF (DSP%ISFIXED(CC,1)) THEN
          SELECT CASE (MESHP%CELLTYPE(CC))
          CASE(0)
             PRINT*, 'Node ', MESHP%NODEIND(CC), ', cell ', CC
          CASE(1)
             PRINT*, 'Edge ', MESHP%EDGEIND(CC,1), ', cell ', CC
          CASE(2)
             PRINT*, 'Reservoir ', MESHP%RESVIND(CC), ', cell ', CC
          CASE DEFAULT
             PRINT*, 'this mesh cell is an unknown case:', CC, MESHP%CELLTYPE(CC)
          END SELECT
       ENDIF
    ENDDO

    DSP%VARRAD = VARRAD
    DSP%USEEDGEFLOW = DOFLOW

    ! Treat concentrations as 3D if explicitly flagged, OR if working
    ! with meshed reservoirs
    DSP%CONC3D = USERESVELEMENTS.OR.CONCENTRATIONS3D

    ! set up global reservoir kinetics
    IF (USEGLOBALRESV) THEN
       DSP%GLOBALRESVKR = GLOBALRESVKR
       DSP%GLOBALRESVKMR = GLOBALRESVKMR
       DSP%GLOBALRESVKOUT = GLOBALRESVKOUT
       DSP%GLOBALRESVKMOUT = GLOBALRESVKMOUT
       DSP%PERMTOGLOBALRESV = PERMTOGLOBALRESV
    ENDIF
  END SUBROUTINE SETPARAMDYNSYS
  
  SUBROUTINE SETUPDYNSYS(DSP,MESHP,NFIELD)
    ! allocate all the arrays for a dynamical system
    ! assuming mesh object has been set up already
    USE KEYS, ONLY : FIXNODES
    IMPLICIT NONE
    TYPE(DYNSYSTEM), POINTER :: DSP
    TYPE(MESH), POINTER :: MESHP
    INTEGER, INTENT(IN) :: NFIELD

    
    IF (.NOT.MESHP%ARRAYSET.OR..NOT.MESHP%CELLSET) THEN
       PRINT*, 'ERROR IN SETUPDYNSYS: mesh object not fully set up'
       STOP 1
    ENDIF
   
    DSP%NFIELD = NFIELD
   
    ALLOCATE(DSP%FIELDS(MESHP%NCELL,NFIELD),DSP%VEL(MESHP%NCELL,MESHP%MAXDEG))
    ALLOCATE(DSP%DCOEFF(NFIELD))   
    ALLOCATE(DSP%ISFIXED(MESHP%NCELL,NFIELD),DSP%FIXVALS(MESHP%NCELL,NFIELD))
    ALLOCATE(DSP%MOBILEFIELD(NFIELD))
    ALLOCATE(DSP%ACTRATE(MESHP%NCELL), DSP%DEPRATE(MESHP%NCELL))
    ALLOCATE(DSP%ISPERM(MESHP%NCELL),DSP%PERM(MESHP%NCELL,NFIELD),DSP%CEXT(NFIELD))
    DSP%PERM = 0D0
    DSP%FIELDS = 0D0
    DSP%ISPERM = .FALSE.
    
    DSP%MESHP=>MESHP
    
    DSP%ARRAYSET = .TRUE.    

  END SUBROUTINE SETUPDYNSYS

  SUBROUTINE CLEANUPDYNSYS(DSP)
    ! clean up dynamical system object
    IMPLICIT NONE
    TYPE(DYNSYSTEM), POINTER :: DSP

    DEALLOCATE(DSP%FIELDS, DSP%VEL, DSP%DCOEFF, DSP%ISFIXED, &
         & DSP%FIXVALS,DSP%MOBILEFIELD, DSP%ACTRATE, DSP%DEPRATE)
  END SUBROUTINE CLEANUPDYNSYS
  
END MODULE DYNSYSUTIL
