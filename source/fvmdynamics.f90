MODULE FVMDYNAMICS
  ! Run dynamic evolution using finite volume method
  USE NETWORKUTIL, ONLY : NETWORK
  USE MESHUTIL, ONLY : MESH
  USE DYNSYSUTIL, ONLY : DYNSYSTEM
  USE KEYS, ONLY : VERBOSE
  USE GENUTIL, ONLY : PI
  
  IMPLICIT NONE
  ! global variable: bad coder, no biscuit. Just for use in debugging
  DOUBLE PRECISION :: GLOBALTEST = 0D0
  
CONTAINS  
  SUBROUTINE RUNDYNAMICS(NETP,DSP,NSTEP,DELT,VELCONTROL,STARTRATE,STOPRATE, &
       & RUNSPEED,DOCLOSEBOUNDRATE, CLOSEBOUNDRATEPERLEN, OPENBOUNDRATE,&
       & DOCLOSEBOUNDONCE, PCLOSEBOUNDPERLEN)
    ! dynamics with random velocities along network edges
    ! VELCONTROL: string setting how velocities on edges are controlled
    
    USE KEYS, ONLY : OUTFILE, SNAPSHOTFILE, OUTPUTEVERY, PRINTEVERY, &
         & SNAPSHOTEVERY, TRACKFIXNODEFLUX, VERBOSE, DOFLOW, CONTFILE, &
         & OUTPUTEVERYSWITCH, OUTPUTEVERYSTART, TRACKFLUXPERM, OUTPUTTOTFLUXONLY, &
         & LOGSNAPSHOT, NSNAPSHOT,OUTPUT1FIELD      
    !USE mt19937, ONLY : GRND
    USE NETWORKUTIL, ONLY : NETWORK, OUTPUTNETWORK
    USE MESHUTIL, ONLY : MESH, OUTPUTMESH, CLOSEBOUNDS
    USE DYNSYSUTIL, ONLY : DYNSYSTEM, OUTPUTFIELDS, INTEGRATEFIELD
    USE NODECONTRACTUTIL, ONLY : OUTPUT_NODECONTRACT, GETVELS_NODECONTRACT
    
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    TYPE(DYNSYSTEM), POINTER :: DSP
    INTEGER, INTENT(IN) :: NSTEP
    DOUBLE PRECISION, INTENT(IN) :: DELT, STARTRATE, STOPRATE, RUNSPEED
    CHARACTER(LEN=*), INTENT(IN) :: VELCONTROL
    LOGICAL, INTENT(IN) :: DOCLOSEBOUNDRATE, DOCLOSEBOUNDONCE
    DOUBLE PRECISION, INTENT(IN) :: CLOSEBOUNDRATEPERLEN, OPENBOUNDRATE, PCLOSEBOUNDPERLEN
    
    TYPE(MESH), POINTER :: MESHP
    DOUBLE PRECISION :: INFO(1)
    INTEGER :: STEP
    LOGICAL :: EDGECLOSED(NETP%NEDGE)
    DOUBLE PRECISION :: FLUX(DSP%MESHP%NCELL,DSP%NFIELD)
    DOUBLE PRECISION :: DFDT(DSP%MESHP%NCELL,DSP%NFIELD)
    INTEGER :: NFIELD, OU
    DOUBLE PRECISION :: RUNVELS(NETP%NEDGE), TIMEHASRUN(NETP%NEDGE)
    LOGICAL :: ISRUN(NETP%NEDGE)
    DOUBLE PRECISION :: CURTIME, INTFIELD(DSP%NFIELD), TOTLEN
    INTEGER :: FIXCELLS(DSP%MESHP%NCELL,DSP%NFIELD), NFIX(DSP%NFIELD)
    INTEGER :: ABSORBERNODES(NETP%NNODE,DSP%NFIELD), NABS(DSP%NFIELD)
    INTEGER :: FC, NC, CC, EC, DEG, SC
    DOUBLE PRECISION :: PSTART, U, PSTOP, U2, DV
    INTEGER :: CELLLIST(DSP%MESHP%NCELL)
    INTEGER :: SNAPSTEPLIST(NSNAPSHOT), NEXTSNAPIND
    LOGICAL :: TAKESNAP, PERMON
    DOUBLE PRECISION :: LOGMAX,DLOG
    DOUBLE PRECISION :: CFIELD(DSP%MESHP%NCELL), BFIELD(DSP%MESHP%NCELL)

    IF (.NOT.NETP%ARRAYSET) THEN
       PRINT*, 'ERROR IN RUNDYNAMICS: network not set up'
       STOP 1
    ELSEIF (.NOT.DSP%ARRAYSET) THEN
       PRINT*, 'ERROR IN RUNDYNAMICS: dynamic system not set up'
       STOP 1
    ENDIF        
    
    
    NFIELD = DSP%NFIELD    
    MESHP=>DSP%MESHP
    ! full list of cell indices
    CELLLIST = (/(CC, CC = 1,MESHP%NCELL)/)

    EDGECLOSED = .FALSE.    
    IF (DOCLOSEBOUNDONCE) THEN
       ! randomly close some edge mesh boundaries once at start of sim
       CALL CLOSEBOUNDS(MESHP,PCLOSEBOUNDPERLEN,0D0,EDGECLOSED)
       PRINT*, 'Fraction of edges initially blocked:', DBLE(COUNT(EDGECLOSED))/NETP%NEDGE
    ELSE
       EDGECLOSED = .FALSE.
    ENDIF

    
    ! Set up list of fixed cells
    NFIX = 0
    DO CC = 1,MESHP%NCELL
       DO FC = 1,NFIELD
          IF (DSP%ISFIXED(CC,FC)) THEN
             NFIX(FC) = NFIX(FC)+1
             FIXCELLS(NFIX(FC),FC) = CC
          ENDIF
       END DO
    ENDDO

    RUNVELS = 0D0 ! speed on each edge
    TIMEHASRUN = 0D0 ! times a run has been going on
    ISRUN = .FALSE. ! which edges are currently running

    CURTIME = 0D0
    DFDT = 0D0

    
    IF (SNAPSHOTEVERY.GT.0.OR.(LOGSNAPSHOT.AND.NSNAPSHOT.GT.0)) THEN
       ! dump original field values
       INFO = CURTIME
       CALL OUTPUTFIELDS(DSP,SNAPSHOTFILE,INFO,.FALSE.,OUTPUT1FIELD)

       IF (VELCONTROL.EQ.'NODECONTRACTIONS') THEN
          ! output node contraction status
          CALL OUTPUT_NODECONTRACT(NETP,CONTFILE,INFO,.FALSE.)
       ENDIF
    ENDIF


    IF (LOGSNAPSHOT) THEN
       ! decide when to take snapshots, logarithmically spaced
       LOGMAX = LOG(DBLE(NSTEP))
       DLOG = LOGMAX/(NSNAPSHOT-1) ! logarithmic separation btwn snapshots
       DO SC = 1,NSNAPSHOT
          SNAPSTEPLIST(SC) = NINT(EXP(DLOG*(SC-1)))
          IF (SC.EQ.1) THEN
             IF (SNAPSTEPLIST(1).EQ.0) THEN
                SNAPSTEPLIST(1) = 1
             ENDIF
          ELSEIF(SNAPSTEPLIST(SC).LE.SNAPSTEPLIST(SC-1)) THEN ! make separations between snapshots at least 1 step
             SNAPSTEPLIST(SC) = SNAPSTEPLIST(SC-1)+1
          ENDIF
       ENDDO

       ! count which is the next snapshot to take
       NEXTSNAPIND = 1
    END IF
    IF (VERBOSE) PRINT*, 'SNAPSHOTS WILL BE TAKEN ON THESE STEPS:', SNAPSTEPLIST
    
    ! output file for tracking total flux out of fixed or absorber nodes
    OU = 55
    OPEN(FILE=OUTFILE,UNIT=OU,STATUS='UNKNOWN')


    CALL INTEGRATEFIELD(DSP,INTFIELD,TOTLEN,.true.)
    PRINT*, 'TESTX1:',  CURTIME, TOTLEN
    PRINT*, 'TESTX1B:', INTFIELD/TOTLEN, NFIELD
    PRINT*, 'TESTX2:', DBLE(COUNT(EDGECLOSED))
    PRINT*, 'TESTX3:', NETP%NEDGE
    PRINT*, 'STEP, TIME, AVG FIELDS, NFIELD, EDGECLOSED: ', &
               & 0, CURTIME, INTFIELD/TOTLEN, NFIELD,&
               & DBLE(COUNT(EDGECLOSED))/NETP%NEDGE

   
    ! DO CC = 1,MESHP%NCELL
    !    PRINT*, CC, MESHP%POS(CC,2),DSP%ISFIXED(CC,1), DSP%FIXVALS(CC,1)
    ! ENDDO

    ! probability of starting a processive run over an edge on any given step
    Pstart = 1D0-EXP(-DELT*STARTRATE)
    PSTOP = 1D0-EXP(-DELT*STOPRATE)

    ! how much should node volume change on each step (for node contractions)
    ! Note volumes are in terms of tube length (V/(pi a^2))
    ! will only give the correct flow speed for degree 1 nodes
    DV = RUNSPEED*DELT

    ISRUN = .FALSE.
    RUNVELS = 0D0
    DSP%VEL = 0D0   

    CALL OUTPUTNETWORK(NETP,'outputnetwork.net',.FALSE.)
    
   
    DO STEP = 1,NSTEP
       ! update run velocities on edges
       IF (DOFLOW) THEN
          SELECT CASE(VELCONTROL)
          CASE('RANDVEL')
             ! randomly generated edge velocities, Poisson on/off
             CALL GETVELS_RANDVEL(DSP,NETP,PSTART,PSTOP,ISRUN,RUNVELS,(STEP.EQ.1),RUNSPEED)
          CASE('RANDVELSWITCH')
             ! randomly generated edge velocities, instant switching
             CALL GETVELS_RANDVELSWITCH(DSP,NETP,PSTOP,RUNVELS,(STEP.EQ.1),RUNSPEED)
          CASE('NODECONTRACTIONS')
             ! node contraction and expansion
             CALL GETVELS_NODECONTRACT(DSP,NETP,PSTART,DV,RUNSPEED)
          CASE('EDGECONTRACTIONS')
             ! contraction and expansion at random mesh points on edges
             PRINT*, 'EDGECONTRACTIONS velocity control not yet set up'
             STOP 1
          CASE('FIXEDVELS')
             RUNVELS = RUNSPEED
             CALL GETVELS_FIXED(DSP,NETP,RUNVELS)

             !print*, DSP%VEL(:,1)
          CASE DEFAULT
             PRINT*, 'ERROR IN RUNDYNAMICS: unknown velocity control: ', VELCONTROL
             STOP 1          
          END SELECT       
       ENDIF       

       ! Update closed and open cell boundaries (blockages along edges)
       IF (DOCLOSEBOUNDRATE) THEN
          CALL CLOSEBOUNDS(MESHP,CLOSEBOUNDRATEPERLEN,OPENBOUNDRATE,EDGECLOSED,DELT)
       ENDIF
       
       ! ---------------
       ! Propagate field forward in time
       ! ---------------

       IF (DSP%BUFFERTYPE.EQ.2) THEN! rapid equilibration
          !          VERBOSE = STEP.GT.80D3
          CALL EULERSTEPEQUIL(DSP,DELT,DFDT,FLUX)          
       ELSE
          CALL EULERSTEP(DSP,DELT,DFDT,FLUX)
       ENDIF

       IF ((STEP.LE.OUTPUTEVERYSWITCH.AND.MOD(STEP,OUTPUTEVERYSTART).EQ.0).OR.&
            & MOD(STEP,OUTPUTEVERY).EQ.0) THEN
          ! output flux to file

          IF (TRACKFLUXPERM) THEN
             ! track flux out of permeable cells
             WRITE(OU,*) CURTIME, NFIELD, DELT, 1, COUNT(DSP%ISPERM)

             IF (OUTPUTTOTFLUXONLY) THEN
                DO FC = 1,DSP%NFIELD
                   WRITE(OU,*) FC, SUM(pack(FLUX(:,1),DSP%ISPERM))
                ENDDO
             ELSE
                WRITE(OU,*) PACK(CELLLIST,DSP%ISPERM)
                DO FC = 1,DSP%NFIELD
                   WRITE(OU,*) FC, pack(FLUX(:,FC),DSP%ISPERM), PACK(DSP%PERM(:,FC),DSP%ISPERM)
                ENDDO
             ENDIF
          ELSE
             WRITE(OU,*) CURTIME, NFIELD, DELT, 0, NFIX(1:NFIELD)

             ! output flux out of fixed cells          
             DO FC = 1,NFIELD
                IF (OUTPUTTOTFLUXONLY) THEN
                   ! flux out of all nodes together
                   WRITE(OU,*) FC, SUM(FLUX(FIXCELLS(1:NFIX(FC),FC),FC))
                ELSE
                   IF (NFIX(FC).GT.0) THEN
                      WRITE(OU,*) FIXCELLS(1:NFIX(FC),FC)
                   ELSE
                      WRITE(OU,*) 0
                   ENDIF
                   ! flux out of those nodes
                   WRITE(OU,*) FC, FLUX(FIXCELLS(1:NFIX(FC),FC),FC)
                ENDIF
             ENDDO
          ENDIF
          CALL FLUSH(OU)
       ENDIF

       ! propagate field
        IF (DSP%TRACKDCDT) THEN       
        ! propagate total ligand directly
           CFIELD = DSP%FIELDS(:,1)*(1 + DSP%FIELDS(:,2)/(DSP%FIELDS(:,1)+DSP%KDEQUIL))
          ! PRINT*, 'TESTX1:', SUM(DFDT(:,1)*MESHP%VOL), SUM(CFIELD*MESHP%VOL)
           CFIELD = CFIELD+DFDT(:,1)*DELT
           DSP%FIELDS(:,2) = DSP%FIELDS(:,2)+DFDT(:,2)*DELT
           
           BFIELD = DSP%KDEQUIL+DSP%FIELDS(:,2) - CFIELD           
           DSP%FIELDS(:,1) = (-BFIELD + SQRT(BFIELD**2 + 4*CFIELD*DSP%KDEQUIL))/2
        ELSE
           ! propagate free ligand, assumes dfdt is already
           ! the converted derivative for the free ligand
           DSP%FIELDS = DSP%FIELDS + DFDT*DELT
        ENDIF
        
       CURTIME = CURTIME+DELT

       ! if doing periodic global permeability, turn permeability off and on
       IF (DSP%PERIODGLOBALPERM.GT.0) THEN
          PERMON = (MOD(CURTIME,DSP%PERIODGLOBALPERM).LT.DSP%DURGLOBALPERM)
          IF (.NOT.ANY(DSP%ISPERM).AND.PERMON) THEN
             PRINT*, 'turn on permeability', curtime
          ELSEIF (ANY(DSP%ISPERM).AND..NOT.PERMON) THEN
             PRINT*, 'turn off permeability', curtime
          ENDIF                    
          DSP%ISPERM = PERMON
       ENDIF
              
       IF (MINVAL(DSP%FIELDS(:,1)) < -1D0) THEN
          PRINT*, 'ERROR: NEGATIVE CONCENTRATION'
          DO CC = 1,MESHP%NCELL
             PRINT*, CC, MESHP%CELLTYPE(CC), DSP%FIELDS(CC,1)
          ENDDO
          STOP 1
       ENDIF

       ! should we save a snapshot?
       IF (LOGSNAPSHOT) THEN
          TAKESNAP= (NEXTSNAPIND.LE.NSNAPSHOT.AND.STEP.EQ.SNAPSTEPLIST(NEXTSNAPIND))
          IF (TAKESNAP) THEN
             NEXTSNAPIND = NEXTSNAPIND+1
          ENDIF
       ELSE
          TAKESNAP = MOD(STEP,SNAPSHOTEVERY).EQ.0
       ENDIF
       
       IF (TAKESNAP) THEN
          ! save snapshot of the fields to a file
          INFO = CURTIME
          !PRINT*, 'TESTX1:', STEP, CURTIME, NEXTSNAPIND, snapsteplist(nextsnapind)
          CALL OUTPUTFIELDS(DSP,SNAPSHOTFILE,INFO,.TRUE.,OUTPUT1FIELD)

          IF (VELCONTROL.EQ.'NODECONTRACTIONS') THEN
             ! output node contraction status
             CALL OUTPUT_NODECONTRACT(NETP,CONTFILE,INFO,.TRUE.)
          ENDIF
       END IF

  !     print*, 'TESTX1', STEP, DSP%FIELDS(1,2), SUM(DSP%FIELDS(:,2)*DSP%MESHP%VOL)
       
       IF (MOD(STEP,PRINTEVERY).EQ.0) THEN

          IF (DSP%BUFFERTYPE.EQ.2) THEN! equilibrated buffers
             CFIELD = DSP%FIELDS(:,1)*(1+DSP%FIELDS(:,2)/(DSP%FIELDS(:,1)+DSP%KDEQUIL))
          ELSE
             CFIELD = DSP%FIELDS(:,DSP%NFIELD)
          ENDIF
          CALL INTEGRATEFIELD(DSP,INTFIELD,TOTLEN,.true.)          
          
          PRINT*, 'STEP, TIME, AVG FIELDS, NFIELD, EDGECLOSED, TOTAVG: ', &
               & STEP, CURTIME, INTFIELD/TOTLEN, NFIELD,&
               & DBLE(COUNT(EDGECLOSED))/NETP%NEDGE, SUM(CFIELD*MESHP%VOL)
       ENDIF
    ENDDO

  END SUBROUTINE RUNDYNAMICS

  !SUBROUTINE GETVELS_NODECONTRACT(DSP,NETP,PSTART,DV,FLOWSPEED)
  !END SUBROUTINE GETVELS_NODECONTRACT

   SUBROUTINE  GETVELS_RANDVEL(DSP,NETP,PSTART,PSTOP,ISRUN,RUNVELS,INITIALIZE,runspeed)
     ! update edge velocities:
     ! turn off with probability PSTOP, restart with probability PSTART
    ! RUNSPEED: flow speed
    ! RUNVELS: current velocity on each edge
    ! updates both RUNVELS and DSP%vel (flow across each cell boundary)
    ! INITIALIZE: if set, then need to initialize runvels for the first time
    USE KEYS, ONLY : VERBOSE
    USE mt19937, ONLY : GRND
    
    IMPLICIT NONE
    TYPE(DYNSYSTEM), POINTER :: DSP
    TYPE(NETWORK), POINTER :: NETP
    DOUBLE PRECISION, INTENT(IN) :: RUNSPEED, PStart,PSTOP
    LOGICAL, INTENT(INOUT) :: ISRUN(NETP%NEDGE)
    DOUBLE PRECISION, INTENT(INOUT) :: RUNVELS(NETP%NEDGE)
    LOGICAL, INTENT(IN) :: INITIALIZE
    INTEGER :: EC, CC, DEG
    DOUBLE PRECISION :: U, U2

    IF (RUNSPEED.GT.TINY(1D0)) THEN

       IF (INITIALIZE) THEN
          DO EC = 1,NETP%NEDGE            
             ! set direction of flow on this edge
             U = GRND()
             IF (U.LT.PSTART/(PSTART+PSTOP)) THEN ! start running in random dir
                ISRUN(EC) = .TRUE.
                U2= GRND()
                IF (U2.LE.0.5D0) THEN
                   RUNVELS(EC) = -RUNSPEED
                ELSE
                   RUNVELS(EC) = RUNSPEED
                ENDIF                
                IF (VERBOSE) PRINT*, 'Edge flow initiated:', EC, RUNVELS(EC)
             ELSE
                RUNVELS(EC) = 0D0
                ISRUN(EC) = .FALSE.
             END IF
          END DO
       ELSE ! update existing velocities
          DO EC = 1,NETP%NEDGE
             IF (ISRUN(EC)) THEN                
                ! edge is currently running, decide whether to stop
                U = GRND()
                IF (U.LT.PSTOP) THEN
                   ! stop flow
                   RUNVELS(EC) = 0D0
                   ISRUN(EC) = .FALSE.
                   IF (VERBOSE) PRINT*, 'Edge flow has stopped:', EC, RUNVELS(EC)
                ENDIF
             ELSE
                ! edge is stopped, decide whether to restart
                U = GRND()
                IF (U.LT.PSTART) THEN
                   ! start flow in random dir
                   ISRUN(EC) = .TRUE.
                   U2= GRND()
                   IF (U2.LE.0.5D0) THEN
                      RUNVELS(EC) = -RUNSPEED
                   ELSE
                      RUNVELS(EC) = RUNSPEED
                   ENDIF
                   IF (VERBOSE) PRINT*, 'Edge flow initiated:', EC, RUNVELS(EC)
                ENDIF
             END IF
          END DO
       ENDIF

       ! Copy run velocities to cell boundaries
       DSP%VEL = 0D0
       DO CC = 1,DSP%MESHP%NCELL
          DEG = DSP%MESHP%DEG(CC)
          DSP%VEL(CC,1:DEG) = RUNVELS(DSP%MESHP%BOUNDEDGE(CC,1:DEG))
       ENDDO
    ELSE
       RUNVELS = 0D0
       DSP%VEL = 0D0
    END IF
  END SUBROUTINE GETVELS_RANDVEL
  
  SUBROUTINE  GETVELS_RANDVELSWITCH(DSP,NETP,PSWITCH,RUNVELS,INITIALIZE,runspeed)
    ! update edge velocities, switching direction instantaneously
    ! RUNSPEED: flow speed
    ! PSWITCH: probability of switching on a given timestep
    ! RUNVELS: current velocity on each edge
    ! updates both RUNVELS and DSP%vel (flow across each cell boundary)
    ! INITIALIZE: if set, then need to initialize runvels for the first time
    USE KEYS, ONLY : VERBOSE
    USE mt19937, ONLY : GRND
    
    IMPLICIT NONE
    TYPE(DYNSYSTEM), POINTER :: DSP
    TYPE(NETWORK), POINTER :: NETP
    DOUBLE PRECISION, INTENT(IN) :: RUNSPEED, PSWITCH
    DOUBLE PRECISION, INTENT(INOUT) :: RUNVELS(NETP%NEDGE)
    LOGICAL, INTENT(IN) :: INITIALIZE
    INTEGER :: EC, CC, DEG
    DOUBLE PRECISION :: U, U2

    IF (RUNSPEED.GT.TINY(1D0).OR.ANY(NETP%FIXEDGEVEL)) THEN

       IF (INITIALIZE) THEN
          DO EC = 1,NETP%NEDGE
             IF (NETP%FIXEDGEVEL(EC)) THEN
                RUNVELS(EC) = NETP%FIXEDGEVELVAL(EC)
             ELSE
                ! set direction of flow on this edge
                U = GRND()
                IF (U.LE.0.5D0) THEN
                   RUNVELS(EC) = -RUNSPEED
                ELSE
                   RUNVELS(EC) = RUNSPEED                   
                END IF
             ENDIF
             IF (VERBOSE) PRINT*, 'Edge flow initiated:', EC, RUNVELS(EC)
          END DO
       ELSE ! Flip existing velocities
          DO EC = 1,NETP%NEDGE             
             ! edge is currently running, decide whether to switch
             IF (.NOT.NETP%FIXEDGEVEL(EC)) THEN
                U = GRND()
                IF (U.LT.PSWITCH) THEN
                   ! switch flow direction on this edge
                   RUNVELS(EC) = -RUNVELS(EC)
                   IF (VERBOSE) PRINT*, 'Edge flow has switched:', EC, RUNVELS(EC)
                END IF
             ENDIF
          END DO
       ENDIF

       !IF (VERBOSE) PRINT*, '# positive flows:', COUNT(RUNVELS>0)

       ! Copy run velocities to cell boundaries
       DSP%VEL = 0D0
       DO CC = 1,DSP%MESHP%NCELL
          DEG = DSP%MESHP%DEG(CC)
          DSP%VEL(CC,1:DEG) = RUNVELS(DSP%MESHP%BOUNDEDGE(CC,1:DEG))        
       ENDDO
    ELSE
       RUNVELS = 0D0
       DSP%VEL = 0D0
    END IF
  END SUBROUTINE GETVELS_RANDVELSWITCH

  SUBROUTINE GETVELS_FIXED(DSP,NETP,EDGEVELS)
    IMPLICIT NONE
    TYPE(DYNSYSTEM), POINTER :: DSP
    TYPE(NETWORK), POINTER :: NETP
    DOUBLE PRECISION, INTENT(IN) :: EDGEVELS(NETP%NEDGE)
    INTEGER :: CC, DEG
    
    ! Fix velocities to specific values on edges
    
    ! Copy run velocities to cell boundaries
       DSP%VEL = 0D0
       DO CC = 1,DSP%MESHP%NCELL
          DEG = DSP%MESHP%DEG(CC)
          DSP%VEL(CC,1:DEG) = EDGEVELS(DSP%MESHP%BOUNDEDGE(CC,1:DEG))        
       ENDDO
    
     END SUBROUTINE GETVELS_FIXED
  
  SUBROUTINE EULERSTEP(DSP,DELT,DFDT,FLUX)
    ! DFDT = time derivative of fields on all cells
    ! FLUX = flux going into fixed or absorber nodes
    IMPLICIT NONE
    TYPE(DYNSYSTEM), POINTER :: DSP
    DOUBLE PRECISION, INTENT(IN) :: DELT
    DOUBLE PRECISION, INTENT(OUT) :: DFDT(:,:)
    DOUBLE PRECISION, INTENT(OUT) :: FLUX(:,:)
    DOUBLE PRECISION :: FLUXDIFF(DSP%NFIELD)
    TYPE(MESH), POINTER :: MESHP
    INTEGER :: CC, DEG, FC, BC, BCt
    DOUBLE PRECISION :: WSHIFT(DSP%NFIELD), WAVG(DSP%NFIELD), FLUXADV(DSP%NFIELD)
    DOUBLE PRECISION :: DFREE

    IF (DSP%MESHP%USEGLOBALRESV) THEN
       PRINT*, 'ERROR: non-equilibrated euler steps with global reservoir are not currently implemented.'
       STOP 1
    ENDIF
    
    FLUX = 0D0
    DFDT = 0D0

    MESHP=>DSP%MESHP
    
    DO CC = 1,MESHP%NCELL
       DEG = MESHP%DEG(CC)       

       ! diffusive flux    
       FLUXDIFF = 0D0
       DO BCT = 1,DEG ! boundary counter
          BC = MESHP%BOUNDS(CC,BCT) ! boundary cell
          IF (BC.GT.0) THEN
             ! flux across each boundary defined as
             ! D*(w_(j+1)-w_j)/h+
             FLUXDIFF = FLUXDIFF + DSP%DCOEFF*(DSP%FIELDS(BC,:) - DSP%FIELDS(CC,:))/MESHP%LENPM(CC,BCt)
          ENDIF
       ENDDO

       ! Advective flux: via Lax-Wendroff discretization
       FLUXADV = 0D0
       DO BCT = 1,DEG
          ! approximate field values on boundary after half a timestep
          ! positive bounddir means + velocities point out from this cell
          BC = MESHP%BOUNDS(CC,BCT) ! boundary cell

          IF (BC.GT.0) THEN
             ! Weighted average of field on boundary
             
             WAVG = (MESHP%LEN(BC)/MESHP%DEG(BC)*DSP%FIELDS(CC,:) &
                  & + MESHP%LEN(CC)/MESHP%DEG(CC)*DSP%FIELDS(BC,:))/MESHP%LENPM(CC,BCT)
             WSHIFT = WAVG - DELT/2/MESHP%LENPM(CC,BCT)*DSP%VEL(CC,BCT)*MESHP%BOUNDDIR(CC,BCT)*(DSP%FIELDS(BC,:) - DSP%FIELDS(CC,:))
             
!             WSHIFT = DSP%FIELDS(CC,:) + &
!                  & 0.5*(1-MESHP%BOUNDDIR(CC,BCT)*DSP%VEL(CC,BCT)&
             !                  & *DELT/MESHP%LENPM(CC,BCT))*(DSP%FIELDS(BC,:) - DSP%FIELDS(CC,:))
             
             FLUXADV = FLUXADV - MESHP%BOUNDDIR(CC,BCT)*DSP%VEL(CC,BCT)*WSHIFT
          ENDIF
       ENDDO

       ! Turn off advective and diffusive flux for non-mobile fields
       DO FC = 1,DSP%NFIELD
          IF (.NOT.DSP%MOBILEFIELD(FC)) THEN
             FLUXDIFF(FC) = 0D0; FLUXADV(FC) = 0D0
          ENDIF
       ENDDO       
       
       DFDT(CC,:) = (FLUXDIFF+FLUXADV)/MESHP%VOL(CC)
       
       !IF (CC.EQ.1) PRINT*, 'TESTX2:', FLUXDIFF, FLUXADV, MESHP%VOL(CC)
       !if (cc.eq.1) print*, 'testx3:', meshp%lenpm(cc,:)       


       IF (DSP%BUFFERTYPE.EQ.1) THEN
          ! binding to buffer proteins with explicit on/off rates
          ! Free ligand
          !IF (CC.EQ.1) print*, 'TESTX1:', CC, DSP%KON, DSP%KOFF, DSP%FIELDS(CC,:)
          DFREE = - DSP%KON*DSP%FIELDS(CC,1)*DSP%FIELDS(CC,3) + DSP%KOFF*DSP%FIELDS(CC,2)
          DFDT(CC,1) = DFDT(CC,1) + DFREE
          ! Bound ligand
          DFDT(CC,2) = DFDT(CC,2) - DFREE
          ! free protein
          DFDT(CC,3) = DFDT(CC,3) + DFREE
       ELSEIF (DSP%BUFFERTYPE.EQ.2) THEN
          PRINT*, 'ERROR IN EULER STEP: rapid equilibration not yet set up'
          STOP 1
       ENDIF
       IF (DSP%DOACTIVATION) THEN
          ! two fields: activated and inactivated molecule
          ! change in active molecules
          DFDT(CC,1) = DFDT(CC,1) + DSP%ACTRATE(CC)*DSP%FIELDS(CC,2) - DSP%DEPRATE(CC)*DSP%FIELDS(CC,1)
          ! DEBUGGING
!          DFDT(CC,1) = DFDT(CC,1) + DSP%ACTRATE(CC)*(1-DSP%FIELDS(CC,1))
          ! change in inactive molecules
          DFDT(CC,2) = DFDT(CC,2) - DSP%ACTRATE(CC)*DSP%FIELDS(CC,2)
       ENDIF
      
       IF (DSP%ISPERM(CC)) THEN
          ! connection with external concentrations at permeable node                    
          ! can have different permeability for each field
          ! save flux to external environment only
          FLUX(CC,:) = -DSP%PERM(CC,:)*(DSP%CEXT - DSP%FIELDS(CC,:))
          DFDT(CC,:) = DFDT(CC,:) - FLUX(CC,:)/MESHP%VOL(CC)
       ELSE
          ! track overall change in conc at this node
          FLUX(CC,:) = DFDT(CC,:)*MESHP%VOL(CC)
       ENDIF
      
    ENDDO

    ! no change in fixed cells
    DO FC = 1,DSP%NFIELD
       DO CC = 1,MESHP%NCELL         
          ! IF (DSP%ACTRATE(CC,FC).GT.0) THEN
          !    PRINT*, 'TESTX2:', FC, CC, DSP%ISFIXED(CC,FC), DSP%ACTRATE(CC,FC), DSP%FIELDS(CC,FC)  
          ! ENDIF
          IF (DSP%ISFIXED(CC,FC)) DFDT(CC,FC) = 0D0
       ENDDO
    ENDDO
  END SUBROUTINE EULERSTEP


  SUBROUTINE EULERSTEPEQUIL(DSP,DELT,DFDT,FLUX)
    ! Forward Euler step for explicit rapid equilibration
    ! Assumes dynamic system has 2 fields: (1) free ligand L and (2) total protein T
    
    ! DFDT = time derivative of fields on all cells
    ! FLUX = flux going into fixed or absorber nodes
    
    IMPLICIT NONE
    TYPE(DYNSYSTEM), POINTER :: DSP
    DOUBLE PRECISION, INTENT(IN) :: DELT
    DOUBLE PRECISION, INTENT(OUT) :: DFDT(:,:)
    DOUBLE PRECISION, INTENT(OUT) :: FLUX(:,:)
    DOUBLE PRECISION :: FLUXDIFF(DSP%NFIELD)
    DOUBLE PRECISION :: FLUXPUMP(DSP%NFIELD), TOTFLUXPUMP
    TYPE(MESH), POINTER :: MESHP
    INTEGER :: CC, DEG, FC, BC, BCt
    DOUBLE PRECISION :: WSHIFT(DSP%NFIELD), WAVG(DSP%NFIELD), FLUXADV(DSP%NFIELD)
    DOUBLE PRECISION :: DFREE
    DOUBLE PRECISION :: BFIELD(DSP%MESHP%NCELL), CFIELD(DSP%MESHP%NCELL), LKD
    DOUBLE PRECISION :: DR, DSCL, A1, A2, ABOUND,tmp
    DOUBLE PRECISION :: GLOBALRESVFIELD(DSP%NFIELD), TOTPERMFLUX

    
    IF (DSP%NFIELD.NE.2.or.DSP%BUFFERTYPE.NE.2.OR.DSP%DOACTIVATION) THEN
       PRINT*, 'ERROR: EULERSTEPEQUIL only works with 2 fields and Buffer Type 2. No activation'
       STOP 1
    ENDIF

    FLUXPUMP = 0D0    
    FLUX = 0D0
    DFDT = 0D0

    MESHP=>DSP%MESHP
    
    IF (MESHP%USEGLOBALRESV) THEN
       ! initial field value in global reservoir
       GLOBALRESVFIELD = 0D0
       GLOBALRESVFIELD(1) = DSP%FIELDS(MESHP%GLOBALRESVIND,1)
    END IF
    
    ! get field of bound ligand
    BFIELD = DSP%FIELDS(:,1)*DSP%FIELDS(:,2)/(DSP%FIELDS(:,1)+DSP%KDEQUIL)
    ! Get field of total ligand
    CFIELD = DSP%FIELDS(:,1) + BFIELD

    TOTFLUXPUMP = 0D0
    DO CC = 1,MESHP%NCELL
       DEG = MESHP%DEG(CC)       

       ! IF (ANY(DSP%FIELDS(CC,:).LT.-1D-10)) THEN
       !    PRINT*, 'Negative field!'
       !    PRINT*, CC, MESHP%CELLTYPE(CC), MESHP%NODEIND(CC), MESHP%EDGEIND(CC,1), MESHP%RESVIND(CC)
       !    PRINT*, DSP%FIELDS(CC,:)
       !    STOP 1
       ! ENDIF       
       
       ! diffusive flux (of total ligand and total protein)
       FLUXDIFF = 0D0; 
       DO BCT = 1,DEG ! boundary counter
          BC = MESHP%BOUNDS(CC,BCT) ! boundary cell
          IF (BC.GT.0) THEN
             IF (MESHP%BOUNDCLOSED(CC,BCT)) CYCLE ! this boundary is closed off

             IF (DSP%CONC3D) THEN
                ! using meshed reservoir elements, work entirely with 3D concentrations
                ! scaling factor for diffusivity
                ! Only for boundaries along an edge, only if using varrad
                IF (DSP%VARRAD.AND.MESHP%CELLTYPE(CC).EQ.1.AND.MESHP%CELLTYPE(BC).EQ.1) THEN
                   DR = MESHP%RAD(BC)-MESHP%RAD(CC)
                   ! This is formula 1.8 in Berezhkovskii, 2007 (originally from Reguerra, 2001)
                   DSCL = SQRT(1 + (DR/MESHP%LENPM(CC,BCT))**2)
                ELSE
                   DSCL = 1D0
                ENDIF
                
                ! area of cell boundary
                ABOUND = MESHP%BOUNDAREA(CC,BCT)

                ! flux for total ligand
                ! This is formula 1.6 from Berezhkovskii, 2007 (where FIELDS is c/A, the 3D concentration)
                FLUXDIFF(1) = FLUXDIFF(1) &
                     & + DSP%DCOEFF(1)/DSCL*ABOUND*(DSP%FIELDS(BC,1) &
                     & - DSP%FIELDS(CC,1))/MESHP%LENPM(CC,BCt) &
                     & + DSP%DCOEFF(2)/DSCL*ABOUND*(BFIELD(BC) - BFIELD(CC))/MESHP%LENPM(CC,BCT)
                
                IF (.NOT.DSP%UNIFORMBUFFER) THEN ! spatially varying buffer sites
                   ! flux for total protein
                   FLUXDIFF(2) = FLUXDIFF(2) + &
                        & DSP%DCOEFF(2)/DSCL*ABOUND*&
                        & (DSP%FIELDS(BC,2) - DSP%FIELDS(CC,2))/MESHP%LENPM(CC,BCT)
                ENDIF
             ELSEIF (DSP%VARRAD) THEN
                ! mesh cells can have varying radii
                ! but still working with 1D concentrations
                ! flux across boundary = D(x) A(x) d/dx[C/A]
                ! D(x) = D0/sqrt(1+R'^2)

                ! scaling factor for diffusivity
                DR = MESHP%RAD(BC)-MESHP%RAD(CC)
                DSCL = SQRT(1 + (DR/MESHP%LENPM(CC,BCT))**2)
                ! areas in the two membrane cell, and average at boundary
                A1 = PI*MESHP%RAD(BC)**2; A2 = PI*MESHP%RAD(CC)**2
                ABOUND = (A1 + A2)/2
                !ABOUND = MESHP%BOUNDAREA(CC,BCT)
                
                ! flux for total ligand
                FLUXDIFF(1) = FLUXDIFF(1) &
                     & + DSP%DCOEFF(1)/DSCL*ABOUND*(DSP%FIELDS(BC,1)/A1 &
                     & - DSP%FIELDS(CC,1)/A2)/MESHP%LENPM(CC,BCt) &
                     & + DSP%DCOEFF(2)/DSCL*ABOUND*(BFIELD(BC)/A1 - BFIELD(CC)/A2)/MESHP%LENPM(CC,BCT)
                IF (.NOT.DSP%UNIFORMBUFFER) THEN ! spatially varying buffer conc
                   ! flux for total protein
                   FLUXDIFF(2) = FLUXDIFF(2) + &
                        & DSP%DCOEFF(2)/DSCL*ABOUND*&
                        & (DSP%FIELDS(BC,2)/A1 - DSP%FIELDS(CC,2)/A2)/MESHP%LENPM(CC,BCT)
                ENDIF
             ELSE
                ! flux across each boundary defined as
                ! D*(w_(j+1)-w_j)/h+
                ! get diffusive flux of TOTAL LIGAND

                FLUXDIFF(1) = FLUXDIFF(1) &
                     & + DSP%DCOEFF(1)*(DSP%FIELDS(BC,1) - DSP%FIELDS(CC,1))/MESHP%LENPM(CC,BCt) &
                     & + DSP%DCOEFF(2)*(BFIELD(BC) - BFIELD(CC))/MESHP%LENPM(CC,BCt)               
                IF (.NOT.DSP%UNIFORMBUFFER) THEN  ! spatially varying buffer conc
                   ! get diffusive flux of  total protein 
                   FLUXDIFF(2) = FLUXDIFF(2) + DSP%DCOEFF(2)*(DSP%FIELDS(BC,2) - DSP%FIELDS(CC,2))/MESHP%LENPM(CC,BCt)
                ENDIF
             END IF
          ENDIF
       ENDDO

       ! Advective flux: via Lax-Wendroff discretization
       ! WARNING: flows are not set up with reservoir elements or with varrad
       FLUXADV = 0D0
       IF (DSP%USEEDGEFLOW) THEN
          IF (DSP%VARRAD.OR.DSP%CONC3D) THEN
             PRINT*, 'USEVARRAD AND CONC3D not set up together with edge flows'
             STOP 1
          ENDIF
          
          ! WARNING: this has not been thoroughly thought through for the case where there are flows
          ! going across boundaries to large reservoirs!       
          DO BCT = 1,DEG
             ! approximate field values on boundary after half a timestep
             ! positive bounddir means + velocities point out from this cell
             BC = MESHP%BOUNDS(CC,BCT) ! boundary cell

             IF (BC.GT.0) THEN
                IF (MESHP%BOUNDCLOSED(CC,BCT)) CYCLE ! this boundary is closed off

                ! TOTAL ligand
                ! Weighted average of field on boundary
               
                IF (.NOT.DSP%MOBILEFIELD(2)) THEN
                   ! immobile proteins, only free ligand feels flow
                   WAVG(1) = (MESHP%LEN(BC)/MESHP%DEG(BC)*DSP%FIELDS(CC,1) &
                        & + MESHP%LEN(CC)/MESHP%DEG(CC)*DSP%FIELDS(BC,1))/MESHP%LENPM(CC,BCT)
                   WSHIFT(1) = WAVG(1) - DELT/2/MESHP%LENPM(CC,BCT)*DSP%VEL(CC,BCT)*MESHP%BOUNDDIR(CC,BCT)&
                        & *(DSP%FIELDS(BC,1) - DSP%FIELDS(CC,1))
                   WAVG(2) =0D0; WSHIFT(2) = 0D0
                ELSE
                   WAVG(1) = (MESHP%LEN(BC)/MESHP%DEG(BC)*CFIELD(CC) &
                        & + MESHP%LEN(CC)/MESHP%DEG(CC)*CFIELD(BC))/MESHP%LENPM(CC,BCT)
                   WSHIFT(1) = WAVG(1) - DELT/2/MESHP%LENPM(CC,BCT)*DSP%VEL(CC,BCT)*MESHP%BOUNDDIR(CC,BCT)*(CFIELD(BC) - CFIELD(CC))

                   ! total protein
                   WAVG(2) = (MESHP%LEN(BC)/MESHP%DEG(BC)*DSP%FIELDS(CC,2) &
                        & + MESHP%LEN(CC)/MESHP%DEG(CC)*DSP%FIELDS(BC,2))/MESHP%LENPM(CC,BCT)
                   WSHIFT(2) = WAVG(2) - &
                        & DELT/2/MESHP%LENPM(CC,BCT)*DSP%VEL(CC,BCT)*MESHP%BOUNDDIR(CC,BCT)&
                        & *(DSP%FIELDS(BC,2) - DSP%FIELDS(CC,2))                     
                ENDIF
                ! advective flux for total ligand and total protein
                FLUXADV = FLUXADV - MESHP%BOUNDDIR(CC,BCT)*DSP%VEL(CC,BCT)*WSHIFT
             ENDIF
          ENDDO
       ENDIF

       ! Flux into this mesh cell due to pumping from the global reservoir
       IF (MESHP%USEGLOBALRESV.and.MESHP%CELLTYPE(CC).NE.3) THEN
          FLUXPUMP(1) = MESHP%SA(CC)*DSP%GLOBALRESVKR*GLOBALRESVFIELD(1)/(GLOBALRESVFIELD(1)+DSP%GLOBALRESVKMR)
          TOTFLUXPUMP = TOTFLUXPUMP + FLUXPUMP(1)
       END IF
       
       DFDT(CC,:) = (FLUXDIFF+FLUXADV+FLUXPUMP)/MESHP%VOL(CC)
     
       ! DO FC = 1,DSP%NFIELD
       !    IF (.NOT.DSP%MOBILEFIELD(FC)) THEN
       !       ! field not allowed to move
       !       DFDT(CC,FC) = 0D0
       !    ENDIF
       ! END DO

       ! Returned FLUX array contains flux for each permeable cell
       ! in units of particles per second
       FLUX(CC,:) = DFDT(CC,:)*MESHP%VOL(CC)

       IF (DSP%ISPERM(CC)) THEN
          ! connection with external concentrations at permeable node                    
          ! can have different permeability for each field
          ! save flux to external environment only
          IF (MESHP%USEGLOBALRESV.AND.DSP%PERMTOGLOBALRESV) THEN
             ! permeability flux depends on concentration in global reservoir
             FLUX(CC,:) = -DSP%PERM(CC,:)*(GLOBALRESVFIELD - DSP%FIELDS(CC,:))
          ELSE
             FLUX(CC,:) = -DSP%PERM(CC,:)*(DSP%CEXT - DSP%FIELDS(CC,:))
          ENDIF
          DFDT(CC,:) = DFDT(CC,:) - FLUX(CC,:)/MESHP%VOL(CC)
       ENDIF
       
       IF (.NOT.DSP%TRACKDCDT) THEN
          ! get change in free ligand from delta total lig and delta total prot
          LKD = DSP%FIELDS(CC,1) + DSP%KDEQUIL;
          IF (DSP%UNIFORMBUFFER) THEN ! spatially constant buffer
             DFDT(CC,1) = DFDT(CC,1)/ &            
                  & (1 + DSP%FIELDS(CC,2)*DSP%KDEQUIL/LKD**2)
          ELSE
             DFDT(CC,1) = (DFDT(CC,1) - DSP%FIELDS(CC,1)*DFDT(CC,2)/LKD)/ &            
                  & (1 + DSP%FIELDS(CC,2)*DSP%KDEQUIL/LKD**2)
          ENDIF
       ENDIF
    ENDDO
    
    IF (MESHP%USEGLOBALRESV) THEN             
       CC = MESHP%GLOBALRESVIND
       FLUX(CC,1) = -(TOTFLUXPUMP &
            & + DSP%GLOBALRESVKOUT*GLOBALRESVFIELD(1)/(GLOBALRESVFIELD(1)+DSP%GLOBALRESVKMOUT) &
            & )
       IF (DSP%PERMTOGLOBALRESV) THEN
          ! total flux out of permeable nodes   
          TOTPERMFLUX = SUM(pack(FLUX(:,1),DSP%ISPERM))
          FLUX(CC,1) = FLUX(CC,1) + TOTPERMFLUX/MESHP%VOL(CC)
       ENDIF
       DFDT(CC,1) = FLUX(CC,1)/MESHP%VOL(CC)
    END IF

   
    
    ! check that total mass is not changing
    !GLOBALTEST = GLOBALTEST+SUM(FLUX(:,1))
    !PRINT*, 'TESTX1:', SUM(FLUX(:,1)), GLOBALTEST, SUM(CFIELD*DSP%MESHP%VOL)   
    
    
    ! no change in fixed cells
    DO FC = 1,DSP%NFIELD
       DO CC = 1,MESHP%NCELL         
          ! IF (DSP%ACTRATE(CC,FC).GT.0) THEN
          !    PRINT*, 'TESTX2:', FC, CC, DSP%ISFIXED(CC,FC), DSP%ACTRATE(CC,FC), DSP%FIELDS(CC,FC)  
          ! ENDIF
          IF (DSP%ISFIXED(CC,FC)) DFDT(CC,FC) = 0D0
       ENDDO
    ENDDO
   
  END SUBROUTINE EULERSTEPEQUIL
END MODULE FVMDYNAMICS
