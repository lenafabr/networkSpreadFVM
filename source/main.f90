PROGRAM MAIN
  USE NETWORKUTIL, ONLY : NETWORK, CLEANUPNETWORK
  USE KEYS, ONLY : ACTION
  IMPLICIT NONE

  TYPE(NETWORK), TARGET :: NET
  TYPE(NETWORK), POINTER :: NETP

  NETP=>NET

  CALL READKEY

  SELECT CASE(ACTION)
  CASE('RUNDYNAMICS')
     CALL RUNDYNAMICSDRIVER(NETP,ACTION)
  CASE DEFAULT
     PRINT*, 'Not a valid action:', ACTION
  END SELECT

  CALL  CLEANUPNETWORK(NETP)
  
CONTAINS  
  
  SUBROUTINE RUNDYNAMICSDRIVER(NETP,ACTION)
    ! run dynamic simulation of field evolution over a network
    USE KEYS, ONLY : NSTEP, DELT, PRINTEVERY, SNAPSHOTEVERY, SNAPSHOTFILE, &
         & NETFILE, MINMESHPT, MAXMESHSIZE, NSTARTEDGE, STARTEDGES, &
         & NSTARTNODE, STARTNODES,OUTPUTEVERY, MESHFILE,&
         & MAXSTARTEDGE, STARTNODERAD, VELCONTROL,&
         & RESVVOL, RESVMIX, STARTCONC, startequil, NFIELD, NCONT, &
         & STARTRATE, STOPRATE, RUNSPEED, RANDFIXNODES, FIXVALS, FIXNODES, NFIX,&
         & NSTARTPOS,STARTPOS
    USE NETWORKUTIL, ONLY : NETWORKFROMFILE!, SETUPFIXNODES
    USE MESHUTIL, ONLY : MESH, SETUPNETWORKMESH, CLEANUPMESH, OUTPUTMESH
    USE DYNSYSUTIL, ONLY : DYNSYSTEM, SETUPDYNSYS, SETPARAMDYNSYS, &
         & CLEANUPDYNSYS, INITIALIZEFIELDNODES, INITIALIZEFIELDNEARPOS, &
         & INITIALIZEFIELDEDGES, GETINITCONC
    USE FVMDYNAMICS, ONLY : RUNDYNAMICS
    USE NODECONTRACTUTIL, ONLY : INITIALIZE_NODECONTRACT
    
    IMPLICIT NONE
    TYPE(MESH), TARGET :: MESHOBJ
    TYPE(MESH), POINTER :: MESHP
    TYPE(DYNSYSTEM), TARGET :: DYNSYS
    TYPE(DYNSYSTEM), POINTER :: DSP
    TYPE(NETWORK), POINTER :: NETP
    CHARACTER(LEN=*) :: ACTION
    TYPE(NETWORK), TARGET :: NETAUX
    INTEGER :: CC, RC, NC
    DOUBLE PRECISION :: INITCONC(NFIELD)
    DOUBLE PRECISION, ALLOCATABLE :: FIELDSVELS(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: STARTPOSNODES(:,:)
    
    MESHP=>MESHOBJ
    DSP=>DYNSYS
    
    ! set up network structure
    CALL NETWORKFROMFILE(NETP,NETFILE)    

     ! set up reservoir volumes, external concentration
     DO RC = 1,NETP%NRESV
       NETP%RESVVOL(RC) = RESVVOL(RC)
       NETP%RESVMIX(RC) = RESVMIX(RC)
    ENDDO
     
     ! set up mesh on network
     CALL SETUPNETWORKMESH(MESHP,NETP,MAXMESHSIZE,MINMESHPT)     

     
     ! set up and parameterize dynamical system    
     CALL SETUPDYNSYS(DSP,MESHP,NFIELD)
     CALL SETPARAMDYNSYS(DSP,NETP)


     
     ! Set up contraction information
     IF (VELCONTROL.EQ.'NODECONTRACTIONS') THEN
        CALL INITIALIZE_NODECONTRACT(NETP,NCONT)
        PRINT*, 'INITIAL CONTRACTIONS:'
        DO CC = 1,NETP%NCONT
           PRINT*, CC, NETP%CONTNODES(CC,:), NETP%NODEVOLS(NETP%CONTNODES(CC,:))
        ENDDO
     ENDIF
    
    IF (ACTION.EQ.'RUNDYNAMICS_CONTRACTION') THEN
       ! initialize contraction positions
       PRINT*, 'ERROR IN DRIVER: contractions not yet set up for FVM simulations'
       STOP 1
    ENDIF
   
    ! Dump out mesh structure
    PRINT*, 'outputting mesh to: ', MESHFILE
    CALL OUTPUTMESH(MESHP,MESHFILE)

    ! initialize the fields
    CALL GETINITCONC(DSP%NFIELD,STARTCONC,STARTEQUIL,DSP%KDEQUIL,INITCONC)

    IF (NSTARTNODE.GT.0) THEN
       IF (STARTNODERAD.GT.0) THEN         
          ALLOCATE(STARTPOSNODES(NSTARTNODE,NETP%DIM))
          DO NC = 1,NSTARTNODE
             STARTPOSNODES(NC,:) = NETP%NODEPOS(STARTNODES(NC),:)
          ENDDO
          CALL INITIALIZEFIELDNEARPOS(NETP,DSP,STARTPOSNODES,INITCONC,STARTNODERAD)
          DEALLOCATE(STARTPOSNODES)          
       ELSE
          CALL INITIALIZEFIELDNODES(NETP,DSP,STARTNODES(1:NSTARTNODE),INITCONC)
       ENDIF
    ELSEIF(NSTARTPOS.GT.0) THEN ! Initialize near specific positions
       CALL INITIALIZEFIELDNEARPOS(NETP,DSP,STARTPOS(1:NSTARTPOS,1:NETP%DIM),INITCONC,STARTNODERAD)
    ELSEIF(NSTARTEDGE.GT.0) THEN              
       CALL INITIALIZEFIELDEDGES(NETP,DSP,STARTEDGES(1:NSTARTEDGE),INITCONC)
    ELSE
       PRINT*, 'ERROR: nstartnode and nstartedge both zero. Dont know how to initialize.'
       STOP 1
    ENDIF

    CALL RUNDYNAMICS(NETP,DSP,NSTEP, DELT, VELCONTROL,STARTRATE,STOPRATE,RUNSPEED)
    
    CALL CLEANUPMESH(MESHP)
    CALL CLEANUPDYNSYS(DSP)
    
  END SUBROUTINE RUNDYNAMICSDRIVER
  
END PROGRAM MAIN
