MODULE FIELDDYNAMICS
  ! dynamic evolution of fields on a mesh
  IMPLICIT NONE

CONTAINS
  SUBROUTINE RUNDIFFDYNAMICS(NETP, DELT, MAXSTEPS)
    ! diffusive dynamics on a meshed networks
    ! DELT: timestep
    ! MAXSTEPS: max number of timesteps to run for

    ! ------ THIS SUBROUTINE HAS NOT BEEN COMPILED OR TESTED YET -------
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    DOUBLE PRECISION, INTENT(IN) :: DELT
    INTEGER, INTENT(IN) :: MAXSTEPS

    IF (.NOT.(NETP%ARRAYSET.AND.NETP%STRUCTURESET)) THEN
       PRINT*, 'ERROR IN RUNDIFFDYNAMICS: network structure is not set up'
       STOP 1
    ENDIF
    IF (.NOT.NETP%MESHSET) THEN
       PRINT*, 'ERROR IN RUNDIFFUSIONDYNAMICS: mesh on network is not set up'
       STOP 1
    ENDIF

    CURTIME = 0D0
    FIELD0 = FIELD
    
    DO STEP = 1,MAXSTEPS
       ! take a single euler step
       CALL EULERSTEPDIFF(NETP, CONC, NEWCONC)
       
       FIELD = NEWFIELD
       CURTIME = CURTIME + DELT

       IF (MOD(STEP,SNAPSHOTEVERY).EQ.0) THEN
          ! dump out a snapshot of the mesh
          CALL OUTPUTMESHFIELD(NETP,FIELD,SNAPSHOTFILE,.TRUE.)
       ENDIF
    ENDDO
    
  END SUBROUTINE RUNDIFFDYNAMICS

  SUBROUTINE GETTIMEDERIV(NETP, CONC, D, DCDT)
    ! calculate time derivative in concentration field
    ! CONC: concentration field
    ! D: diffusivity
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    DOUBLE PRECISION, INTENT(IN) :: CONC(:,:), D
    DOUBLE PRECISION, INTENT(OUT) :: DCDT(NETP%MESHDIM1, NETP%NEDGE)
    

    ! -------------
    ! get diffusive terms
    ! ---------------
    
    ! for each mesh point along an edge (no end-points)
    DO EC = 1, NETP%NEDGE
       NPT = NETP%NMESH(EC)
       !second derivative approximation
       D2C(2:NPT-1,EC) =  (CONC(1:NPT-2,EC) + CONC(3:NPT,EC) - 2*CONC(2:NPT-1,EC))/NETP%DX(EC)**2
       DCDT(2:NPT-1,EC) = D*D2C(2:NPT-1,EC)
    ENDDO
    
    ! get total flux into each node
    FLUX = 0D0
    DO NC = 1,NETP%NNODE
       FLUX(NC) = 0
       DO CC = 1,NETP%NODEDEG(NC)
          EC = NETP%NODEEDGE(NC,CC)
          IF (NETP%EDGENODE(EC,1).EQ.NC) THEN
             ! this is an outgoing edge
             FLUX(NC) = FLUX(NC) + (CONC(2,EC)-CONC(1,EC))/NETP%DX(EC)**2
          ELSEIF (NETP%EDGENODE(EC,2).EQ.NC) THEN
             ! this is an incoming edge
             NPT = NETP%NMESH(EC)
             FLUX(NC) = FLUX(NC) - (CONC(NPT,EC)-CONC(NPT-1,EC))/NETP%DX(EC)**2
          ELSE
             PRINT*, 'ERROR IN GETTIMEDERIV: &
                  & something is wrong in node / edge connectivity', &
                  & NC, EC, NETP%EDGENODE(EC,:)
             STOP 1
          ENDIF             
       ENDDO
       ! derivative for the endpoints of each relevant edge
       ! WARNING: this is redundant (tracks the node points multiple times)
       DO CC = 1,NETP%NODEDEG(NC)          
          DCDT()
       ENDDO
    ENDDO
  END SUBROUTINE GETTIMEDERIV
END MODULE FIELDDYNAMICS
