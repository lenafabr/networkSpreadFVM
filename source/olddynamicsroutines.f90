   SUBROUTINE EULERSTEP_ADVDIFF_UPWIND(NETP,FIELD,Dcoeff,VELMESH,DFDT,FLUXSAVE)
    ! calculate a forward euler propagation step (dF/dt)
     ! assuming the field is spreading by advection-diffusion
     ! DCOEFF is the diffusion coefficient
     ! VELMESH gives the flow velocity on each mesh point
     ! Uses upwind differencing scheme for advection component of the time derivative
    USE MESHUTIL, ONLY : UPDATEFIELDNODEVALS
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    DOUBLE PRECISION, INTENT(IN) :: FIELD(NETP%MESHDIM1,NETP%NEDGE), VELMESH(NETP%MESHDIM1,NETP%NEDGE,2), DCOEFF
    DOUBLE PRECISION, INTENT(OUT) :: DFDT(NETP%MESHDIM1,NETP%NEDGE),FLUXSAVE(NETP%NNODE)
    INTEGER :: EC, NC, N1, N2, PC, NPT
    DOUBLE PRECISION :: DX,  FLUX(NETP%NNODE), VF(3)

    DFDT = 0D0
    
    FLUX = 0D0
    
    DO EC = 1,NETP%NEDGE ! for each edge
      ! print*, 'testx1:', ec, field(1,ec)
       
       NPT = NETP%NMESH(EC) ! numer of meshpoints along edge
       DX = NETP%DX(EC) ! spacing between mesh points
                    
       DO PC = 2,NPT-1
          ! central approximation for 2nd derivative for intermediate points
       !   print*, 'TESTX1:', DX, FIELD(PC-1:PC+1,EC), DCOEFF
          DFDT(PC,EC) = DCOEFF*(FIELD(PC-1,EC)-2*FIELD(PC,EC)+FIELD(PC+1,EC))/DX**2
       ENDDO

       ! advective terms      
       DO PC = 1,NPT
           ! incoming flow
          IF (PC.GT.1) THEN
           !  IF (NETP%MESHISCONT(PC-1,EC)) THEN ! prior point is a contraction
                
            ! ELSE
                IF (VELMESH(PC-1,EC,2).GE.0) THEN
                   DFDT(PC,EC) = DFDT(PC,EC) + VELMESH(PC-1,EC,2)*FIELD(PC-1,EC)/DX
                ENDIF
             !ENDIF
          ENDIF
          IF (PC.LT.NPT) THEN
             IF (VELMESH(PC+1,EC,1).LT.0) THEN
                DFDT(PC,EC) = DFDT(PC,EC) - VELMESH(PC+1,EC,1)*FIELD(PC+1,EC)/DX
             ENDIF
          ENDIF

          ! outgoing flow
          IF (VELMESH(PC,EC,1).LT.0) THEN
             DFDT(PC,EC) = DFDT(PC,EC) + VELMESH(PC,EC,1)*FIELD(PC,EC)/DX
          ENDIF
          IF (VELMESH(PC,EC,2).GT.0) THEN
             DFDT(PC,EC) = DFDT(PC,EC) - VELMESH(PC,EC,2)*FIELD(PC,EC)/DX
          ENDIF
          
       ENDDO
      
       ! add overall flux from this edge into end-points
       N1 = NETP%EDGENODE(EC,1)
       ! diffusive flux and advective flux
       FLUX(N1) = FLUX(N1) + DCOEFF*(FIELD(2,EC)-FIELD(1,EC))/DX !- VELMESH(2,EC)*FIELD(2,EC)
       IF (VELMESH(1,EC,2).LT.0) THEN ! inflow from this edge
          FLUX(N1) = FLUX(N1) - VELMESH(1,EC,2)*FIELD(2,EC)
       ELSE ! outflow to this edge
          FLUX(N1) = FLUX(N1) - VELMESH(1,EC,2)*FIELD(1,EC)
       ENDIF
       N2 = NETP%EDGENODE(EC,2)
       FLUX(N2) = FLUX(N2) - DCOEFF*(FIELD(NPT,EC)-FIELD(NPT-1,EC))/DX !+ VELMESH(NPT-1,EC)*FIELD(NPT-1,EC)
       IF (VELMESH(NPT,EC,1).GT.0) THEN ! inflow from this edge
          FLUX(N2) = FLUX(N2) + VELMESH(NPT,EC,1)*FIELD(NPT-1,EC)
       ELSE !outflow to this edge
          FLUX(N2) = FLUX(N2) + VELMESH(NPT,EC,1)*FIELD(NPT,EC)
       ENDIF
    END DO

    ! save the flux before scaling and setting to zero at absorbing nodes
    FLUXSAVE = FLUX
    
    ! scale flux at each node by the overall mesh length associated with that node

    DO NC = 1,NETP%NNODE
       IF (NETP%NODEABS(NC).OR.NETP%NODEFIX(NC)) THEN
          ! absorbing node: fixed concentration
          FLUX(NC) = 0D0
       ELSE             
          FLUX(NC) = FLUX(NC)/NETP%NODEDX(NC)
       ENDIF
    ENDDO

    
    ! redistribute the scaled flux onto the endpoints in the mesh array
    CALL UPDATEFIELDNODEVALS(NETP,FLUX,DFDT)
   
    
  END SUBROUTINE EULERSTEP_ADVDIFF_UPWIND
  
  SUBROUTINE EULERSTEP_DIFFUSION(NETP,FIELD, DFDT,FLUXSAVE)
    ! calculate a forward euler propagation step (dF/dt)
    ! assuming purely diffusive spread of particles
    USE MESHUTIL, ONLY : UPDATEFIELDNODEVALS
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    DOUBLE PRECISION, INTENT(IN) :: FIELD(NETP%MESHDIM1,NETP%NEDGE)
    DOUBLE PRECISION, INTENT(OUT) :: DFDT(NETP%MESHDIM1,NETP%NEDGE),FLUXSAVE(NETP%NNODE)
    INTEGER :: EC, NC, N1, N2, PC, NPT
    DOUBLE PRECISION :: DX, AVGDX, FLUX(NETP%NNODE)

    FLUX = 0D0
    
    DO EC = 1,NETP%NEDGE ! for each edge
      ! print*, 'testx1:', ec, field(1,ec)
       
       NPT = NETP%NMESH(EC) ! numer of meshpoints along edge
       DX = NETP%DX(EC) ! spacing between mesh points

              
       ! central approximation for 2nd derivative for intermediate points
       DO PC = 2,NPT-1
          DFDT(PC,EC) = (FIELD(PC-1,EC)-2*FIELD(PC,EC)+FIELD(PC+1,EC))/DX**2
       END DO

       ! add overall flux from this edge into end-points
       N1 = NETP%EDGENODE(EC,1)
      ! print*, 'testx4', ec, n1, field(1:2,ec), dx
       FLUX(N1) = FLUX(N1) + (FIELD(2,EC)-FIELD(1,EC))/DX
       N2 = NETP%EDGENODE(EC,2)
       FLUX(N2) = FLUX(N2) - (FIELD(NPT,EC)-FIELD(NPT-1,EC))/DX
    END DO

    ! save the flux before scaling and setting to zero at absorbing nodes
    FLUXSAVE = FLUX
    
    ! scale flux at each node by the overall mesh length associated with that node

    DO NC = 1,NETP%NNODE
       IF (NETP%NODEABS(NC)) THEN
          ! absorbing node: fixed concentration
          FLUX(NC) = 0D0
       ELSE
          FLUX(NC) = FLUX(NC)/NETP%NODEDX(NC)
       ENDIF
    ENDDO

    
    ! redistribute the scaled flux onto the endpoints in the mesh array
    CALL UPDATEFIELDNODEVALS(NETP,FLUX,DFDT)
   
    
  END SUBROUTINE EULERSTEP_DIFFUSION

SUBROUTINE EULERDYNAMICS(NETP,FIELDS, NSTEP, DELT, PRINTEVERY,SNAPSHOTEVERY,SNAPSHOTFILE,OUTPUTEVERY,OUTFILE)
    ! run forward euler dynamics, evolving the field with time
    ! NSTEP = total number of steps to run    
    ! DELT = timestep
    ! SNAPSHOTEVERY = how often to dump snapshots
    ! SNAPSHOTFILE = file in which to dump snapshots
    USE MESHUTIL, ONLY : OUTPUTMESHFIELD, INTEGRATEMESHFIELD
    
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    DOUBLE PRECISION, INTENT(INOUT) :: FIELDS(:,:,:)
    INTEGER, INTENT(IN) :: NSTEP, PRINTEVERY, SNAPSHOTEVERY,OUTPUTEVERY
    DOUBLE PRECISION, INTENT(IN) :: DELT
    CHARACTER*100, INTENT(IN) :: SNAPSHOTFILE, OUTFILE
    DOUBLE PRECISION :: FIELD0(NETP%MESHDIM1,NETP%NEDGE), DFDT(NETP%MESHDIM1,NETP%NEDGE), FLUX(NETP%NNODE)
    INTEGER :: STEP, N1, N2, NPT, EC, ABSORBERNODES(NETP%NNODE), NABS, OU, NC
    DOUBLE PRECISION :: DX, TOTINT, TOTFLUX
    DOUBLE PRECISION :: CURTIME

    
    FIELD0 = FIELDS(:,:,1)
    DFDT = 0D0
    
    CURTIME = 0D0
    ! save original field to file
    CALL OUTPUTMESHFIELD(NETP,FIELDS,CURTIME,SNAPSHOTFILE,.FALSE.)

    NABS = 0
    DO NC = 1,NETP%NNODE
       IF (NETP%NODEABS(NC)) THEN
          NABS=NABS+1
          ABSORBERNODES(NABS) = NC
       ENDIF
    ENDDO

    !output flux to file
    OU=55
    OPEN(FILE=OUTFILE,UNIT=OU,STATUS='UNKNOWN')
    
    DO STEP = 1,NSTEP
       ! calculate time derivative
       CALL EULERSTEP_DIFFUSION(NETP,FIELDS(:,:,1),DFDT,FLUX)

       IF (MOD(STEP,OUTPUTEVERY).EQ.0) THEN
          ! output flux to file

          ! set of absorbing nodes
          WRITE(OU,*) CURTIME
          WRITE(OU,*) ABSORBERNODES(1:NABS)
          ! flux out of those nodes
          WRITE(OU,*) FLUX(ABSORBERNODES(1:NABS))
          CALL FLUSH(OU)
       END IF
       
       ! propagate forward in time
       FIELDS(:,:,1) = FIELDS(:,:,1)+DFDT*DELT
      
       CURTIME = CURTIME + DELT
       
       IF (MOD(STEP,SNAPSHOTEVERY).EQ.0) THEN
          ! save a snapshot of the field to file
          ! appending to end of file
          CALL OUTPUTMESHFIELD(NETP,FIELDS,CURTIME,SNAPSHOTFILE,.TRUE.)
       ENDIF
       
       IF (MOD(STEP,PRINTEVERY).EQ.0) THEN
           ! Calculate total integral of the field
          ! and total flux out of network
          CALL INTEGRATEMESHFIELD(NETP,FIELDS(:,:,1),TOTINT)
!         CALL INTEGRATEMESHFIELD(NETP,DFDT,TOTFLUX)
!         TOTFLUX = -TOTFLUX

          IF (MOD(STEP,PRINTEVERY).EQ.0) THEN
             PRINT*, 'STEP ', STEP, TOTINT
          ENDIF
!         IF (MOD(STEP,OUTPUTEVERY).EQ.0) THEN
!            WRITE(OU,*) STEP, CURTIME, TOTINT, TOTFLUX
!         ENDIF
          
       ENDIF
    ENDDO

    CLOSE(OU)
  END SUBROUTINE EULERDYNAMICS
