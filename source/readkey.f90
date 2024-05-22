SUBROUTINE READKEY
  ! this subroutine reads in keywords from a parameter file
  ! it sets the various global variables defined in KEYS module
  ! name of the parameter file is param.* where * is a keyword argument
  ! if no keyword argument is supplied, the default is just a file called param
  ! The EXTRAPARAMFILES keyword will allow extra parameter files to be 
  ! read in as well

  USE KEYS
  USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
  USE GENUTIL

  IMPLICIT NONE

  ! ---- stuff for inputing the parameter file in free format --------
  CHARACTER*100 :: ARG ! command line argument
  INTEGER :: NUMARG ! number of command line arguments
  INTEGER :: NITEMS ! number of items on the line in the parameter file
  INTEGER :: PF ! input file unit
  LOGICAL :: FILEEND=.FALSE. ! done reading file?
  CHARACTER*100 :: WORD ! keyword
  ! -------------- for reading multiple parameter files --------  
  INTEGER, PARAMETER :: MAXNFILES = 10
  CHARACTER*100 :: PARAMFILES(MAXNFILES)
  INTEGER :: NPARAMFILES, NPARAMREAD
  ! ------ for initializing random number generator
  INTEGER :: TIMEVAL(8), SEED
  ! ---------------- temporary variables ---------------
  INTEGER :: DUMI, I, TMPI, DUMI1, DUMI2, DUMI3, rc, FC
  CHARACTER*100 :: DUMSTR
  LOGICAL :: LDUM, CONTSEPSET
  DOUBLE PRECISION :: TMP
  LOGICAL :: DOACTIVATION

  ! ------------------------
  ! set variable defaults
  ! ------------------------
  ACTION = 'NONE'
  RNGSEED = 0
  VERBOSE = .FALSE.

  DOFLOW = .TRUE. ! include fluid flow
  
  ! input/output  
  OUTFILE = '*.out'
  DUMPSNAPSHOTS = .false. ! periodically dump chain snapshots
  SNAPSHOTEVERY = 1 ! how often to dump snapshots
  SNAPSHOTFILE = '*.snap.out' ! snapshot file
  APPENDSNAPSHOTS = .FALSE. ! append snapshots to file rather than replacing
  PRINTEVERY = 1000
  NETFILE = '*.net' ! File containing network structure
  NETOUTFILE='*.netout.out' ! File for outputting intermediate network structure
  CONTFILE = '*.cont.out' ! File containing info on contractions

  ! snapshots always include one at 0 time
  LOGSNAPSHOT = .FALSE. ! logarithmic spacing of snapshots?
  NSNAPSHOT = 0 ! number of snapshots to use with logarithmic spacing
  
  OUTPUTEVERY = 1D3 ! output flux every so many steps
  ! output flux more often at the start
  OUTPUTEVERYSTART = 1D1
  OUTPUTEVERYSWITCH = 0

  ! Only output 1 field to snapshots
  ! default (negative) means output all
  OUTPUT1FIELD = -1
  ! assume spatially uniform buffer
  UNIFORMBUFFER = .FALSE.
  
  ! network parameters
  MAXBRANCH = 12 ! max number of branches per node

  NSTEP = 1E6 ! Number of steps to propagate for
  DELT = 1D-4 ! time step for BD propagation
  DCOEFF = 1D0 ! diffusion coefficient

  ! contractions
  NCONT = 0 ! number of contractions in network
  CONTSEP = 0D0 ! minimum separation btwn contractions
  CONTSEPSET = .FALSE. ! separation explicitly set?
  CONTSPEED = 1D0 ! speed of contraction opening/closing
  CONTFLOW = 1D0 ! fluid flow associated with each contraction opening/closing
  OPENRATE = 1D0 ! rate at which fully closed contractions start to open
  FIXCONTRACTIONS = .FALSE.

  NODEVOL = 1D0
  
  ! Random runs
  STARTRATE = 0D0
  STOPRATE = 1D0
  RUNSPEED = 0D0
  SWITCHVELS = .FALSE.
  
  ! Mesh on network
  MAXMESHSIZE = 0.1D0 ! maximum size of mesh interval
  MINMESHPT = 2 ! minimum interior cells per edge
  MESHFILE = '*.mesh.txt' ! file for output of mesh info

  ! list of absorber nodes (ie: calcium exit sites)
  ABSORBERS = 0
  NABS = 0 ! number of absorber nodes
  NFIX = 0 ! number of fixed nodes
  NFIXCELL = 0 ! number of fixed cells
  NFIXPT = 0 ! Number of fixed points
  FIXNODEFROMNETFILE = .false. ! determine fixed nodes based on network file
  ! if fixing nodes from netfile, pick random subset of the ones marked
  ! to actually get fixed
  ! default to none
  FIXSUBSETNODES = -1
  
  ! Where to pick points in a circle
  FIXPTCENT = 0D0; FIXPTRAD = 1D0

  ! maximum allowed distance from picked point to a node
  FIXPTMAXDIST = 1D5
  
  ! randomly pick some number of network nodes to fix
  RANDFIXNODES = .FALSE.
  ! randomly pick reservoirs to fix
  RANDFIXRESV = .FALSE.
  ! randomly pick some number of mesh cells to fix
  RANDFIXCELLS = .FALSE.
  ! Randomly pick points in a circle and fix mesh cells nearest them
  RANDFIXPTS = .FALSE.

  ! By default do not close any mesh boundaries
  ! probability of closing some edge boundaries permanently  
  PCLOSEBOUNDPERLEN = -1D0
  DOCLOSEBOUNDONCE = .FALSE.
  ! rates of opening and closing edge boundaries  
  CLOSEBOUNDRATEPERLEN = -1D0
  OPENBOUNDRATE = -1D0
  DOCLOSEBOUNDRATE = .FALSE.
  
  ! edges on which the field starts (overwrites nodes)
  NSTARTEDGE = 0
  STARTEDGES = 0
  ! nodes on which field starts 
  NSTARTNODE = 0
  STARTNODES = 0
  NSTARTPOS = 0
  STARTPOS = 0D0
  ! radius around nodes for starting bolus
  STARTNODERAD = -1
  ! starting concentration
  STARTCONC = 1D0
  ! Start with equilibrated concentrations?
  STARTEQUIL = 0
  ! Background concentration specified?
  SETBACKGROUNDCONC = .FALSE.
  BACKGROUNDCONC = 0D0
  ! start from permeable nodes?
  STARTFROMPERMEABLE = .FALSE.
  
  ! reservoir volumes and surface areas
  RESVVOL= 1D0
  RESVSA = 1D0
  RESVLEN = 1D0
  ! are reservoirs well-mixed?
  RESVMIX = .TRUE.
  DORESERVOIRS = .FALSE.
  ! file to output reservoir data
  RESVFILE = '*.resv.out'
  ! allow fixing of nodes attached to reservoirs
  ALLOWFIXEDRESV = .FALSE.
  ! some specific reservoirs are not allowed to be fixed
  ! this setting is temporary
  ALLOWRESVFIX = .TRUE.
  
  ! permeability of permeable nodes (units of length per time, unless usepermprefactor is true)
  ! if usepermprefactor: this parameter is multiplied
  ! by mesh cell length and diffusivity to give your permeability
  PERMEABILITY =0D0
  NPERM = 0 ! number of permeable nodes

  ! periodically turn on global permeability
  ! this sets the period (do not do this if negative)
  PERIODGLOBALPERM = -1D0
  ! this sets the duration of the permeable stretch
  DURGLOBALPERM = 0D0
  
  NPERMPOS = 0 ! number of permeable center positions
  ! permeability around each position
  POSPERMEABILITY = 0D0
  ! position of each permeability center
  PERMPOS = 0D0
  
  RANDPERMNODES = .FALSE. 
  ! external concentration (*pi*a^2 to get units of per length)
  CEXT = 0D0
  ! permeability info read in from network file
  PERMFROMFILE = .FALSE.
  TRACKFLUXPERM = .TRUE.
  ! make all cells within a certain radius of a permeable node also permeable
  ! Also used to make permeable all cells within a radius of a permeable position
  PERMNEARNODEDIST = -1D0
  ! permeability supplied is actually prefactor to be multiplied by cell length and diffusivity
  USEPERMPREFACTOR = .TRUE.
  
  ! number of fields
  NFIELD = 1
  MOBILEFIELD = .TRUE.

  ! evolve buffer protein concentations
  DOBUFFER = .FALSE.
  KON = 100D0
  KOFF = 100D0
  KDEQUIL = -1D0
  FASTEQUIL = .FALSE. ! rapid buffer equilibration

  ! warn when a field hits zero?
  WARNZEROFIELD=.FALSE.

  ! output flux out of fixed nodes rather than absorber nodes
  TRACKFIXNODEFLUX = .FALSE.

  ! fix all cells near the fixed nodes
  ! default is to fix the specified node cell only
  FIXNEARNODEDIST = -1

  ! Fix mesh cells within a rectangle
  FIXRECTANGLE = -1D0
  
  ! include velocities in snapshots?
  SNAPSHOTVEL = 0

  ! How to control velocities
  VELCONTROL = 'RANDVEL'

  ! Fixing edge velocities
  NFIXEDGEVEL = 0
  FIXEDGEVEL = 0
  FIXEDGEVELVAL = 0D0

  ! Activation rate near certain nodes
  DOACTIVATION = .FALSE.
  ACTRATE = 0D0
  ACTNEARNODEDIST = 1D0
  NACT = 0
  ACTNODES = 0
  ! global depletion rate
  DEPRATE = 0D0

  ! allow for variable radii of mesh elements
  USEVARRAD = .FALSE.
  ! parameters for randomizing edge radii
  EDGERADRANDTYPE = 'NONE'
  EDGERADRANDPARAMS = 0
  EDGERADBASE = SQRT(1/PI) ! default edge radius
  EDGERADVARTYPE = 'CONSTANT'
  EDGERADVARPARAMS = 0
  
  ! if positive, predefines a network dimension
  ! otherwise, dimension set to # items in NODE row of network file - 2
  NETWORKDIM = 0

  ! only output the total flux (not individual nodes)
  OUTPUTTOTFLUXONLY = .TRUE.

  ! input reservoir elements from file
  USERESVELEMENTS = .FALSE.
  RESVELEMENTFILE = '*.resv.txt'

  ! concentrations are expressed as 3D rather than 1D quantities
  CONCENTRATIONS3D = .FALSE.

  ! global reservoir (ie: cytoplasm) with pumping kinetics
  USEGLOBALRESV = .FALSE.
  ! Permeable cells release to global reservoir?
  PERMTOGLOBALRESV = .FALSE.
  ! volume of global reservoir (scaled by pi*a^2 if working with 1D concs)
  GLOBALRESVOL = 1D0
  ! rate constants and saturation concentrations for
  ! recovery from global reservoir
  GLOBALRESVKR = 1D0
  GLOBALRESVKMR = 1D0
  ! and pumping out of global reservoir
  GLOBALRESVKOUT = 0D0
  GLOBALRESVKMOUT = 1D0
  ! initial concentration for global reservoir
  GLOBALRESVSTART = 0D0

  ! track changes in total calcium directly
  TRACKDCDT = .FALSE.
  
  ! -------------------------
  ! Read in all parameter files, starting with the ones specified on command line
  ! --------------------------

  PF = 55 ! i/o unit number to be used for parameter files

  ! get input parameter files from command line
  NPARAMFILES = 0
  NUMARG = COMMAND_ARGUMENT_COUNT()  
  IF (NUMARG==0) THEN
     NPARAMFILES = 1
     PARAMFILES(1) = 'param'
     ARG = ''
  ELSE
     DO I = 1,NUMARG
        CALL GETARG(I, ARG)
        NPARAMFILES = NPARAMFILES + 1
        WRITE(DUMSTR,'(A)') 'param.' //TRIM(ADJUSTL(ARG))!//'.txt'
        PARAMFILES(NPARAMFILES) = DUMSTR
     ENDDO
     ! reset arg to its original value
     IF (NUMARG.GT.1) CALL GETARG(1,ARG)
  ENDIF

  NPARAMREAD = 0 ! keep track of how many files have been read
  DO WHILE (NPARAMREAD.LT.NPARAMFILES)
     NPARAMREAD = NPARAMREAD + 1

     PRINT*, 'Reading parameter file: ', PARAMFILES(NPARAMREAD)
     INQUIRE(FILE=PARAMFILES(NPARAMREAD),EXIST=LDUM)
     IF (.NOT.LDUM) THEN
        PRINT*, 'ERROR in READKEY: Parameter file ', TRIM(ADJUSTL(PARAMFILES(NPARAMREAD))), ' does not exist.'
        STOP 1
     ENDIF
     OPEN(UNIT=PF, FILE=PARAMFILES(NPARAMREAD), STATUS='OLD')

     ! read in the keywords one line at a time
     DO 
        CALL READLINE(PF,FILEEND,NITEMS)
        IF (FILEEND.and.nitems.eq.0) EXIT

        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE

        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)

        ! Skip any empty lines or any comment lines
        IF (WORD(1:1).EQ.'#') CYCLE

        SELECT CASE(WORD) ! pick which keyword
        CASE('ABSORBERS')
           print*, 'ABSORBERS keyword is currently buggy. Do not use.'
           STOP 1
           
           CALL READI(FC) ! which field is this for
           IF (FC.GT.MAXNFIELD) THEN
              PRINT*, 'ERROR: absorbers has an invalid field number', FC
              STOP 1
           ENDIF
           DO I = 1,NITEMS-2
              NABS(FC) = NABS(FC) + 1
              IF (NABS(FC).GT.MAXNABSORBER) THEN
                 PRINT*, 'ERROR IN READKEY: too many absorber nodes', FC, NABS(FC), MAXNABSORBER
                 STOP 1
              ENDIF
              CALL READI(ABSORBERS(NABS(FC),FC))
           ENDDO        
        CASE('ACTION')
           CALL READA(ACTION, CASESET=1)
        CASE('ACTNODE')
           DOACTIVATION = .TRUE.
           NACT = NACT + 1
           IF (NACT.GT.MAXNACT) THEN
              PRINT*, 'ERROR IN READKEY: too many activation nodes', NACT, MAXNACT
              STOP 1
           ENDIF
           CALL READI(ACTNODES(NACT))
           CALL READF(ACTRATE(NACT))
           CALL READF(ACTNEARNODEDIST(nact))
        CASE('ALLOWFIXEDRESV')
           IF (NITEMS.GT.1) THEN
              CALL READO(ALLOWFIXEDRESV)
           ELSE
              ALLOWFIXEDRESV = .TRUE.
           ENDIF
        CASE('BACKGROUNDCONC')
           SETBACKGROUNDCONC = .TRUE.
           DO I = 1,NITEMS-1
              CALL READF(BACKGROUNDCONC(I))
           ENDDO          
        CASE('CEXT')
           DO I = 1,MIN(MAXNFIELD,NITEMS-1)
              CALL READF(CEXT(I))
           ENDDO
        CASE('CLOSEBOUNDONCE')
           ! close some set of boundaries once at start of sim
           CALL READF(PCLOSEBOUNDPERLEN)
           IF (PCLOSEBOUNDPERLEN.GE.0) DOCLOSEBOUNDONCE = .TRUE.
        CASE('CLOSEBOUNDRATE')
           ! close and reopen boundaries at some rate
           CALL READF(CLOSEBOUNDRATEPERLEN)
           CALL READF(OPENBOUNDRATE)
           IF (CLOSEBOUNDRATEPERLEN.GE.0.AND.OPENBOUNDRATE.GE.0) &
                & DOCLOSEBOUNDRATE = .TRUE.
        CASE('CONCENTRATIONS3D')
           IF (NITEMS.GT.1) THEN
              CALL READO(CONCENTRATIONS3D)
           ELSE
              CONCENTRATIONS3D = .TRUE.
           ENDIF
        CASE('CONTFILE')
           CALL READA(CONTFILE)
        CASE('CONTSEP')
           CALL READF(CONTSEP)
           CONTSEPSET = .TRUE.
        CASE('CONTSPEED')
           CALL READF(CONTSPEED)
           CALL READF(CONTFLOW)
        CASE('DEPRATE')
           CALL READI(FC) ! which field is this for
           CALL READF(DEPRATE)
        CASE('DOBUFFER')
           IF (NITEMS.GT.1) THEN
              CALL READO(DOBUFFER)
           ELSE
              DOBUFFER = .TRUE.
           ENDIF
        CASE('DOFLOW')
           IF (NITEMS.GT.1) THEN
              CALL READO(DOFLOW)
           ELSE
              DOFLOW = .TRUE.
           ENDIF
        CASE('DORESERVOIRS')
           IF (NITEMS.GT.1) THEN
              CALL READO(DORESERVOIRS)
           ELSE
              DORESERVOIRS = .TRUE.
           ENDIF           
        CASE('OPENRATE')
           CALL READF(OPENRATE)
        CASE('DCOEFF')
           DO FC = 1,MIN(NITEMS-1,MAXNFIELD)
              CALL READF(DCOEFF(FC))
           ENDDO
        CASE('DELT')
           CALL READF(DELT)
        CASE('EDGERADBASE') 
           CALL READF(EDGERADBASE)
        CASE('EDGERADRAND')
           CALL READA(EDGERADRANDTYPE,1)
           DO I = 1,MIN(NITEMS-2,10)
              CALL READF(EDGERADRANDPARAMS(I))
           ENDDO
        CASE('EDGERADVAR')
           CALL READA(EDGERADVARTYPE,1)
           DO I = 1,MIN(NITEMS-2,10)
              CALL READF(EDGERADVARPARAMS(I))
           ENDDO
        CASE('FASTEQUIL')
           IF (NITEMS.GT.1) THEN
              CALL READO(FASTEQUIL)
           ELSE
              FASTEQUIL = .TRUE.
           ENDIF         
        CASE('FIXCONTRACTIONS')
           IF (NITEMS.GT.1) THEN
              CALL READO(FIXCONTRACTIONS)
           ELSE
              FIXCONTRACTIONS = .TRUE.
           ENDIF
        CASE('FIXEDGEVEL')
           NFIXEDGEVEL = NFIXEDGEVEL + 1
           CALL READI(FIXEDGEVEL(NFIXEDGEVEL)) ! which edge to fix
           CALL READF(FIXEDGEVELVAL(NFIXEDGEVEL)) ! value to fix velocity
        CASE('FIXNEARNODEDIST')
           CALL READF(FIXNEARNODEDIST)
        CASE('FIXNODE')
           CALL READI(FC) ! which field is this for
           NFIX(FC) = NFIX(FC) + 1
           IF (NFIX(FC).GT.MAXNABSORBER) THEN
              PRINT*, 'ERROR IN READKEY: too many fixed nodes', NFIX, MAXNABSORBER
              STOP 1
           ENDIF
           CALL READI(FIXNODES(NFIX(FC),FC)) ! which node is being fixed
           CALL READF(FIXVALS(NFIX(FC),FC)) ! what value is it fixed to
        CASE('FIXNODEFROMNETFILE')
           IF (NITEMS.GT.1) THEN
              CALL READO(FIXNODEFROMNETFILE)
           ELSE
              FIXNODEFROMNETFILE = .TRUE.
           ENDIF
        CASE('FIXRECTANGLE')
           ! get field, then minx miny, dx, dy, value
           CALL READI(FC)
           IF (FC.GT.MAXNFIELD) THEN
              PRINT*, 'ERROR: impossible field for fixrectangle', FC
              STOP 1
           ENDIF
           DO I = 1,5
              CALL READF(FIXRECTANGLE(FC,I))
           ENDDO
        CASE('FIXSUBSETNODES')
           ! if fixing nodes from the network file, pick a random subset of
           ! the nodes labeled FN to actually get fixed
           FIXNODEFROMNETFILE = .TRUE.
           CALL READI(FIXSUBSETNODES)
        CASE('GLOBALRESERVOIR')
           ! recovery of particles from global reservoir
           USEGLOBALRESV = .TRUE.
           ! volume and kinetic parameters for pumping from global reservoir
           CALL READF(GLOBALRESVOL)
           CALL READF(GLOBALRESVKR)
           CALL READF(GLOBALRESVKMR)
           IF (NITEMS.GT.4) THEN
              CALL READF(GLOBALRESVKOUT)
              CALL READF(GLOBALRESVKMOUT)
           ENDIF
           IF (NITEMS.GT.6) THEN
              CALL READO(PERMTOGLOBALRESV)
           END IF
        CASE('GLOBALRESERVOIRSTART')
           CALL READF(GLOBALRESVSTART)
        CASE('KDEQUIL')
           CALL READF(KDEQUIL)
        CASE('KOFF')
           CALL READF(KOFF)
        CASE('KON')
           CALL READF(KON)
        CASE('MAXBRANCH')
           CALL READI(MAXBRANCH)
        CASE('MESHFILE')           
           CALL READA(MESHFILE)
        CASE('MESHSIZE')
           CALL READF(MAXMESHSIZE)
           IF (NITEMS.GT.2) CALL READI(MINMESHPT)
        CASE('MOBILEFIELD')
           DO FC = 1,MIN(NITEMS-1,MAXNFIELD)
              CALL READO(MOBILEFIELD(FC))
           ENDDO
        CASE('NETFILE')
           CALL READA(NETFILE)
           IF (NITEMS.GT.2) CALL READA(NETOUTFILE)
        CASE('NETWORKDIM')
           CALL READI(NETWORKDIM)
        CASE('NCONT')
           CALL READI(NCONT)
        CASE('NOFIXRESV')
           DO FC = 1,NITEMS-1
              CALL READI(DUMI)
              IF (DUMI<MAXNRESV) ALLOWRESVFIX(DUMI) = .FALSE.              
           ENDDO    
        CASE('NODEVOL')
           CALL READF(NODEVOL)
        CASE('NSTEP')
           CALL READI(NSTEP)
        CASE('OUTFILE')
           CALL READA(OUTFILE)
           IF (NITEMS.GT.2) CALL READI(OUTPUTEVERY)
        CASE('OUTPUT1FIELD')
           CALL READI(OUTPUT1FIELD)
        CASE('OUTPUTEVERY')
           CALL READI(OUTPUTEVERY)
        CASE('OUTPUTEVERYSTART')
           CALL READI(OUTPUTEVERYSWITCH)
           CALL READI(OUTPUTEVERYSTART)
        CASE('OUTPUTTOTFLUXONLY')
           IF (NITEMS.EQ.1) THEN
              OUTPUTTOTFLUXONLY = .TRUE.
           ELSE
              CALL READO(OUTPUTTOTFLUXONLY)
           ENDIF        
        CASE('PERIODICGLOBALPERM')
           CALL READF(PERIODGLOBALPERM)
           CALL READF(DURGLOBALPERM)
           DO I = 1,min(NITEMS-3,MAXNFIELD)
              ! permeability for each field
              CALL READF(PERMEABILITY(1,I))
           ENDDO           
        CASE('PERMFROMFILE')
           IF (NITEMS.GT.1) THEN
              CALL READO(PERMFROMFILE)
           ELSE
              PERMFROMFILE = .TRUE.
           ENDIF
        CASE('PERMNEARNODEDIST')
           CALL READF(PERMNEARNODEDIST)
        CASE('PERMNODE')
           NPERM = NPERM +1
           IF (NPERM.GT.MAXNABSORBER) THEN
              PRINT*, 'ERROR IN PERMNODE:', NPERM, MAXNABSORBER
              STOP 1
           ENDIF
           CALL READI(PERMNODES(NPERM)) ! which node is permeable
           DO I = 1,min(NITEMS-2,MAXNFIELD)
              ! permeability for each field
              CALL READF(PERMEABILITY(NPERM,I))
           ENDDO
        CASE('PERMPOS')
           NPERMPOS = NPERMPOS +1
           IF (NPERMPOS.GT.MAXNABSORBER) THEN
              PRINT*, 'ERROR IN PERMPOS:', NPERMPOS, MAXNABSORBER
              STOP 1
           ENDIF
           ! read in coordinates
           ! WARNING: nondefault NETWORKDIM needs to be specified before this line
           DO I = 1,NETWORKDIM
              CALL READF(PERMPOS(NPERMPOS,I))
           ENDDO           
           DO I = 1,min(NITEMS-3,MAXNFIELD)
              ! permeability for each field
              CALL READF(POSPERMEABILITY(NPERMPOS,I))
           ENDDO
        CASE('PRINTEVERY')
           CALL READI(PRINTEVERY)
        CASE('RANDFIXCELLS')
           RANDFIXCELLS = .TRUE.
           CALL READI(FC) ! which field is this for
           CALL READI(NFIXCELL(FC)) ! how many cells to fix?
           CALL READF(TMP) ! What value to fix to
           FIXVALS(:,FC) = TMP
        CASE('RANDFIXNODES')
           RANDFIXNODES = .TRUE.
           RANDFIXRESV = .FALSE.
           CALL READI(FC) ! which field is this for
           CALL READI(NFIX(FC)) ! how many nodes to fix?
           CALL READF(TMP) ! What value to fix to
           FIXVALS(:,FC) = TMP
        CASE('RANDFIXPTS')
           RANDFIXPTS = .TRUE.
           CALL READI(FC) ! which field is this for
           CALL READI(NFIXPT(FC)) ! how many points to fix
           CALL READF(TMP) ! what value to fix to
           FIXVALS(:,FC) = TMP
           CALL READF(FIXPTRAD) ! radius of circle in which to pick points
           IF (NITEMS.GT.5) THEN ! center of circle
              CALL READF(FIXPTCENT(1)); CALL READF(FIXPTCENT(2))
           ENDIF
           IF (NITEMS.GT.7) THEN
              ! if selected point is more than this distance from a fixable
              ! mesh cell, then ignore and pick another point
              CALL READF(FIXPTMAXDIST)
           ENDIF
        CASE('RANDFIXPTSEXCLUDE')
           ! exclude a circular region from randomly selected fixed points
           CALL READF(FIXPTEXRAD)
           DO I = 1,MIN(NITEMS-2,3)
              CALL READF(FIXPTEXCENT(I))
           ENDDO
        CASE('RANDFIXRESV')
           RANDFIXRESV = .TRUE.
           RANDFIXNODES = .FALSE.
           CALL READI(FC) ! which field is this for
           CALL READI(NFIX(FC)) ! how many reservoirs to fix?
           CALL READF(TMP) ! What value to fix to
           FIXVALS(:,FC) = TMP
        CASE('RANDPERMNODES')
           RANDPERMNODES = .TRUE.
           CALL READI(NPERM)
           IF (NPERM.GT.MAXNABSORBER) THEN
              PRINT*, 'ERROR IN RANDPERMNODES:', NPERM, MAXNABSORBER
              STOP 1
           ENDIF
           DO I = 1,NITEMS-2
              CALL READF(PERMEABILITY(1,I))
           ENDDO
        CASE('RESVELEMENTS')
           USERESVELEMENTS = .TRUE.
           CALL READA(RESVELEMENTFILE)
        CASE('RESVFILE')
           CALL READA(RESVFILE)
        CASE('RESVINFO')
           CALL READI(RC)
           IF (RC.GT.MAXNRESV) THEN
              PRINT*, 'ERROR: reservoir index is too large in readkey. Skipping.'
              STOP 1
           ELSEIF (RC.LE.0) THEN
              ! assume all reservoirs have same parameters
              CALL READF(RESVVOL(1))
              CALL READF(RESVSA(1))
              CALL READF(RESVLEN(1))
              IF (NITEMS.GT.5) CALL READO(RESVMIX(1))
              RESVVOL = RESVVOL(1)
              RESVSA = RESVSA(1)
              RESVLEN = RESVLEN(1)
              RESVMIX = RESVMIX(1)
           ELSE
              CALL READF(RESVVOL(RC))
              CALL READF(RESVSA(RC))
              CALL READF(RESVLEN(RC))
              IF (NITEMS.GT.5) CALL READO(RESVMIX(RC))
           ENDIF            
        CASE('RNGSEED')
           CALL READI(RNGSEED)
        CASE('RUNSPEED')
           CALL READF(RUNSPEED)       
        CASE('SNAPSHOTFILE')
           CALL READA(SNAPSHOTFILE)
        CASE('SNAPSHOTLOG')
           ! logarithmic spacing of snapshots
           DUMPSNAPSHOTS = .TRUE.
           LOGSNAPSHOT = .TRUE.
           CALL READI(NSNAPSHOT)
        CASE('SNAPSHOTVEL')
           IF (NITEMS.GT.1) THEN
              CALL READI(SNAPSHOTVEL)
           ELSE
              SNAPSHOTVEL = 1
           ENDIF
        CASE('SNAPSHOTS')
           DUMPSNAPSHOTS = .TRUE.
           IF (NITEMS.GT.1) CALL READI(SNAPSHOTEVERY)
           IF (NITEMS.GT.2) CALL READA(SNAPSHOTFILE)
           IF (NITEMS.GT.3) CALL READO(APPENDSNAPSHOTS)        
        CASE('STARTEDGES')
           DO I = 1,NITEMS-1
              NSTARTEDGE = NSTARTEDGE + 1
              IF (NSTARTEDGE.GT.MAXSTARTEDGE) THEN
                 PRINT*, 'ERROR IN READKEY: too many startedge', NSTARTEDGE, MAXSTARTEDGE
                 STOP 1
              ENDIF
              CALL READI(STARTEDGES(NSTARTEDGE))
           ENDDO
        CASE('STARTCONC')
           DO FC = 1,MIN(NITEMS-1,MAXNFIELD)
              CALL READF(STARTCONC(FC))
           ENDDO
        CASE('STARTFROMPERMEABLE')
           IF (NITEMS.GT.1) THEN
              CALL READO(STARTFROMPERMEABLE)
           ELSE           
              STARTFROMPERMEABLE = .TRUE.
           ENDIF
        CASE('STARTEQUIL')
           CALL READI(STARTEQUIL)
        CASE('STARTNODERAD')
           CALL READF(STARTNODERAD)
        CASE('STARTNODES')
           DO I = 1,NITEMS-1
              NSTARTNODE = NSTARTNODE + 1
              IF (NSTARTNODE.GT.MAXSTARTNODE) THEN
                 PRINT*, 'ERROR IN READKEY: too many startnode', NSTARTNODE, MAXSTARTNODE
                 STOP 1
              ENDIF
              CALL READI(STARTNODES(NSTARTNODE))
           ENDDO
        CASE('STARTPOS')
           NSTARTPOS = NSTARTPOS+1
           IF (NSTARTPOS.GT.MAXSTARTNODE) THEN
              PRINT*, 'ERROR IN READKEY: too many startPOS', NSTARTPOS, MAXSTARTNODE
              STOP 1
           ENDIF
           DO I = 1,NITEMS-1
              CALL READF(STARTPOS(NSTARTPOS,I))
           ENDDO
        CASE('STARTRATE')
           CALL READF(STARTRATE)
        CASE('STOPRATE')
           CALL READF(STOPRATE)
        CASE('SWITCHVELS')
           IF (NITEMS.GT.1) THEN
              CALL READO(SWITCHVELS)
           ELSE
              SWITCHVELS = .TRUE.
           ENDIF
        CASE('TRACKDCDT')           
           CALL READO(TRACKDCDT)
        CASE('TRACKFLUXPERM')
           IF (NITEMS.GT.1) THEN
              CALL READO(TRACKFLUXPERM)
           ELSE
              TRACKFLUXPERM = .TRUE.
           ENDIF
        CASE('UNIFORMBUFFER')
           IF (NITEMS.GT.1) THEN
              CALL READO(UNIFORMBUFFER)
           ELSE
              UNIFORMBUFFER = .TRUE.
           ENDIF        
        CASE('USEPERMPREFACTOR')
           IF (NITEMS.GT.1) THEN
              CALL READO(USEPERMPREFACTOR)
           ELSE
              USEPERMPREFACTOR = .TRUE.
           ENDIF
        CASE('USEVARRAD')
           IF (NITEMS.GT.1) THEN
              CALL READO(USEVARRAD)
           ELSE
              USEVARRAD = .TRUE.
           ENDIF
        CASE('VELCONTROL')
           CALL READA(VELCONTROL, CASESET=1)
        CASE('VERBOSE')
           IF (NITEMS.GT.1) THEN
              CALL READO(VERBOSE)
           ELSE
              VERBOSE =.TRUE.
           ENDIF
        CASE('WARNZEROFIELD')
           IF (NITEMS.GT.1) THEN
              CALL READO(WARNZEROFIELD)
           ELSE
              WARNZEROFIELD = .TRUE.
           ENDIF
        CASE DEFAULT
           print*, 'ERROR: unidentified keyword ', TRIM(WORD), " Will ignore."
        END SELECT
     ENDDO
     CLOSE(PF)
  ENDDO

  ! do not allow fixing of any reservoirs
  IF (.NOT.ALLOWFIXEDRESV) ALLOWRESVFIX = .FALSE.

  ! -----------------
  ! check validity of some values, raise errors or adjust as necessary
  ! -----------------  

  IF (DOBUFFER) THEN
     IF (FASTEQUIL) THEN
        NFIELD = 2 ! fields: free ligand, total protein
     ELSE
        NFIELD = 3 ! fields: free ligand, bound ligand, free protein
     ENDIF
  ELSE
     NFIELD = 1
  ENDIF

  IF (DOACTIVATION) THEN
     IF (DOBUFFER) THEN
        PRINT*, 'Activation and buffer currently not implemented together.'
        STOP 1
     ENDIF
     NFIELD = 2
  ENDIF
  
  ! ----------- fix file names -----------
  CALL REPLACESUBSTR(OUTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(SNAPSHOTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(NETFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(CONTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(MESHFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(NETOUTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(resvFILE,'*',TRIM(ADJUSTL(ARG)))
  ! ---------------------------

  ! default minimal separation is mesh size
  IF (.NOT.CONTSEPSET) CONTSEP = MAXMESHSIZE
  
  ! Initiate random number generator 
  IF (RNGSEED.EQ.0) THEN
     ! use the current time of day in milliseconds
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = TIMEVAL(5)*3600*1000 + TIMEVAL(6)*60*1000 + TIMEVAL(7)*1000 + TIMEVAL(8)
  ELSEIF (RNGSEED.EQ.-1) THEN
     ! use the last 5 characters in the command-line argument
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)))    
  ELSEIF (RNGSEED.EQ.-2) THEN
     ! use the last 4 characters in the command-line argument 
     ! and additionally the millisecond time 
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)),TIMEVAL(8))
  ELSE
     ! use this seed directly
     SEED = RNGSEED
  ENDIF

  print*, 'Initiating Mersenne twister random number generator with seed:', SEED
  CALL SGRND(SEED)

  print*, '------------Parameter values : -------------------'
  print*, 'ACTION: ', TRIM(ADJUSTL(ACTION))
  print*, 'Output file: ', TRIM(OUTFILE)
  print*, 'Network file: ', TRIM(NETFILE)
  PRINT*, 'Starting edges:', STARTEDGES(1:NSTARTEDGE)
  PRINT*, 'Start on nodes:', STARTNODES(1:NSTARTNODE)
  DO FC = 1,NFIELD
     print*, 'Absorber nodes:', FC, ABSORBERS(1:NABS(FC),FC)
     if (.not.RANDFIXNODES.AND..NOT.FIXNODEFROMNETFILE) PRINT*, 'Fixed nodes:', FC, nfix(fc), FIXNODES(1:NFIX(FC),FC)
  ENDDO
  IF (DOACTIVATION) THEN
     DO I = 1,NACT
        PRINT*, 'Activation node:', I, ACTNODES(I), ACTRATE(I), ACTNEARNODEDIST(I)
     ENDDO
  ENDIF

  PRINT*, 'Buffer proteins? Reservoirs?:', DOBUFFER, DORESERVOIRS
  IF (DOBUFFER) THEN
     IF (FASTEQUIL) THEN
        PRINT*, 'Dissociation constant:', KDEQUIL
     ELSE
        PRINT*, 'Kon, Koff:', KON, KOFF
     ENDIF     
  ENDIF

  print*, 'Work with 3D concentrations? Or meshed reservoirs?:', CONCENTRATIONS3D, USERESVELEMENTS
  
  print*, 'Flows along edges?:', DOFLOW, RUNSPEED  
  
 !print*, 'CONTSPEED, CONTFLOW, OPENRATE:', CONTSPEED, CONTFLOW, OPENRATE
  print*, 'DCOEFF:', dcoeff
  print*, 'MOBILEFIELD:', MOBILEFIELD
  IF (DUMPSNAPSHOTS) THEN
     IF (LOGSNAPSHOT) THEN
        PRINT*, 'Dumping snapshots logarithmically. Total number ', NSNAPSHOT,'.  In file: ', TRIM(ADJUSTL(SNAPSHOTFILE))
     ELSE
        PRINT*, 'Dumping snapshot every', SNAPSHOTEVERY,'steps. In file:', TRIM(ADJUSTL(SNAPSHOTFILE))
     ENDIF
  ENDIF
  IF (USEVARRAD) THEN
     PRINT*, 'Allowing for variable mesh cell radii'
  ENDIF
     PRINT*, 'EDGERADRAND: ', EDGERADRANDTYPE, EDGERADRANDPARAMS
     PRINT*, 'EDGERADVAR: ', EDGERADVARTYPE, EDGERADVARPARAMS
  print*, '----------------------------------------------------'

END SUBROUTINE READKEY
