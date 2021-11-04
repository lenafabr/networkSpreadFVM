MODULE KEYS
  ! keyword parameters that are globally used in many different places in the code
  IMPLICIT NONE

   ! -------- General program control ---------------
  CHARACTER*100 :: ACTION
  INTEGER :: RNGSEED
  LOGICAL :: VERBOSE  
  
  ! ----------------------
  ! Output / input
  ! -----------------------
  CHARACTER*100 :: OUTFILE, SNAPSHOTFILE, NETFILE, CONTFILE, NETOUTFILE
  LOGICAL :: DUMPSNAPSHOTS, RESTART, APPENDSNAPSHOTS
  INTEGER :: SNAPSHOTEVERY
  INTEGER :: PRINTEVERY, OUTPUTEVERY, OUTPUTEVERYSWITCH, OUTPUTEVERYSTART
  ! logarithmic spacing of snapshots?
  LOGICAL :: LOGSNAPSHOT
  INTEGER :: NSNAPSHOT ! how many snapshots to use with logarithmic spacing
  
  ! ------------
  ! network geometry and setup
  ! ------------
  INTEGER :: NETWORKDIM ! predefined network dimension
  INTEGER :: MAXBRANCH ! max number of branches (per node) that can be input in the network
  INTEGER, PARAMETER :: MAXSTARTEDGE = 10000
  INTEGER :: NSTARTEDGE, STARTEDGES(MAXSTARTEDGE) ! starting edges for propagating field
  INTEGER, PARAMETER :: MAXSTARTnode= 10000
  ! Maximum cycle length
  INTEGER, PARAMETER :: MAXLOOPLEN = 100
  INTEGER :: NSTARTNODE, STARTnodes(MAXSTARTNODE) ! starting nodes for propagating field
  INTEGER :: NSTARTPOS
  DOUBLE PRECISION :: STARTPOS(MAXSTARTNODE,3)
  ! radius for starting around specific nodes
  DOUBLE PRECISION :: STARTNODERAD
  
  INTEGER, PARAMETER :: MAXNFIELD = 5 ! maximum allowed concentration fields
  LOGICAL :: MOBILEFIELD(MAXNFIELD)
  
  ! -----------------
  ! solving diffusion equations and propagating particles
  ! ---------------- 
  INTEGER :: NPART ! number of particles to propagate
  INTEGER :: NSTEP ! number of steps to propagate for
  DOUBLE PRECISION :: DELT ! time-step for BD sims
  DOUBLE PRECISION :: DCOEFF(MAXNFIELD) ! sets time units

  ! ------------------
  ! edge contractions
  ! ------------------
  INTEGER :: NCONT ! number of contractions existing at a given time
  DOUBLE PRECISION :: CONTSEP ! minimal separation btwn contractions
  DOUBLE PRECISION :: CONTSPEED ! speed of contraction open/closing
  DOUBLE PRECISION :: CONTFLOW ! fluid flow from each active contraction
  DOUBLE PRECISION :: OPENRATE ! rate of contraction beginning to open
  LOGICAL :: FIXCONTRACTIONS ! do not allow the contractions to jump
  DOUBLE PRECISION :: NODEVOL ! default node volume

  CHARACTER(LEN=100) :: VELCONTROL ! how are velocities controlled?
  LOGICAL :: DOFLOW ! include flow?
  
  ! Random runs
  DOUBLE PRECISION :: STARTRATE, STOPRATE ! rate at which to start new runs, and stop them
  DOUBLE PRECISION :: RUNSPEED ! speed of a run
  ! instantaneously switch velocity direction rather than stopping runs
  LOGICAL :: SWITCHVELS
  
  ! ----------
  ! meshing on edges
  ! ----------
  DOUBLE PRECISION :: MAXMESHSIZE ! maximum allowed mesh interval
  INTEGER :: MINMESHPT ! minimum mesh points on each edge
  CHARACTER(LEN=100) :: MESHFILE ! file for output of mesh info
  
  ! --------------
  ! dynamics of fields over mesh
  ! --------------
  INTEGER, PARAMETER :: MAXNABSORBER = 100  
  INTEGER :: NABS(MAXNFIELD), NFIX(MAXNFIELD), NFIELD
  INTEGER :: ABSORBERS(MAXNABSORBER,MAXNFIELD), FIXNODES(MAXNABSORBER,MAXNFIELD)
  DOUBLE PRECISION :: FIXVALS(MAXNABSORBER,MAXNFIELD)
  DOUBLE PRECISION :: FIXNEARNODEDIST
  DOUBLE PRECISION :: FIXRECTANGLE(MAXNFIELD,5)
  LOGICAL :: RANDFIXNODES
  LOGICAL :: TRACKFIXNODEFLUX ! track flux out of fixed nodes rather than absorber nodes
  ! determine fixed nodes from the network file
  LOGICAL :: FIXNODEFROMNETFILE

  ! Activation and depletion
  INTEGER, PARAMETER :: MAXNACT = 100
  DOUBLE PRECISION :: ACTRATE(MAXNACT), DEPRATE
  DOUBLE PRECISION :: ACTNEARNODEDIST(MAXNACT)
  INTEGER :: NACT, ACTNODES(MAXNACT)
  
  
  ! Reservoirs
  INTEGER, PARAMETER :: MAXNRESV = 1000
  DOUBLE PRECISION :: RESVVOL(MAXNRESV), RESVSA(MAXNRESV)
  LOGICAL :: RESVMIX(MAXNRESV), DORESERVOIRS
  CHARACTER(LEN=100) :: RESVFILE
  LOGICAL :: ALLOWFIXEDRESV

  ! permeability
  ! permeability for each field and external concentrations
  DOUBLE PRECISION :: PERMEABILITY(MAXNABSORBER,MAXNFIELD), CEXT(MAXNFIELD)
  INTEGER :: NPERM ! number of permeable nodes
  INTEGER :: PERMNODES(MAXNABSORBER) ! which nodes are permeable
  LOGICAL :: TRACKFLUXPERM ! track flux out of permeable nodes
  LOGICAL :: OUTPUTTOTFLUXONLY ! output tot flux, not individual nodes
  ! permeable nodes selected randomly, or read in from file
  LOGICAL :: RANDPERMNODES, PERMFROMFILE
  DOUBLE PRECISION:: PERMNEARNODEDIST
  ! permeability parameter supplied is actually a prefactor multiplied by mesh length and species diffusivity?
  LOGICAL :: USEPERMPREFACTOR
  
  ! start with fixed conc on network (normalize to 1 if <0)
  ! optionally, set a background concentration (separate from starting concs around a node)
  DOUBLE PRECISION :: STARTCONC(maxnfield), BACKGROUNDCONC(MAXNFIELD)
  LOGICAL :: SETBACKGROUNDCONC ! will we be using the background conc?
  INTEGER :: STARTEQUIL

  ! Evolve buffer proteins? Rapidly equilibrating proteins?
  LOGICAL :: DOBUFFER, FASTEQUIL
  DOUBLE PRECISION :: KON, KOFF !(on/off rates for buffer prots)
  DOUBLE PRECISION :: KDEQUIL ! (dissociation constant for rapidly equilibrating proteins)

  ! warn when field hits zero?
  LOGICAL :: WARNZEROFIELD

  ! include boundary velocities in snapshots
  INTEGER :: SNAPSHOTVEL

  ! edges with fixed velocities
  INTEGER :: NFIXEDGEVEL
  INTEGER, PARAMETER :: MAXFIXEDGEVEL = 100
  INTEGER :: FIXEDGEVEL(MAXFIXEDGEVEL)
  DOUBLE PRECISION :: FIXEDGEVELVAL(MAXFIXEDGEVEL)
END MODULE KEYS
