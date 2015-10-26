!
! ################# MPI variables module ####################################
! ###########################################################################
MODULE MPI_ENV
IMPLICIT NONE
INTEGER, PARAMETER:: Npe = 1
!
! MYID2   = ID of processor for COMM%COMM1D
! LID     = Local ID for COMM%LOCAL
! COMM1D  = COMMUNICATOR for 64 communication
! LOCAL   = COMMUNICATOR for 8 communication(cluster)
! NUMPROC = Number of available processors
! MOVING  = Data set for moving particles
! GTREE   = Data set for ghost TREE communications
! POS     = Data set for particle positions
! DEST    = Destiny node for COMM%COMM1D
! SRC     = Source node for COMM%COMM1D
! rc      = initializer/finalizer of COMMUNICATOR
TYPE MP
   INTEGER:: MYID2, LID, COMM1D, LOCAL, NUMPROC, MOVING, GTREE, POS, &
        DEST(7), SRC(7), rc
END TYPE MP
!
END MODULE MPI_ENV
!
! ################# General purpose variables module ########################
! ###########################################################################
MODULE DATASTR
IMPLICIT NONE
!
! Nparam  = Number of particle kind (2 = ion/electron)
! Nref    = Reference number for thermostat: 3 for PBC, 6 for vacuum
! Ndist   = Interval for radial distribution
! eps     = Coefficient for Coulomb force 
!           ( in eV*micrometer = 1/8.8542E-12C^2/N/m^2/4/pi)
! Basic unit - 1micrometer/ 1eV/ 1amu
! 
INTEGER, PARAMETER:: Nparam = 2,  Nref = 6, Ndist = 10000, Nradius = 50
REAL*8,  PARAMETER:: tps = 1.0180505305371923D-1 ! in nanoseconds
REAL*8,  PARAMETER:: eps = 0.0014399835661145824D0
REAL*8,  PARAMETER:: F_multi = 0.2 ! Factor for multipole expansion
!
! Particle data type
! xx      = current position
! xv      = current velocity
! ff      = current force by interaction of electrons
TYPE PT
   REAL*8 :: xx(3), xv(3), ff(3)
   INTEGER:: id
END TYPE PT
!
! Particle information of neighboring nodes
! xx      = current position
TYPE GS
   REAL*8 :: xx(3)
END TYPE GS
!
! Copy set for TREE/ghost TREE
! xx      = Center position of TREE box
! q       = Effective charge of twig/leaf
! lx      = Size of TREE box
! dx      = distance between box center and charge center
TYPE GH
   REAL*8 :: xx(3), qeff, lx, dx, D(3), Q(3), QQ(3)
END TYPE GH   
!
! Time data type
! Nloop   = Number of loops
! Nrest   = Number frequency for restart 
! Ndump   = Number frequency for dump
! Nsamp   = Number frequency for sampling
! tmax    = Maximum simulation time
! trest   = Restart file dump interval
! tdump   = Output file dump interval
! tsamp   = Sampling interval
! tnow    = Current time
TYPE TM
   INTEGER:: Nloop, Nrest, Ndump, Nsamp
   REAL*8 :: tmax, trest, tdump, tsamp, tnow, self, comm, other
END TYPE TM
!
! Parameter data type
! q       = Charge of each particle kind
! Te      = Initial electron temperature
! Ti      = Initial ion temperature
! xm      = Mass of particle
! rs      = Initial screening parameter for Kelbg potential 
! sd_ratio= Opening criterion for TREE algorithm
! rest    = Restart option(on/off)
TYPE PM
   REAL*8:: q(Nparam),  Te, Ti, xm(Nparam), rs, sd_ratio
   CHARACTER*3:: rest
END TYPE PM
!
! System variable data type
! temp    = ambient temperature
! Te      = ambient electron temperature
! Ti      = ambient ion temperature
! mv2     = Twice of kinetic energy
! mv2el   = Twice of kinetic energy of electrons
! mv2ion  = Twice of kinetic energy of ions
! Epot    = Potential energy of system
! Epot_e  = Potential energy by interaction of electrons
! Epot_i  = Potential energy by interaction of ions
! kbm     = total mass and kb parameter to measure instant temperature
! kbm_e   = kbm for electron temperature
!    _rad =     with respect to radial temperature
!    _cir =     with respect to circumferential temperature
! kbm_i   = kbm for ion temperature
!    _rad =     with respect to radial temperature
!    _cir =     with respect to circumferential temperature
! rs      = ambient parameter for Kelbg potential
! Te_rad  = ambient electron temperature along radial direction
! Ti_rad  = ambient ion temperature along radial direction
! kelbg   = How to use kelbg parameter
TYPE ST
   REAL*8:: temp, Te, Ti, mv2, mv2el, mv2ion, Epot, Epot_e, Epot_i,  &
        kbm, kbm_e, kbm_i, kbm_e_rad, kbm_e_cir, kbm_i_rad, kbm_i_cir, &
        rs, Te_rad, Ti_rad
   CHARACTER*5::kelbg
END TYPE ST
!
! Various number type
! Nion    = Number of ions of single node
! Nel     = Number of electrons
! NTion   = Number of total ions of total nodes
!         = During the loop,  number of total ions of single cluster
! NTel    = Number of total electrons of total nodes
!         = During the loop,  number of total electrons of single cluster
! NTpt    = Number of total particles of total nodes
! Nnode   = Maximum size of nodal information of copied TREE
! Ncrit   = Criterion for copying TREE
! dx_ee   = Radial distribution interval for electron-electron interactions
! dx_ei   = Radial distribution interval for electron-ion interactions
! dx_ii   = Radial distribution interval for ion-ion interactions
! RDF_ei  = RDF of electron-ion interactions
TYPE NM
   INTEGER:: Nion, Nel, NTion, NTel, NTpt, Npt, Ncrit(7)
   INTEGER:: Nnode, RDF_ei(Ndist+1), RDF_ii(Ndist+1)
   REAL*8 :: dx_ee, dx_ei, dx_ii, dr
END TYPE NM
!
CONTAINS
  !
  ! ############### Particle set reconfiguration functions ###################
  ! ##########################################################################
  !
  ! Resize pointer p into n length
  FUNCTION REALLOC_PT(pc, nn)
    TYPE(PT), POINTER:: pc(:), REALLOC_PT(:)
    INTEGER, intent(in):: nn
    INTEGER:: nold, ierr
    ALLOCATE(REALLOC_PT(1:nn), STAT = ierr)
    IF(ierr /=0) STOP "allocate error"
    IF(.NOT.ASSOCIATED(pc)) RETURN
    nold = MIN(size(pc),nn)
    REALLOC_PT(1:nold) = pc(1:nold)
    DEALLOCATE(pc)
  END FUNCTION REALLOC_PT
  !
  ! Remove some particles from pointer p
  FUNCTION REMOVE_PT(pk, l, n, m)
    TYPE(PT), POINTER:: pk(:), REMOVE_PT(:)
    INTEGER, intent(in):: l, n, m(n)
    INTEGER:: ierr, ii, kk, ni, nf
    ALLOCATE(REMOVE_PT(1:l-n), STAT = ierr)
    IF(ierr /=0) STOP "allocate error"
    ni = 1
    nf = 0
    kk = 1
    DO ii=1, n
       ni = nf + 1
       nf = m(ii) - ii
       REMOVE_PT(ni:nf) = pk(kk:m(ii)-1)
       kk = m(ii)+1
    END DO
    IF (m(n) < l) REMOVE_PT(nf+1:l-n ) = pk(kk:l)
    DEALLOCATE(pk)
  END FUNCTION REMOVE_PT
  !
  ! Add some particles to pointer p
  FUNCTION ADD_PT(px, ll, nl, qx)
    TYPE(PT), POINTER:: px(:), ADD_PT(:)
    INTEGER, intent(in):: ll, nl
    TYPE(PT)           :: qx(nl)
    INTEGER:: ierr
    ALLOCATE(ADD_PT(1:ll+nl), STAT = ierr)
    IF(ierr /=0) STOP "allocate error"
    ADD_PT(1:ll) = px(1:ll)
    ADD_PT(ll+1:ll+nl) = qx(1:nl)
    DEALLOCATE(px)
  END FUNCTION ADD_PT
!
! ################# electron migration between processors ####################
! ############################################################################
SUBROUTINE MIGRATIONE(NS, qel, COMM)
USE MPI_ENV
IMPLICIT NONE
INCLUDE 'mpif.h'
INTERFACE
   FUNCTION id_FINDER(x, y, z)
     IMPLICIT NONE
     REAL*8, INTENT(IN) :: x, y, z
     INTEGER            :: id_FINDER
     INTEGER            :: i, j, k
   END FUNCTION id_FINDER
END INTERFACE
!
INTEGER,  PARAMETER    :: Nsize = 100
TYPE(NM), INTENT(INOUT):: NS
TYPE(MP), INTENT(IN)   :: COMM
TYPE(PT), POINTER      :: qel(:)
TYPE(PT)  :: q(Nsize,7), qsend(Nsize), qrec(Nsize,7)
INTEGER :: Nsend(7), Nout(Npe), Nsum, Nel(Npe), Nrec_e(7)
INTEGER :: i, j, k, k_e, id, Nhead, DEST(7), SRC(7), &
     Nprune, Dprune(Nsize), displ(Npe), ISTATUS(MPI_STATUS_SIZE), IERR
!
! Electron migration #########################################################
! Cluster -> HEAD -> HEAD -> Cluster
Nsend = 0
Nhead = INT(COMM%MYID2/Npe)
Nprune = 0
DO i=1, NS%Nel
   id = id_FINDER(qel(i)%xx(1), qel(i)%xx(2), qel(i)%xx(3))
   IF (id /= Nhead) THEN
      k = MOD(id-Nhead+8,8)
      Nsend(k) = Nsend(k) + 1
      Nprune = Nprune + 1
      q(Nsend(k), k)%xx(:) = qel(i)%xx(:)
      q(Nsend(k), k)%xv(:) = qel(i)%xv(:)
      q(Nsend(k), k)%ff(:) = qel(i)%ff(:)
      q(Nsend(k), k)%id    = qel(i)%id
      Dprune(Nprune) = i
      IF (Nprune > Nsize) STOP "====== Not enough e_migration array ======="
   END IF
END DO
!
! Prune electron data sets
IF (Nprune > 0) THEN      
   qel =>  REMOVE_PT(qel, NS%Nel, Nprune, Dprune(1:Nprune) )
   NS%Nel = NS%Nel - Nprune
END IF
!
DO i=1, 7   
   ! Gather pruned data into HEAD
   CALL MPI_ALLGATHER(Nsend(i), 1, MPI_INTEGER, Nout, 1, MPI_INTEGER, &
        COMM%LOCAL, IERR)
   Nsum = 0
   DO j=1, Npe
      Nsum = Nsum + Nout(j)
   END DO
   displ(1) = 0   
   DO j=1, Npe-1
      displ(j+1) = Nout(j)+displ(j)
   END DO
   CALL MPI_GATHERV(q(:,i), Nsend(i), COMM%MOVING, qsend, Nout, displ, &
        COMM%MOVING, 0, COMM%LOCAL, IERR)
   ! Send/receive between HEADs
   IF (COMM%MYID2 == Nhead*Npe) THEN
      DEST(i) = MOD((Nhead + i)*Npe, Npe*8)
      SRC(i)  = MOD((Nhead - i)*Npe+Npe*8, Npe*8)
      CALL MPI_SENDRECV(Nsum, 1, MPI_INTEGER, DEST(i), 0, &
           Nrec_e(i), 1, MPI_INTEGER, SRC(i), 0, COMM%COMM1D, ISTATUS, IERR)
      CALL MPI_SENDRECV(qsend, Nsum, COMM%MOVING, DEST(i), 0, qrec(:,i), &
           Nrec_e(i), COMM%MOVING, SRC(i), 0, COMM%COMM1D, ISTATUS, IERR)
   END IF
END DO
!
! Check particle number for each node in cluster
CALL MPI_ALLGATHER(NS%Nel, 1,  MPI_INTEGER, Nel, 1, MPI_INTEGER, &
     COMM%LOCAL, IERR)
CALL MPI_BCAST(Nrec_e, 7, MPI_INTEGER, 0, COMM%LOCAL, IERR)
DO i=1, 7
   CALL MPI_BCAST(qrec(:,i), Nrec_e(i), COMM%MOVING, 0, COMM%LOCAL, IERR)
END DO
!
! Distribute incoming particle information into minimum particle node
DO i=1, 7
   !
   ! Find the node whose particles are the least
   k_e = MINLOC(Nel,1)
   !
   IF (Nrec_e(i) > 0 .AND. COMM%LID == (k_e-1)) THEN
      qel => ADD_PT(qel, NS%Nel, Nrec_e(i), qrec(1:Nrec_e(i),i))
      NS%Nel = NS%Nel + Nrec_e(i)
   END IF
   Nel(k_e) = Nel(k_e) + Nrec_e(i)
END DO
CALL MPI_BARRIER(COMM%COMM1D, IERR)
!
END SUBROUTINE MIGRATIONE
!
! ################# ion migration between processors #########################
! ############################################################################
SUBROUTINE MIGRATIONI(NS, qion, COMM)
USE MPI_ENV
IMPLICIT NONE
INCLUDE 'mpif.h'
INTERFACE
   FUNCTION id_FINDER(x, y, z)
     IMPLICIT NONE
     REAL*8, INTENT(IN) :: x, y, z
     INTEGER            :: id_FINDER
     INTEGER            :: i, j, k
   END FUNCTION id_FINDER
END INTERFACE
!
INTEGER,  PARAMETER    :: Nsize = 100
TYPE(NM), INTENT(INOUT):: NS
TYPE(MP), INTENT(IN)   :: COMM
TYPE(PT), POINTER      :: qion(:)
TYPE(PT):: q(Nsize,7), qsend(Nsize), qrec(Nsize,7)
INTEGER :: Nsend(7), Nout(Npe), Nsum, Nrec_i(7), Nion(Npe)
INTEGER :: Nhead, i, j,  k, k_i, id, displ(Npe), DEST(7), SRC(7), &
     Nprune, Dprune(Nsize), ISTATUS(MPI_STATUS_SIZE), IERR
!
! Ion migration #############################################################
! Cluster -> HEAD -> HEAD -> Cluster
Nsend = 0
Nhead = INT(COMM%MYID2/Npe)
Nprune = 0
DO i=1, NS%Nion
   id = id_FINDER(qion(i)%xx(1), qion(i)%xx(2), qion(i)%xx(3))
   IF (id /= Nhead) THEN
      k = MOD(id-Nhead+8,8)
      Nsend(k) = Nsend(k) + 1
      Nprune = Nprune + 1
      q(Nsend(k), k)%xx(:) = qion(i)%xx(:)
      q(Nsend(k), k)%xv(:) = qion(i)%xv(:)
      q(Nsend(k), k)%ff(:) = qion(i)%ff(:)
      q(Nsend(k), k)%id    = qion(i)%id
      Dprune(Nprune) = i
      IF (Nprune > Nsize) STOP "====== Not enough i_migration array ======="
   END IF
END DO
!
! Prune ion data sets
IF (Nprune > 0) THEN
   qion =>  REMOVE_PT(qion, NS%Nion, Nprune, Dprune(1:Nprune) )
   NS%Nion = NS%Nion - Nprune
END IF
!
DO i=1, 7   
   ! Gather pruned data into HEAD
   CALL MPI_ALLGATHER(Nsend(i), 1, MPI_INTEGER, Nout, 1, MPI_INTEGER, &
        COMM%LOCAL, IERR)
   Nsum = 0
   DO j=1, Npe
	Nsum = Nsum + Nout(j)
   END DO
   displ(1) = 0   
   DO j=1, Npe-1
      displ(j+1) = Nout(j)+displ(j)
   END DO
   CALL MPI_GATHERV(q(:,i), Nsend(i), COMM%MOVING, qsend, Nout, displ, &
        COMM%MOVING, 0, COMM%LOCAL, IERR)
   ! Send/receive between HEADs
   IF (COMM%MYID2 == Nhead*Npe) THEN
      DEST(i) = MOD((Nhead + i)*Npe, Npe*8)
      SRC(i)  = MOD((Nhead - i)*Npe+Npe*8, Npe*8)
      CALL MPI_SENDRECV(Nsum, 1, MPI_INTEGER, DEST(i), 0, &
           Nrec_i(i), 1, MPI_INTEGER, SRC(i), 0, COMM%COMM1D, ISTATUS, IERR)
      CALL MPI_SENDRECV(qsend, Nsum, COMM%MOVING, DEST(i), 0, qrec(:,i), &
           Nrec_i(i), COMM%MOVING, SRC(i), 0, COMM%COMM1D, ISTATUS, IERR)
   END IF
END DO
!
! Check particle number for each node in cluster
CALL MPI_ALLGATHER(NS%Nion, 1,  MPI_INTEGER, Nion, 1, MPI_INTEGER, &
     COMM%LOCAL, IERR)
CALL MPI_BCAST(Nrec_i, 7, MPI_INTEGER, 0, COMM%LOCAL, IERR)
DO i=1,7
   CALL MPI_BCAST(qrec(:,i), Nrec_i(i), COMM%MOVING, 0, COMM%LOCAL, IERR)
END DO
!
! Distribute incoming particle information into minimum particle node
DO i=1, 7
   !
   ! Find the node whose particles are the least
   k_i = MINLOC(Nion,1)
   !
   IF (Nrec_i(i) > 0 .AND. COMM%LID == (k_i-1)) THEN
      qion => ADD_PT(qion, NS%Nion, Nrec_i(i), qrec(1:Nrec_i(i),i))
      NS%Nion = NS%Nion + Nrec_i(i)
   END IF
   Nion(k_i) = Nion(k_i) + Nrec_i(i)
END DO
!
END SUBROUTINE MIGRATIONI
!
END MODULE DATASTR
