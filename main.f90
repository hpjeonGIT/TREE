!
! ######### Program for effective MD of cold plasma analysis ###############
! ################ #########################################################
! 
! Managing two kind of TREEs for electrons and ions, multiple time step
! scheme has been implemented. Basic loops are:
! 1. Data parsing
! 1-1. Distribute along cluster.
! 1-2. Redistribute along the nodes of each cluster
! 2. Electron TREE build
! 2-1. Copy neighboring electrons inside of same cluster
! 3. Ion TREE Build
! 3-1. Copy neighboring ions inside of same cluster
! 3-2. Ion TREE is built every ten iterations
! 4. Coulomb charge and force estimation
! 5. Copy TREE into array form
! 6. SEND/RECV between clusters for electron ghost TREE
! 6-1. Communicate between heads of the clusters
! 6-2. Recopy received TREE into all nodes of each cluster
! 6-3. Rebuild ghost(received) TREEs
! 6-4. Coulomg charge and force
! 7. SEND/RECV between clusters for ion ghost TREE (done at every ten loops)
! 7-1. Communicate between heads of the clusters
! 7-2. Recopy received TREE into all nodes of each cluster
! 7-3. Rebuild ghost(received) TREEs
! 7-4. Coulomb charge and force
! 8. Migration/restart/sampling
! 9. Time update
!
! Time integration is done by velocity Verlet. Coulomb force is calculated 
! using Kelbg potential method.
!
! ############################################################################
!
! Program has been developed by:
! Byoungseon Jeon
! Graduate student, Department of Applied Science, UC.Davis
! GRA, Theoretical division, Los Alamos National Laboratory.
! July 08, 2006
!
! Multipole Acceptance Criteria modified
! Nov. 10. 2006
!
PROGRAM OCTREE_HYBRID
USE DATASTR
USE BINTREE
USE MPI_ENV
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(PT), POINTER:: qel(:), qion(:)
TYPE(TM):: ttime
TYPE(PM):: param
TYPE(ST):: sys
TYPE(NM):: NS
TYPE(MP):: COMM
REAL*8  :: dt
!
! INTERNAL VARIABLES
REAL     :: time0, time1, time2, secnds
REAL*8   :: a_time
INTEGER  :: i, k, IERR, Nloop_max, Nfreq, OMP_GET_NUM_THREADS, &
     OMP_GET_THREAD_NUM, TID, NTD

CHARACTER(LEN=20):: ENERFILE
!
time0 = 0
time1 = secnds(time0)
ENERFILE = "ener000.dat"
!
! ###########################################################################
! MPI Initialization
CALL MPI_INITIALIZE(COMM, NS)
!$OMP PARALLEL PRIVATE(TID, NTD)
TID = OMP_GET_THREAD_NUM()
NTD = OMP_GET_NUM_THREADS()
IF ( COMM%MYID2 == 0 .AND. TID == 0) THEN
   PRINT *, "Number of available threads = ", NTD
END IF
!$OMP END PARALLEL
!
! ###########################################################################
! File open for sampling data 
IF ( COMM%MYID2 == 0) THEN
   OPEN(UNIT=25,file=ENERFILE)
   WRITE(25,200)      
END IF
200 FORMAT("# time(ns)-   T_el.  -  T_ion  -  T_el_radial  -  T_ion_radial &
         & -  T_el_circumf  -  T_ion_circumf  -  Potential Energy  - &
         & Total energy  -  Number of Rydberg atoms")
!
! ###########################################################################
! Parsing input data and allocate pointer variables
!
CALL Init(NS, qel, qion, ttime, param, sys, dt, COMM)
ttime%self = 0.D0
ttime%comm = 0.D0
ttime%other= 0.D0
sys%rs = param%rs
IF (param%rest == 'OFF') THEN
   ! Initial velocity from Normal distribution
   CALL Vinit(NS, qel, qion, param, ttime, sys, dt)
ELSE
   ! Initial velocity from restart file
   CALL Vinit_REST(NS, qel, qion, param, ttime, sys, dt)
END IF
!
! Time parameters
ttime%Ndump = NINT(ttime%tdump/dt)
ttime%Nrest = NINT(ttime%trest/dt)
ttime%Nsamp = NINT(ttime%tsamp/dt)
ttime%Nloop = 0
ttime%tnow = 0.0
ttime%tmax = ttime%tmax+0.5D0*dt
Nloop_max = NINT(ttime%tmax/dt)-1
!
! RDF
NS%RDF_ei = 0.D0
NS%RDF_ii = 0.D0
!
! Time integration loop starts ###############################################
DO WHILE (ttime%Nloop <= Nloop_max)
   a_time = MPI_Wtime()
   ttime%tnow  = ttime%tnow  + dt
   ttime%Nloop = ttime%Nloop + 1
   ! Velocity verlet 1st step
   CALL VVerletNotemp1(NS, qel, qion, param, dt, ttime)
   ttime%self = ttime%self + MPI_Wtime() - a_time   
   CALL Force (NS, qel, qion, param, ttime)
   a_time = MPI_Wtime()
   CALL VVerletNotemp2(NS, qel, qion, param, sys, dt, ttime)
   !
   ! Electron migration between processors ##################################
   IF (MOD(ttime%Nloop,10) == 1) THEN
      CALL MIGRATIONE(NS, qel, COMM)
   END IF
   !
   ! Ion migration between processors #######################################
   IF (MOD(ttime%Nloop,100) == 1) THEN
      CALL MIGRATIONI(NS, qion, COMM)
   END IF
   !
   ! RDF print
   IF (MOD(ttime%Nloop,ttime%Ndump) == 0) THEN
      CALL RDF(NS, ttime, COMM)
   END IF
   ! 
   ! Print sampled data
   sys%Epot = sys%Epot_e + sys%Epot_i
   IF (MOD(ttime%Nloop,ttime%Nsamp) == 0) THEN
      CALL TEMPORAL(NS, sys, ttime, COMM)
   END IF
   ttime%self = ttime%self + MPI_Wtime() - a_time
   !
   ! Restart file print
   IF (MOD(ttime%Nloop,ttime%Nrest) == 0) THEN
      CALL Restart(NS, qel, qion, ttime, COMM)
      IF (COMM%MYID2 == 0) THEN
         CLOSE(25)
         Nfreq = INT(ttime%Nloop/ttime%Nrest)
         WRITE(ENERFILE,480) Nfreq
         OPEN(UNIT=25, file=ENERFILE)
         WRITE(25,200)
      END IF
   END IF
END DO
480 FORMAT("ener",I3.3,".dat")
!
! Close XYZ file for output
IF (COMM%MYID2 == 0) THEN
   CLOSE(25)
END IF
!
! Wall clock
time2 = secnds(time1)
!time0 = MPI_Wtime()
PRINT '("Wall time is" , ES11.3, " at CPU ", I3)', time2, COMM%MYID2
PRINT '(3(1X, ES11.3))', ttime%self, ttime%comm, ttime%other
!PRINT *, COMM%MYID2, Nsum, Nrsum
!
! MPI END
CALL MPI_Finalize(COMM%rc)
DEALLOCATE(qel, qion)
!
! Main routine ends here. Below are other subroutines
CONTAINS
!
! ####### Initialization and input data parsing routine ######################
! Unit normalization
! Length: 1. means 1micrometer = 1E-6 m
! Mass: 1. means 1.6605E-27 kg (a.m.u.)
! Energy: 1. means 1 eV = 1.6022E-19 J
! Time: 1. means .1018 ns = 1.018E-10 sec
!
SUBROUTINE Init(NS, qel, qion, ttime, param, sys, dt, COMM)
INTERFACE
   FUNCTION id_FINDER(x, y, z)
     IMPLICIT NONE
     REAL*8, INTENT(IN) :: x, y, z
     INTEGER            :: id_FINDER
     INTEGER            :: i, j, k
   END FUNCTION id_FINDER
END INTERFACE
!
TYPE(PT), POINTER       :: qel(:), qion(:)
TYPE(NM), INTENT(INOUT) :: NS
TYPE(TM), INTENT(INOUT) :: ttime
TYPE(PM), INTENT(INOUT) :: param
TYPE(ST), INTENT(INOUT) :: sys
TYPE(MP), INTENT(IN)    :: COMM
REAL*8,   INTENT(INOUT) :: dt
!
! INTERNAL VARIABLES
INTEGER  :: i, j, k, IERR, id, Ntag(8), nid
INTEGER:: Ntemp
REAL*8   :: x, y, z, xv, yv, zv, Tel, Tion, xm1, xm2
CHARACTER(len=5):: dummy
!
!
![[[[[[[[[[[[[[[[ "control.prm" parsing and data arrangement ]]]]]]]]]]]]]]]
OPEN(UNIT=15, file="control.prm")
READ(15,*) dummy
!
! charge q
READ(15,*) dummy
DO i = 1, Nparam
   READ(15,*) param%q(i)
END DO
!
! Mass for each particle kind
READ(15,*) dummy
DO i = 1, Nparam
   READ(15,*) param%xm(i)
END DO
!
! Number of particles
READ(15,*) dummy
READ(15,*) NS%NTel ! Number of electrons
READ(15,*) NS%NTion ! Number of ions
NS%NTpt = NS%NTel + NS%NTion
!
! Initial temperatures
READ(15,*) dummy
READ(15,*) param%Te
READ(15,*) param%Ti
!
! Time parameters
! ############################################################################
! Input data time unit = ns
! Normalized time unit in code = 0.1ns
READ(15,*) dummy
READ(15,*) ttime%tmax, ttime%trest, ttime%tdump, ttime%tsamp, dt
ttime%tmax  = ttime%tmax  / tps
ttime%trest = ttime%trest / tps
ttime%tdump = ttime%tdump / tps
ttime%tsamp = ttime%tsamp / tps
dt = dt / tps
!
! Screening paramters
READ(15,*) dummy
READ(15,*) sys%kelbg, param%rs
!
! Restart parameter
READ(15,*) dummy
READ(15,*) param%rest
!
! RDF parameter
READ(15,*) dummy
READ(15,*) NS%dx_ee, NS%dx_ei, NS%dx_ii
!
! TREE opening criterion
READ(15,*) dummy
READ(15,*) param%sd_ratio
!
! Close "control.prm"
CLOSE(15)
!
!
![[[[[[[[[[[[[[[[[ "input.dat" parsing and data arrangement ]]]]]]]]]]]]]]]]]
OPEN (UNIT=11, file="input.xyz", STATUS = "OLD")
!
! Number of electrons/ions per each node
NS%Npt = INT(NS%NTpt/COMM%NUMPROC)
!
! Allocate number of electrons/ions
NS%Npt  = INT(NS%Npt*3.0)
!
! Allocate electrons/ions
ALLOCATE(qel(NS%Npt), qion(NS%Npt))
!
! ############################################################################
! Read position data
! 
NS%Nel = 0
NS%Nion = 0
Ntag = 0
READ(11,*) dummy
READ(11,*) dummy
IF (param%rest == 'OFF') THEN
   !
   ! Restart is turned off - read only positions
   j = 0
   DO i = 1, NS%NTel
      READ(11,*) dummy, x, y, z, nid
      k = id_FINDER(x,y,z)+1
      id = Npe*(k - 1) + Ntag(k)
      Ntag(k) = MOD(Ntag(k)+1,Npe)
      IF (id == COMM%MYID2) THEN
         j = j + 1
         qel(j)%xx(1) = x
         qel(j)%xx(2) = y
         qel(j)%xx(3) = z
         qel(j)%id    = nid
         NS%Nel = NS%Nel + 1
      END IF
   END DO
   j = 0
   DO i = 1, NS%NTion
      READ(11,*) dummy, x, y, z, nid
      k = id_FINDER(x,y,z)+1
      id = Npe*(k - 1) + Ntag(k)
      Ntag(k) = MOD(Ntag(k)+1,Npe)
      IF (id == COMM%MYID2) THEN         
         j = j + 1
         qion(j)%xx(1) = x
         qion(j)%xx(2) = y
         qion(j)%xx(3) = z
         qion(j)%id    = nid
         NS%Nion = NS%Nion + 1
      END IF
   END DO
ELSE
   !
   ! Read from restart file - velocity loaded
   j = 0
   sys%mv2el = 0.0
   xm1 = param%xm(1)
   DO i = 1, NS%NTel
      READ(11,*) dummy, x, y, z, xv, yv, zv, nid
      !
      k = id_FINDER(x,y,z)+1
      id = Npe*(k - 1) + Ntag(k)
      Ntag(k) = MOD(Ntag(k)+1,Npe)
      IF (id == COMM%MYID2) THEN         
         j = j + 1
         qel(j)%xx(1) = x
         qel(j)%xx(2) = y
         qel(j)%xx(3) = z
         qel(j)%xv(1) = xv
         qel(j)%xv(2) = yv
         qel(j)%xv(3) = zv
         qel(j)%id    = nid
         sys%mv2el = sys%mv2el + xm1*(xv**2+yv**2+zv**2)
         NS%Nel = NS%Nel + 1
      END IF
   END DO
   j = 0
   sys%mv2ion = 0.0
   xm2 = param%xm(2)
   DO i = 1, NS%NTion
      READ(11,*) dummy, x, y, z, xv, yv, zv, nid
      !
      k = id_FINDER(x,y,z)+1
      id = Npe*(k - 1) + Ntag(k)
      Ntag(k) = MOD(Ntag(k)+1,Npe)
      IF (id == COMM%MYID2) THEN         
         j = j + 1
         qion(j)%xx(1) = x
         qion(j)%xx(2) = y
         qion(j)%xx(3) = z
         qion(j)%xv(1) = xv
         qion(j)%xv(2) = yv
         qion(j)%xv(3) = zv
         qion(j)%id    = nid
         NS%Nion = NS%Nion + 1
         sys%mv2ion = sys%mv2ion + xm2*(xv**2+yv**2+zv**2)  
      END IF
   END DO
END IF
NS%Npt = NS%Nel + NS%Nion
qel  => REALLOC_PT(qel,  NS%Nel)
qion => REALLOC_PT(qion, NS%Nion)
!
! Marginal size of buffer
NS%Nnode = 100*NS%Npt
!
! Close "input.bin"
CLOSE(11)
!
! ############################################################################
! Parameter for instantaneous temperature calculation
sys%kbm = 8.617343D-5*REAL(3*NS%NTel+3*NS%NTion-Nref)
sys%kbm_e = 8.617343D-5*REAL(3*NS%NTel-Nref)
sys%kbm_i = 8.617343D-5*REAL(3*NS%NTion-Nref)
sys%kbm_e_rad = 8.617343D-5*REAL(NS%NTel-Nref)
sys%kbm_e_cir = 8.617343D-5*REAL(2*NS%NTel-Nref)
sys%kbm_i_rad = 8.617343D-5*REAL(NS%NTion-Nref)
sys%kbm_i_cir = 8.617343D-5*REAL(2*NS%NTion-Nref)
Ntemp = 0
Tel = 0.0D0
Tion = 0.0D0
!
! temperature estimation for restart
IF (param%rest /= 'OFF') THEN
   CALL MPI_REDUCE(sys%mv2el, Tel, 1, MPI_REAL8, MPI_SUM, 0, &
        COMM%COMM1D, IERR)
   CALL MPI_REDUCE(sys%mv2ion, Tion, 1, MPI_REAL8, MPI_SUM, 0, &
        COMM%COMM1D, IERR)
   sys%temp = (Tel+Tion)/sys%kbm
   CALL MPI_BCAST(sys%temp, 1, MPI_REAL8, 0, COMM%COMM1D, IERR)
END IF
!
! Number of particle check
CALL MPI_REDUCE(NS%Npt, Ntemp, 1, MPI_INTEGER, MPI_SUM, 0, &
     COMM%COMM1D, IERR)
IF ( COMM%MYID2 == 0) THEN
   IF (Ntemp /= NS%NTpt) THEN
      PRINT *, "parsing error", Ntemp, "while given number is", NS%NTpt
      STOP
   END IF
END IF
!
RETURN
!
END SUBROUTINE Init
!
! ####################### particle gathering for electrons - 2
! Using gatherv, electron informations are gathered
SUBROUTINE PTEgathering2(NS, qel, gel, COMM)
USE DATASTR
USE MPI_ENV
USE BINTREE
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(NM):: NS
TYPE(MP):: COMM
TYPE(PT):: qel(NS%Nel)
TYPE(GS), POINTER:: gel(:)
!
INTEGER:: Nel(Npe), displ(Npe)
INTEGER  :: i, IERR
TYPE(GS) :: el_buffer(NS%Nel)
!
el_buffer(:)%xx(1) = qel(:)%xx(1)
el_buffer(:)%xx(2) = qel(:)%xx(2)
el_buffer(:)%xx(3) = qel(:)%xx(3)
!
CALL MPI_Allgather(NS%Nel, 1, MPI_INTEGER, Nel, 1, MPI_INTEGER, &
     COMM%LOCAL, IERR)
NS%NTel = 0
DO i=1,Npe
   NS%NTel = NS%NTel + Nel(i)
END DO
ALLOCATE(gel(NS%NTel))
displ(1) = 0
DO i = 1, Npe-1
   displ(i+1) = Nel(i)+ displ(i)
END DO
CALL MPI_AllGatherv(el_buffer, NS%Nel, COMM%POS, gel, Nel, displ, &
     COMM%POS, COMM%LOCAL, IERR)
!
RETURN
END SUBROUTINE PTEgathering2
!
! ####################### particle gathering for ion
! Using gatherv, ion informations are gathered
SUBROUTINE PTIgathering2(NS, qion, gion, COMM)
USE DATASTR
USE MPI_ENV
USE BINTREE
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(NM):: NS
TYPE(MP):: COMM
TYPE(PT):: qion(NS%Nion)
TYPE(GS), POINTER:: gion(:)
!
INTEGER:: Nion(Npe), displ(Npe)
INTEGER  :: i, IERR
TYPE(GS) :: ion_buffer(NS%Nion)
!
ion_buffer(:)%xx(1) = qion(:)%xx(1)
ion_buffer(:)%xx(2) = qion(:)%xx(2)
ion_buffer(:)%xx(3) = qion(:)%xx(3)
!
CALL MPI_Allgather(NS%Nion, 1, MPI_INTEGER, Nion, 1, MPI_INTEGER, &
     COMM%LOCAL, IERR)
NS%NTion = 0
DO i=1, Npe
   NS%NTion = NS%NTion + Nion(i)
END DO
ALLOCATE(gion(NS%NTion))
displ(1) = 0
DO i = 1, Npe-1
   displ(i+1) = Nion(i)+ displ(i)
END DO
CALL MPI_AllGatherv(ion_buffer, NS%Nion, COMM%POS, gion, Nion, displ, &
     COMM%POS, COMM%LOCAL, IERR)
!
RETURN
END SUBROUTINE PTIgathering2
!
!
!
SUBROUTINE Force (NS, qel, qion, param, ttime)
USE DATASTR
USE MPI_ENV
USE BINTREE
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(NM), INTENT(INOUT)::NS
TYPE(PT), INTENT(INOUT)::qel(NS%Nel), qion(NS%Nion)
TYPE(PM), INTENT(IN)::param
TYPE(TM), INTENT(INOUT)::ttime
!
!
TYPE(BT), POINTER:: Roote, Rooti
TYPE(GT), POINTER:: gRoote, gRooti, temp
TYPE(GS), POINTER:: gel(:), gion(:)
TYPE(GH):: node(NS%Nnode,7), rnode(NS%Nnode,7)
INTEGER :: arch(NS%Nnode,7), rarch(NS%Nnode,7)
INTEGER :: Nnode(7), Narch(7), Nrnod(7), Nrarc(7)
INTEGER  :: i, k, ISTATUS(MPI_STATUS_SIZE), IERR, Ncrit, ireq(28), &
     STATUSARR(MPI_STATUS_SIZE,28), irq(14), SARR(MPI_STATUS_SIZE,14)
LOGICAL:: Log_head
REAL*8:: a_time, b_time
!   
!
! Build local tree and sending tree
a_time = MPI_Wtime()
NULLIFY(temp)
ALLOCATE(temp)
NULLIFY(Roote)
ALLOCATE(Roote)
CALL PTEgathering2(NS, qel, gel, COMM)
NULLIFY(Rooti)
ALLOCATE(Rooti)
CALL PTIgathering2(NS, qion, gion, COMM)
CALL BuildLocalTree(NS, gel, gion, Roote, Rooti, ttime)
CALL LocalCoulombSum(NS, Roote, Rooti, qel, qion, sys, param, ttime)
IF (MOD(COMM%MYID2,Npe) == 0) THEN
   Log_head = .TRUE.
ELSE
   Log_head = .FALSE.
END IF!
Nnode = 0
Narch = 0
ttime%self = ttime%self + MPI_Wtime() - a_time
a_time = MPI_Wtime()
!
! Copy TREE
IF (Log_head) THEN
   DO i = 1, 7
      Ncrit = NS%Ncrit(i)
      CALL CopyTree(Roote, arch(:,i), node(:,i), Nnode(i), Narch(i), &
           NS, param, Ncrit)
   END DO
END IF
CALL DEALLOC_TREE(Roote)
DEALLOCATE(gel)
ttime%other = ttime%other + MPI_Wtime() - a_time
!
a_time = MPI_Wtime()
CALL MPI_BARRIER(COMM%COMM1D, IERR)
!
! Size of TREE send/receive
IF (Log_head) THEN
   DO i = 1, 7
      CALL MPI_ISEND(Nnode(i), 1, MPI_INTEGER, COMM%DEST(i), i*2-1, &
           COMM%COMM1D, irq(2*i-1), IERR)
      CALL MPI_RECV(Nrnod(i), 1, MPI_INTEGER, COMM%SRC(i), i*2-1, &
           COMM%COMM1D,  ISTATUS, IERR)
      CALL MPI_ISEND(Narch(i), 1, MPI_INTEGER, COMM%DEST(i), i*2, &
           COMM%COMM1D, irq(2*i), IERR)
      CALL MPI_RECV(Nrarc(i), 1, MPI_INTEGER, COMM%SRC(i), i*2, &
           COMM%COMM1D, ISTATUS, IERR)
   END DO
   CALL MPI_Waitall(14, irq, SARR, IERR)
END IF
!
! TREE array send/receive
IF (Log_head) THEN
   DO i = 1, 7
      CALL MPI_ISEND(node(1:Nnode(i),i), Nnode(i), COMM%GTREE, &
           COMM%DEST(i), i*2-1, COMM%COMM1D, ireq(4*i-3), IERR)
      CALL MPI_IRECV(rnode(1:Nrnod(i),i), Nrnod(i), COMM%GTREE, &
           COMM%SRC(i), i*2-1, COMM%COMM1D, ireq(4*i-2), IERR)
      CALL MPI_ISEND(arch(1:Narch(i),i), Narch(i), MPI_INTEGER, &
           COMM%DEST(i), i*2, COMM%COMM1D, ireq(4*i-1), IERR)
      CALL MPI_IRECV(rarch(1:Nrarc(i),i), Nrarc(i), MPI_INTEGER, &
           COMM%SRC(i), i*2, COMM%COMM1D, ireq(4*i),   IERR)
   END DO
   CALL MPI_Waitall(28, ireq, STATUSARR, IERR)
END IF
!
! Copy size of TREE onto neighboring nodes of the cluster
CALL MPI_BCAST(Nrarc, 7, MPI_INTEGER, 0, COMM%LOCAL, IERR)
CALL MPI_BCAST(Nrnod, 7, MPI_INTEGER, 0, COMM%LOCAL, IERR)
!
! Copy TREE array onto neighboring nodes of the cluster
DO i = 1,7
   CALL MPI_BCAST(rarch(:,i),Nrarc(i), MPI_INTEGER, 0, COMM%LOCAL, IERR)
   CALL MPI_BCAST(rnode(:,i),Nrnod(i), COMM%GTREE,  0, COMM%LOCAL, IERR)
END DO
ttime%comm = ttime%comm + MPI_Wtime() - a_time
!
! Ghost TREE build for electrons
a_time = MPI_Wtime()
DO i = 1, 7      
   NULLIFY(gRoote)
   ALLOCATE(gRoote)
   gRoote%Root => temp
   !
   ! Rebuild ghost tree using received array data
   Nnode(i) = 0
   CALL RebuildGhostTree(gRoote, temp, rarch(1:Nrarc(i),i), &
        rnode(1:Nrnod(i),i), Nnode(i), Nrnod(i), Nrarc(i))
   !
   ! Estimate Coulomb sum
   CALL GhostCoulombSum_e(NS, qel, qion, sys, param, gRoote)
   CALL DEALLOC_gTREE(gRoote)
END DO
ttime%other = ttime%other + MPI_Wtime() - a_time
a_time = MPI_Wtime()
!
! ##################### CPU SWAP for ION #################################
Nnode = 0
Narch = 0
!
! Copy TREE
IF (Log_head) THEN
   DO i = 1, 7
      Ncrit = NS%Ncrit(i)
      CALL CopyTree(Rooti, arch(:,i), node(:,i), Nnode(i), Narch(i), &
           NS, param, Ncrit)
   END DO
END IF
CALL DEALLOC_TREE(Rooti)
DEALLOCATE(gion)
CALL MPI_BARRIER(COMM%COMM1D, IERR)
ttime%other = ttime%other + MPI_Wtime() - a_time
!
a_time = MPI_Wtime()
!
! Size of TREE send/receive
IF (Log_head) THEN
   DO i = 1, 7
      CALL MPI_ISEND(Nnode(i), 1, MPI_INTEGER, COMM%DEST(i), i*2-1, &
           COMM%COMM1D, irq(2*i-1), IERR)
      CALL MPI_RECV(Nrnod(i), 1, MPI_INTEGER, COMM%SRC(i), i*2-1, &
           COMM%COMM1D,  ISTATUS, IERR)
      CALL MPI_ISEND(Narch(i), 1, MPI_INTEGER, COMM%DEST(i), i*2, &
           COMM%COMM1D, irq(2*i), IERR)
      CALL MPI_RECV(Nrarc(i), 1, MPI_INTEGER, COMM%SRC(i), i*2, &
           COMM%COMM1D, ISTATUS, IERR)
   END DO
   CALL MPI_Waitall(14, irq, SARR, IERR)
END IF
!
! TREE array send/receive
IF (Log_head) THEN
   DO i = 1, 7
      CALL MPI_ISEND(node(1:Nnode(i),i), Nnode(i), COMM%GTREE, &
           COMM%DEST(i), i*2-1, COMM%COMM1D, ireq(4*i-3), IERR)
      CALL MPI_IRECV(rnode(1:Nrnod(i),i), Nrnod(i), COMM%GTREE, &
           COMM%SRC(i), i*2-1, COMM%COMM1D, ireq(4*i-2), IERR)
      CALL MPI_ISEND(arch(1:Narch(i),i), Narch(i), MPI_INTEGER, &
           COMM%DEST(i), i*2, COMM%COMM1D, ireq(4*i-1), IERR)
      CALL MPI_IRECV(rarch(1:Nrarc(i),i), Nrarc(i), MPI_INTEGER, &
           COMM%SRC(i), i*2, COMM%COMM1D, ireq(4*i),   IERR)
   END DO
   CALL MPI_Waitall(28, ireq, STATUSARR, IERR)
END IF
!
! Copy size of TREE onto neighboring nodes of the cluster
CALL MPI_BCAST(Nrarc, 7, MPI_INTEGER, 0, COMM%LOCAL, IERR)
CALL MPI_BCAST(Nrnod, 7, MPI_INTEGER, 0, COMM%LOCAL, IERR)
!
! Copy TREE array onto neighboring nodes of the cluster
DO i = 1, 7
   CALL MPI_BCAST(rarch(:,i),Nrarc(i), MPI_INTEGER, 0, COMM%LOCAL, IERR)
   CALL MPI_BCAST(rnode(:,i),Nrnod(i), COMM%GTREE,  0, COMM%LOCAL, IERR)
END DO
ttime%comm = ttime%comm + MPI_Wtime() - a_time
a_time = MPI_Wtime()
!
! Ghost TREE build for ions
DO i = 1, 7
   NULLIFY(gRooti)
   ALLOCATE(gRooti)
   Nnode(i) = 0
   CALL RebuildGhostTree(gRooti, temp, rarch(1:Nrarc(i),i), &
        rnode(1:Nrnod(i),i), Nnode(i), Nrnod(i), Nrarc(i))
   CALL GhostCoulombSum_i(NS, qel, qion, sys, param, gRooti)
   CALL DEALLOC_gTREE(gRooti)
END DO
DEALLOCATE(temp)
ttime%other = ttime%other + MPI_Wtime() - a_time
!
END SUBROUTINE FORCE
!
!
! ###################### Veolcity Initialization ##############################
!
! With given temperature, initial velocity is given to each particle.
! But distrubtion follows Gaussian(Normalized) random distribution
!
SUBROUTINE Vinit(NS, qel, qion, param, ttime, sys, dt)
USE DataStr
IMPLICIT NONE
INTERFACE
   FUNCTION fluct(x)
     IMPLICIT NONE
     REAL*8, INTENT(IN):: x
     REAL*8:: fluct, r, v1, v2
     REAL*8:: rand1, rand2, ran2
   END FUNCTION fluct
END INTERFACE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::qel(NS%Nel), qion(NS%Nion)
TYPE(PM)::param
TYPE(TM)::ttime
TYPE(ST)::sys
REAL*8  ::dt
!
INTEGER:: i, j, k, Nr
REAL*8 :: xm, lambda, T0, kb, xv(NS%Npt,3)
TYPE(PT)::gel(NS%Nel), gion(NS%Nion)
!
kb = 8.617343D-5
sys%mv2 = 0.0D0
Nr = Nref
!
! Electron velocity initialization
xm = param%xm(1)
DO i = 1, NS%Nel
   DO j=1, 3
      xv(i,j) = fluct(SQRT(param%Te))
      sys%mv2 = sys%mv2 + xm*xv(i,j)**2 
   END DO
END DO
T0 = sys%mv2/kb/REAL(3*NS%Nel - Nr)
lambda = SQRT(param%Te/T0)
DO i = 1, NS%Nel
   DO j = 1,3
      qel(i)%xv(j) = xv(i,j)*lambda      
      gel(i)%xx(j) = qel(i)%xx(j)
   END DO
END DO
!
! Ion velocity initialization
sys%mv2 = 0.0
xm = param%xm(2)
DO i = 1, NS%Nion
   k = i+NS%Nel
   DO j=1, 3
      xv(k,j) = fluct(SQRT(param%Ti))
      sys%mv2 = sys%mv2 + xm*xv(k,j)**2 
   END DO
END DO
T0 = sys%mv2/kb/REAL(3*NS%Nion - Nr)
lambda = SQRT(param%Ti/T0)
DO i = 1, NS%Nion
   k = i+NS%Nel
   DO j = 1,3
      qion(i)%xv(j) = xv(k,j)*lambda
      gion(i)%xx(j) = qion(i)%xx(j)
   END DO
END DO
sys%Te = param%Te
sys%Ti = param%Ti
sys%temp = (sys%kbm_e*param%Te+sys%kbm_i*param%Ti)/sys%kbm
!
CALL Force (NS, gel, gion, param, ttime)
DO i = 1, NS%Nel
   DO k = 1,3
      qel(i)%ff(k) = gel(i)%ff(k)
   END DO
END DO
DO i = 1, NS%Nion
   DO k = 1,3
      qion(i)%ff(k) = gion(i)%ff(k)
   END DO
END DO
!
RETURN
END SUBROUTINE Vinit
!
! ################### Initializing routine for restart #####################
! Initialize particle force. Velocity terms come from restart file and
! don't need to be initialized
! 
SUBROUTINE Vinit_REST(NS, qel, qion, param, ttime, sys, dt)
USE DATASTR
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM)::NS
TYPE(PT)::qel(NS%Nel), qion(NS%Nion), gel(NS%Nel), gion(NS%Nion)
TYPE(PM)::param
TYPE(TM)::ttime
TYPE(ST)::sys
REAL*8  ::dt
!
INTEGER:: i, k
!
!
DO i = 1, NS%Nel
   DO k = 1,3
      gel(i)%xx(k) = qel(i)%xx(k)
   END DO
END DO
DO i = 1, NS%Nion
   DO k = 1,3
      gion(i)%xx(k) = qion(i)%xx(k)
   END DO
END DO
CALL Force (NS, gel, gion, param, ttime)
DO i = 1, NS%Nel
   DO k = 1,3
      qel(i)%ff(k) = gel(i)%ff(k)
   END DO
END DO
DO i = 1, NS%Nion
   DO k = 1,3
      qion(i)%ff(k) = gion(i)%ff(k)
   END DO
END DO
!
RETURN
END SUBROUTINE Vinit_REST
!
END PROGRAM OCTREE_HYBRID
!
! ID - finder function ######################################################
! ALLOCATE processor id for particle position ###############################
FUNCTION id_FINDER(x, y, z)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: x, y, z
  INTEGER            :: id_FINDER
  INTEGER            :: i, j, k
  IF ( z < 0.0) THEN
     k = 1
  ELSE
     k = 2
  END IF
  !
  IF (y < 0.0) THEN
     j = 1
  ELSE
     j = 2
  END IF
  IF (x < 0.0) THEN
     i = 1
  ELSE
     i = 2
  END IF
  id_FINDER = 4*(k-1) + 2*(j-1) + i - 1
END FUNCTION id_FINDER
!
! Criterion for copying local TREE ###########################################
SUBROUTINE CriterionFind(id1, id2, Ncrit)
  USE DATASTR
  USE MPI_ENV
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: id1, id2
  INTEGER, INTENT(OUT) :: Ncrit
  INTEGER              :: i, N, Loctn(8,3)
  Loctn = RESHAPE ((/1,2, 1,2, 1,2, 1,2, 1,1, 2,2, 1,1, 2,2, &
       1,1,1,1, 2,2,2,2/), (/8,3/))
  N = INT(id1/Npe) + 1
  i = INT(id2/Npe) + 1
  !
  IF (Loctn(N,1) == Loctn(i,1)) THEN
     IF (Loctn(N,2) == Loctn(i,2)) THEN
        !
        ! FACE contact for x-y plane
        ! (i,j,:)
        Ncrit = 3
     ELSE
        IF (Loctn(N,3) == Loctn(i,3)) THEN
           !
           ! FACE contact for x-z plane
           ! (i,:,k)
           Ncrit = 2
        ELSE
           !
           ! LINE contact for x-axis           
           ! (i,:,:)
           Ncrit = 4
        END IF
     END IF
  ELSE 
     IF (Loctn(N,2) == Loctn(i,2)) THEN
        IF (Loctn(N,3) == Loctn(i,3)) THEN
           !
           ! FACE contact for y-z plane
           ! (:,j,k)
           Ncrit = 1
        ELSE
           !
           ! LINE contact for y-axis
           ! (:,j,:)
           Ncrit = 5
        END IF
     ELSE
        IF (Loctn(N,3) == Loctn(i,3)) THEN
           !
           ! LINE contact for z-axis
           ! (:,:,k)
           Ncrit = 6
        ELSE
           !
           ! POINT contact
           ! (:,:,:)
           Ncrit = 7
        END IF
     END IF
  END IF
  !
  RETURN
END SUBROUTINE CriterionFind
