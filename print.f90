!
! ##################### Binary restart file generator  ########################
!
!
SUBROUTINE Restart(NS, qel, qion, ttime, COMM)
USE DATASTR
USE MPI_ENV
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(TM):: ttime
TYPE(PT):: qel(NS%Nel), qion(NS%Nion)
TYPE(MP):: COMM
!
!
INTEGER:: i, j, k, Nfreq, ccnt(8*Npe), displ(8*Npe), IERR, N, NTel, NTion, &
     id_ion(NS%NTpt), id_el(NS%NTpt), id_ion_l(NS%Nion), id_el_l(NS%Nel), &
     idspl(8*Npe)
REAL*8 :: ixx(NS%Nion*3), ixv(NS%Nion*3), exx(NS%Nel*3), exv(NS%Nel*3), &
     xx(NS%NTpt*3), xv(NS%NTpt*3)
CHARACTER(LEN=20):: FILENAME
!
! Dump position and velocity into local variables
k = 0
DO i=1, NS%Nel
   DO j=1, 3
      k = k + 1
      exx(k) = qel(i)%xx(j)
      exv(k) = qel(i)%xv(j)
   END DO
   id_el_l(i) = qel(i)%id
END DO
k = 0
DO i=1, NS%Nion
   DO j=1, 3
      k = k + 1
      ixx(k) = qion(i)%xx(j)
      ixv(k) = qion(i)%xv(j)
   END DO
   id_ion_l(i) = qion(i)%id
END DO
!
! Electron data dump
N = NS%Nel*3
CALL MPI_Allgather(N, 1, MPI_INTEGER, ccnt, 1, MPI_INTEGER, COMM%COMM1D, IERR)
displ(1) = 0
idspl(1) = 0
NTel = ccnt(8*Npe)
DO i=1, 8*Npe-1
   displ(i+1) = ccnt(i) + displ(i)
   idspl(i+1) = ccnt(i)/3 + idspl(i)
   NTel = NTel + ccnt(i)
END DO
NTel = NTel/3
CALL MPI_GatherV(id_el_l, NS%Nel, MPI_INTEGER, id_el, ccnt/3, idspl, &
     MPI_INTEGER, 0, COMM%COMM1D, IERR)
CALL MPI_GatherV(exx, N, MPI_REAL8, xx, ccnt, displ, MPI_REAL8, 0, &
     COMM%COMM1D, IERR)
CALL MPI_GatherV(exv, N, MPI_REAL8, xv, ccnt, displ, MPI_REAL8, 0, &
     COMM%COMM1D, IERR)
IF (COMM%MYID2 == 0) THEN
   !
   ! File name decision
   Nfreq = INT(ttime%Nloop/ttime%Nrest)
   WRITE(FILENAME,55) Nfreq
   OPEN(UNIT=30,file=FILENAME)
   WRITE(30,*) NS%NTpt
   WRITE(30,*) "frame = ", Nfreq, "  energy = 0"
   k = 1
   DO i=1, NTel
      WRITE(30, 80) xx(k), xx(k+1), xx(k+2), xv(k), xv(k+1), xv(k+2), id_el(i)
      k = k + 3
   END DO
END IF
55 FORMAT("res", I3.3, ".xyz")
80 FORMAT("P  ", 6(ES14.6, 1X), I8)
!
! Ion data dump
N = NS%Nion*3
CALL MPI_Allgather(N, 1, MPI_INTEGER, ccnt, 1, MPI_INTEGER, COMM%COMM1D, IERR)
displ(1) = 0
idspl(1) = 0
NTion = ccnt(Npe*8)
DO i=1, 8*Npe-1
   displ(i+1) = ccnt(i) + displ(i)
   idspl(i+1) = ccnt(i)/3 + idspl(i)
   NTion = NTion + ccnt(i)
END DO
NTion = NTion/3
CALL MPI_GatherV(id_ion_l, NS%Nion, MPI_INTEGER, id_ion, ccnt/3, idspl, &
     MPI_INTEGER, 0, COMM%COMM1D, IERR)
CALL MPI_GatherV(ixx, N, MPI_REAL8, xx, ccnt, displ, MPI_REAL8, 0, &
     COMM%COMM1D, IERR)
CALL MPI_GatherV(ixv, N, MPI_REAL8, xv, ccnt, displ, MPI_REAL8, 0, &
     COMM%COMM1D, IERR)
IF (COMM%MYID2 == 0) THEN
   k = 1
   DO i=1, NTion
      WRITE(30, 90) xx(k), xx(k+1), xx(k+2), xv(k), xv(k+1), xv(k+2), id_ion(i)
      k = k + 3
   END DO
   CLOSE(30)
END IF
90 FORMAT("O  ", 6(ES14.6, 1X), I8)
!
100 RETURN
END SUBROUTINE Restart
