!
! ##################### Binary restart file generator  ########################
!
!
SUBROUTINE RDF(NS, ttime, COMM)
USE DATASTR
USE MPI_ENV
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM), INTENT(INOUT):: NS
TYPE(TM), INTENT(IN):: ttime
TYPE(MP), INTENT(IN):: COMM
!
!
INTEGER:: i, j, Nfreq, IERR
REAL*8 :: dv, drr, dr, dxx, dx, Ntau, Nei(Ndist+1), Nii(Ndist+1), dx_temp
CHARACTER(LEN=20):: FILENAME
!
! RDF sum for e-i
CALL MPI_REDUCE(NS%RDF_ei, Nei, Ndist+1, MPI_REAL8, MPI_SUM, 0, &
     COMM%COMM1D, IERR)
CALL MPI_REDUCE(NS%RDF_ii, Nii, Ndist+1, MPI_REAL8, MPI_SUM, 0, &
     COMM%COMM1D, IERR)
!
Ntau = REAL(ttime%Ndump)
IF (COMM%MYID2 == 0) THEN
   !
   ! File name decision
   Nfreq = INT(ttime%Nloop/ttime%Ndump)
   WRITE(FILENAME,50) Nfreq
   OPEN(UNIT=30,file=FILENAME)
   WRITE(30,100)
   DO i=1, Ndist+1
      WRITE(30,200) NS%dx_ei*REAL(i), Nei(i)/Ntau, &
           NS%dx_ii*REAL(i), Nii(i)/Ntau
   END DO
   CLOSE(30)
END IF
50  FORMAT("rdf",I3.3,".dat")
100 FORMAT("# R_ei,   RDF_ei       R_ii     RDF_ii     scaled RDF_ii")
200 FORMAT(5(ES14.6, 1X))
!
! Initialize
NS%RDF_ii = 0.D0
NS%RDF_ei = 0.D0
!
RETURN
END SUBROUTINE RDF
!
! ##################### Binary restart file generator  ########################
!
!
SUBROUTINE TEMPORAL(NS, sys, ttime, COMM)
USE DATASTR
USE MPI_ENV
IMPLICIT NONE
INCLUDE 'mpif.h'
!
! COMMUNICATION VARIABLES
TYPE(NM), INTENT(INOUT):: NS
TYPE(ST), INTENT(INOUT):: sys
TYPE(TM), INTENT(IN):: ttime
TYPE(MP), INTENT(IN):: COMM
!
!
INTEGER:: IERR
REAL*8 :: Epot, Tion, Tel, Te_rad, Ti_rad, Te_cir, Ti_cir
!   
Epot = 0.0D0
Tel = 0.0D0
Tion = 0.0D0
Te_rad = 0.0D0
Ti_rad = 0.0D0
!
sys%Epot = sys%Epot_e + sys%Epot_i
CALL MPI_REDUCE(sys%Epot, Epot, 1, MPI_REAL8, MPI_SUM, 0, COMM%COMM1D, IERR)
CALL MPI_REDUCE(sys%mv2el, Tel, 1, MPI_REAL8, MPI_SUM, 0, COMM%COMM1D, IERR)
CALL MPI_REDUCE(sys%mv2ion, Tion, 1, MPI_REAL8, MPI_SUM, 0, COMM%COMM1D, IERR)
CALL MPI_REDUCE(sys%Te_rad,Te_rad, 1, MPI_REAL8, MPI_SUM, 0, COMM%COMM1D, IERR)
CALL MPI_REDUCE(sys%Ti_rad,Ti_rad, 1, MPI_REAL8, MPI_SUM, 0, COMM%COMM1D, IERR)
sys%Te = Tel/sys%kbm_e
sys%Ti = Tion/sys%kbm_i
Te_cir = (Tel - Te_rad)/sys%kbm_e_cir
Ti_cir = (Tion - Ti_rad)/sys%kbm_i_cir
Te_rad = Te_rad/sys%kbm_e_rad
Ti_rad = Ti_rad/sys%kbm_i_rad
IF ( COMM%MYID2 == 0) THEN
   WRITE(25,500)  ttime%tnow*tps, sys%Te, sys%Ti, &
        Te_rad, Ti_rad, Te_cir, Ti_cir, Epot, 0.5D0*(Tel+Tion) + Epot
END IF
500 FORMAT(9(ES14.6, 1x))
!
RETURN
END SUBROUTINE TEMPORAL
