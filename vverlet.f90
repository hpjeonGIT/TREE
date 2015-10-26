!
! ####### Pure Verlet time integration routine with No Thermostat #########
! 
! Purely integration routine with no thermostat. Basically, this routine is
! applied when new potential is implemented - Any potential should not increase
! the given temperature or increase temperature as iteration goes. 
! Also for NVE (microcanonical ensemble) simulation.
!
SUBROUTINE VVerletNotemp1(NS, qel, qion, param, dt, ttime)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM), INTENT(IN)::NS
TYPE(PT), INTENT(INOUT)::qel(NS%Nel), qion(NS%Nion)
TYPE(PM), INTENT(IN)::param
TYPE(TM), INTENT(IN)::ttime
REAL*8,   INTENT(IN)::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k
REAL*8 :: xm
!
! Electron update
xm = param%xm(1)
DO i = 1, NS%Nel
   DO k = 1, 3
      qel(i)%xv(k) = qel(i)%xv(k) + 0.5D0*dt*qel(i)%ff(k)/xm
      qel(i)%xx(k) = qel(i)%xx(k) + dt*qel(i)%xv(k)
   END DO
END DO
!
! Ion update
xm = param%xm(2)
DO i = 1, NS%Nion
   DO k = 1, 3
      qion(i)%xv(k) = qion(i)%xv(k) + 0.5D0*dt*qion(i)%ff(k)/xm
      qion(i)%xx(k) = qion(i)%xx(k) + dt*qion(i)%xv(k)
   END DO
END DO
!
RETURN
!
END SUBROUTINE VVerletNotemp1
!
! ####### Pure Verlet time integration routine for second stage #########
! #######################################################################
!
SUBROUTINE VVerletNotemp2(NS, qel, qion, param, sys, dt, ttime)
USE DataStr
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM), INTENT(IN)::NS
TYPE(PT), INTENT(INOUT)::qel(NS%Nel), qion(Ns%Nion)
TYPE(PM), INTENT(IN)::param
TYPE(ST), INTENT(INOUT)::sys
TYPE(TM), INTENT(IN)::ttime
REAL*8,   INTENT(IN)::dt
!
! INTERNAL VARIABLES
INTEGER:: i, k
REAL*8 :: xm
!
! Electron update
xm = param%xm(1)
DO i = 1, NS%Nel
   DO k = 1, 3
      qel(i)%xv(k) = qel(i)%xv(k) + 0.5D0*dt*qel(i)%ff(k)/xm
   END DO
END DO
!
! Ion update
! Be careful on MOD number criterion - not same to verlet1
xm = param%xm(2)
DO i = 1, NS%Nion
   DO k = 1, 3
      qion(i)%xv(k) = qion(i)%xv(k) + 0.5D0*dt*qion(i)%ff(k)/xm
   END DO
END DO
IF (MOD(ttime%Nloop,ttime%Nsamp) == 0) THEN
   sys%mv2el  = 0.0D0
   sys%mv2ion = 0.0D0
   sys%Te_rad = 0.0D0
   sys%Ti_rad = 0.0D0
   !  
   ! Electron update
   xm = param%xm(1)
   DO i=1, NS%Nel
      sys%mv2el = sys%mv2el + &
           xm*(qel(i)%xv(1)**2 + qel(i)%xv(2)**2 + qel(i)%xv(3)**2)
      sys%Te_rad = sys%Te_rad + xm*(qel(i)%xx(1)*qel(i)%xv(1) + &
           qel(i)%xx(2)*qel(i)%xv(2) + qel(i)%xx(3)*qel(i)%xv(3))**2/ &
           (qel(i)%xx(1)**2 + qel(i)%xx(2)**2 + qel(i)%xx(3)**2)
   END DO
   !
   ! Ion update
   xm = param%xm(2)
   DO i = 1, NS%Nion
      sys%mv2ion = sys%mv2ion + &
           xm*(qion(i)%xv(1)**2 + qion(i)%xv(2)**2 + qion(i)%xv(3)**2)
      sys%Ti_rad = sys%Ti_rad + xm*(qion(i)%xx(1)*qion(i)%xv(1) + &
           qion(i)%xx(2)*qion(i)%xv(2) + qion(i)%xx(3)*qion(i)%xv(3))**2/ &
           (qion(i)%xx(1)**2 + qion(i)%xx(2)**2 + qion(i)%xx(3)**2)
   END DO
END IF
!
sys%mv2 = sys%mv2ion+sys%mv2el
!
RETURN
!
END SUBROUTINE VVerletNotemp2
!
! ############## Random Gaussian(Normal) Distribution Function ################
! 
! For stochastic thermostat, fluctuation dissipation theorem is implemented.
! Basically, random number generator which follows Gaussian distribution is
! needed to implement this thermal noise force.
! Random number is generated using FORTRAN intrinsic fucntion - RANDOM_SEED
! and RANDOM_NUMBER. But during the implementation, it is found that those
! intrinsic functions may not work well under false seeding - to use this
! routine on new machine or new compiler, please make sure that the 
! distribution follows zero-mean with suitable deviation.
!
! This function provides a random number with zero-mean and deviation x 
! along Gaussian distribution.
! <p> = 0.0
! <p**2> = x**2
! 
FUNCTION fluct(x)
  IMPLICIT NONE
  REAL*8, INTENT(IN):: x
  REAL*8:: fluct, r, v1, v2
  REAL*8:: rand1, rand2, ran2
  !
  ! Initialization
  r=1.
  DO WHILE (r.ge.1.D0)
     CALL RANDOM_NUMBER(rand1)
     CALL RANDOM_NUMBER(rand2)
     v1 = 2.D0*rand1 - 1.D0
     v2 = 2.D0*rand2 - 1.D0
     r = v1*v1+v2*v2
  END DO
  fluct = v1*SQRT(-2.D0*log(r)/r)*x
  RETURN
END FUNCTION fluct
