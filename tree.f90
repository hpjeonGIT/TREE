MODULE BINTREE
IMPLICIT NONE
!
TYPE BT
   !
   ! Binary Tree (local)
   ! 
   ! ileaf = tag of leaf (1: ion, 0: twig, -1: electron)
   ! pt = tag for particle number
   ! icut = how to do binary search
   ! xx = coordinate of search box
   ! lx = length of search box
   ! qeff = effective charge of node
   ! dx = distance between box center and charge center
   INTEGER:: ileaf, pt
   REAL*8 :: xx(3), lx(3), qeff, dx, D(3), Q(3), QQ(3)
   TYPE (BT), POINTER:: Child1=> null(), Child2 => null(), &
        Child3 => null(), Child4 => null(), Child5 => null(), &
        Child6 => null(), Child7 => null(), Child8 => null(), &
        Root => null()
END TYPE BT
!
TYPE GT
   !
   ! Ghost Tree 
   !
   INTEGER:: ileaf
   REAL*8 :: xx(3), lx, qeff, dx, D(3), Q(3), QQ(3)
   TYPE (GT), POINTER:: Child1=> null(), Child2 => null(), &
        Child3 => null(), Child4 => null(), Child5 => null(), &
        Child6 => null(), Child7 => null(), Child8 => null(), &
        Root => null()
END TYPE GT
!
CONTAINS
!
!
! ###########################################################
SUBROUTINE BuildLocalTree(NS, gel, gion, Roote, Rooti, ttime)
  USE DATASTR
  IMPLICIT NONE
  TYPE(NM), INTENT(IN):: NS
  TYPE(GS), INTENT(IN):: gel(NS%NTel), gion(NS%NTion)
  TYPE(TM), INTENT(IN):: ttime
  TYPE(BT), Pointer   :: Roote, Rooti, Loce, Loci
  INTEGER :: i, bias
  REAL*8  :: minx, miny, minz, maxx, maxy, maxz
  !
  NULLIFY(Loce)
  ALLOCATE(Loce)
  Roote%Root => Loce
  minx = MINVAL(gel(:)%xx(1))
  maxx = MAXVAL(gel(:)%xx(1))
  miny = MINVAL(gel(:)%xx(2))
  maxy = MAXVAL(gel(:)%xx(2))
  minz = MINVAL(gel(:)%xx(3))
  maxz = MAXVAL(gel(:)%xx(3))
  Roote%xx(1) = (maxx+minx)/2.D0
  Roote%xx(2) = (maxy+miny)/2.D0
  Roote%xx(3) = (maxz+minz)/2.D0
  Roote%lx(1) = maxx-minx
  Roote%lx(2) = maxy-miny
  Roote%lx(3) = maxz-minz
  Roote%ileaf = 0
  Roote%dx    = 0.0D0
  bias = 1
  !
  ! Electron Tree build
  DO i=1, NS%NTel
     CALL Insert_Tree_e(Roote, gel, Loce, i, bias, NS)
  END DO
  CALL PSeudoE(Roote, gel, NS)
  DEALLOCATE(Loce)
  !
  ! Ion Tree build
  NULLIFY(Loci)
  ALLOCATE(Loci)
  Rooti%Root => Loci
  minx = MINVAL(gion(:)%xx(1))
  maxx = MAXVAL(gion(:)%xx(1))
  miny = MINVAL(gion(:)%xx(2))
  maxy = MAXVAL(gion(:)%xx(2))
  minz = MINVAL(gion(:)%xx(3))
  maxz = MAXVAL(gion(:)%xx(3))
  Rooti%xx(1) = (maxx+minx)/2.D0
  Rooti%xx(2) = (maxy+miny)/2.D0
  Rooti%xx(3) = (maxz+minz)/2.D0
  Rooti%lx(1) = maxx-minx
  Rooti%lx(2) = maxy-miny
  Rooti%lx(3) = maxz-minz
  Rooti%ileaf = 0
  Rooti%dx    = 0.0D0
  bias = 1
  DO i=1, NS%NTion
     CALL Insert_Tree_i(Rooti, gion, Loci, i, bias, NS)
  END DO
  CALL PSeudoI(Rooti, gion, NS)
  DEALLOCATE(Loci)
  RETURN
END SUBROUTINE BuildLocalTree
!
! ################ e- TREE build #########################################
RECURSIVE SUBROUTINE Insert_Tree_e(P, gel, R, n, bias, NS)
  USE DATASTR
  TYPE(NM), INTENT(IN):: NS
  TYPE(GS), INTENT(IN):: gel(NS%NTel)
  INTEGER, INTENT(IN) :: n
  TYPE(BT), POINTER:: P, R
  INTEGER             :: m, bias
  REAL*8 :: ijk(3,8)
  ! 
  IF (.NOT. ASSOCIATED(P)) THEN     
     ALLOCATE(P)
     ijk = RESHAPE((/ -.25D0,-.25D0,-.25D0,  +.25D0,-.25D0,-.25D0, &
          & -.25D0,+.25D0,-.25D0,  +.25D0,+.25D0,-.25D0, &
          & -.25D0,-.25D0,+.25D0,  +.25D0,-.25D0,+.25D0, &
          & -.25D0,+.25D0,+.25D0, +.25D0, +.25D0, +.25D0 /), (/3,8/))
     P%Root => R
     P%pt = n
     P%ileaf = -1
     P%qeff = -1.D0
     P%dx = 0.0D0
     P%xx(1) = R%xx(1) + ijk(1,bias)*R%lx(1)
     P%xx(2) = R%xx(2) + ijk(2,bias)*R%lx(2)
     P%xx(3) = R%xx(3) + ijk(3,bias)*R%lx(3)
     P%lx(:) = R%lx(:)*0.5D0
  ELSE
     IF (P%ileaf == 0) THEN
        !
        ! twig
        P%qeff = P%qeff - 1.D0
        IF ( gel(n)%xx(3) < P%xx(3) ) THEN
           ! 1-4
           IF ( gel(n)%xx(2) < P%xx(2) ) THEN
              !
              ! 1-2
              IF ( gel(n)%xx(1) < P%xx(1) ) THEN
                 bias = 1
                 CALL Insert_Tree_e(P%Child1, gel, P, n, bias, NS)
              ELSE
                 bias = 2
                 CALL Insert_Tree_e(P%Child2, gel, P, n, bias, NS)
              END IF
           ELSE
              IF ( gel(n)%xx(1) < P%xx(1) ) THEN
                 bias = 3
                 CALL Insert_Tree_e(P%Child3, gel, P, n, bias, NS)
              ELSE
                 bias = 4
                 CALL Insert_Tree_e(P%Child4, gel, P, n, bias, NS)
              END IF
           END IF
        ELSE
           ! 5-8
           IF ( gel(n)%xx(2) < P%xx(2) ) THEN
              !
              ! 5-6
              IF ( gel(n)%xx(1) < P%xx(1) ) THEN
                 bias = 5
                 CALL Insert_Tree_e(P%Child5, gel, P, n, bias, NS)
              ELSE
                 bias = 6
                 CALL Insert_Tree_e(P%Child6, gel, P, n, bias, NS)
              END IF
           ELSE
              IF ( gel(n)%xx(1) < P%xx(1) ) THEN
                 bias = 7
                 CALL Insert_Tree_e(P%Child7, gel, P, n, bias, NS)
              ELSE
                 bias = 8
                 CALL Insert_Tree_e(P%Child8, gel, P, n, bias, NS)
              END IF
           END IF
        END IF
     ELSE 
        !
        ! has been leaf of electron
        P%ileaf = 0
        m = P%pt
        P%pt = 0
        P%qeff = P%qeff - 1.D0
        !
        ! Former resident particle
        IF ( gel(m)%xx(3) < P%xx(3) ) THEN
           ! 1-4
           IF ( gel(m)%xx(2) < P%xx(2) ) THEN
              !
              ! 1-2
              IF ( gel(m)%xx(1) < P%xx(1) ) THEN
                 bias = 1
                 CALL Insert_Tree_e(P%Child1, gel, P, m, bias, NS)
              ELSE
                 bias = 2
                 CALL Insert_Tree_e(P%Child2, gel, P, m, bias, NS)
              END IF
           ELSE
              IF ( gel(m)%xx(1) < P%xx(1) ) THEN
                 bias = 3
                 CALL Insert_Tree_e(P%Child3, gel, P, m, bias, NS)
              ELSE
                 bias = 4
                 CALL Insert_Tree_e(P%Child4, gel, P, m, bias, NS)
              END IF
           END IF
        ELSE
           ! 5-8
           IF ( gel(m)%xx(2) < P%xx(2) ) THEN
              !
              ! 5-6
              IF ( gel(m)%xx(1) < P%xx(1) ) THEN
                 bias = 5
                 CALL Insert_Tree_e(P%Child5, gel, P, m, bias, NS)
              ELSE
                 bias = 6
                 CALL Insert_Tree_e(P%Child6, gel, P, m, bias, NS)
              END IF
           ELSE
              IF ( gel(m)%xx(1) < P%xx(1) ) THEN
                 bias = 7
                 CALL Insert_Tree_e(P%Child7, gel, P, m, bias, NS)
              ELSE
                 bias = 8
                 CALL Insert_Tree_e(P%Child8, gel, P, m, bias, NS)
              END IF
           END IF
        END IF
        !
        ! Newly coming particle
        IF ( gel(n)%xx(3) < P%xx(3) ) THEN
           ! 1-4
           IF ( gel(n)%xx(2) < P%xx(2) ) THEN
              !
              ! 1-2
              IF ( gel(n)%xx(1) < P%xx(1) ) THEN
                 bias = 1
                 CALL Insert_Tree_e(P%Child1, gel, P, n, bias, NS)
              ELSE
                 bias = 2
                 CALL Insert_Tree_e(P%Child2, gel, P, n, bias, NS)
              END IF
           ELSE
              IF ( gel(n)%xx(1) < P%xx(1) ) THEN
                 bias = 3
                 CALL Insert_Tree_e(P%Child3, gel, P, n, bias, NS)
              ELSE
                 bias = 4
                 CALL Insert_Tree_e(P%Child4, gel, P, n, bias, NS)
              END IF
           END IF
        ELSE
           ! 5-8
           IF ( gel(n)%xx(2) < P%xx(2) ) THEN
              !
              ! 5-6
              IF ( gel(n)%xx(1) < P%xx(1) ) THEN
                 bias = 5
                 CALL Insert_Tree_e(P%Child5, gel, P, n, bias, NS)
              ELSE
                 bias = 6
                 CALL Insert_Tree_e(P%Child6, gel, P, n, bias, NS)
              END IF
           ELSE
              IF ( gel(n)%xx(1) < P%xx(1) ) THEN
                 bias = 7
                 CALL Insert_Tree_e(P%Child7, gel, P, n, bias, NS)
              ELSE
                 bias = 8
                 CALL Insert_Tree_e(P%Child8, gel, P, n, bias, NS)
              END IF
           END IF
        END IF
     END IF
  END IF
END SUBROUTINE Insert_Tree_e
!
! ################ ion TREE build ########################################
RECURSIVE SUBROUTINE Insert_Tree_i(P, gion, R, n, bias, NS)
  USE DATASTR
  TYPE(NM), INTENT(IN):: NS
  TYPE(GS), INTENT(IN):: gion(NS%NTion)
  INTEGER,  INTENT(IN):: n
  TYPE(BT), POINTER   :: P, R
  INTEGER             :: m, bias
  REAL*8 :: ijk(3,8)
  IF (.NOT. ASSOCIATED(P)) THEN
     ALLOCATE(P)
     ijk = RESHAPE((/ -.25D0,-.25D0,-.25D0,  +.25D0,-.25D0,-.25D0, &
          & -.25D0,+.25D0,-.25D0, +.25D0,+.25D0,-.25D0,  &
          & -.25D0,-.25D0,+.25D0,  +.25D0,-.25D0,+.25D0, &
          & -.25D0,+.25D0,+.25D0, +.25D0, +.25D0, +.25D0 /), (/3,8/))
     P%Root => R
     P%pt = n
     P%ileaf = 1
     P%qeff = +1.D0    
     P%dx = 0.0D0
     P%xx(1) = R%xx(1) + ijk(1,bias)*R%lx(1)
     P%xx(2) = R%xx(2) + ijk(2,bias)*R%lx(2)
     P%xx(3) = R%xx(3) + ijk(3,bias)*R%lx(3)
     P%lx(:) = R%lx(:)*0.5D0
  ELSE
     IF (P%ileaf == 0) THEN
        !
        ! twig
        P%qeff = P%qeff + 1.D0
        IF ( gion(n)%xx(3) < P%xx(3)) THEN
           ! 1-4
           IF ( gion(n)%xx(2) < P%xx(2) ) THEN
              !
              ! 1-2
              IF ( gion(n)%xx(1) < P%xx(1) ) THEN
                 bias = 1
                 CALL Insert_Tree_i(P%Child1, gion, P, n, bias, NS)
              ELSE
                 bias = 2
                 CALL Insert_Tree_i(P%Child2, gion, P, n, bias, NS)
              END IF
           ELSE
              IF ( gion(n)%xx(1) < P%xx(1) ) THEN
                 bias = 3
                 CALL Insert_Tree_i(P%Child3, gion, P, n, bias, NS)
              ELSE
                 bias = 4
                 CALL Insert_Tree_i(P%Child4, gion, P, n, bias, NS)
              END IF
           END IF
        ELSE
           ! 5-8
           IF ( gion(n)%xx(2) < P%xx(2) ) THEN
              !
              ! 5-6
              IF ( gion(n)%xx(1) < P%xx(1) ) THEN
                 bias = 5
                 CALL Insert_Tree_i(P%Child5, gion, P, n, bias, NS)
              ELSE
                 bias = 6
                 CALL Insert_Tree_i(P%Child6, gion, P, n, bias, NS)
              END IF
           ELSE
              IF ( gion(n)%xx(1) < P%xx(1) ) THEN
                 bias = 7
                 CALL Insert_Tree_i(P%Child7, gion, P, n, bias, NS)
              ELSE
                 bias = 8
                 CALL Insert_Tree_i(P%Child8, gion, P, n, bias, NS)
              END IF
           END IF
        END IF
     ELSE
        !
        ! has been leaf of ion
        P%ileaf = 0
        m = P%pt
        P%pt = 0
        P%qeff = P%qeff + 1.D0
        !
        ! Former resident particle(ion)
        IF ( gion(m)%xx(3) < P%xx(3)) THEN
           ! 1-4
           IF ( gion(m)%xx(2) < P%xx(2) ) THEN
              !
              ! 1-2
              IF ( gion(m)%xx(1) < P%xx(1) ) THEN
                 bias = 1
                 CALL Insert_Tree_i(P%Child1, gion, P, m, bias, NS)
              ELSE
                 bias = 2
                 CALL Insert_Tree_i(P%Child2, gion, P, m, bias, NS)
              END IF
           ELSE
              IF ( gion(m)%xx(1) < P%xx(1) ) THEN
                 bias = 3
                 CALL Insert_Tree_i(P%Child3, gion, P, m, bias, NS)
              ELSE
                 bias = 4
                 CALL Insert_Tree_i(P%Child4, gion, P, m, bias, NS)
              END IF
           END IF
        ELSE
           ! 5-8
           IF ( gion(m)%xx(2) < P%xx(2) ) THEN
              !
              ! 5-6
              IF ( gion(m)%xx(1) < P%xx(1) ) THEN
                 bias = 5
                 CALL Insert_Tree_i(P%Child5, gion, P, m, bias, NS)
              ELSE
                 bias = 6
                 CALL Insert_Tree_i(P%Child6, gion, P, m, bias, NS)
              END IF
           ELSE
              IF ( gion(m)%xx(1) < P%xx(1) ) THEN
                 bias = 7
                 CALL Insert_Tree_i(P%Child7, gion, P, m, bias, NS)
              ELSE
                 bias = 8
                 CALL Insert_Tree_i(P%Child8, gion, P, m, bias, NS)
              END IF
           END IF
        END IF
        !
        ! Newly coming particle
        IF ( gion(n)%xx(3) < P%xx(3)) THEN
           ! 1-4
           IF ( gion(n)%xx(2) < P%xx(2) ) THEN
              !
              ! 1-2
              IF ( gion(n)%xx(1) < P%xx(1) ) THEN
                 bias = 1
                 CALL Insert_Tree_i(P%Child1, gion, P, n, bias, NS)
              ELSE
                 bias = 2
                 CALL Insert_Tree_i(P%Child2, gion, P, n, bias, NS)
              END IF
           ELSE
              IF ( gion(n)%xx(1) < P%xx(1) ) THEN
                 bias = 3
                 CALL Insert_Tree_i(P%Child3, gion, P, n, bias, NS)
              ELSE
                 bias =4
                 CALL Insert_Tree_i(P%Child4, gion, P, n, bias, NS)
              END IF
           END IF
        ELSE
           ! 5-8
           IF ( gion(n)%xx(2) < P%xx(2) ) THEN
              !
              ! 5-6
              IF ( gion(n)%xx(1) < P%xx(1) ) THEN
                 bias = 5
                 CALL Insert_Tree_i(P%Child5, gion, P, n, bias, NS)
              ELSE
                 bias = 6
                 CALL Insert_Tree_i(P%Child6, gion, P, n, bias, NS)
              END IF
           ELSE
              IF ( gion(n)%xx(1) < P%xx(1) ) THEN
                 bias = 7
                 CALL Insert_Tree_i(P%Child7, gion, P, n, bias, NS)
              ELSE
                 bias = 8
                 CALL Insert_Tree_i(P%Child8, gion, P, n, bias, NS)
              END IF
           END IF
        END IF
     END IF
  END IF
END SUBROUTINE Insert_Tree_i
!
! #########################
RECURSIVE SUBROUTINE PSeudoE(P, gel, NS)
  USE DATASTR
  TYPE(NM), INTENT(IN):: NS
  TYPE(GS), INTENT(IN):: gel(NS%NTel)
  TYPE(BT), POINTER   :: P
  REAL*8  :: qsum, xx(3,8), ox(3), q(8)
  LOGICAL :: tag(8)
  INTEGER :: i
  !
  IF  (P%ileaf < 0) THEN 
     ! electron
     P%xx(:) = gel(P%pt)%xx(:)
     P%D(:)  = 0.D0
     P%Q(:)  = 0.D0
     P%QQ(:) = 0.D0
  ELSE 
     tag = .FALSE.
     P%D(:) = 0.D0
     P%Q(:) = 0.D0
     P%QQ(:) = 0.D0
     qsum = 0.D0
     ox(:) = P%xx(:)
     P%xx(:) = 0.0D0
     IF (ASSOCIATED(P%Child1)) THEN
        CALL PSeudoE(P%Child1, gel, NS)
        qsum    = qsum    + P%Child1%qeff
        xx(:,1) = P%Child1%xx(:)
        P%xx(:) = P%xx(:) + P%Child1%qeff*P%Child1%xx(:)
        tag(1)  = .TRUE.
        q(1)    = P%Child1%qeff
     END IF
     IF (ASSOCIATED(P%Child2)) THEN
        CALL PSeudoE(P%Child2, gel, NS)
        qsum    = qsum    + P%Child2%qeff
        xx(:,2) = P%Child2%xx(:)
        P%xx(:) = P%xx(:) + P%Child2%qeff*P%Child2%xx(:)
        tag(2)  = .TRUE.
        q(2)    = P%Child2%qeff
     END IF
     IF (ASSOCIATED(P%Child3)) THEN
        CALL PSeudoE(P%Child3, gel, NS)
        qsum    = qsum    + P%Child3%qeff
        xx(:,3) = P%Child3%xx(:)
        P%xx(:) = P%xx(:) + P%Child3%qeff*P%Child3%xx(:)
        tag(3)  = .TRUE.
        q(3)    = P%Child3%qeff
     END IF
     IF (ASSOCIATED(P%Child4)) THEN
        CALL PSeudoE(P%Child4, gel, NS)
        qsum    = qsum    + P%Child4%qeff
        xx(:,4) = P%Child4%xx(:)
        P%xx(:) = P%xx(:) + P%Child4%qeff*P%Child4%xx(:)
        tag(4)  = .TRUE.
        q(4)    = P%Child4%qeff
     END IF
     IF (ASSOCIATED(P%Child5)) THEN
        CALL PSeudoE(P%Child5, gel, NS)
        qsum    = qsum    + P%Child5%qeff
        xx(:,5) = P%Child5%xx(:)
        P%xx(:) = P%xx(:) + P%Child5%qeff*P%Child5%xx(:)
        tag(5)  = .TRUE.
        q(5)    = P%Child5%qeff
     END IF
     IF (ASSOCIATED(P%Child6)) THEN
        CALL PSeudoE(P%Child6, gel, NS)
        qsum    = qsum    + P%Child6%qeff
        xx(:,6) = P%Child6%xx(:)
        P%xx(:) = P%xx(:) + P%Child6%qeff*P%Child6%xx(:)
        tag(6)  = .TRUE.
        q(6)    = P%Child6%qeff
     END IF
     IF (ASSOCIATED(P%Child7)) THEN
        CALL PSeudoE(P%Child7, gel, NS)
        qsum    = qsum    + P%Child7%qeff
        xx(:,7) = P%Child7%xx(:)
        P%xx(:) = P%xx(:) + P%Child7%qeff*P%Child7%xx(:)
        tag(7)  = .TRUE.
        q(7)    = P%Child7%qeff
     END IF
     IF (ASSOCIATED(P%Child8)) THEN
        CALL PSeudoE(P%Child8, gel, NS)
        qsum    = qsum    + P%Child8%qeff
        xx(:,8) = P%Child8%xx(:)
        P%xx(:) = P%xx(:) + P%Child8%qeff*P%Child8%xx(:)
        tag(8)  = .TRUE.
        q(8)    = P%Child8%qeff
     END IF
     P%xx(:) = P%xx(:)/qsum
     P%dx = DSQRT((P%xx(1)-ox(1))**2 + (P%xx(2)-ox(2))**2 + (P%xx(3)-ox(3))**2)
     DO i=1,8
        IF (tag(i)) THEN
           P%D(1)  = P%D(1)  + q(i)*(xx(1,i) - P%xx(1))
           P%D(2)  = P%D(2)  + q(i)*(xx(2,i) - P%xx(2))
           P%D(3)  = P%D(3)  + q(i)*(xx(3,i) - P%xx(3))
           P%Q(1)  = P%Q(1)  + q(i)*(xx(1,i) - P%xx(1))**2/2.D0
           P%Q(2)  = P%Q(2)  + q(i)*(xx(2,i) - P%xx(2))**2/2.D0
           P%Q(3)  = P%Q(3)  + q(i)*(xx(3,i) - P%xx(3))**2/2.D0
           P%QQ(1) = P%QQ(1) + q(i)*(xx(1,i) - P%xx(1))*(xx(2,i) - P%xx(2))
           P%QQ(2) = P%QQ(2) + q(i)*(xx(2,i) - P%xx(2))*(xx(3,i) - P%xx(3))
           P%QQ(3) = P%QQ(3) + q(i)*(xx(3,i) - P%xx(3))*(xx(1,i) - P%xx(1))
        END IF
     END DO
  END IF
END SUBROUTINE PSeudoE
!
! #########################
RECURSIVE SUBROUTINE PSeudoI(P, gion, NS)
  USE DATASTR
  TYPE(NM), INTENT(IN):: NS
  TYPE(GS), INTENT(IN):: gion(NS%NTion)
  TYPE(BT), POINTER   :: P
  REAL*8  :: qsum, xx(3,8), ox(3), q(8)
  LOGICAL :: tag(8)
  INTEGER :: i
  !
  IF (P%ileaf > 0) THEN 
     ! ion
     P%xx(:) = gion(P%pt)%xx(:)
     P%D(:)  = 0.D0
     P%Q(:)  = 0.D0
     P%QQ(:) = 0.D0
  ELSE 
     tag = .FALSE.
     P%D(:) = 0.D0
     P%Q(:) = 0.D0
     P%QQ(:) = 0.D0
     qsum = 0.D0
     ox(:) = P%xx(:)
     P%xx(:) = 0.0D0
     IF (ASSOCIATED(P%Child1)) THEN
        CALL PSeudoI(P%Child1, gion, NS)
        qsum    = qsum    + P%Child1%qeff
        xx(:,1) = P%Child1%xx(:)
        P%xx(:) = P%xx(:) + P%Child1%qeff*P%Child1%xx(:)
        tag(1)  = .TRUE.
        q(1)    = P%Child1%qeff
     END IF
     IF (ASSOCIATED(P%Child2)) THEN
        CALL PSeudoI(P%Child2, gion, NS)
        qsum    = qsum    + P%Child2%qeff
        xx(:,2) = P%Child2%xx(:)
        P%xx(:) = P%xx(:) + P%Child2%qeff*P%Child2%xx(:)
        tag(2)  = .TRUE.
        q(2)    = P%Child2%qeff
     END IF
     IF (ASSOCIATED(P%Child3)) THEN
        CALL PSeudoI(P%Child3, gion, NS)
        qsum    = qsum    + P%Child3%qeff
        xx(:,3) = P%Child3%xx(:)
        P%xx(:) = P%xx(:) + P%Child3%qeff*P%Child3%xx(:)
        tag(3)  = .TRUE.
        q(3)    = P%Child3%qeff
     END IF
     IF (ASSOCIATED(P%Child4)) THEN
        CALL PSeudoI(P%Child4, gion, NS)
        qsum    = qsum    + P%Child4%qeff
        xx(:,4) = P%Child4%xx(:)
        P%xx(:) = P%xx(:) + P%Child4%qeff*P%Child4%xx(:)
        tag(4)  = .TRUE.
        q(4)    = P%Child4%qeff
     END IF
     IF (ASSOCIATED(P%Child5)) THEN
        CALL PSeudoI(P%Child5, gion, NS)
        qsum    = qsum    + P%Child5%qeff
        xx(:,5) = P%Child5%xx(:)
        P%xx(:) = P%xx(:) + P%Child5%qeff*P%Child5%xx(:)
        tag(5)  = .TRUE.
        q(5)    = P%Child5%qeff
     END IF
     IF (ASSOCIATED(P%Child6)) THEN
        CALL PSeudoI(P%Child6, gion, NS)
        qsum    = qsum    + P%Child6%qeff
        xx(:,6) = P%Child6%xx(:)
        P%xx(:) = P%xx(:) + P%Child6%qeff*P%Child6%xx(:)
        tag(6)  = .TRUE.
        q(6)    = P%Child6%qeff
     END IF
     IF (ASSOCIATED(P%Child7)) THEN
        CALL PSeudoI(P%Child7, gion, NS)
        qsum    = qsum    + P%Child7%qeff
        xx(:,7) = P%Child7%xx(:)
        P%xx(:) = P%xx(:) + P%Child7%qeff*P%Child7%xx(:)
        tag(7)  = .TRUE.
        q(7)    = P%Child7%qeff
     END IF
     IF (ASSOCIATED(P%Child8)) THEN
        CALL PSeudoI(P%Child8, gion, NS)
        qsum    = qsum    + P%Child8%qeff
        xx(:,8) = P%Child8%xx(:)
        P%xx(:) = P%xx(:) + P%Child8%qeff*P%Child8%xx(:)
        tag(8)  = .TRUE.
        q(8)    = P%Child8%qeff
     END IF
     P%xx(:) = P%xx(:)/qsum
     P%dx = DSQRT((P%xx(1)-ox(1))**2 + (P%xx(2)-ox(2))**2 + (P%xx(3)-ox(3))**2)
     DO i=1,8
        IF (tag(i)) THEN
           P%D(1)  = P%D(1)  + q(i)*(xx(1,i) - P%xx(1))
           P%D(2)  = P%D(2)  + q(i)*(xx(2,i) - P%xx(2))
           P%D(3)  = P%D(3)  + q(i)*(xx(3,i) - P%xx(3))
           P%Q(1)  = P%Q(1)  + q(i)*(xx(1,i) - P%xx(1))**2/2.
           P%Q(2)  = P%Q(2)  + q(i)*(xx(2,i) - P%xx(2))**2/2.
           P%Q(3)  = P%Q(3)  + q(i)*(xx(3,i) - P%xx(3))**2/2.
           P%QQ(1) = P%QQ(1) + q(i)*(xx(1,i) - P%xx(1))*(xx(2,i) - P%xx(2))
           P%QQ(2) = P%QQ(2) + q(i)*(xx(2,i) - P%xx(2))*(xx(3,i) - P%xx(3))
           P%QQ(3) = P%QQ(3) + q(i)*(xx(3,i) - P%xx(3))*(xx(1,i) - P%xx(1))
        END IF
     END DO
  END IF
END SUBROUTINE PSeudoI
!
! #########################
SUBROUTINE LocalCoulombSum(NS, Roote, Rooti, qel, qion, sys, param, ttime)
  USE DATASTR
  IMPLICIT NONE
  !
  TYPE(NM), INTENT(INOUT):: NS
  TYPE(PT), INTENT(INOUT):: qel(NS%Nel), qion(NS%Nion)
  TYPE(ST), INTENT(INOUT):: sys
  TYPE(PM), INTENT(IN):: param
  TYPE(TM), INTENT(IN):: ttime
  TYPE(BT), Pointer:: Roote, Rooti
  INTEGER :: i, Rei(Ndist+1), Rii(Ndist+1)
  REAL*8  :: Epot
  sys%Epot_e = 0.0D0
  qel(:)%ff(1) = 0.0D0
  qel(:)%ff(2) = 0.0D0
  qel(:)%ff(3) = 0.0D0
  qion(:)%ff(1) = 0.0D0
  qion(:)%ff(2) = 0.0D0
  qion(:)%ff(3) = 0.0D0
  Epot = 0.0D0
  !$OMP PARALLEL DO SHARED(NS, Roote, qel, param) PRIVATE(i) REDUCTION(+:Epot)
  DO i=1, NS%Nel
     CALL Search_Tree_ee(NS, Roote, qel, i, param, Epot)
  END DO
  !$OMP END PARALLEL DO
  sys%Epot_e = sys%Epot_e + Epot
  Epot = 0.0D0
  Rei = 0
  !$OMP PARALLEL DO SHARED(NS, Roote, qion, param) PRIVATE(i) REDUCTION(+:Epot,Rei)
  DO i=1, NS%Nion
     CALL Search_Tree_ie(NS, Roote, qion, i, param, Epot, Rei)
  END DO
  !$OMP END PARALLEL DO
  sys%Epot_e = sys%Epot_e + Epot
  Epot = 0.0D0  
  !$OMP PARALLEL DO SHARED(NS, Rooti, qel, param) PRIVATE(i) REDUCTION(+:Epot,Rei)
  DO i=1, NS%Nel
     CALL Search_Tree_ei(NS, Rooti, qel, i, param, Epot, Rei)
  END DO
  !$OMP END PARALLEL DO
  sys%Epot_e = sys%Epot_e + Epot
  NS%RDF_ei(:) = NS%RDF_ei(:) + Rei(:)
  !
  sys%Epot_i = 0.0D0
  Epot = 0.0D0
  Rii = 0.0D0
  !$OMP PARALLEL DO SHARED(NS, Rooti, qion, param) PRIVATE(i) REDUCTION(+:Epot,Rii)
  DO i=1, NS%Nion
     CALL Search_Tree_ii(NS, Rooti, qion, i, param, Epot, Rii)
  END DO
  !$OMP END PARALLEL DO
  sys%Epot_i = sys%Epot_i + Epot     
  NS%RDF_ii(:) = NS%RDF_ii(:) + Rii(:)
  !
  RETURN
END SUBROUTINE LocalCoulombSum
!
!
! #########################
RECURSIVE SUBROUTINE Search_Tree_ee(NS, P, qel, i, param, Epot)
  USE DATASTR
  IMPLICIT NONE
  !
  TYPE(NM), INTENT(INOUT):: NS
  TYPE(PT), INTENT(INOUT):: qel(NS%Nel)
  TYPE(PM), INTENT(IN):: param
  INTEGER,  INTENT(IN):: i
  TYPE(BT), Pointer:: P
  INTEGER :: k
  REAL*8  :: s, r, ff, xr(3), r2, R3, R5, R7, D(6), Q(10), Epot
  !
  IF (P%ileaf == 0) THEN
     !
     ! PP - twig
     s = MAXVAL(P%lx(:))
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qel(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     r = DSQRT(r2)
     IF ( (s/param%sd_ratio + P%dx) < r) THEN
        ff = -1.0D0*eps
        R3 = r2*r
        R5 = R3*r2
        R7 = R5*r2
        D(1) = -1.D0/R3 + 3.D0*xr(1)**2/R5
        D(2) = -1.D0/R3 + 3.D0*xr(2)**2/R5
        D(3) = -1.D0/R3 + 3.D0*xr(3)**2/R5
        D(4) = 3.D0*xr(1)*xr(2)/R5
        D(5) = 3.D0*xr(2)*xr(3)/R5
        D(6) = 3.D0*xr(3)*xr(1)/R5
        Q(1) = 15.D0*xr(1)**3/R7 - 9.D0*xr(1)/R5
        Q(2) = 15.D0*xr(2)**3/R7 - 9.D0*xr(2)/R5
        Q(3) = 15.D0*xr(3)**3/R7 - 9.D0*xr(3)/R5
        Q(4) = 15.D0*xr(1)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(5) = 15.D0*xr(1)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(6) = 15.D0*xr(2)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(7) = 15.D0*xr(2)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(8) = 15.D0*xr(3)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(9) = 15.D0*xr(3)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(10)= 15.D0*xr(1)*xr(2)*xr(3)/R7
        Epot = Epot + 0.5D0*ff*(P%qeff/r + &
             (P%D(1)*xr(1)+P%D(2)*xr(2)+P%D(3)*xr(3))/R3 + &
             P%Q(1)*D(1) + P%Q(2)*D(2)+P%Q(3)*D(3) + &
             P%QQ(1)*D(4) + P%QQ(2)*D(5)+P%QQ(3)*D(6))
        qel(i)%ff(1) = qel(i)%ff(1) + ff*(P%qeff*xr(1)/R3 + &
             D(1)*P%D(1) + D(4)*P%D(2) + D(6)*P%D(3) + &
             Q(1)*P%Q(1) + Q(6)*P%Q(2) + Q(8)*P%Q(3) + &
             Q(4)*P%QQ(1) + Q(5)*P%QQ(3) + Q(10)*P%QQ(2))
        qel(i)%ff(2) = qel(i)%ff(2) + ff*(P%qeff*xr(2)/R3 + &
             D(4)*P%D(1) + D(2)*P%D(2) + D(5)*P%D(3) + &
             Q(4)*P%Q(1) + Q(2)*P%Q(2) + Q(9)*P%Q(3) + &
             Q(6)*P%QQ(1) + Q(10)*P%QQ(3) + Q(7)*P%QQ(2))
        qel(i)%ff(3) = qel(i)%ff(3) + ff*(P%qeff*xr(3)/R3 + &
             D(6)*P%D(1) + D(5)*P%D(2) + D(3)*P%D(3) + &
             Q(5)*P%Q(1) + Q(7)*P%Q(2) + Q(3)*P%Q(3) + &
             Q(10)*P%QQ(1) + Q(8)*P%QQ(3) + Q(9)*P%QQ(2))
     ELSE
        IF (ASSOCIATED(P%Child1)) &
             CALL Search_Tree_ee(NS, P%Child1, qel, i, param, Epot)
        IF (ASSOCIATED(P%Child2)) &
             CALL Search_Tree_ee(NS, P%Child2, qel, i, param, Epot)
        IF (ASSOCIATED(P%Child3)) &
             CALL Search_Tree_ee(NS, P%Child3, qel, i, param, Epot)
        IF (ASSOCIATED(P%Child4)) &
             CALL Search_Tree_ee(NS, P%Child4, qel, i, param, Epot)
        IF (ASSOCIATED(P%Child5)) &
             CALL Search_Tree_ee(NS, P%Child5, qel, i, param, Epot)
        IF (ASSOCIATED(P%Child6)) &
             CALL Search_Tree_ee(NS, P%Child6, qel, i, param, Epot)
        IF (ASSOCIATED(P%Child7)) &
             CALL Search_Tree_ee(NS, P%Child7, qel, i, param, Epot)
        IF (ASSOCIATED(P%Child8)) &
             CALL Search_Tree_ee(NS, P%Child8, qel, i, param, Epot)
     END IF
  ELSE
     !
     ! leaf of electron
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qel(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     IF ( r2 > 0.01D0 ) THEN
        r = DSQRT(r2)
        ff = eps
        !
        ! Generic Coulomb
        Epot = Epot + 0.5D0*ff/r
        DO k=1,3
           qel(i)%ff(k) = qel(i)%ff(k) + ff*xr(k)/r/r2
        END DO
     END IF
  END IF
END SUBROUTINE Search_Tree_ee
!
! #########################
RECURSIVE SUBROUTINE Search_Tree_ie(NS, P, qion, i, param, Epot, Rei)
  USE DATASTR
  IMPLICIT NONE
  !
  TYPE(NM), INTENT(INOUT):: NS
  TYPE(PT), INTENT(INOUT):: qion(NS%Nion)
  TYPE(PM), INTENT(IN):: param
  INTEGER,  INTENT(IN):: i
  TYPE(BT), Pointer:: P
  INTEGER :: k, id, Rei(Ndist+1)
  REAL*8  :: s, r, ff, xr(3), r2, R3, R5, R7, D(6), Q(10), Epot
  !
  IF (P%ileaf == 0) THEN
     !
     ! PP - twig
     s = MAXVAL(P%lx(:))
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qion(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     r = DSQRT(r2)
     IF ( (s/param%sd_ratio + P%dx) < r) THEN
        ff = eps
        R3 = r2*r
        R5 = R3*r2
        R7 = R5*r2
        D(1) = -1.D0/R3 + 3.D0*xr(1)**2/R5
        D(2) = -1.D0/R3 + 3.D0*xr(2)**2/R5
        D(3) = -1.D0/R3 + 3.D0*xr(3)**2/R5
        D(4) = 3.D0*xr(1)*xr(2)/R5
        D(5) = 3.D0*xr(2)*xr(3)/R5
        D(6) = 3.D0*xr(3)*xr(1)/R5
        Q(1) = 15.D0*xr(1)**3/R7 - 9.D0*xr(1)/R5
        Q(2) = 15.D0*xr(2)**3/R7 - 9.D0*xr(2)/R5
        Q(3) = 15.D0*xr(3)**3/R7 - 9.D0*xr(3)/R5
        Q(4) = 15.D0*xr(1)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(5) = 15.D0*xr(1)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(6) = 15.D0*xr(2)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(7) = 15.D0*xr(2)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(8) = 15.D0*xr(3)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(9) = 15.D0*xr(3)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(10)= 15.D0*xr(1)*xr(2)*xr(3)/R7
        Epot = Epot + 0.5D0*ff*(P%qeff/r + &
             (P%D(1)*xr(1)+P%D(2)*xr(2)+P%D(3)*xr(3))/R3 + &
             P%Q(1)*D(1) + P%Q(2)*D(2)+P%Q(3)*D(3) + &
             P%QQ(1)*D(4) + P%QQ(2)*D(5)+P%QQ(3)*D(6))
        qion(i)%ff(1) = qion(i)%ff(1) + ff*(P%qeff*xr(1)/R3 + &
             D(1)*P%D(1) + D(4)*P%D(2) + D(6)*P%D(3) + &
             Q(1)*P%Q(1) + Q(6)*P%Q(2) + Q(8)*P%Q(3) + &
             Q(4)*P%QQ(1) + Q(5)*P%QQ(3) + Q(10)*P%QQ(2))
        qion(i)%ff(2) = qion(i)%ff(2) + ff*(P%qeff*xr(2)/R3 + &
             D(4)*P%D(1) + D(2)*P%D(2) + D(5)*P%D(3) + &
             Q(4)*P%Q(1) + Q(2)*P%Q(2) + Q(9)*P%Q(3) + &
             Q(6)*P%QQ(1) + Q(10)*P%QQ(3) + Q(7)*P%QQ(2))
        qion(i)%ff(3) = qion(i)%ff(3) + ff*(P%qeff*xr(3)/R3 + &
             D(6)*P%D(1) + D(5)*P%D(2) + D(3)*P%D(3) + &
             Q(5)*P%Q(1) + Q(7)*P%Q(2) + Q(3)*P%Q(3) + &
             Q(10)*P%QQ(1) + Q(8)*P%QQ(3) + Q(9)*P%QQ(2))
     ELSE
        IF (ASSOCIATED(P%Child1)) &
             CALL Search_Tree_ie(NS, P%Child1, qion, i, param, Epot, Rei)
        IF (ASSOCIATED(P%Child2)) &
             CALL Search_Tree_ie(NS, P%Child2, qion, i, param, Epot, Rei)
        IF (ASSOCIATED(P%Child3)) &
             CALL Search_Tree_ie(NS, P%Child3, qion, i, param, Epot, Rei)
        IF (ASSOCIATED(P%Child4)) &
             CALL Search_Tree_ie(NS, P%Child4, qion, i, param, Epot, Rei)
        IF (ASSOCIATED(P%Child5)) &
             CALL Search_Tree_ie(NS, P%Child5, qion, i, param, Epot, Rei)
        IF (ASSOCIATED(P%Child6)) &
             CALL Search_Tree_ie(NS, P%Child6, qion, i, param, Epot, Rei)
        IF (ASSOCIATED(P%Child7)) &
             CALL Search_Tree_ie(NS, P%Child7, qion, i, param, Epot, Rei)
        IF (ASSOCIATED(P%Child8)) &
             CALL Search_Tree_ie(NS, P%Child8, qion, i, param, Epot, Rei)
     END IF
  ELSE
     !
     ! leaf of electron
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qion(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     r = DSQRT(r2)
     ff = -1.D0*eps
     !
     ! Kelbg potential
     Epot = Epot + 0.5D0*ff*(1.-EXP(-r/param%rs))/r
     DO k=1,3
        qion(i)%ff(k) = qion(i)%ff(k) + &
             & ff*xr(k)*(1.D0/r - EXP(-r/param%rs)*(1.D0/r+1.D0/param%rs))/r2
     END DO
     !
     ! RDF
     id = INT(r/NS%dx_ei)+1
     IF (id <= Ndist + 1 .AND. id > 0 ) THEN
        NS%RDF_ei(id) = NS%RDF_ei(id) + 1
     END IF
  END IF
END SUBROUTINE Search_Tree_ie
!
!
! #########################
RECURSIVE SUBROUTINE Search_Tree_ei(NS, P, qel, i, param, Epot, Rei)
  USE DATASTR
  IMPLICIT NONE
  !
  TYPE(NM), INTENT(INOUT):: NS
  TYPE(PT), INTENT(INOUT):: qel(NS%Nel)
  TYPE(PM), INTENT(IN):: param
  INTEGER,  INTENT(IN):: i
  TYPE(BT), Pointer:: P
  INTEGER :: k, id, Rei(Ndist+1)
  REAL*8  :: s, r, ff, xr(3), r2, R3, R5, R7, D(6), Q(10), Epot
  !
  IF (P%ileaf == 0) THEN
     !
     ! PP - twig
     s = MAXVAL(P%lx(:))
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qel(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     r = DSQRT(r2)
     IF ( (s/param%sd_ratio + P%dx) < r) THEN
        ff = -1.0D0*eps
        R3 = r2*r
        R5 = R3*r2
        R7 = R5*r2
        D(1) = -1.D0/R3 + 3.D0*xr(1)**2/R5
        D(2) = -1.D0/R3 + 3.D0*xr(2)**2/R5
        D(3) = -1.D0/R3 + 3.D0*xr(3)**2/R5
        D(4) = 3.D0*xr(1)*xr(2)/R5
        D(5) = 3.D0*xr(2)*xr(3)/R5
        D(6) = 3.D0*xr(3)*xr(1)/R5
        Q(1) = 15.D0*xr(1)**3/R7 - 9.D0*xr(1)/R5
        Q(2) = 15.D0*xr(2)**3/R7 - 9.D0*xr(2)/R5
        Q(3) = 15.D0*xr(3)**3/R7 - 9.D0*xr(3)/R5
        Q(4) = 15.D0*xr(1)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(5) = 15.D0*xr(1)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(6) = 15.D0*xr(2)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(7) = 15.D0*xr(2)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(8) = 15.D0*xr(3)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(9) = 15.D0*xr(3)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(10)= 15.D0*xr(1)*xr(2)*xr(3)/R7
        Epot = Epot + 0.5D0*ff*(P%qeff/r + &
             (P%D(1)*xr(1)+P%D(2)*xr(2)+P%D(3)*xr(3))/R3 + &
             P%Q(1)*D(1) + P%Q(2)*D(2)+P%Q(3)*D(3) + &
             P%QQ(1)*D(4) + P%QQ(2)*D(5)+P%QQ(3)*D(6))
        qel(i)%ff(1) = qel(i)%ff(1) + ff*(P%qeff*xr(1)/R3 + &
             D(1)*P%D(1) + D(4)*P%D(2) + D(6)*P%D(3) + &
             Q(1)*P%Q(1) + Q(6)*P%Q(2) + Q(8)*P%Q(3) + &
             Q(4)*P%QQ(1) + Q(5)*P%QQ(3) + Q(10)*P%QQ(2))
        qel(i)%ff(2) = qel(i)%ff(2) + ff*(P%qeff*xr(2)/R3 + &
             D(4)*P%D(1) + D(2)*P%D(2) + D(5)*P%D(3) + &
             Q(4)*P%Q(1) + Q(2)*P%Q(2) + Q(9)*P%Q(3) + &
             Q(6)*P%QQ(1) + Q(10)*P%QQ(3) + Q(7)*P%QQ(2))
        qel(i)%ff(3) = qel(i)%ff(3) + ff*(P%qeff*xr(3)/R3 + &
             D(6)*P%D(1) + D(5)*P%D(2) + D(3)*P%D(3) + &
             Q(5)*P%Q(1) + Q(7)*P%Q(2) + Q(3)*P%Q(3) + &
             Q(10)*P%QQ(1) + Q(8)*P%QQ(3) + Q(9)*P%QQ(2))
     ELSE
        IF (ASSOCIATED(P%Child1)) &
             CALL Search_Tree_ei(NS, P%Child1, qel, i, param, Epot, Rei)
        IF (ASSOCIATED(P%Child2)) &
             CALL Search_Tree_ei(NS, P%Child2, qel, i, param, Epot, Rei)
        IF (ASSOCIATED(P%Child3)) &
             CALL Search_Tree_ei(NS, P%Child3, qel, i, param, Epot, Rei)
        IF (ASSOCIATED(P%Child4)) &
             CALL Search_Tree_ei(NS, P%Child4, qel, i, param, Epot, Rei)
        IF (ASSOCIATED(P%Child5)) &
             CALL Search_Tree_ei(NS, P%Child5, qel, i, param, Epot, Rei)
        IF (ASSOCIATED(P%Child6)) &
             CALL Search_Tree_ei(NS, P%Child6, qel, i, param, Epot, Rei)
        IF (ASSOCIATED(P%Child7)) &
             CALL Search_Tree_ei(NS, P%Child7, qel, i, param, Epot, Rei)
        IF (ASSOCIATED(P%Child8)) &
             CALL Search_Tree_ei(NS, P%Child8, qel, i, param, Epot, Rei)
     END IF
  ELSE 
     !
     ! leaf of ion
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qel(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     r = DSQRT(r2)
     ff = -1.D0*eps
     !
     ! Kelbg potential
     Epot = Epot + 0.5D0*ff*(1.-EXP(-r/param%rs))/r
     DO k=1,3
        qel(i)%ff(k) = qel(i)%ff(k) + &
             ff*xr(k)*(1.D0/r - EXP(-r/param%rs)*(1.D0/r+1.D0/param%rs))/r2
     END DO
     !
     ! RDF
     id = INT(r/NS%dx_ei)+1
     IF (id <= Ndist + 1 .AND. id > 0 ) THEN
        Rei(id) = Rei(id) + 1
     END IF
     !
  END IF
END SUBROUTINE Search_Tree_ei
!
! #########################
RECURSIVE SUBROUTINE Search_Tree_ii(NS, P, qion, i, param, Epot, Rii)
  USE DATASTR
  IMPLICIT NONE
  !
  TYPE(NM), INTENT(INOUT):: NS
  TYPE(PT), INTENT(INOUT):: qion(NS%Nion)
  TYPE(PM), INTENT(IN):: param
  INTEGER,  INTENT(IN):: i
  TYPE(BT), Pointer:: P
  INTEGER :: k, id, Rii(Ndist+1)
  REAL*8  :: s, r, ff, xr(3), r2, R3, R5, R7, D(6), Q(10), Epot
  !
  IF (P%ileaf == 0) THEN
     !
     ! PP - twig
     s = MAXVAL(P%lx(:))
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qion(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     r = DSQRT(r2)
     IF ( (s/param%sd_ratio + P%dx) < r) THEN
        ff = eps
        R3 = r2*r
        R5 = R3*r2
        R7 = R5*r2
        D(1) = -1.D0/R3 + 3.D0*xr(1)**2/R5
        D(2) = -1.D0/R3 + 3.D0*xr(2)**2/R5
        D(3) = -1.D0/R3 + 3.D0*xr(3)**2/R5
        D(4) = 3.D0*xr(1)*xr(2)/R5
        D(5) = 3.D0*xr(2)*xr(3)/R5
        D(6) = 3.D0*xr(3)*xr(1)/R5
        Q(1) = 15.D0*xr(1)**3/R7 - 9.D0*xr(1)/R5
        Q(2) = 15.D0*xr(2)**3/R7 - 9.D0*xr(2)/R5
        Q(3) = 15.D0*xr(3)**3/R7 - 9.D0*xr(3)/R5
        Q(4) = 15.D0*xr(1)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(5) = 15.D0*xr(1)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(6) = 15.D0*xr(2)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(7) = 15.D0*xr(2)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(8) = 15.D0*xr(3)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(9) = 15.D0*xr(3)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(10)= 15.D0*xr(1)*xr(2)*xr(3)/R7
        Epot = Epot + 0.5D0*ff*(P%qeff/r + &
             (P%D(1)*xr(1)+P%D(2)*xr(2)+P%D(3)*xr(3))/R3 + &
             P%Q(1)*D(1) + P%Q(2)*D(2)+P%Q(3)*D(3) + &
             P%QQ(1)*D(4) + P%QQ(2)*D(5)+P%QQ(3)*D(6))
        qion(i)%ff(1) = qion(i)%ff(1) + ff*(P%qeff*xr(1)/R3 + &
             D(1)*P%D(1) + D(4)*P%D(2) + D(6)*P%D(3) + &
             Q(1)*P%Q(1) + Q(6)*P%Q(2) + Q(8)*P%Q(3) + &
             Q(4)*P%QQ(1) + Q(5)*P%QQ(3) + Q(10)*P%QQ(2))
        qion(i)%ff(2) = qion(i)%ff(2) + ff*(P%qeff*xr(2)/R3 + &
             D(4)*P%D(1) + D(2)*P%D(2) + D(5)*P%D(3) + &
             Q(4)*P%Q(1) + Q(2)*P%Q(2) + Q(9)*P%Q(3) + &
             Q(6)*P%QQ(1) + Q(10)*P%QQ(3) + Q(7)*P%QQ(2))
        qion(i)%ff(3) = qion(i)%ff(3) + ff*(P%qeff*xr(3)/R3 + &
             D(6)*P%D(1) + D(5)*P%D(2) + D(3)*P%D(3) + &
             Q(5)*P%Q(1) + Q(7)*P%Q(2) + Q(3)*P%Q(3) + &
             Q(10)*P%QQ(1) + Q(8)*P%QQ(3) + Q(9)*P%QQ(2))
     ELSE
        IF (ASSOCIATED(P%Child1)) &
             CALL Search_Tree_ii(NS, P%Child1, qion, i, param, Epot, Rii)
        IF (ASSOCIATED(P%Child2)) &
             CALL Search_Tree_ii(NS, P%Child2, qion, i, param, Epot, Rii)
        IF (ASSOCIATED(P%Child3)) &
             CALL Search_Tree_ii(NS, P%Child3, qion, i, param, Epot, Rii)
        IF (ASSOCIATED(P%Child4)) &
             CALL Search_Tree_ii(NS, P%Child4, qion, i, param, Epot, Rii)
        IF (ASSOCIATED(P%Child5)) &
             CALL Search_Tree_ii(NS, P%Child5, qion, i, param, Epot, Rii)
        IF (ASSOCIATED(P%Child6)) &
             CALL Search_Tree_ii(NS, P%Child6, qion, i, param, Epot, Rii)
        IF (ASSOCIATED(P%Child7)) &
             CALL Search_Tree_ii(NS, P%Child7, qion, i, param, Epot, Rii)
        IF (ASSOCIATED(P%Child8)) &
             CALL Search_Tree_ii(NS, P%Child8, qion, i, param, Epot, Rii)
     END IF
  ELSE
     !
     ! leaf of ion
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qion(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     IF ( r2 > 0.01D0 ) THEN
        r = DSQRT(r2)
        ff = eps
        !
        ! Generic Coulomb
        Epot = Epot + 0.5D0*ff/r
        DO k=1,3
           qion(i)%ff(k) = qion(i)%ff(k) + ff*xr(k)/r/r2
        END DO
        !
        ! RDF
        id = INT(r/NS%dx_ii)+1
        IF (id <= Ndist + 1 .AND. id > 0 ) THEN
           Rii(id) = Rii(id) + 1
        END IF
        !
     END IF
  END IF
END SUBROUTINE Search_Tree_ii
!
! ############## DEALLOCATE all linked pointer #####################
RECURSIVE SUBROUTINE DEALLOC_TREE(P)
  IMPLICIT NONE
  TYPE(BT), POINTER:: P
  IF (ASSOCIATED(P%Child1)) CALL DEALLOC_TREE(P%Child1)
  IF (ASSOCIATED(P%Child2)) CALL DEALLOC_TREE(P%Child2)
  IF (ASSOCIATED(P%Child3)) CALL DEALLOC_TREE(P%Child3)
  IF (ASSOCIATED(P%Child4)) CALL DEALLOC_TREE(P%Child4)
  IF (ASSOCIATED(P%Child5)) CALL DEALLOC_TREE(P%Child5)
  IF (ASSOCIATED(P%Child6)) CALL DEALLOC_TREE(P%Child6)
  IF (ASSOCIATED(P%Child7)) CALL DEALLOC_TREE(P%Child7)
  IF (ASSOCIATED(P%Child8)) CALL DEALLOC_TREE(P%Child8)
  DEALLOCATE(P)
  !
END SUBROUTINE DEALLOC_TREE
RECURSIVE SUBROUTINE DEALLOC_gTREE(P)
  IMPLICIT NONE
  TYPE(GT), POINTER:: P
  IF (ASSOCIATED(P%Child1)) CALL DEALLOC_gTREE(P%Child1)
  IF (ASSOCIATED(P%Child2)) CALL DEALLOC_gTREE(P%Child2)
  IF (ASSOCIATED(P%Child3)) CALL DEALLOC_gTREE(P%Child3)
  IF (ASSOCIATED(P%Child4)) CALL DEALLOC_gTREE(P%Child4)
  IF (ASSOCIATED(P%Child5)) CALL DEALLOC_gTREE(P%Child5)
  IF (ASSOCIATED(P%Child6)) CALL DEALLOC_gTREE(P%Child6)
  IF (ASSOCIATED(P%Child7)) CALL DEALLOC_gTREE(P%Child7)
  IF (ASSOCIATED(P%Child8)) CALL DEALLOC_gTREE(P%Child8)
  DEALLOCATE(P)
  !
END SUBROUTINE DEALLOC_gTREE
!
! ############### Copy local tree into array form #######################
! #######################################################################
SUBROUTINE CopyTree(P, arch, node, Nnode, Narch, NS, param, N)
  USE DATASTR
  IMPLICIT NONE
  !  
  TYPE(NM), INTENT(INOUT):: NS
  TYPE(GH), INTENT(INOUT):: node(NS%Nnode)
  TYPE(PM), INTENT(IN):: param
  INTEGER,  INTENT(INOUT):: arch(NS%Nnode), Narch, Nnode
  TYPE(BT), Pointer:: P
  INTEGER:: N
  REAL*8  :: x, y, z
  x = 0.0D0
  y = 0.0D0
  z = 0.0D0
  SELECT CASE(N)
  CASE (1) ! FACE contact along y-z plane
     CALL CopyTree1(P, arch, node, Nnode, Narch, NS, param, x)
  CASE (2) ! FACE contact along x-z plane
     CALL CopyTree2(P, arch, node, Nnode, Narch, NS, param, y)
  CASE (3) ! FACE contact along x-y plane
     CALL CopyTree3(P, arch, node, Nnode, Narch, NS, param, z)
  CASE (4) ! LINE contact along x axis
     CALL CopyTree4(P, arch, node, Nnode, Narch, NS, param, y, z)     
  CASE (5) ! LINE contact along y axis
     CALL CopyTree5(P, arch, node, Nnode, Narch, NS, param, x, z)
  CASE (6) ! LINE contact along z axis
     CALL CopyTree6(P, arch, node, Nnode, Narch, NS, param, x, y)
  CASE (7) ! POINT contact along origin
     CALL CopyTree7(P, arch, node, Nnode, Narch, NS, param, x, y, z)
  CASE DEFAULT
     STOP " CRITERION FOR COPYING TREE ERROR"
  END SELECT
  Narch = Narch + 1
  arch(Narch) = -1
  RETURN
END SUBROUTINE CopyTree
!
RECURSIVE SUBROUTINE CopyTree1(P, arch, node, Nnode, Narch, NS, param, x)
  USE DATASTR
  IMPLICIT NONE
  !  
  TYPE(NM), INTENT(IN):: NS
  TYPE(GH), INTENT(INOUT):: node(NS%Nnode)
  TYPE(PM), INTENT(IN):: param
  INTEGER, INTENT(INOUT):: arch(NS%Nnode), Narch, Nnode
  REAL*8,  INTENT(IN):: x
  TYPE(BT), Pointer   :: P
  Nnode = Nnode + 1 ! position tag
  IF (P%ileaf == 0) THEN
     node(Nnode)%qeff  = P%qeff
     node(Nnode)%xx(:) = P%xx(:)
     node(Nnode)%lx    = MAXVAL(P%lx(:))
     node(Nnode)%dx    = P%dx
     node(Nnode)%D(:)  = P%D(:)
     node(Nnode)%Q(:)  = P%Q(:)
     node(Nnode)%QQ(:) = P%QQ(:)
     Narch = Narch + 1
     arch(Narch) = Nnode
     !
     IF ( (node(Nnode)%lx/param%sd_ratio + node(Nnode)%dx) > &
          ABS(node(Nnode)%xx(1) - x) ) THEN
        !
        ! PP - twig
        IF (ASSOCIATED(P%Child1)) THEN
           CALL CopyTree1(P%Child1, arch, node, Nnode, Narch, NS, param, x)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child2)) THEN
           CALL CopyTree1(P%Child2, arch, node, Nnode, Narch, NS, param, x)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child3)) THEN
           CALL CopyTree1(P%Child3, arch, node, Nnode, Narch, NS, param, x)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child4)) THEN
           CALL CopyTree1(P%Child4, arch, node, Nnode, Narch, NS, param, x)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child5)) THEN
           CALL CopyTree1(P%Child5, arch, node, Nnode, Narch, NS, param, x)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child6)) THEN
           CALL CopyTree1(P%Child6, arch, node, Nnode, Narch, NS, param, x)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child7)) THEN
           CALL CopyTree1(P%Child7, arch, node, Nnode, Narch, NS, param, x)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child8)) THEN
           CALL CopyTree1(P%Child8, arch, node, Nnode, Narch, NS, param, x)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
     END IF
  ELSE
     !
     ! leaf
     node(Nnode)%qeff  = P%qeff
     node(Nnode)%xx(:) = P%xx(:)
     node(Nnode)%lx    = MAXVAL(P%lx(:))
     node(Nnode)%dx    = 0.0D0
     node(Nnode)%D(:)  = 0.D0
     node(Nnode)%Q(:)  = 0.D0
     node(Nnode)%QQ(:) = 0.D0
     Narch = Narch + 1
     arch(Narch) = Nnode
  END IF
END SUBROUTINE CopyTree1
!
RECURSIVE SUBROUTINE CopyTree2(P, arch, node, Nnode, Narch, NS, param, y)
  USE DATASTR
  IMPLICIT NONE
  !  
  TYPE(NM), INTENT(IN):: NS
  TYPE(GH), INTENT(INOUT):: node(NS%Nnode)
  TYPE(PM), INTENT(IN):: param
  INTEGER, INTENT(INOUT):: arch(NS%Nnode), Narch, Nnode
  REAL*8,  INTENT(IN):: y
  TYPE(BT), Pointer   :: P
  Nnode = Nnode + 1 ! position tag
  IF (P%ileaf == 0) THEN
     node(Nnode)%qeff  = P%qeff
     node(Nnode)%xx(:) = P%xx(:)
     node(Nnode)%lx    = MAXVAL(P%lx(:))
     node(Nnode)%dx    = P%dx
     node(Nnode)%D(:)  = P%D(:)
     node(Nnode)%Q(:)  = P%Q(:)
     node(Nnode)%QQ(:) = P%QQ(:)
     Narch = Narch + 1
     arch(Narch) = Nnode
     !
     IF ( (node(Nnode)%lx/param%sd_ratio + node(Nnode)%dx) > &
          ABS(node(Nnode)%xx(2) - y) ) THEN
        !
        ! PP - twig
        IF (ASSOCIATED(P%Child1)) THEN
           CALL CopyTree2(P%Child1, arch, node, Nnode, Narch, NS, param, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child2)) THEN
           CALL CopyTree2(P%Child2, arch, node, Nnode, Narch, NS, param, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child3)) THEN
           CALL CopyTree2(P%Child3, arch, node, Nnode, Narch, NS, param, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child4)) THEN
           CALL CopyTree2(P%Child4, arch, node, Nnode, Narch, NS, param, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child5)) THEN
           CALL CopyTree2(P%Child5, arch, node, Nnode, Narch, NS, param, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child6)) THEN
           CALL CopyTree2(P%Child6, arch, node, Nnode, Narch, NS, param, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child7)) THEN
           CALL CopyTree2(P%Child7, arch, node, Nnode, Narch, NS, param, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child8)) THEN
           CALL CopyTree2(P%Child8, arch, node, Nnode, Narch, NS, param, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
     END IF
  ELSE
     !
     ! leaf
     node(Nnode)%qeff  = P%qeff
     node(Nnode)%xx(:) = P%xx(:)
     node(Nnode)%lx    = MAXVAL(P%lx(:))
     node(Nnode)%dx    = 0.0D0
     node(Nnode)%D(:)  = 0.D0
     node(Nnode)%Q(:)  = 0.D0
     node(Nnode)%QQ(:) = 0.D0
     Narch = Narch + 1
     arch(Narch) = Nnode
  END IF
END SUBROUTINE CopyTree2
!
RECURSIVE SUBROUTINE CopyTree3(P, arch, node, Nnode, Narch, NS, param, z)
  USE DATASTR
  IMPLICIT NONE
  !  
  TYPE(NM), INTENT(IN):: NS
  TYPE(GH), INTENT(INOUT):: node(NS%Nnode)
  TYPE(PM), INTENT(IN):: param
  INTEGER, INTENT(INOUT):: arch(NS%Nnode), Narch, Nnode
  REAL*8,  INTENT(IN):: z
  TYPE(BT), Pointer   :: P
  Nnode = Nnode + 1 ! position tag
  IF (P%ileaf == 0) THEN
     node(Nnode)%qeff  = P%qeff
     node(Nnode)%xx(:) = P%xx(:)
     node(Nnode)%lx    = MAXVAL(P%lx(:))
     node(Nnode)%dx    = P%dx
     node(Nnode)%D(:)  = P%D(:)
     node(Nnode)%Q(:)  = P%Q(:)
     node(Nnode)%QQ(:) = P%QQ(:)
     Narch = Narch + 1
     arch(Narch) = Nnode
     !
     IF ( (node(Nnode)%lx/param%sd_ratio + node(Nnode)%dx) > &
          ABS(node(Nnode)%xx(3) - z) ) THEN
        !
        ! PP - twig
        IF (ASSOCIATED(P%Child1)) THEN
           CALL CopyTree3(P%Child1, arch, node, Nnode, Narch, NS, param, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child2)) THEN
           CALL CopyTree3(P%Child2, arch, node, Nnode, Narch, NS, param, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child3)) THEN
           CALL CopyTree3(P%Child3, arch, node, Nnode, Narch, NS, param, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child4)) THEN
           CALL CopyTree3(P%Child4, arch, node, Nnode, Narch, NS, param, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child5)) THEN
           CALL CopyTree3(P%Child5, arch, node, Nnode, Narch, NS, param, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child6)) THEN
           CALL CopyTree3(P%Child6, arch, node, Nnode, Narch, NS, param, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child7)) THEN
           CALL CopyTree3(P%Child7, arch, node, Nnode, Narch, NS, param, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child8)) THEN
           CALL CopyTree3(P%Child8, arch, node, Nnode, Narch, NS, param, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
     END IF
  ELSE
     !
     ! leaf
     node(Nnode)%qeff  = P%qeff
     node(Nnode)%xx(:) = P%xx(:)
     node(Nnode)%lx    = MAXVAL(P%lx(:))
     node(Nnode)%dx    = 0.0D0
     node(Nnode)%D(:)  = 0.D0
     node(Nnode)%Q(:)  = 0.D0
     node(Nnode)%QQ(:) = 0.D0
     Narch = Narch + 1
     arch(Narch) = Nnode
  END IF
END SUBROUTINE CopyTree3
!
RECURSIVE SUBROUTINE CopyTree4(P, arch, node, Nnode, Narch, NS, param, y, z)
  USE DATASTR
  IMPLICIT NONE
  !  
  TYPE(NM), INTENT(IN):: NS
  TYPE(GH), INTENT(INOUT):: node(NS%Nnode)
  TYPE(PM), INTENT(IN):: param
  INTEGER, INTENT(INOUT):: arch(NS%Nnode), Narch, Nnode
  REAL*8,  INTENT(IN):: y, z
  TYPE(BT), Pointer   :: P
  Nnode = Nnode + 1 ! position tag
  IF (P%ileaf == 0) THEN
     node(Nnode)%qeff  = P%qeff
     node(Nnode)%xx(:) = P%xx(:)
     node(Nnode)%lx    = MAXVAL(P%lx(:))
     node(Nnode)%dx    = P%dx
     node(Nnode)%D(:)  = P%D(:)
     node(Nnode)%Q(:)  = P%Q(:)
     node(Nnode)%QQ(:) = P%QQ(:)
     Narch = Narch + 1
     arch(Narch) = Nnode
     !
     IF ( (node(Nnode)%lx/param%sd_ratio + node(Nnode)%dx) >  &
          DSQRT((node(Nnode)%xx(2)-y)**2 + (node(Nnode)%xx(3)-z)**2) ) THEN
        !
        ! PP - twig
        IF (ASSOCIATED(P%Child1)) THEN
           CALL CopyTree4(P%Child1, arch, node, Nnode, Narch, NS, param, y, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child2)) THEN
           CALL CopyTree4(P%Child2, arch, node, Nnode, Narch, NS, param, y, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child3)) THEN
           CALL CopyTree4(P%Child3, arch, node, Nnode, Narch, NS, param, y, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child4)) THEN
           CALL CopyTree4(P%Child4, arch, node, Nnode, Narch, NS, param, y, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child5)) THEN
           CALL CopyTree4(P%Child5, arch, node, Nnode, Narch, NS, param, y, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child6)) THEN
           CALL CopyTree4(P%Child6, arch, node, Nnode, Narch, NS, param, y, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child7)) THEN
           CALL CopyTree4(P%Child7, arch, node, Nnode, Narch, NS, param, y, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child8)) THEN
           CALL CopyTree4(P%Child8, arch, node, Nnode, Narch, NS, param, y, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
     END IF
  ELSE
     !
     ! leaf
     node(Nnode)%qeff  = P%qeff
     node(Nnode)%xx(:) = P%xx(:)
     node(Nnode)%lx    = MAXVAL(P%lx(:))
     node(Nnode)%dx    = 0.0D0
     node(Nnode)%D(:)  = 0.D0
     node(Nnode)%Q(:)  = 0.D0
     node(Nnode)%QQ(:) = 0.D0
     Narch = Narch + 1
     arch(Narch) = Nnode
  END IF
END SUBROUTINE CopyTree4
!
RECURSIVE SUBROUTINE CopyTree5(P, arch, node, Nnode, Narch, NS, param, x, z)
  USE DATASTR
  IMPLICIT NONE
  !  
  TYPE(NM), INTENT(IN):: NS
  TYPE(GH), INTENT(INOUT):: node(NS%Nnode)
  TYPE(PM), INTENT(IN):: param
  INTEGER, INTENT(INOUT):: arch(NS%Nnode), Narch, Nnode
  REAL*8,  INTENT(IN):: x, z
  TYPE(BT), Pointer   :: P
  Nnode = Nnode + 1 ! position tag
  IF (P%ileaf == 0) THEN
     node(Nnode)%qeff  = P%qeff
     node(Nnode)%xx(:) = P%xx(:)
     node(Nnode)%lx    = MAXVAL(P%lx(:))
     node(Nnode)%dx    = P%dx
     node(Nnode)%D(:)  = P%D(:)
     node(Nnode)%Q(:)  = P%Q(:)
     node(Nnode)%QQ(:) = P%QQ(:)
     Narch = Narch + 1
     arch(Narch) = Nnode
     !
     IF ( (node(Nnode)%lx/param%sd_ratio + node(Nnode)%dx) >  &
          DSQRT((node(Nnode)%xx(1)-x)**2 + (node(Nnode)%xx(3)-z)**2) ) THEN
        !
        ! PP - twig
        IF (ASSOCIATED(P%Child1)) THEN
           CALL CopyTree5(P%Child1, arch, node, Nnode, Narch, NS, param, x, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child2)) THEN
           CALL CopyTree5(P%Child2, arch, node, Nnode, Narch, NS, param, x, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child3)) THEN
           CALL CopyTree5(P%Child3, arch, node, Nnode, Narch, NS, param, x, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child4)) THEN
           CALL CopyTree5(P%Child4, arch, node, Nnode, Narch, NS, param, x, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child5)) THEN
           CALL CopyTree5(P%Child5, arch, node, Nnode, Narch, NS, param, x, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child6)) THEN
           CALL CopyTree5(P%Child6, arch, node, Nnode, Narch, NS, param, x, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child7)) THEN
           CALL CopyTree5(P%Child7, arch, node, Nnode, Narch, NS, param, x, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child8)) THEN
           CALL CopyTree5(P%Child8, arch, node, Nnode, Narch, NS, param, x, z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
     END IF
  ELSE
     !
     ! leaf
     node(Nnode)%qeff  = P%qeff
     node(Nnode)%xx(:) = P%xx(:)
     node(Nnode)%lx    = MAXVAL(P%lx(:))
     node(Nnode)%dx    = 0.0D0
     node(Nnode)%D(:)  = 0.D0
     node(Nnode)%Q(:)  = 0.D0
     node(Nnode)%QQ(:) = 0.D0
     Narch = Narch + 1
     arch(Narch) = Nnode
  END IF
END SUBROUTINE CopyTree5
!
RECURSIVE SUBROUTINE CopyTree6(P, arch, node, Nnode, Narch, NS, param, x, y)
  USE DATASTR
  IMPLICIT NONE
  !  
  TYPE(NM), INTENT(IN):: NS
  TYPE(GH), INTENT(INOUT):: node(NS%Nnode)
  TYPE(PM), INTENT(IN):: param
  INTEGER, INTENT(INOUT):: arch(NS%Nnode), Narch, Nnode
  REAL*8,  INTENT(IN):: x, y
  TYPE(BT), Pointer   :: P
  Nnode = Nnode + 1 ! position tag
  IF (P%ileaf == 0) THEN
     node(Nnode)%qeff  = P%qeff
     node(Nnode)%xx(:) = P%xx(:)
     node(Nnode)%lx    = MAXVAL(P%lx(:))
     node(Nnode)%dx    = P%dx
     node(Nnode)%D(:)  = P%D(:)
     node(Nnode)%Q(:)  = P%Q(:)
     node(Nnode)%QQ(:) = P%QQ(:)
     Narch = Narch + 1
     arch(Narch) = Nnode
     !
     IF ( (node(Nnode)%lx/param%sd_ratio + node(Nnode)%dx) >  &
          DSQRT((node(Nnode)%xx(1)-x)**2 + (node(Nnode)%xx(2)-y)**2) ) THEN
        !
        ! PP - twig
        IF (ASSOCIATED(P%Child1)) THEN
           CALL CopyTree6(P%Child1, arch, node, Nnode, Narch, NS, param, x, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child2)) THEN
           CALL CopyTree6(P%Child2, arch, node, Nnode, Narch, NS, param, x, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child3)) THEN
           CALL CopyTree6(P%Child3, arch, node, Nnode, Narch, NS, param, x, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child4)) THEN
           CALL CopyTree6(P%Child4, arch, node, Nnode, Narch, NS, param, x, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child5)) THEN
           CALL CopyTree6(P%Child5, arch, node, Nnode, Narch, NS, param, x, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child6)) THEN
           CALL CopyTree6(P%Child6, arch, node, Nnode, Narch, NS, param, x, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child7)) THEN
           CALL CopyTree6(P%Child7, arch, node, Nnode, Narch, NS, param, x, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child8)) THEN
           CALL CopyTree6(P%Child8, arch, node, Nnode, Narch, NS, param, x, y)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
     END IF
  ELSE
     !
     ! leaf
     node(Nnode)%qeff  = P%qeff
     node(Nnode)%xx(:) = P%xx(:)
     node(Nnode)%lx    = MAXVAL(P%lx(:))
     node(Nnode)%dx    = 0.0D0
     node(Nnode)%D(:)  = 0.D0
     node(Nnode)%Q(:)  = 0.D0
     node(Nnode)%QQ(:) = 0.D0
     Narch = Narch + 1
     arch(Narch) = Nnode
  END IF
END SUBROUTINE CopyTree6
!
RECURSIVE SUBROUTINE CopyTree7(P, arch, node, Nnode, Narch, NS, param,x,y,z)
  USE DATASTR
  IMPLICIT NONE
  !  
  TYPE(NM), INTENT(IN):: NS
  TYPE(GH), INTENT(INOUT):: node(NS%Nnode)
  TYPE(PM), INTENT(IN):: param
  INTEGER, INTENT(INOUT):: arch(NS%Nnode), Narch, Nnode
  REAL*8,  INTENT(IN):: x, y, z
  TYPE(BT), Pointer   :: P
  Nnode = Nnode + 1 ! position tag
  IF (P%ileaf == 0) THEN
     node(Nnode)%qeff  = P%qeff
     node(Nnode)%xx(:) = P%xx(:)
     node(Nnode)%lx    = MAXVAL(P%lx(:))
     node(Nnode)%dx    = P%dx
     node(Nnode)%D(:)  = P%D(:)
     node(Nnode)%Q(:)  = P%Q(:)
     node(Nnode)%QQ(:) = P%QQ(:)
     Narch = Narch + 1
     arch(Narch) = Nnode
     !
     IF ( (node(Nnode)%lx/param%sd_ratio + node(Nnode)%dx) >  &
          DSQRT((node(Nnode)%xx(1)-x)**2 +(node(Nnode)%xx(2)-y)**2 + &
          (node(Nnode)%xx(3)-z)**2) ) THEN
        !
        ! PP - twig
        IF (ASSOCIATED(P%Child1)) THEN
           CALL CopyTree7(P%Child1, arch, node, Nnode, Narch, NS, param, x,y,z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child2)) THEN
           CALL CopyTree7(P%Child2, arch, node, Nnode, Narch, NS, param, x,y,z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child3)) THEN
           CALL CopyTree7(P%Child3, arch, node, Nnode, Narch, NS, param, x,y,z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child4)) THEN
           CALL CopyTree7(P%Child4, arch, node, Nnode, Narch, NS, param, x,y,z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child5)) THEN
           CALL CopyTree7(P%Child5, arch, node, Nnode, Narch, NS, param, x,y,z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child6)) THEN
           CALL CopyTree7(P%Child6, arch, node, Nnode, Narch, NS, param, x,y,z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child7)) THEN
           CALL CopyTree7(P%Child7, arch, node, Nnode, Narch, NS, param, x,y,z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
        IF (ASSOCIATED(P%Child8)) THEN
           CALL CopyTree7(P%Child8, arch, node, Nnode, Narch, NS, param, x,y,z)
           Narch = Narch + 1
           arch(Narch) = -1
        END IF
     END IF
  ELSE
     !
     ! leaf
     node(Nnode)%qeff  = P%qeff
     node(Nnode)%xx(:) = P%xx(:)
     node(Nnode)%lx    = MAXVAL(P%lx(:))
     node(Nnode)%dx    = 0.0D0
     node(Nnode)%D(:)  = 0.D0
     node(Nnode)%Q(:)  = 0.D0
     node(Nnode)%QQ(:) = 0.D0
     Narch = Narch + 1
     arch(Narch) = Nnode
  END IF
END SUBROUTINE CopyTree7
!
! ############## Rebuild ghost tree using array data ####################
! #######################################################################
RECURSIVE SUBROUTINE RebuildGhostTree(P, R, rarch, rnode, Nnode, Nrnod, Nrarc)
  USE DATASTR
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)   :: Nrnod, Nrarc, rarch(Nrarc)
  INTEGER, INTENT(INOUT):: Nnode
  TYPE(GH),INTENT(IN)     :: rnode(Nrnod)
  TYPE(GT), Pointer:: P, R
  Nnode = Nnode + 1
  IF (.NOT.ASSOCIATED(P)) THEN
     ALLOCATE(P)
     P%Root => R
     P%qeff  = rnode(rarch(Nnode))%qeff
     P%xx(:) = rnode(rarch(Nnode))%xx(:)
     P%lx    = rnode(rarch(Nnode))%lx
     P%dx    = rnode(rarch(Nnode))%dx
     P%D(:)  = rnode(rarch(Nnode))%D(:)
     P%Q(:)  = rnode(rarch(Nnode))%Q(:)
     P%QQ(:) = rnode(rarch(Nnode))%QQ(:)
  ELSE
     P%Root => R
     P%qeff  = rnode(rarch(Nnode))%qeff
     P%xx(:) = rnode(rarch(Nnode))%xx(:)
     P%lx    = rnode(rarch(Nnode))%lx
     P%dx    = rnode(rarch(Nnode))%dx
     P%D(:)  = rnode(rarch(Nnode))%D(:)
     P%Q(:)  = rnode(rarch(Nnode))%Q(:)
     P%QQ(:) = rnode(rarch(Nnode))%QQ(:)
  END IF
  
  IF (Nnode < Nrarc) THEN
     IF (rarch(Nnode+1) > 0) THEN
        P%ileaf = 0
        CALL RebuildGhostTree(P%Child1, P, rarch, rnode, Nnode, Nrnod, Nrarc)
        Nnode = Nnode + 1
        IF (rarch(Nnode+1) > 0) THEN
           CALL RebuildGhostTree(P%Child2, P, rarch, rnode, Nnode, Nrnod,Nrarc)
           Nnode = Nnode + 1
        END IF
        IF (rarch(Nnode+1) > 0) THEN
           CALL RebuildGhostTree(P%Child3, P, rarch, rnode, Nnode, Nrnod,Nrarc)
           Nnode = Nnode + 1
        END IF
        IF (rarch(Nnode+1) > 0) THEN
           CALL RebuildGhostTree(P%Child4, P, rarch, rnode, Nnode, Nrnod,Nrarc)
           Nnode = Nnode + 1
        END IF
        IF (rarch(Nnode+1) > 0) THEN
           CALL RebuildGhostTree(P%Child5, P, rarch, rnode, Nnode, Nrnod,Nrarc)
           Nnode = Nnode + 1
        END IF
        IF (rarch(Nnode+1) > 0) THEN
           CALL RebuildGhostTree(P%Child6, P, rarch, rnode, Nnode, Nrnod,Nrarc)
           Nnode = Nnode + 1
        END IF
        IF (rarch(Nnode+1) > 0) THEN
           CALL RebuildGhostTree(P%Child7, P, rarch, rnode, Nnode, Nrnod,Nrarc)
           Nnode = Nnode + 1
        END IF
        IF (rarch(Nnode+1) > 0) THEN
           CALL RebuildGhostTree(P%Child8, P, rarch, rnode, Nnode, Nrnod,Nrarc)
           Nnode = Nnode + 1
        END IF
     ELSE
        IF (P%qeff > 0.) THEN
           P%ileaf = 1
        ELSE
           P%ileaf = -1
        END IF
     END IF
  END IF
  !
END SUBROUTINE RebuildGhostTree
! Coulomb sum for Ghost Tree
! #########################
SUBROUTINE GhostCoulombSum_e(NS, qel, qion, sys, param, P)
  USE DATASTR
  IMPLICIT NONE
  !
  TYPE(NM), INTENT(INOUT):: NS
  TYPE(PT), INTENT(INOUT):: qel(NS%Nel), qion(NS%Nion)
  TYPE(ST), INTENT(INOUT):: sys
  TYPE(PM), INTENT(IN):: param
  TYPE(GT), Pointer:: P
  INTEGER :: i, Rei(Ndist+1)
  REAL*8  :: Epot
  Epot = 0.0D0
  Rei(:) = 0
  !$OMP PARALLEL DO SHARED(NS, P, qel, param) PRIVATE(i) REDUCTION(+:Epot)
  DO i=1, NS%Nel
     CALL Ghost_Tree_ee(NS, qel, i, param, P, Epot)
  END DO
  !$OMP END PARALLEL DO
  sys%Epot_e = sys%Epot_e + Epot
  Epot = 0.0
  !$OMP PARALLEL DO SHARED(NS, P, qion, param) PRIVATE(i) REDUCTION(+:Epot, Rei)
  DO i=1, NS%Nion
     CALL Ghost_Tree_ie(NS, qion, i, param, P, Epot, Rei)
  END DO
  !$OMP END PARALLEL DO
  sys%Epot_e = sys%Epot_e + Epot
  NS%RDF_ei(:) = NS%RDF_ei(:) + Rei(:)
  !
  RETURN
END SUBROUTINE GhostCoulombSum_e
!
! #########################
RECURSIVE SUBROUTINE Ghost_Tree_ee(NS, qel, i, param, P, Epot)
  USE DATASTR
  IMPLICIT NONE
  !
  TYPE(NM), INTENT(INOUT):: NS
  TYPE(PT), INTENT(INOUT):: qel(NS%Nel)
  TYPE(PM), INTENT(IN):: param
  INTEGER,  INTENT(IN) :: i
  TYPE(GT), Pointer:: P
  INTEGER :: k, id, Rei(Ndist+1)
  REAL*8  :: s, r, ff, xr(3), r2, R3, R5, R7, D(6), Q(10), Epot
  !
  IF (P%ileaf == 0) THEN
     !
     ! PP - twig
     s = P%lx
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qel(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     r = DSQRT(r2)
     IF ( (s/param%sd_ratio + P%dx) < r) THEN
        ff = -1.0D0*eps
        R3 = r2*r
        R5 = R3*r2 
        R7 = R5*r2
        D(1) = -1.D0/R3 + 3.D0*xr(1)**2/R5
        D(2) = -1.D0/R3 + 3.D0*xr(2)**2/R5
        D(3) = -1.D0/R3 + 3.D0*xr(3)**2/R5
        D(4) = 3.D0*xr(1)*xr(2)/R5
        D(5) = 3.D0*xr(2)*xr(3)/R5
        D(6) = 3.D0*xr(3)*xr(1)/R5
        Q(1) = 15.D0*xr(1)**3/R7 - 9.D0*xr(1)/R5
        Q(2) = 15.D0*xr(2)**3/R7 - 9.D0*xr(2)/R5
        Q(3) = 15.D0*xr(3)**3/R7 - 9.D0*xr(3)/R5
        Q(4) = 15.D0*xr(1)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(5) = 15.D0*xr(1)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(6) = 15.D0*xr(2)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(7) = 15.D0*xr(2)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(8) = 15.D0*xr(3)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(9) = 15.D0*xr(3)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(10)= 15.D0*xr(1)*xr(2)*xr(3)/R7
        Epot = Epot + 0.5D0*ff*(P%qeff/r + &
             (P%D(1)*xr(1)+P%D(2)*xr(2)+P%D(3)*xr(3))/R3 + &
             P%Q(1)*D(1) + P%Q(2)*D(2)+P%Q(3)*D(3) + &
             P%QQ(1)*D(4) + P%QQ(2)*D(5)+P%QQ(3)*D(6))
        qel(i)%ff(1) = qel(i)%ff(1) + ff*(P%qeff*xr(1)/R3 + &
             D(1)*P%D(1) + D(4)*P%D(2) + D(6)*P%D(3) + &
             Q(1)*P%Q(1) + Q(6)*P%Q(2) + Q(8)*P%Q(3) + &
             Q(4)*P%QQ(1) + Q(5)*P%QQ(3) + Q(10)*P%QQ(2))
        qel(i)%ff(2) = qel(i)%ff(2) + ff*(P%qeff*xr(2)/R3 + &
             D(4)*P%D(1) + D(2)*P%D(2) + D(5)*P%D(3) + &
             Q(4)*P%Q(1) + Q(2)*P%Q(2) + Q(9)*P%Q(3) + &
             Q(6)*P%QQ(1) + Q(10)*P%QQ(3) + Q(7)*P%QQ(2))
        qel(i)%ff(3) = qel(i)%ff(3) + ff*(P%qeff*xr(3)/R3 + &
             D(6)*P%D(1) + D(5)*P%D(2) + D(3)*P%D(3) + &
             Q(5)*P%Q(1) + Q(7)*P%Q(2) + Q(3)*P%Q(3) + &
             Q(10)*P%QQ(1) + Q(8)*P%QQ(3) + Q(9)*P%QQ(2))
     ELSE
        IF (ASSOCIATED(P%Child1)) &
             CALL Ghost_Tree_ee(NS, qel, i, param, P%Child1, Epot)
        IF (ASSOCIATED(P%Child2)) &
             CALL Ghost_Tree_ee(NS, qel, i, param, P%Child2, Epot)
        IF (ASSOCIATED(P%Child3)) &
             CALL Ghost_Tree_ee(NS, qel, i, param, P%Child3, Epot)
        IF (ASSOCIATED(P%Child4)) &
             CALL Ghost_Tree_ee(NS, qel, i, param, P%Child4, Epot)
        IF (ASSOCIATED(P%Child5)) &
             CALL Ghost_Tree_ee(NS, qel, i, param, P%Child5, Epot)
        IF (ASSOCIATED(P%Child6)) &
             CALL Ghost_Tree_ee(NS, qel, i, param, P%Child6, Epot)
        IF (ASSOCIATED(P%Child7)) &
             CALL Ghost_Tree_ee(NS, qel, i, param, P%Child7, Epot)
        IF (ASSOCIATED(P%Child8)) &
             CALL Ghost_Tree_ee(NS, qel, i, param, P%Child8, Epot)

     END IF
  ELSE
     !
     ! leaf of electron
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qel(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     r = DSQRT(r2)
     ff = -1.D0*eps*P%qeff
     !
     ! Generic Coulomb
     Epot = Epot + 0.5D0*ff/r
     DO k=1,3
        qel(i)%ff(k) = qel(i)%ff(k) + ff*xr(k)/r/r2
     END DO
  END IF
END SUBROUTINE Ghost_Tree_ee
!
! #########################
RECURSIVE SUBROUTINE Ghost_Tree_ie(NS, qion, i, param, P, Epot, Rei)
  USE DATASTR
  IMPLICIT NONE
  !
  TYPE(NM), INTENT(INOUT):: NS
  TYPE(PT), INTENT(INOUT):: qion(NS%Nion)
  TYPE(PM), INTENT(IN):: param
  INTEGER,  INTENT(IN):: i
  TYPE(GT), Pointer:: P
  INTEGER :: k, id, Rei(Ndist+1)
  REAL*8  :: s, r, ff, xr(3), r2, R3, R5, R7, D(6), Q(10), Epot
  !
  IF (P%ileaf == 0) THEN
     !
     ! PP - twig
     s = P%lx
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qion(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     r = DSQRT(r2)
     IF ( (s/param%sd_ratio + P%dx) < r) THEN
        ff = 1.D0*eps
        R3 = r2*r
        R5 = R3*r2
        R7 = R5*r2
        D(1) = -1.D0/R3 + 3.D0*xr(1)**2/R5
        D(2) = -1.D0/R3 + 3.D0*xr(2)**2/R5
        D(3) = -1.D0/R3 + 3.D0*xr(3)**2/R5
        D(4) = 3.D0*xr(1)*xr(2)/R5
        D(5) = 3.D0*xr(2)*xr(3)/R5
        D(6) = 3.D0*xr(3)*xr(1)/R5
        Q(1) = 15.D0*xr(1)**3/R7 - 9.D0*xr(1)/R5
        Q(2) = 15.D0*xr(2)**3/R7 - 9.D0*xr(2)/R5
        Q(3) = 15.D0*xr(3)**3/R7 - 9.D0*xr(3)/R5
        Q(4) = 15.D0*xr(1)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(5) = 15.D0*xr(1)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(6) = 15.D0*xr(2)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(7) = 15.D0*xr(2)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(8) = 15.D0*xr(3)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(9) = 15.D0*xr(3)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(10)= 15.D0*xr(1)*xr(2)*xr(3)/R7
        Epot = Epot + 0.5D0*ff*(P%qeff/r + &
             (P%D(1)*xr(1)+P%D(2)*xr(2)+P%D(3)*xr(3))/R3 + &
             P%Q(1)*D(1) + P%Q(2)*D(2)+P%Q(3)*D(3) + &
             P%QQ(1)*D(4) + P%QQ(2)*D(5)+P%QQ(3)*D(6))
        qion(i)%ff(1) = qion(i)%ff(1) + ff*(P%qeff*xr(1)/R3 + &
             D(1)*P%D(1) + D(4)*P%D(2) + D(6)*P%D(3) + &
             Q(1)*P%Q(1) + Q(6)*P%Q(2) + Q(8)*P%Q(3) + &
             Q(4)*P%QQ(1) + Q(5)*P%QQ(3) + Q(10)*P%QQ(2))
        qion(i)%ff(2) = qion(i)%ff(2) + ff*(P%qeff*xr(2)/R3 + &
             D(4)*P%D(1) + D(2)*P%D(2) + D(5)*P%D(3) + &
             Q(4)*P%Q(1) + Q(2)*P%Q(2) + Q(9)*P%Q(3) + &
             Q(6)*P%QQ(1) + Q(10)*P%QQ(3) + Q(7)*P%QQ(2))
        qion(i)%ff(3) = qion(i)%ff(3) + ff*(P%qeff*xr(3)/R3 + &
             D(6)*P%D(1) + D(5)*P%D(2) + D(3)*P%D(3) + &
             Q(5)*P%Q(1) + Q(7)*P%Q(2) + Q(3)*P%Q(3) + &
             Q(10)*P%QQ(1) + Q(8)*P%QQ(3) + Q(9)*P%QQ(2))
     ELSE
        IF (ASSOCIATED(P%Child1)) &
             CALL Ghost_Tree_ie(NS, qion, i, param, P%Child1, Epot, Rei)
        IF (ASSOCIATED(P%Child2)) &
             CALL Ghost_Tree_ie(NS, qion, i, param, P%Child2, Epot, Rei)
        IF (ASSOCIATED(P%Child3)) &
             CALL Ghost_Tree_ie(NS, qion, i, param, P%Child3, Epot, Rei)
        IF (ASSOCIATED(P%Child4)) &
             CALL Ghost_Tree_ie(NS, qion, i, param, P%Child4, Epot, Rei)
        IF (ASSOCIATED(P%Child5)) &
             CALL Ghost_Tree_ie(NS, qion, i, param, P%Child5, Epot, Rei)
        IF (ASSOCIATED(P%Child6)) &
             CALL Ghost_Tree_ie(NS, qion, i, param, P%Child6, Epot, Rei)
        IF (ASSOCIATED(P%Child7)) &
             CALL Ghost_Tree_ie(NS, qion, i, param, P%Child7, Epot, Rei)
        IF (ASSOCIATED(P%Child8)) &
             CALL Ghost_Tree_ie(NS, qion, i, param, P%Child8, Epot, Rei)
     END IF
  ELSE
     !
     ! leaf of electron
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qion(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     r = DSQRT(r2)
     ff = 1.D0*eps*P%qeff
     !
     ! Kelbg potential
     Epot = Epot + 0.5D0*ff*(1.-EXP(-r/param%rs))/r
     DO k=1,3
        qion(i)%ff(k) = qion(i)%ff(k) + &
             ff*xr(k)*(1.D0/r - EXP(-r/param%rs)*(1.D0/r+1.D0/param%rs))/r2
     END DO
     !
     ! RDF
     id = INT(r/NS%dx_ei)+1
     IF (id <= Ndist + 1 .AND. id > 0 ) THEN
        Rei(id) = Rei(id) + 1
     END IF
  END IF
END SUBROUTINE Ghost_Tree_ie
! Coulomb sum for Ghost Tree
! #########################
SUBROUTINE GhostCoulombSum_i(NS, qel, qion, sys, param, P)
  USE DATASTR
  IMPLICIT NONE
  !
  TYPE(NM), INTENT(INOUT):: NS
  TYPE(PT), INTENT(INOUT):: qel(NS%Nel), qion(NS%Nion)
  TYPE(ST), INTENT(INOUT):: sys
  TYPE(PM), INTENT(IN):: param
  TYPE(GT), Pointer:: P
  INTEGER :: i, Rei(Ndist+1), Rii(Ndist+1)
  REAL*8  :: Epot
  Epot = 0.0D0
  Rei = 0
  !$OMP PARALLEL DO SHARED(NS, P, qel, param) PRIVATE(i) REDUCTION(+:Epot,Rei)
  DO i=1, NS%Nel
     CALL Ghost_Tree_ei(NS, qel, i, param, P, Epot, Rei)
  END DO
  !$OMP END PARALLEL DO
  sys%Epot_e = sys%Epot_e + Epot
  NS%RDF_ei(:) = NS%RDF_ei(:) + Rei(:)
  Epot = 0.0D0
  Rii = 0
  !$OMP PARALLEL DO SHARED(NS, P, qion, param) PRIVATE(i) REDUCTION(+:Epot,Rii)
  DO i=1, NS%Nion
     CALL Ghost_Tree_ii(NS, qion, i, param, P, Epot, Rii)
  END DO
  !$OMP END PARALLEL DO
  sys%Epot_i = sys%Epot_i + Epot
  NS%RDF_ii(:) = NS%RDF_ii(:) + Rii(:)
  !
  RETURN
END SUBROUTINE GhostCoulombSum_i
!
! #########################
RECURSIVE SUBROUTINE Ghost_Tree_ei(NS, qel, i, param, P, Epot, Rei)
  USE DATASTR
  IMPLICIT NONE
  !
  TYPE(NM), INTENT(INOUT):: NS
  TYPE(PT), INTENT(INOUT):: qel(NS%Nel)
  TYPE(PM), INTENT(IN):: param
  INTEGER,  INTENT(IN) :: i
  TYPE(GT), Pointer:: P
  INTEGER :: k, id, Rei(Ndist+1)
  REAL*8  :: s, r, ff, xr(3), r2, R3, R5, R7, D(6), Q(10), Epot
  !
  IF (P%ileaf == 0) THEN
     !
     ! PP - twig
     s = P%lx
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qel(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     r = DSQRT(r2)
     IF ( (s/param%sd_ratio + P%dx) < r) THEN
        ff = -1.0D0*eps
        R3 = r2*r
        R5 = R3*r2
        R7 = R5*r2
        D(1) = -1.D0/R3 + 3.D0*xr(1)**2/R5
        D(2) = -1.D0/R3 + 3.D0*xr(2)**2/R5
        D(3) = -1.D0/R3 + 3.D0*xr(3)**2/R5
        D(4) = 3.D0*xr(1)*xr(2)/R5
        D(5) = 3.D0*xr(2)*xr(3)/R5
        D(6) = 3.D0*xr(3)*xr(1)/R5
        Q(1) = 15.D0*xr(1)**3/R7 - 9.D0*xr(1)/R5
        Q(2) = 15.D0*xr(2)**3/R7 - 9.D0*xr(2)/R5
        Q(3) = 15.D0*xr(3)**3/R7 - 9.D0*xr(3)/R5
        Q(4) = 15.D0*xr(1)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(5) = 15.D0*xr(1)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(6) = 15.D0*xr(2)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(7) = 15.D0*xr(2)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(8) = 15.D0*xr(3)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(9) = 15.D0*xr(3)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(10)= 15.D0*xr(1)*xr(2)*xr(3)/R7
        Epot = Epot + 0.5D0*ff*(P%qeff/r + &
             (P%D(1)*xr(1)+P%D(2)*xr(2)+P%D(3)*xr(3))/R3 + &
             P%Q(1)*D(1) + P%Q(2)*D(2)+P%Q(3)*D(3) + &
             P%QQ(1)*D(4) + P%QQ(2)*D(5)+P%QQ(3)*D(6))
        qel(i)%ff(1) = qel(i)%ff(1) + ff*(P%qeff*xr(1)/R3 + &
             D(1)*P%D(1) + D(4)*P%D(2) + D(6)*P%D(3) + &
             Q(1)*P%Q(1) + Q(6)*P%Q(2) + Q(8)*P%Q(3) + &
             Q(4)*P%QQ(1) + Q(5)*P%QQ(3) + Q(10)*P%QQ(2))
        qel(i)%ff(2) = qel(i)%ff(2) + ff*(P%qeff*xr(2)/R3 + &
             D(4)*P%D(1) + D(2)*P%D(2) + D(5)*P%D(3) + &
             Q(4)*P%Q(1) + Q(2)*P%Q(2) + Q(9)*P%Q(3) + &
             Q(6)*P%QQ(1) + Q(10)*P%QQ(3) + Q(7)*P%QQ(2))
        qel(i)%ff(3) = qel(i)%ff(3) + ff*(P%qeff*xr(3)/R3 + &
             D(6)*P%D(1) + D(5)*P%D(2) + D(3)*P%D(3) + &
             Q(5)*P%Q(1) + Q(7)*P%Q(2) + Q(3)*P%Q(3) + &
             Q(10)*P%QQ(1) + Q(8)*P%QQ(3) + Q(9)*P%QQ(2))
     ELSE
        IF (ASSOCIATED(P%Child1)) &
             CALL Ghost_Tree_ei(NS, qel, i, param, P%Child1, Epot, Rei)
        IF (ASSOCIATED(P%Child2)) &
             CALL Ghost_Tree_ei(NS, qel, i, param, P%Child2, Epot, Rei)
        IF (ASSOCIATED(P%Child3)) &
             CALL Ghost_Tree_ei(NS, qel, i, param, P%Child3, Epot, Rei)
        IF (ASSOCIATED(P%Child4)) &
             CALL Ghost_Tree_ei(NS, qel, i, param, P%Child4, Epot, Rei)
        IF (ASSOCIATED(P%Child5)) &
             CALL Ghost_Tree_ei(NS, qel, i, param, P%Child5, Epot, Rei)
        IF (ASSOCIATED(P%Child6)) &
             CALL Ghost_Tree_ei(NS, qel, i, param, P%Child6, Epot, Rei)
        IF (ASSOCIATED(P%Child7)) &
             CALL Ghost_Tree_ei(NS, qel, i, param, P%Child7, Epot, Rei)
        IF (ASSOCIATED(P%Child8)) &
             CALL Ghost_Tree_ei(NS, qel, i, param, P%Child8, Epot, Rei)
     END IF
  ELSE 
     !
     ! leaf of ion
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qel(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     r = DSQRT(r2)
     ff = -1.D0*eps*P%qeff
     !
     ! Kelbg potential
     Epot = Epot +0.5D0*ff*(1.-EXP(-r/param%rs))/r
     DO k=1,3
        qel(i)%ff(k) = qel(i)%ff(k) + &
             ff*xr(k)*(1.D0/r - EXP(-r/param%rs)*(1.D0/r+1.D0/param%rs))/r2
     END DO
     !
     ! RDF
     id = INT(r/NS%dx_ei)+1
     IF (id <= Ndist + 1 .AND. id > 0 ) THEN
        Rei(id) = Rei(id) + 1
     END IF
     !
  END IF
END SUBROUTINE Ghost_Tree_ei
!
! #########################
RECURSIVE SUBROUTINE Ghost_Tree_ii(NS, qion, i, param, P, Epot, Rii)
  USE DATASTR
  IMPLICIT NONE
  !
  TYPE(NM), INTENT(INOUT):: NS
  TYPE(PT), INTENT(INOUT):: qion(NS%Nion)
  TYPE(PM), INTENT(IN):: param
  INTEGER,  INTENT(IN):: i
  TYPE(GT), Pointer:: P
  INTEGER :: k, id, Rii(Ndist+1)
  REAL*8  :: s, r, ff, xr(3), r2, R3, R5, R7, D(6), Q(10), Epot
  !
  IF (P%ileaf == 0) THEN
     !
     ! PP - twig
     s = P%lx
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qion(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     r = DSQRT(r2)
     IF ( (s/param%sd_ratio + P%dx) < r) THEN
        ff = 1.D0*eps
        R3 = r2*r
        R5 = R3*r2
        R7 = R5*r2
        D(1) = -1.D0/R3 + 3.D0*xr(1)**2/R5
        D(2) = -1.D0/R3 + 3.D0*xr(2)**2/R5
        D(3) = -1.D0/R3 + 3.D0*xr(3)**2/R5
        D(4) = 3.D0*xr(1)*xr(2)/R5
        D(5) = 3.D0*xr(2)*xr(3)/R5
        D(6) = 3.D0*xr(3)*xr(1)/R5
        Q(1) = 15.D0*xr(1)**3/R7 - 9.D0*xr(1)/R5
        Q(2) = 15.D0*xr(2)**3/R7 - 9.D0*xr(2)/R5
        Q(3) = 15.D0*xr(3)**3/R7 - 9.D0*xr(3)/R5
        Q(4) = 15.D0*xr(1)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(5) = 15.D0*xr(1)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(6) = 15.D0*xr(2)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(7) = 15.D0*xr(2)**2*xr(3)/R7 - 3.D0*xr(3)/R5
        Q(8) = 15.D0*xr(3)**2*xr(1)/R7 - 3.D0*xr(1)/R5
        Q(9) = 15.D0*xr(3)**2*xr(2)/R7 - 3.D0*xr(2)/R5
        Q(10)= 15.D0*xr(1)*xr(2)*xr(3)/R7
        Epot = Epot + 0.5D0*ff*(P%qeff/r + &
             (P%D(1)*xr(1)+P%D(2)*xr(2)+P%D(3)*xr(3))/R3 + &
             P%Q(1)*D(1) + P%Q(2)*D(2)+P%Q(3)*D(3) + &
             P%QQ(1)*D(4) + P%QQ(2)*D(5)+P%QQ(3)*D(6))
        qion(i)%ff(1) = qion(i)%ff(1) + ff*(P%qeff*xr(1)/R3 + &
             D(1)*P%D(1) + D(4)*P%D(2) + D(6)*P%D(3) + &
             Q(1)*P%Q(1) + Q(6)*P%Q(2) + Q(8)*P%Q(3) + &
             Q(4)*P%QQ(1) + Q(5)*P%QQ(3) + Q(10)*P%QQ(2))
        qion(i)%ff(2) = qion(i)%ff(2) + ff*(P%qeff*xr(2)/R3 + &
             D(4)*P%D(1) + D(2)*P%D(2) + D(5)*P%D(3) + &
             Q(4)*P%Q(1) + Q(2)*P%Q(2) + Q(9)*P%Q(3) + &
             Q(6)*P%QQ(1) + Q(10)*P%QQ(3) + Q(7)*P%QQ(2))
        qion(i)%ff(3) = qion(i)%ff(3) + ff*(P%qeff*xr(3)/R3 + &
             D(6)*P%D(1) + D(5)*P%D(2) + D(3)*P%D(3) + &
             Q(5)*P%Q(1) + Q(7)*P%Q(2) + Q(3)*P%Q(3) + &
             Q(10)*P%QQ(1) + Q(8)*P%QQ(3) + Q(9)*P%QQ(2))
     ELSE
        IF (ASSOCIATED(P%Child1)) &
             CALL Ghost_Tree_ii(NS, qion, i, param, P%Child1, Epot, Rii)
        IF (ASSOCIATED(P%Child2)) &
             CALL Ghost_Tree_ii(NS, qion, i, param, P%Child2, Epot, Rii)
        IF (ASSOCIATED(P%Child3)) &
             CALL Ghost_Tree_ii(NS, qion, i, param, P%Child3, Epot, Rii)
        IF (ASSOCIATED(P%Child4)) &
             CALL Ghost_Tree_ii(NS, qion, i, param, P%Child4, Epot, Rii)
        IF (ASSOCIATED(P%Child5)) &
             CALL Ghost_Tree_ii(NS, qion, i, param, P%Child5, Epot, Rii)
        IF (ASSOCIATED(P%Child6)) &
             CALL Ghost_Tree_ii(NS, qion, i, param, P%Child6, Epot, Rii)
        IF (ASSOCIATED(P%Child7)) &
             CALL Ghost_Tree_ii(NS, qion, i, param, P%Child7, Epot, Rii)
        IF (ASSOCIATED(P%Child8)) &
             CALL Ghost_Tree_ii(NS, qion, i, param, P%Child8, Epot, Rii)
     END IF
  ELSE
     !
     ! leaf of ion
     r2 = 0.0D0
     DO k=1, 3
        xr(k) = qion(i)%xx(k)-P%xx(k)
        r2 = r2 + xr(k)**2
     END DO
     r = DSQRT(r2)
     ff = 1.D0*eps*P%qeff
     !
     ! Generic Coulomb
     Epot = Epot + 0.5D0*ff/r
     DO k=1,3
        qion(i)%ff(k) = qion(i)%ff(k) + ff*xr(k)/r/r2
     END DO
     !
     ! RDF
     id = INT(r/NS%dx_ii)+1
     IF (id <= Ndist + 1 .AND. id > 0 ) THEN
        Rii(id) = Rii(id) + 1
     END IF
     !
  END IF
END SUBROUTINE Ghost_Tree_ii
END MODULE BINTREE
