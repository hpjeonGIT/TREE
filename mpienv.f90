!
! #############################################################################
! MPI coordinate/data type/communicating pair configuration ###################
SUBROUTINE MPI_initialize(COMM, NS)
USE DATASTR
USE MPI_ENV
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE(NM):: NS
TYPE(MP):: COMM
CHARACTER(len=256):: hostname
LOGICAL:: PERIOD, REORDER
INTEGER:: IERR, MYID, offset(2), Oldtype(2), blockcnt(2), i, j, extent, &
     Nhead, COLORR, KEYY, Ncrit
!
CALL MPI_INIT(IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, COMM%NUMPROC, IERR)
IF (COMM%NUMPROC /= Npe*8) STOP "========= Available processor error ========="
CALL MPI_GET_PROCESSOR_NAME(hostname, i,COMM%rc)
PERIOD  = .TRUE.
REORDER = .TRUE.
CALL MPI_CART_CREATE(MPI_COMM_WORLD, 1, COMM%NUMPROC, PERIOD, REORDER, &
     COMM%COMM1D, IERR)
CALL MPI_COMM_RANK(COMM%COMM1D, COMM%MYID2, IERR)
!
! New communicator for nodes inside of cluster
COLORR = INT(COMM%MYID2/Npe)
KEYY = COLORR
CALL MPI_COMM_SPLIT(COMM%COMM1D, COLORR, KEYY, COMM%LOCAL, IERR)
CALL MPI_COMM_RANK(COMM%LOCAL, COMM%LID, IERR)
!
! New data type for moving particles
offset(:) = 0
Oldtype(1) = MPI_REAL8
blockcnt(1) = 9
Oldtype(2) = MPI_INTEGER
blockcnt(2) = 1
CALL MPI_TYPE_EXTENT(MPI_INTEGER, extent, IERR)
offset(2) = extent
CALL MPI_TYPE_STRUCT(2, blockcnt, offset, Oldtype, COMM%MOVING, IERR)
CALL MPI_TYPE_COMMIT(COMM%MOVING, IERR)
!
! New data type for GHOST TREE information
blockcnt(1) = 15
!blockcnt(2) = 1
!Oldtype(2) = MPI_LOGICAL
CALL MPI_TYPE_EXTENT(MPI_LOGICAL, extent, IERR)
!offset(2) = extent
CALL MPI_TYPE_STRUCT(1, blockcnt(1), offset(1), Oldtype(1), COMM%GTREE, IERR)
CALL MPI_TYPE_COMMIT(COMM%GTREE, IERR)
!
! New data type for POS particle information
blockcnt(1) = 3
CALL MPI_TYPE_STRUCT(1, blockcnt(1), offset(1), Oldtype(1), COMM%POS, IERR)
CALL MPI_TYPE_COMMIT(COMM%POS, IERR)
!
Nhead = INT(COMM%MYID2/Npe)*Npe
DO i=1, 7
   COMM%DEST(i) = MOD(Nhead + Npe*i, 8*Npe)
   COMM%SRC(i)  = MOD(Nhead - Npe*i + 8*Npe, 8*Npe)
   j = COMM%DEST(i)
   CALL CriterionFind(j, Nhead, Ncrit)
   NS%Ncrit(i) = Ncrit
END DO
RETURN
END SUBROUTINE MPI_INITIALIZE
