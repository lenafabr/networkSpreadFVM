PROGRAM MAIN
  ! minimalist fortran syntax testing
  USE GENUTIL
  USE MT19937, only : randunifcircle
  IMPLICIT NONE

  CHARACTER(LEN=100) :: STR
  INTEGER :: INT
  INTEGER :: LIST(10), ARRAY(3,3)
  DOUBLE PRECISION :: LIST2(10)
  INTEGER :: INDS(2)
  INTEGER :: VAL, IND, c
  DOUBLE PRECISION :: VALS(10000,2)
  
  LIST = (/(INT,INT=1,10)/)
  LIST2 = LIST**2

  !LIST(LIST2>8) = 5
  
  !PRINT*, LIST

  ! Test weighted random selection
!  DO C = 1,10000
!     CALL RANDSELECT1_INT(LIST,VAL,IND,LIST2)
!     PRINT*, C, IND, VAL
!  ENDDO

  ! Test weighted random selection without replacement
  ! DO  C = 1,1000
  !    CALL RANDSELECT_INT(LIST,2,.FALSE.,VALS,INDS,LIST2)
  !    PRINT*, C, INDS
  ! ENDDO

  ! Test sampling in a circle
  VALS = RANDUNIFCIRCLE(10000,(/1D0,2D0/),3D0)
  DO C = 1,10000
     PRINT*, VALS(C,:)
  ENDDO
END PROGRAM MAIN
