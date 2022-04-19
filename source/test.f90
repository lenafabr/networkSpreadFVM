PROGRAM MAIN
  ! minimalist fortran syntax testing
  USE GENUTIL
  IMPLICIT NONE

  CHARACTER(LEN=100) :: STR
  INTEGER :: INT
  INTEGER :: LIST(10), ARRAY(3,3)
  DOUBLE PRECISION :: LIST2(10)
  INTEGER :: INDS(2),VALS(2)
  INTEGER :: VAL, IND, c
  
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
  DO  C = 1,1000
     CALL RANDSELECT_INT(LIST,2,.FALSE.,VALS,INDS,LIST2)
     PRINT*, C, INDS
  ENDDO
END PROGRAM MAIN
