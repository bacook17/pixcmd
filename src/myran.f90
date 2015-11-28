FUNCTION MYRAN()

  !turns ran1 into a function, rather than a subroutine

  USE nr, ONLY : ran1; USE nrtype
  IMPLICIT NONE

  REAL(SP) :: myran

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  CALL ran1(myran)

END FUNCTION MYRAN
