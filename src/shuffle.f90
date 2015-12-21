SUBROUTINE SHUFFLE(arr)

  USE pixcmd_utils, ONLY : myran

  IMPLICIT NONE
  INTEGER, DIMENSION(:), INTENT(inout) :: arr
  INTEGER :: i,nn,j,tmp

 !------------------------------------------------------------!

  nn = SIZE(arr)

  DO i=nn,1,-1
     j = INT(myran()*(i-1))+1
     tmp = arr(i)
     arr(i) = arr(j)
     arr(j) = tmp
  ENDDO


END SUBROUTINE SHUFFLE
