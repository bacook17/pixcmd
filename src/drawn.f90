FUNCTION DRAWN(nn)

  !Routine to draw from a distribution, either Poisson if N<N1
  !or Gaussian if N>N1

  USE nr, ONLY : poidev,gasdev; USE nrtype

  IMPLICIT NONE
  
  REAL(SP), INTENT(in) :: nn
  REAL(SP) :: drawn,gdev
  INTEGER :: n1=10

  !----------------------------------------------------------------!

  IF (nn.LE.n1) THEN
     drawn = poidev(nn)
  ELSE
     CALL GASDEV(gdev)
     drawn = gdev*SQRT(nn)+nn
  ENDIF


END FUNCTION DRAWN
