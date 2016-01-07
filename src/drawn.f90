FUNCTION DRAWN(nn)

  !Routine to draw from a distribution, either Poisson if N<N1
  !or Gaussian if N>N1

  USE nr, ONLY : gasdev; USE nrtype
  USE pixcmd_vars; USE pixcmd_utils, ONLY : mypoidev
  IMPLICIT NONE
  
  REAL(SP), INTENT(in) :: nn
  REAL(SP), DIMENSION(npix,npix) :: drawn
  REAL(SP), DIMENSION(npix) :: gdev
  INTEGER :: i

  !----------------------------------------------------------------!

  IF (nn.LE.maxpoidev) THEN
     drawn = mypoidev(nn)
  ELSE
     DO i=1,npix
        CALL GASDEV(gdev)
        drawn(i,:) = gdev*SQRT(nn)+nn
     ENDDO
  ENDIF


END FUNCTION DRAWN
