FUNCTION HIST_2D(xx,yy,xarr,yarr,nx,ny,npix)

  USE nr, ONLY : locate; USE nrtype
  IMPLICIT NONE

  INTEGER, INTENT(in) :: nx,ny,npix
  REAL(SP), DIMENSION(npix,npix), INTENT(in) :: xx,yy
  REAL(SP), DIMENSION(nx), INTENT(in) :: xarr
  REAL(SP), DIMENSION(ny), INTENT(in) :: yarr
  REAL(SP), DIMENSION(nx,ny) :: hist_2d
  INTEGER :: i,j,i1,i2

  !------------------------------------------------------------!

  hist_2d = 0.0
  DO i=1,npix
     DO j=1,npix
        IF (xx(i,j).LT.xarr(1).OR.xx(i,j).GT.xarr(nx).OR.&
             yy(i,j).LT.yarr(1).OR.yy(i,j).GT.yarr(ny)) CYCLE
        i1 = locate(xarr,xx(i,j))
        i2 = locate(yarr,yy(i,j))
        hist_2d(i1,i2) = hist_2d(i1,i2)+1.
     ENDDO
  ENDDO


END FUNCTION HIST_2D
