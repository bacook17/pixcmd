FUNCTION GETMODEL(inpos)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : convolve,hist_2d,add_obs_err
  USE nr, ONLY : locate
  IMPLICIT NONE

  REAL(SP), DIMENSION(npar) :: inpos
  REAL(SP), DIMENSION(nx,ny) :: getmodel
  REAL(SP), DIMENSION(npix,npix,nfil) :: f1,cf1,of1
  INTEGER :: ilo,i,j,k
  REAL(SP) :: di

  !------------------------------------------------------------!

  !linearly combine the model components
  f1 = 0.0
  k=1
  DO j=1,nage
     DO i=1,nz
        f1 = f1+10**inpos(k)*model(:,:,:,i,j)
        k=k+1
     ENDDO
  ENDDO

  !convolve with PSF
  cf1 = -2.5*LOG10(convolve(f1))

  !add obs errors
  of1 = add_obs_err(cf1)

  !compute Hess diagram, normalize to unity
  getmodel = hist_2d(of1(:,:,1)-of1(:,:,2),of1(:,:,2),&
       xhess,yhess,nx,ny,npix) / npix**2


END FUNCTION GETMODEL
