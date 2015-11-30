FUNCTION ADD_OBS_ERR(flux,dm,exptime,zpt)

  ! Routine to generate a realization of observational
  ! data given an input exposure time.  Converts from
  ! mags to counts and draws from a Poisson distribution

  USE pixcmd_vars; USE nrtype
  USE nr, ONLY : poidev
  IMPLICIT NONE

  REAL(SP), DIMENSION(npix,npix,nfil), INTENT(in) :: flux
  REAL(SP), INTENT(in) :: dm
  REAL(SP), DIMENSION(nfil), INTENT(in) :: exptime,zpt
  REAL(SP), DIMENSION(npix,npix,nfil) :: add_obs_err
  REAL(SP), DIMENSION(npix,npix) :: cts,cti
  INTEGER :: i,j,k

  !------------------------------------------------------------!

  !convert to counts
  DO k=1,nfil

     !compute total counts 
     cts = 10**(-2./5*(flux(:,:,k)+dm-zpt(k)))*exptime(k)
     
     DO j=1,npix
        DO i=1,npix
           !draw from a Poisson distribution
           cti(i,j) = poidev(cts(i,j))
        ENDDO
     ENDDO
     
     !convert counts back to abs mags
     add_obs_err(:,:,k) = -2.5*LOG10(cti/exptime(k)) + zpt(k) - dm

  ENDDO
  

END FUNCTION ADD_OBS_ERR
