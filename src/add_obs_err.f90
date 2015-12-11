FUNCTION ADD_OBS_ERR(flux)

  ! Routine to generate a realization of observational
  ! data given an input exposure time.  Converts from
  ! mags to counts and draws from a Poisson/Gaussian distribution

  USE pixcmd_vars; USE nrtype
  USE nr, ONLY : poidev, gasdev
  IMPLICIT NONE

  REAL(SP), DIMENSION(npix,npix,nfil), INTENT(in) :: flux
  REAL(SP), DIMENSION(npix,npix,nfil) :: add_obs_err
  REAL(SP), DIMENSION(npix,npix) :: cts,cti
  REAL(SP), DIMENSION(npix) :: gdev
  INTEGER :: i,j,k

  !------------------------------------------------------------!

  !convert to counts
  DO k=1,nfil

     !compute total counts 
     cts = flux(:,:,k) * 10**(-2./5*(dm-zpt(k)))*exptime(k)
     
     DO j=1,npix
        CALL GASDEV(gdev)
        DO i=1,npix
           IF (cts(i,j).LT.0.) CYCLE
           IF (cts(i,j).LT.100.) THEN
              !draw from a Poisson distribution
              cti(i,j) = poidev(cts(i,j))
           ELSE
              !draw from a Guassian distribution
              !this is much faster and OK for N>>1
              cti(i,j) = cts(i,j) + gdev(i)*SQRT(cts(i,j))
           ENDIF
        ENDDO
     ENDDO
     
     !convert counts back to abs mags
     add_obs_err(:,:,k) = -2.5*LOG10(cti/exptime(k)) + zpt(k) - dm

  ENDDO
  

END FUNCTION ADD_OBS_ERR
