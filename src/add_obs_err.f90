FUNCTION ADD_OBS_ERR(flux)

  ! Routine to generate a realization of observational
  ! data given an input exposure time.  Converts from
  ! mags to counts and draws from a Poisson/Gaussian distribution

  USE pixcmd_vars; USE nrtype
  USE nr, ONLY : poidev
  IMPLICIT NONE

  REAL(SP), DIMENSION(npix,npix,nfil), INTENT(in) :: flux
  REAL(SP), DIMENSION(npix,npix,nfil) :: add_obs_err
  REAL(SP), DIMENSION(npix,npix) :: cts,cti,tgd
  INTEGER :: i,j,k

  !------------------------------------------------------------!

  add_obs_err = 0.0

  !convert to counts
  DO k=1,nfil
     
     !we need to separate Gaussian deviate arrays
     !so that the colors get a proper treatment
     IF (k.EQ.1) tgd = gdev
     IF (k.EQ.2) tgd = gdev2

     !compute total counts 
     cts = flux(:,:,k) * 10**(-2./5*(dm-zpt(k)))*exptime(k)
     cti = 0.0

     DO j=1,npix
        DO i=1,npix
           IF (cts(i,j).LT.0.) CYCLE
           IF (cts(i,j).LT.100.) THEN
              !draw from a Poisson distribution
              !cts is almost always >100, unless
              !Mpix<1.  Note that this poidev is a 
              !source of model stochasticity
              cti(i,j) = poidev(cts(i,j))
           ELSE
              !draw from a Gaussian distribution
              !this is much faster and OK for N>>1
              cti(i,j) = cts(i,j) + tgd(i,j)*SQRT(cts(i,j))
           ENDIF
        ENDDO
     ENDDO
     
     !convert to abs mags
     add_obs_err(:,:,k) = -2.5*LOG10(cti/exptime(k)) + zpt(k) - dm

  ENDDO
  

END FUNCTION ADD_OBS_ERR
