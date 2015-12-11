FUNCTION CONVOLVE(arr)

  ! Routine to convolve a 2D image with a PSF
  
  USE pixcmd_vars
  USE nrtype;  USE nr, ONLY : four2, rlft2
  IMPLICIT NONE

  REAL(SP), DIMENSION(npix,npix,nfil), INTENT(inout) :: arr
  REAL(SP), DIMENSION(npix+2*npsf,npix+2*npsf,nfil) :: padarr
  REAL(SP), DIMENSION(npix,npix,nfil) :: convolve
  INTEGER :: i,j,k
  
  REAL(SP), DIMENSION(npix,npix) :: padarr2=0.,padarr3=0.
  COMPLEX(SPC), DIMENSION(npix/2,npix) :: spec_d,spec_out
  COMPLEX(SPC), DIMENSION(npix) :: speq_d,speq_out
  COMPLEX(SPC), DIMENSION(npsf2/2,npsf2) :: spec_p
  COMPLEX(SPC), DIMENSION(npsf2)      :: speq_p

  !------------------------------------------------------------!

  IF (fft.EQ.0) THEN

     !this is the brute force do loop convolution

     padarr(npsf+1:npsf+npix,npsf+1:npsf+npix,:) = arr

     !first pad Y
     padarr(npsf+1:npsf+npix,1:npsf,:) = arr(:,npix-npsf+1:npix,:)
     padarr(npsf+1:npsf+npix,npix+npsf+1:npix+2*npsf,:) = arr(:,1:npsf,:)
     !then pad X
     padarr(1:npsf,:,:) = padarr(npix+npsf+1:npix+2*npsf,:,:)
     padarr(npix+npsf+1:npix+2*npsf,:,:) = padarr(npsf+1:2*npsf,:,:)

     DO j=1,npix
        DO i=1,npix
           DO k=1,nfil
              convolve(i,j,k) = &
                   SUM(padarr(i+npsf-npsf/2:i+npsf+npsf/2,&
                   j+npsf-npsf/2:j+npsf+npsf/2,k)*psf)
           ENDDO
        ENDDO
     ENDDO

  ELSE

     !this is the fancy FFT version (which is much much faster)

     DO k=1,nfil

        spec_out = 0.
        speq_out = 0.
        !padarr2  = arr(:,:,k)

        CALL RLFT2(arr(:,:,k),spec_d,speq_d,1)
        CALL RLFT2(psf2,spec_p,speq_p,1)

        DO i=1,npsf2/2
           speq_out(i) = speq_d(i)*speq_p(i)/(npix**2/2)
           DO j=1,npsf2
              spec_out(i,j) = spec_d(i,j)*spec_p(i,j)/(npix**2/2)
           ENDDO
        ENDDO

        CALL RLFT2(convolve(:,:,k),spec_out,speq_out,-1)
        !convolve(:,:,k) = padarr3

     ENDDO

  ENDIF
  

END FUNCTION CONVOLVE
