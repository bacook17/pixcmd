FUNCTION CONVOLVE(arr)

  ! Routine to convolve a 2D image with a PSF
  
  USE pixcmd_vars
  USE nrtype;  USE nr, ONLY : four2, rlft2
  IMPLICIT NONE

  REAL(SP), DIMENSION(npix,npix,nfil), INTENT(inout) :: arr
  REAL(SP), DIMENSION(npix+2*npsf,npix+2*npsf,nfil) :: padarr
  REAL(SP), DIMENSION(npix,npix,nfil) :: convolve
  INTEGER :: i,j,k,m,n,n2
  
  REAL(SP), DIMENSION(npsf2,npsf2) :: padarr2=0.,padarr3=0.
  COMPLEX(SPC), DIMENSION(npsf2/2,npsf2) :: spec_d,spec_out
  COMPLEX(SPC), DIMENSION(npsf2) :: speq_d,speq_out
  COMPLEX(SPC), DIMENSION(npsf2/2,npsf2) :: spec_p
  COMPLEX(SPC), DIMENSION(npsf2)      :: speq_p

  !------------------------------------------------------------!

  convolve = 0.0

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
                   j+npsf-npsf/2:j+npsf+npsf/2,k)*psfi)
           ENDDO
        ENDDO
     ENDDO

  ELSE

     !this is the fancy FFT version (which is much faster)

     n2 = npix/psf_step

     DO k=1,nfil

        !step=1 -> 0.20s
        !step=2 -> 0.14s
        !step=4 -> 0.12s

        DO m=1,psf_step
           DO n=1,psf_step

              spec_out = 0.
              speq_out = 0.
              padarr2  = arr((m-1)*n2+1:m*n2,(n-1)*n2+1:n*n2,k)

              CALL RLFT2(padarr2,spec_d,speq_d,1)
              CALL RLFT2(psf(:,:,m,n),spec_p,speq_p,1)

              DO i=1,npsf2/2
                 speq_out(i) = speq_d(i)*speq_p(i)/(n2**2/2)
                 DO j=1,npsf2
                    spec_out(i,j) = spec_d(i,j)*spec_p(i,j)/(n2**2/2)
                 ENDDO
              ENDDO

              CALL RLFT2(padarr3,spec_out,speq_out,-1)
              convolve((m-1)*n2+1:m*n2,(n-1)*n2+1:n*n2,k) = padarr3
              
           ENDDO
        ENDDO

     ENDDO

  ENDIF
  

END FUNCTION CONVOLVE
