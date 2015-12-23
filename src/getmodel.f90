FUNCTION GETMODEL(inpos,im)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : convolve,hist_2d,add_obs_err,myran,mypoidev
  USE nr, ONLY : locate,ran1
  IMPLICIT NONE

  REAL(SP), DIMENSION(npar), INTENT(in) :: inpos
  REAL(SP), DIMENSION(npix,npix), OPTIONAL :: im
  REAL(SP), DIMENSION(nx,ny) :: getmodel
  REAL(SP), DIMENSION(npix,npix,nfil) :: f1,cf1,of1
  REAL(SP), DIMENSION(npix,npix) :: narr
  INTEGER :: ilo,i,j,k,m,f,wgti,ii,jj,test_time=0
  REAL(SP) :: di,nnn
  CHARACTER(10) :: time
  REAL(SP), DIMENSION(niso_max) :: imf

  !------------------------------------------------------------!

  IF (test_time.EQ.1) THEN
     CALL DATE_AND_TIME(TIME=time)
     WRITE(*,*) '1 Time '//time(1:2)//':'//time(3:4)//':'//time(5:9)
  ENDIF

  !push the age weights into the IMF weights
  DO i=1,nage
     IF (10**inpos(i+1).LT.1/(npix**2)) THEN
        imf(ageind(i)+1:ageind(i+1)) = 0.0
     ELSE
        imf(ageind(i)+1:ageind(i+1)) = &
             10**inpos(i+1)*iso(ageind(i)+1:ageind(i+1))%imf
     ENDIF
  ENDDO

  !this call takes about a second with 1E7 points
  CALL RAN1(ranarr)

  !compute the model at each pixel
  DO k=1,niso
 
     IF (imf(k).LT.tiny_number) CYCLE

     nnn = 10**inpos(1)*imf(k)

     !treat masses less than minmass as continuously sampled
     IF (iso(k)%mass.LT.minmass) THEN
        DO f=1,2
           f1(:,:,f) = f1(:,:,f)+nnn*iso(k)%bands(f)
        ENDDO
     ELSE
        narr = mypoidev(nnn)
        DO f=1,2
           f1(:,:,f) = f1(:,:,f)+narr*iso(k)%bands(f)
        ENDDO
     ENDIF

  ENDDO

  IF (test_time.EQ.1) THEN
     CALL DATE_AND_TIME(TIME=time)
     WRITE(*,*) '2 Time '//time(1:2)//':'//time(3:4)//':'//time(5:9)
  ENDIF

  !convolve with PSF
  cf1 = convolve(f1)

  IF (test_time.EQ.1) THEN
     CALL DATE_AND_TIME(TIME=time)
     WRITE(*,*) '3 Time '//time(1:2)//':'//time(3:4)//':'//time(5:9)
  ENDIF

  !add obs errors
  of1 = add_obs_err(cf1)

  IF (test_time.EQ.1) THEN
     CALL DATE_AND_TIME(TIME=time)
     WRITE(*,*) '4 Time '//time(1:2)//':'//time(3:4)//':'//time(5:9)
  ENDIF

  IF (PRESENT(im)) im=of1(:,:,2)

  !compute Hess diagram, normalize to unity
  getmodel = hist_2d(of1(:,:,1)-of1(:,:,2),of1(:,:,2),&
       xhess,yhess,nx,ny,npix) / npix**2


END FUNCTION GETMODEL
