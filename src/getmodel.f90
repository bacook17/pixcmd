FUNCTION GETMODEL(inpos,im)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : convolve,hist_2d,add_obs_err,myran,mypoidev,drawn
  USE nr, ONLY : locate,ran1
  IMPLICIT NONE

  INTEGER, PARAMETER :: test_time=0

  REAL(SP), DIMENSION(npar), INTENT(in) :: inpos
  REAL(SP), DIMENSION(npix,npix), OPTIONAL :: im
  REAL(SP), DIMENSION(nx,ny) :: getmodel
  REAL(SP), DIMENSION(npix,npix,nfil) :: f1,cf1,of1
  REAL(SP), DIMENSION(npix,npix) :: narr
  INTEGER  :: i,k,f
  REAL(SP) :: nnn
  CHARACTER(10) :: time
  REAL(SP), DIMENSION(niso_max) :: imf

  !------------------------------------------------------------!

  IF (test_time.EQ.1) THEN
     CALL DATE_AND_TIME(TIME=time)
     WRITE(*,*) '1 Time '//time(1:2)//':'//time(3:4)//':'//time(5:9)
  ENDIF

  !push the age weights into the IMF weights
  imf = 0.0
  DO i=1,nage
     IF (10**inpos(i+nxpar).LT.1./(npix**2).OR.&
          inpos(i+nxpar).LE.prlo) THEN
        imf(ageind(i)+1:ageind(i+1)) = 0.0
     ELSE
        imf(ageind(i)+1:ageind(i+1)) = &
             10**inpos(i+1)*iso(ageind(i)+1:ageind(i+1))%imf
     ENDIF
  ENDDO

  !compute the model at each pixel
  f1 = 0.0
  DO k=1,niso

     IF (imf(k).LT.tiny_number) CYCLE

     nnn = 10**inpos(1)*imf(k)

     !treat masses less than minmass as continuously sampled
     IF (iso(k)%mass.LT.minmass.OR.nnn.GT.minnum) THEN
        DO f=1,2
           f1(:,:,f) = f1(:,:,f)+nnn*iso(k)%bands(f)
        ENDDO
     ELSE
        narr = drawn(nnn)
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

  !add obs errors
  of1 = add_obs_err(cf1)

  !return the I-band image if included in the function call
  IF (PRESENT(im)) im=of1(:,:,2)

  !compute Hess diagram, normalize to unity
  getmodel = hist_2d(of1(:,:,1)-of1(:,:,2),of1(:,:,2),&
       xhess,yhess,nx,ny,npix) / npix**2

  !throw an error if the Hess diagram has no entries
  IF (SUM(getmodel).LT.tiny_number.AND.1.EQ.0) THEN
     WRITE(*,*) 'getmodel=0.0!'
     WRITE(*,*) MINVAL(of1),MAXVAL(of1)
     WRITE(*,*) MINVAL(cf1),MAXVAL(cf1)
     WRITE(*,*) MINVAL(f1),MAXVAL(f1)
     WRITE(*,*) inpos
     STOP
  ENDIF

  IF (test_time.EQ.1) THEN
     CALL DATE_AND_TIME(TIME=time)
     WRITE(*,*) '3 Time '//time(1:2)//':'//time(3:4)//':'//time(5:9)
  ENDIF


END FUNCTION GETMODEL
