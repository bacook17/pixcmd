FUNCTION GETMODEL(inpos,im)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : convolve,hist_2d,add_obs_err,&
       myran,mypoidev,drawn,interp_iso
  USE nr, ONLY : locate,ran1
  IMPLICIT NONE

  INTEGER, PARAMETER :: test_time=1

  REAL(SP), DIMENSION(npar), INTENT(in) :: inpos
  REAL(SP), DIMENSION(npix,npix), OPTIONAL :: im
  REAL(SP), DIMENSION(nx,ny) :: getmodel
  REAL(SP), DIMENSION(npix,npix,nfil) :: f1,cf1,of1
  REAL(SP), DIMENSION(npix,npix) :: narr
  REAL(SP), DIMENSION(nz) :: mdf
  INTEGER  :: i,k,f,z,zniso
  REAL(SP) :: nnn,sfh,red
  CHARACTER(10) :: time
  TYPE(TISO), DIMENSION(niso_max) :: ziso

  !------------------------------------------------------------!

  IF (test_time.EQ.1) THEN
     CALL DATE_AND_TIME(TIME=time)
     WRITE(*,*) '1 Time '//time(1:2)//':'//time(3:4)//':'//time(5:9)
  ENDIF

  !interpolate the isochrones to the input metallicity
  CALL INTERP_ISO(inpos(npar),ziso,zniso)

 
  !compute the model at each pixel
  f1 = 0.0
  DO k=1,zniso

     sfh = inpos(nxpar+ziso(k)%aind)

     IF (sfh.LE.prlo_sfh.OR.10**sfh.LT.0.1/(npix**2)) CYCLE
     
     nnn = 10**sfh*ziso(k)%imf
     
     !treat masses less than minmass as continuously sampled
     IF (ziso(k)%mass.LT.minmass.OR.nnn.GT.minnum) THEN
        DO f=1,2
           red = 10**(-2./5*(red_per_ebv(f)*10**inpos(1)))
           f1(:,:,f) = f1(:,:,f)+nnn*ziso(k)%bands(f)*red
        ENDDO
     ELSE
        IF (nnn.LE.maxpoidev) THEN
           narr = mypoidev(nnn,k)
           
           !remove stars that would likely count as "resolved"
       !    IF (immed.LE.tiny_number) WRITE(*,*) 'GETMODEL ERROR, immed=0.0!'
       !    IF (SUM(narr).LE.10.AND.ziso(k)%bands(1).GT.100.*immed) THEN
       !       narr=0.0
       !    ENDIF
           
        ELSE
           narr = gdev*SQRT(nnn)+nnn
        ENDIF
        DO f=1,2
           red = 10**(-2./5*(red_per_ebv(f)*10**inpos(1)))
           f1(:,:,f) = f1(:,:,f)+narr*ziso(k)%bands(f)*red
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
  IF (incl_obs_err.EQ.1) THEN
     of1 = add_obs_err(cf1)
  ELSE
     !no errors!
     of1 = -2.5*LOG10(cf1)
  ENDIF

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
