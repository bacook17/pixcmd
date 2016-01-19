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
  REAL(SP), DIMENSION(nz) :: mdf
  INTEGER  :: i,k,f,z
  REAL(SP) :: nnn,sfh,red
  CHARACTER(10) :: time

  !------------------------------------------------------------!

  IF (test_time.EQ.1) THEN
     CALL DATE_AND_TIME(TIME=time)
     WRITE(*,*) '1 Time '//time(1:2)//':'//time(3:4)//':'//time(5:9)
  ENDIF

  mdf(1:nz-1) = inpos(1+nage+nxpar:npar)
  mdf(nz)     = 1.0-SUM(10**(mdf(1:nz-1)))
  IF (mdf(nz).LT.0.0) THEN
     WRITE(*,*) 'ERROR: MDF(nz)<0.0!',mdf(nz)
     STOP
  ELSE
     mdf(nz) = LOG10(mdf(nz)+tiny_number)
  ENDIF


  !compute the model at each pixel
  f1 = 0.0
  DO z=1,nz
     DO k=1,niso(z)
     
        sfh = inpos(nxpar+iso(z,k)%aind)

        IF (sfh.LE.prlo_sfh.OR.10**sfh.LT.0.1/(npix**2).OR.&
             mdf(z).LT.prlo_zmet) CYCLE
     
        nnn = 10**sfh*10**mdf(z)*iso(z,k)%imf
     
        !treat masses less than minmass as continuously sampled
        IF (iso(z,k)%mass.LT.minmass.OR.nnn.GT.minnum) THEN
           DO f=1,2
              red = 10**(-2./5*(red_per_ebv(f)*10**inpos(1)))
              f1(:,:,f) = f1(:,:,f)+nnn*iso(z,k)%bands(f)*red
           ENDDO
        ELSE
           IF (nnn.LE.maxpoidev) THEN
              narr = mypoidev(nnn,k)
           ELSE
              narr = gdev*SQRT(nnn)+nnn
           ENDIF
           DO f=1,2
              red = 10**(-2./5*(red_per_ebv(f)*10**inpos(1)))
              f1(:,:,f) = f1(:,:,f)+narr*iso(z,k)%bands(f)*red
           ENDDO
        ENDIF

     ENDDO
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
