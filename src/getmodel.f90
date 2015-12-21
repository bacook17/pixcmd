FUNCTION GETMODEL(inpos)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : convolve,hist_2d,add_obs_err,myran
  USE nr, ONLY : locate
  IMPLICIT NONE

  REAL(SP), DIMENSION(npar) :: inpos
  REAL(SP), DIMENSION(nx,ny) :: getmodel
  REAL(SP), DIMENSION(npix,npix,nfil) :: f1,cf1,of1
  INTEGER :: ilo,i,j,k,m,wgti,ii,jj,test_time=0
  REAL(SP) :: di
  CHARACTER(10) :: time

  !------------------------------------------------------------!

  IF (test_time.EQ.1) THEN
     CALL DATE_AND_TIME(TIME=time)
     WRITE(*,*) '1 Time '//time(1:2)//':'//time(3:4)//':'//time(5:9)
  ENDIF

  !linearly combine the model components
  f1 = model(:,:,:,nz,nage)
  k=1
  DO j=1,nage-1
     DO i=1,nz
        wgti = INT(10**inpos(k)*npix**2)
        DO m=1,wgti
           ii = INT(myran()*npix)
           jj = INT(myran()*npix)
           f1(ii,jj,:) = model(ii,jj,:,i,j)
        ENDDO
        k=k+1
     ENDDO
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

  !compute Hess diagram, normalize to unity
  getmodel = hist_2d(of1(:,:,1)-of1(:,:,2),of1(:,:,2),&
       xhess,yhess,nx,ny,npix) / npix**2


END FUNCTION GETMODEL
