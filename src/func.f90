FUNCTION FUNC(inpos)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : getmodel
  IMPLICIT NONE

  REAL(SP), DIMENSION(:) :: inpos
  REAL(SP) :: func
  REAL(SP), DIMENSION(nx,ny) :: imodel,model_err
  INTEGER :: i,j,nn
 
  !------------------------------------------------------------!

  !get the model
  imodel = getmodel(inpos(1:npar))

  !Poisson uncertainty on the model
  model_err = SQRT(imodel*npix**2)/npix**2
  DO i=1,nx
     DO j=1,ny
        IF (model_err(i,j).LE.tiny_number) model_err(i,j)=1./npix
     ENDDO
  ENDDO

  !compute chi^2
  !func = SUM( (hess_data-imodel)**2  / (hess_err**2+model_err**2) )
  func = SUM( (hess_data-imodel)**2 / hess_err**2 )
  !func = SUM( (hess_data-imodel)**2 )
  

  !-------------------specify priors----------------------!

  !don't let the factors get too low or too high
  DO i=1,npar
     IF (inpos(i).LT.-8.0.OR.inpos(i).GT.0.3) func=huge_number
  ENDDO

  !error checking
  IF (ISNAN(func)) THEN
     WRITE(*,'("FUNC ERROR: chi^2 returned a NaN:",10F10.5)') inpos
     DO i=1,npix
        DO j=1,npix
           IF (ISNAN(imodel(i,j))) WRITE(*,*) i,j
        ENDDO
     ENDDO
     STOP
  ENDIF

  !WRITE(*,'(F15.9,100F10.5)') LOG10(func),inpos(1:10)
  

END FUNCTION FUNC
