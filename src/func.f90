FUNCTION FUNC(inpos)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : getmodel
  IMPLICIT NONE

  REAL(SP), DIMENSION(:) :: inpos
  REAL(SP) :: func
  REAL(SP), DIMENSION(nx,ny) :: imodel,model_err
  INTEGER :: i,j,nn
 
  !------------------------------------------------------------!

  func = 0.0
  nn   = SIZE(inpos)

  !------------------- priors----------------------!
  
  !prior on the mass per pixel
  IF (inpos(1).LT.prlo_m.OR.inpos(1).GT.prhi_m) func=huge_number

  !priors on the weights for each component
  DO i=2,npar
     IF (inpos(i).LT.prlo.OR.inpos(i).GT.prhi) func=huge_number
     IF (ISNAN(inpos(i))) THEN
        WRITE(*,'("FUNC ERROR: inpos returned a NaN:",10F18.2)') inpos
        STOP
     ENDIF
  ENDDO

  !actually compute a model and return chi^2 if within priors
  IF (func.LT.huge_number) THEN

     !get the model
     imodel = getmodel(inpos(1:nn))

     !Poisson uncertainty on the model
     model_err = SQRT(imodel*npix**2)/npix**2
     DO i=1,nx
        DO j=1,ny
           IF (model_err(i,j).LE.tiny_number) model_err(i,j)=1./npix**2
        ENDDO
     ENDDO

     !compute chi^2
     !func = SUM( (hess_data-imodel)**2  / (hess_err**2+model_err**2) )
     func = SUM( (hess_data-imodel)**2 / hess_err**2 )
     !func = SUM( (hess_data-imodel)**2 )
  
  ENDIF


  !check for NaNs
  IF (ISNAN(func)) THEN
     WRITE(*,'("FUNC ERROR: chi^2 returned a NaN:",10F10.5)') inpos
     DO i=1,npix
        DO j=1,npix
           IF (ISNAN(imodel(i,j))) WRITE(*,*) i,j
        ENDDO
     ENDDO
     STOP
  ENDIF

  !WRITE(*,'(F10.4,30(F5.2,1x))') LOG10(func),inpos
  

END FUNCTION FUNC
