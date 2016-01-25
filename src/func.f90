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
  
  !priors for each parameter
  DO i=1,npar
     IF (inpos(i).LT.prlo(i).OR.inpos(i).GT.prhi(i)) func=huge_number
     IF (ISNAN(inpos(i))) THEN
        WRITE(*,'("FUNC ERROR: inpos returned a NaN:",10F18.2)') inpos
        STOP
     ENDIF
  ENDDO

  !actually compute a model and return chi^2 if within priors
  IF (func.LT.huge_number) THEN

     !get the model
     imodel = getmodel(inpos(1:nn))

     IF (SUM(imodel).LT.tiny_number) THEN

        !this means that the model is out of bounds of the 
        !Hess diagram.
        func = huge_number / 10.

     ELSE

        !Poisson uncertainty on the model
        model_err = SQRT(imodel*npix**2)/npix**2
        DO i=1,nx
           DO j=1,ny
              IF (imodel(i,j)*npix**2.EQ.1) THEN
                 model_err(i,j)=model_err(i,j)*100
              ENDIF
              IF (model_err(i,j).LE.tiny_number) model_err(i,j)=1./npix**2
           ENDDO
        ENDDO

        !compute chi^2
        func = SUM( (hess_data-imodel)**2  / (hess_err**2+model_err**2) )
  
     ENDIF
        
  ENDIF


  !check for NaNs
  IF (ISNAN(func)) THEN
     WRITE(*,'("FUNC ERROR: chi^2 returned a NaN:",90F10.3)') inpos
     STOP
  ENDIF

  !WRITE(*,'(F10.6,30(F6.3,1x))') LOG10(func),inpos
  

END FUNCTION FUNC
