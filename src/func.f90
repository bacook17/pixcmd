FUNCTION FUNC(inpos)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : getmodel
  IMPLICIT NONE

  REAL(SP), DIMENSION(:) :: inpos
  REAL(SP) :: func
  REAL(SP), DIMENSION(nx,ny) :: imodel,model_err
  INTEGER :: i
 
  !------------------------------------------------------------!

  !get the model
  imodel = getmodel(inpos)
 
  !Poisson uncertainty on the model
  model_err = SQRT(imodel*npix**2)/npix**2

  !compute chi^2
  func = SUM( (hess_data-imodel)**2 /(hess_err**2+model_err**2) )
  

  !-------------------specify priors----------------------!

  !unphysical
  !IF (SUM(10**inpos(1:5)).GT.1.0) func=huge_number
  
  !dont let the factors get too low or too high
  DO i=1,npar
     IF (inpos(i).LT.-9.0.OR.inpos(i).GT.0.0) func=huge_number
  ENDDO

  !set limits on Mpix
  !IF (inpos(1).LT.mpix0.OR.inpos(1).GT.(mpix0+nm*dmpix)) func=huge_number

  !error checking
  IF (ISNAN(func)) THEN
     WRITE(*,'("FUNC ERROR: chi^2 returned a NaN:",10F10.5)') inpos
     STOP
  ENDIF

  !WRITE(*,'(100F10.5)') LOG10(func),inpos(1:10)
  

END FUNCTION FUNC
