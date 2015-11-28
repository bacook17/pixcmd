FUNCTION FUNC(inpos)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : get_model
  IMPLICIT NONE

  REAL(SP), DIMENSION(:) :: inpos
  REAL(SP) :: func
  REAL(SP), DIMENSION(nx,ny) :: imodel,model_err
  INTEGER :: i

  !------------------------------------------------------------!

  imodel = get_model(inpos)

  !Poisson uncertainty on the model
  model_err = SQRT(imodel*npix**2)/npix**2

  !compute chi^2
  func = SUM( (hess_data-imodel)**2 /(hess_err**2+model_err**2) )
  

  !-------------------specify priors----------------------!

  !unphysical
  IF (SUM(10**inpos(1:5)).GT.1.0) func=huge_number
  
  !dont let the factors get too low
  DO i=1,npar-1
     IF (inpos(i).LT.-9.0) func=huge_number
  ENDDO

  !set limits on the age
  IF (inpos(6).LT.7.0.OR.inpos(6).GT.10.3) func=huge_number

  !error checking
  IF (ISNAN(func)) THEN
     WRITE(*,'("FUNC ERROR: chi^2 returned a NaN:",10F10.5)') inpos
     STOP
  ENDIF

  !WRITE(*,'(10F10.5)') LOG10(func),inpos
  

END FUNCTION FUNC
