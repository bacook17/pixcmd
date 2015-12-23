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

  !------------------- priors----------------------!

  IF (inpos(1).LT.prlo_m.OR.inpos(1).GT.prhi_m) func=huge_number

  DO i=2,npar
     IF (inpos(i).LT.prlo.OR.inpos(i).GT.prhi) func=huge_number
     IF (ISNAN(inpos(i))) THEN
        WRITE(*,'("FUNC ERROR: inpos returned a NaN:",10F18.2)') inpos
        STOP
     ENDIF
  ENDDO

  IF (func.LT.huge_number) THEN

     !get the model
     imodel = getmodel(inpos(1:npar))

     !Poisson uncertainty on the model
     !model_err = SQRT(imodel*npix**2)/npix**2
     !DO i=1,nx
     !   DO j=1,ny
     !      IF (model_err(i,j).LE.tiny_number) model_err(i,j)=1./npix
     !   ENDDO
     !ENDDO

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

  !WRITE(*,'(F10.4,100F8.3)') LOG10(func),inpos
  

END FUNCTION FUNC
