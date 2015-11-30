SUBROUTINE SETUP_MODELS()

  USE pixcmd_vars
  IMPLICIT NONE

  CHARACTER(6)  :: zstr
  CHARACTER(4)  :: mstr
  INTEGER :: i,j,m

  !------------------------------------------------------------!

  CALL GETENV('PIXCMD_HOME',PIXCMD_HOME)

  zmetarr = (/0.0010,0.0025,0.0040,0.0080,0.0190,0.0290/)

  !set up model ages array
  DO i=1,nage
     agesarr(i) = age0+(i-1)*dage
  ENDDO

  !set up Mpix array
  DO i=1,nm
     mpixarr(i) = mpix0+(i-1)*dmpix
  ENDDO

  !set up the Hess arrays
  DO i=1,nx
     xhess(i) = xmin+(i-1)*dx
  ENDDO
  DO i=1,ny
     yhess(i) = ymin+(i-1)*dy
  ENDDO

  !--------------read in the model Hess diagrams---------------!

  DO m=1,nm

     WRITE(mstr,'(F4.2)') mpixarr(m)

     DO i=1,nz

        WRITE(zstr,'(F6.4)') zmetarr(i)
        OPEN(11,FILE=TRIM(PIXCMD_HOME)//'/models/hess_M'//mstr//'_Z'//zstr//&
             '.dat',FORM='UNFORMATTED',STATUS='OLD',access='direct',&
             recl=nage*nx*ny*4,ACTION='READ')
        READ(11,rec=1) model(m,i,:,:,:)
        CLOSE(11)
        !normalize each model to unity
        DO j=1,nage
           model(m,i,j,:,:) = model(m,i,j,:,:)/npix**2
        ENDDO

     ENDDO

  ENDDO


END SUBROUTINE SETUP_MODELS
