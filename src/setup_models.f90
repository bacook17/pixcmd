SUBROUTINE SETUP_MODELS()

  USE pixcmd_vars
  IMPLICIT NONE

  CHARACTER(6)  :: zstr
  CHARACTER(4)  :: mstr
  INTEGER :: i,j,m,stat

  !------------------------------------------------------------!

  CALL GETENV('PIXCMD_HOME',PIXCMD_HOME)

 ! zmetarr = (/0.0010,0.0025,0.0040,0.0080,0.0190,0.0290/)
  zmetarr = (/0.0190/)

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

  !read in the ACS F814W PSF (log PSF in the file)
  OPEN(12,file=TRIM(PIXCMD_HOME)//'/psf/f814w.psf',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'PIXCMD ERROR: PSF file not found'
     STOP
  ENDIF
  DO i=1,npsf
     READ(12,*) psf(:,i)
  ENDDO
  psf = 10**psf

  psf2=0.
  !turn the PSF inside out...
  psf2(1:npsf/2+2,1:npsf/2+2)=psf(npsf/2:npsf,npsf/2:npsf)
  psf2(1:npsf/2+2,npsf2-npsf/2+2:npsf2)=psf(npsf/2:npsf,1:npsf/2-1)
  psf2(npsf2-npsf/2+2:npsf2,1:npsf/2+2) = psf(1:npsf/2-1,npsf/2:npsf)
  psf2(npsf2-npsf/2+2:npsf2,npsf2-npsf/2+2:npsf2) = psf(1:npsf/2-1,1:npsf/2-1)


  !--------------read in the model Hess diagrams---------------!

  DO m=1,nm

     WRITE(mstr,'(F4.2)') mpixarr(m)

     DO i=1,nz

        WRITE(zstr,'(F6.4)') zmetarr(i)
        OPEN(11,FILE=TRIM(PIXCMD_HOME)//'/models/M'//mstr//'_Z'//zstr//&
             '.im',FORM='UNFORMATTED',STATUS='OLD',access='direct',&
             recl=nage*npix*npix*nfil*4,ACTION='READ')
        READ(11,rec=1) model(m,i,:,:,:,:)
        CLOSE(11)

     ENDDO

  ENDDO


END SUBROUTINE SETUP_MODELS
