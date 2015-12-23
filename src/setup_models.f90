SUBROUTINE SETUP_MODELS()

  USE pixcmd_vars; USE nrtype
  USE nr, ONLY : ran1
  IMPLICIT NONE

  CHARACTER(5), DIMENSION(nz)  :: zstr
  CHARACTER(10) :: tag=''
  CHARACTER(4)  :: mstr
  CHARACTER(1)  :: is,js
  INTEGER :: i,j,k,m,t,stat
  REAL(SP) :: iage=0.0,age,d2,d3,d4
  REAL(SP), DIMENSION(nage,npix,npix,nfil) :: tmodel

  !------------------------------------------------------------!

  tag = '_x5FEWER'

  CALL GETENV('PIXCMD_HOME',PIXCMD_HOME)

  OPEN(13,FILE=TRIM(PIXCMD_HOME)//'/isoc/zlegend.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SETUP_MODELS ERROR: zlegend.dat file not found'
     STOP
  ENDIF
  DO i=1,nzskip
     READ(13,*)
  ENDDO
  DO i=1,nz
     READ(13,'(A5,3x,F6.4)') zstr(i),zmetarr(i)
  ENDDO
  CLOSE(13)

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

  !-------------------------------------------------------------------!

  psf=0.

  !read in all of the sub-pixel shifted PSFs
  DO i=1,psf_step
     DO j=1,psf_step

        WRITE(is,'(I1)') i-1
        WRITE(js,'(I1)') j-1

        !read in the ACS F814W PSF (log PSF in the file)
        OPEN(12,file=TRIM(PIXCMD_HOME)//'/psf/f814w_'//is//js//'.psf',&
             STATUS='OLD',iostat=stat,ACTION='READ')
        IF (stat.NE.0) THEN
           WRITE(*,*) 'SETUP_MODELS ERROR: PSF file not found'
           STOP
        ENDIF
        DO k=1,npsf
           READ(12,*) psfi(:,k)
        ENDDO
        CLOSE(12)
        psfi = 10**psfi

        !turn the PSF inside out...
        psf(1:npsf/2+2,1:npsf/2+2,i,j) = psfi(npsf/2:npsf,npsf/2:npsf)
        psf(1:npsf/2+2,npsf2-npsf/2+2:npsf2,i,j) = psfi(npsf/2:npsf,1:npsf/2-1)
        psf(npsf2-npsf/2+2:npsf2,1:npsf/2+2,i,j) = psfi(1:npsf/2-1,npsf/2:npsf)
        psf(npsf2-npsf/2+2:npsf2,npsf2-npsf/2+2:npsf2,i,j) = &
             psfi(1:npsf/2-1,1:npsf/2-1)
        
     ENDDO
  ENDDO

  !-------------------------------------------------------------------!

  !open the isochrone file
  OPEN(10,file=TRIM(PIXCMD_HOME)//'/isoc/MIST_v29_Z'//&
       TRIM(zstr(1))//TRIM(tag)//'.dat',STATUS='OLD',iostat=stat,&
       ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SIM_PIXCMD ERROR: isoc file not found'
     STOP
  ENDIF

  i=1
  ageind(1)=0
  DO t=1,nage

     age = age0+(t-1)*dage
     
     !burn the young models
     DO WHILE (iage-age.LT.-0.01) 
        READ(10,*,IOSTAT=stat) iage
     ENDDO
     BACKSPACE(10)
     i=i-1
     DO WHILE (ABS(iage-age).LT.1E-2) 
        IF (i.GT.niso_max) THEN
           WRITE(*,*) 'SIM_PIXCMD ERROR niso_max reached',i
           STOP
        ENDIF
        i=i+1
        READ(10,*,IOSTAT=stat) iage,iso(i)%mass,d2,d3,iso(i)%imf,&
             iso(i)%bands(1),d4,iso(i)%bands(2)
        iso(i)%age = iage
        IF (stat.NE.0) GOTO 20
     ENDDO
20   CONTINUE

     ageind(t+1) = i-1

  ENDDO

  niso = i-1

  iso(1:niso)%imf      = 10**iso(1:niso)%imf
  !convert to flux
  iso(1:niso)%bands(1) = 10**(-2./5*iso(1:niso)%bands(1))
  iso(1:niso)%bands(2) = 10**(-2./5*iso(1:niso)%bands(2))


END SUBROUTINE SETUP_MODELS
