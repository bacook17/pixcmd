SUBROUTINE SETUP_MODELS()

  USE pixcmd_vars; USE nrtype
  USE nr, ONLY : ran1
  IMPLICIT NONE

  CHARACTER(5), DIMENSION(nz)  :: zstr
  CHARACTER(1)  :: is,js
  INTEGER  :: i,j,k=1,m,n,t,f,stat
  REAL(SP) :: iage=0.0,age,d2,d3,d4,mass,imf,b1,b2,tage

  !------------------------------------------------------------!

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

  !set up the Hess arrays
  DO i=1,nx
     xhess(i) = xmin+(i-1)*dx
  ENDDO
  DO i=1,ny
     yhess(i) = ymin+(i-1)*dy
  ENDDO


  !----------------------------Set up the PSF--------------------------------!

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

  !------------------------Set up the isochrones-----------------------------!

  !open the isochrone file
  OPEN(10,file=TRIM(PIXCMD_HOME)//'/isoc/MIST_v29_Z'//&
       TRIM(zstr(1))//TRIM(iso_tag)//'.dat',STATUS='OLD',iostat=stat,&
       ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SIM_PIXCMD ERROR: isoc file not found'
     STOP
  ENDIF

  !--> This needs to be cleaned up at some point...

  ! nage=7
  agesarr2 = (/6.0,7.0,8.0,8.5,9.0,9.5,10.0,10.2/)
  agesarr  = (/6.5,7.5,8.25,8.75,9.25,9.75,10.1/)

  ! nage=22
  !agesarr2 = (/6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,&
  !     8.6,8.8,9.0,9.2,9.4,9.6,9.8,10.0,10.2,10.4/)
  !agesarr  = (/6.1,6.3,6.5,6.7,6.9,7.1,7.2,7.3,7.5,7.7,7.9,8.1,8.3,&
  !     8.5,8.7,8.9,9.1,9.3,9.5,9.7,9.9,10.1/)

  i       = 1
  nisoage = 0
  tage    = 0.

  DO j=1,1000000

     READ(10,'(F6.2,7F10.4)',IOSTAT=stat) iage,mass,d2,d3,imf,b1,d4,b2
     IF (stat.NE.0) GOTO 21

     DO f=1,22

        age = (f-1)*0.2+age0
        !n=0

        IF ( ABS(iage-age).LT.0.01.AND.iage.GE.age0 ) THEN
           iso(i)%age  = iage
           iso(i)%mass = mass
           iso(i)%imf  = imf
           iso(i)%bands(1) = b1
           iso(i)%bands(2) = b2
           DO m=1,nage
              IF (iage.GE.agesarr2(m).AND.iage.LT.agesarr2(m+1)) k=m
              IF (iage.EQ.agesarr2(nage+1)) k=m
           ENDDO
           !WRITE(9,'(F5.2,1x,I2)') iage,k
           iso(i)%aind = k
           i=i+1
           IF (iage.NE.tage) THEN
              nisoage(k) = nisoage(k)+1
              tage=iage
           ENDIF
        ENDIF

     ENDDO

  ENDDO

  WRITE(*,*) 'SETUP ERROR: did not finish reading in the isochrone file'
  STOP

21 CONTINUE

  niso = i-1

  iso(1:niso)%imf      = 10**iso(1:niso)%imf
  !convert to flux
  iso(1:niso)%bands(1) = 10**(-2./5*iso(1:niso)%bands(1))
  iso(1:niso)%bands(2) = 10**(-2./5*iso(1:niso)%bands(2))

  !need to renormalize the weights since we now lump things togehter
  DO i=1,niso
     iso(i)%imf = iso(i)%imf / nisoage(iso(i)%aind)
  ENDDO



END SUBROUTINE SETUP_MODELS
