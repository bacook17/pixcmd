PROGRAM SIM_PIXCMD

  ! Program to compute a pixel CMD and from it a Hess diagram
  ! syntax: sim_pixcmd.exe log(Mpix) N_Mpix dMpix Z

  ! Note that the model Hess diagrams are likely not completely
  ! converged at Npix=1E3.  Should at some point try 1E4
  
  USE pixcmd_vars; USE pixcmd_utils
  USE nr, ONLY : poidev; USE nrtype 
  IMPLICIT NONE


  
  INTEGER :: i,j,k,m,t,f,iseed,nn,niso,stat,nmpix,ct
  REAL(SP) :: mpix,age,iage=0.,mpx0,dmpx,d1,d2,d3,d4,nnn
  REAL(SP), DIMENSION(nage,npix,npix,nfil)  :: flux=0.
  REAL(SP), DIMENSION(npix,npix,nfil)  :: cflux=0.,oflux=0.
  REAL(SP), DIMENSION(npix,npix) :: narr
  TYPE(TISO), DIMENSION(niso_max) :: iso
  CHARACTER(10) :: time
  CHARACTER(50) :: file='',tag=''
  CHARACTER(15)  :: zstr
  CHARACTER(4)  :: mstr


  !----------------------------------------------------------------!

  IF (IARGC().LT.1) THEN
     mpx0  = 2.0
     nmpix = 1
     dmpx  = 0.1
     zstr  = 'p0.00_x5FEWER'
  ELSE IF (IARGC().NE.4.AND.IARGC().GT.0) THEN
     WRITE(*,*) 'syntax: sim_pixcmd.exe Mpix0 N_Mpix dMpix zmet'
     WRITE(*,*) '  e.g., sim_pixcmd.exe  2.0    10    0.1  p0.00'
     STOP
  ELSE
     CALL GETARG(1,file)
     READ(file,*) mpx0
     CALL GETARG(2,file)
     READ(file,*) nmpix
     CALL GETARG(3,file)
     READ(file,*) dmpx
     CALL GETARG(4,file)
     READ(file,'(A5)') zstr
  ENDIF

  WRITE(*,*) 
  IF (zstr(1:1).EQ.'p') THEN
     WRITE(*,'("[Z/H]=+",A4)') TRIM(zstr(2:))
  ELSE
     WRITE(*,'("[Z/H]=-",A4)') TRIM(zstr(2:))
  ENDIF

  !set a background flux level
  flux = bkgnd

  !initialize the random number generator
  CALL INIT_RANDOM_SEED()

  !setup the model arrays, PSF, etc
  CALL SETUP_MODELS(0)

  CALL DATE_AND_TIME(TIME=time)
  WRITE(*,*) 'Start Time '//time(1:2)//':'//time(3:4)//':'//time(5:6)

  !------------------------Loop on Mpix---------------------------!

  DO m=1,nmpix

     hess = 0.0
     iage = 0.0

     mpix = 10**(mpx0+(m-1)*dmpx)
     WRITE(mstr,'(F4.2)') LOG10(mpix)
     WRITE(*,'("   log Mpix  =",F6.2)') LOG10(mpix)

     !READ(10,*)

     !--------------------Loop on population age--------------------!
     
     DO t=1,nage

        age = age0+(t-1)*dage

        !----------------Read in the input isochrone info----------------!
        
        !burn the young models
        DO WHILE (iage-age.LT.-0.01) 
           READ(10,*,IOSTAT=stat) iage
        ENDDO

        BACKSPACE(10)
        i=0
        iso%imf = -99.
        DO WHILE (ABS(iage-age).LT.1E-2) 
           IF (i.GT.niso_max) THEN
              WRITE(*,*) 'SIM_PIXCMD ERROR niso_max reached',i
              STOP
           ENDIF
           i=i+1
           READ(10,*,IOSTAT=stat) iage,iso(i)%mass,d2,d3,iso(i)%imf,&
                iso(i)%bands(1),d4,iso(i)%bands(2)
           IF (stat.NE.0) GOTO 20
        ENDDO
20      CONTINUE
        
        niso = i-1

   !     WRITE(*,'(F6.2,1x,I4,3F7.2)') age,niso!,iso(niso)%bands,iso(niso)%imf

        iso(1:niso)%imf   = 10**iso(1:niso)%imf
        !convert to flux
        iso(1:niso)%bands(1) = 10**(-2./5*iso(1:niso)%bands(1))
        iso(1:niso)%bands(2) = 10**(-2./5*iso(1:niso)%bands(2))
        
        !----------------Compute the model at each pixel-----------------!
        
  !      CALL DATE_AND_TIME(TIME=time)
  !      WRITE(*,*) 'Start Time '//time(1:2)//':'//time(3:4)//':'//time(5:6)

        !compute the model over an NxN image
        DO k=1,niso

           nnn = mpix*iso(k)%imf

           IF (nnn.GT.1.) THEN
              DO f=1,2
                 flux(t,:,:,f) = flux(t,:,:,f)+nnn*iso(k)%bands(f)
              ENDDO
           ELSE
              narr = mypoidev(nnn)
              DO f=1,2
                 flux(t,:,:,f) = flux(t,:,:,f)+narr*iso(k)%bands(f)
              ENDDO
           ENDIF

        ENDDO

   !     CALL DATE_AND_TIME(TIME=time)
   !     WRITE(*,*) '  End Time '//time(1:2)//':'//time(3:4)//':'//time(5:6)
        
        !-------------------Convolve with the ACS PSF--------------------!
        
        cflux = convolve(flux(t,:,:,:))
       
        !--------------------Add observational errors--------------------!
        
        oflux = add_obs_err(cflux)

        !---------------Compute a Hess diagram in B-I vs. I--------------!

        hess(t,:,:) = hist_2d(oflux(:,:,1)-oflux(:,:,2),oflux(:,:,2),&
             xhess,yhess,nx,ny,npix)

     ENDDO
     
     CLOSE(10)

     !save the PSF-convolved, obs err-included Hess diagram
     OPEN(11,FILE=TRIM(PIXCMD_HOME)//'/models/M'//mstr//'_Z'//TRIM(zstr)//&
          '.hess',FORM='UNFORMATTED',STATUS='REPLACE',access='direct',&
          recl=nage*nx*ny*4)
     WRITE(11,rec=1) hess
     CLOSE(11)

     !save the noise-free model image
     OPEN(11,FILE=TRIM(PIXCMD_HOME)//'/models/M'//mstr//'_Z'//TRIM(zstr)//&
          '.im',FORM='UNFORMATTED',STATUS='REPLACE',access='direct',&
          recl=nage*npix*npix*nfil*4)
     WRITE(11,rec=1) flux
     CLOSE(11)

  ENDDO

  CALL DATE_AND_TIME(TIME=time)
  WRITE(*,*) 'End Time '//time(1:2)//':'//time(3:4)//':'//time(5:6)


END PROGRAM SIM_PIXCMD
