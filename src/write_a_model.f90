PROGRAM WRITE_A_MODEL

  ! Routine to create one model Hess diagram convolved
  ! with an MDF, a SFH, and a P(Mpix)

  USE nrtype; USE pixcmd_utils; USE pixcmd_vars
  USE nr, ONLY : locate,ran1
  IMPLICIT NONE
  
  INTEGER  :: i,j
  REAL(SP) :: mpix=2.0,zmet0,zmets,tau,dt,twgt=0.
  CHARACTER(10) :: input
  CHARACTER(50) :: outfile
  REAL(SP), DIMENSION(npar)      :: pos
  REAL(SP), DIMENSION(nx,ny)     :: hess
  REAL(SP), DIMENSION(npix,npix) :: im
  REAL(SP), DIMENSION(nage)      :: sfh
  REAL(SP), DIMENSION(nz,nage)   :: wgt
  REAL(SP), DIMENSION(nz)        :: mdf,zz

  !------------------------------------------------------------!

  IF (IARGC().EQ.1) THEN
     CALL GETARG(1,input)
     READ(input,'(F4.1)') mpix
  ENDIF

  !output file name
  outfile = 'model_M'//TRIM(input(1:3))//'_cSFH'

  WRITE(*,'("  Output File = ",A50)') outfile
  WRITE(*,'("  log(Mpix)   =",F4.1)') mpix

  !initialize the random number generator
  CALL INIT_RANDOM_SEED()
  CALL RAN1(ranarr)

  !setup the model grid
  CALL SETUP_MODELS()
 
  !-----------------------MDF-----------------------!
  !zz  = LOG10(zmetarr/0.0190)
  !mdf = EXP(-(zz-zmet0)**2/2/zmets)
  mdf = 1.0
  
  !-----------------------SFH-----------------------!
  !tau model
  !sfh = EXP(-(10**10.201-10**agesarr)/1E9/tau)
  !constant SFH
  sfh = 1/10**agesarr(nage)

 
  DO j=1,nage
     DO i=1,nz
        IF (j.EQ.1) THEN
           dt = (10**agesarr(j)-10**(agesarr(j)-dage))
        ELSE
           dt = (10**agesarr(j)-10**agesarr(j-1))
        ENDIF
        wgt(i,j) = mdf(i)*sfh(j)*dt
        twgt = twgt+wgt(i,j)
     ENDDO
  ENDDO

  !transfer the parameters to the parameter array
  pos(1)      = mpix
  pos(2:npar) = LOG10(wgt(1,:))

  !get the model hess diagram and I-band image
  !store as the actual counts, not normalized
  hess = getmodel(pos,im) * npix**2

  !save the PSF-convolved, obs err-included I-band image
  OPEN(11,FILE=TRIM(PIXCMD_HOME)//'/data/'//TRIM(outfile)//'.im',&
       FORM='UNFORMATTED',STATUS='REPLACE',access='direct',&
          recl=npix*npix*4)
  WRITE(11,rec=1) im
  CLOSE(11)

  !save the Hess diagram to file
  OPEN(1,FILE=TRIM(PIXCMD_HOME)//'data/'//TRIM(outfile)//'.hess',&
       FORM='UNFORMATTED',STATUS='REPLACE',access='direct',&
       recl=nx*ny*4)
  WRITE(1,rec=1) hess
  CLOSE(1)


END PROGRAM WRITE_A_MODEL
