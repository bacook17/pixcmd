PROGRAM WRITE_A_MODEL

  ! Routine to create one model Hess diagram convolved
  ! with an MDF, a SFH, and a P(Mpix)

  USE nrtype; USE pixcmd_utils; USE pixcmd_vars
  USE nr, ONLY : locate
  IMPLICIT NONE
  
  INTEGER  :: i,j,k,wgti
  REAL(SP) :: mpix,zmet0,zmets,tau,tot=1.,dt,twgt=0.,wgtmax=0.
  CHARACTER(50) :: outfile
  REAL(SP), DIMENSION(npar)      :: pos
  REAL(SP), DIMENSION(nx,ny)     :: hess
  REAL(SP), DIMENSION(npix,npix) :: im
  REAL(SP), DIMENSION(nage)      :: sfh
  REAL(SP), DIMENSION(nz,nage)   :: wgt
  REAL(SP), DIMENSION(nz)        :: mdf,zz
  CHARACTER(10) :: time

  !------------------------------------------------------------!

  !output file name
  outfile = 'model_M2.0_cSFH'

  !initialize the random number generator
  CALL INIT_RANDOM_SEED()

  !setup the model grid
  CALL SETUP_MODELS()


  !-----------------mass per pixel------------------!
  mpix = 2.0
 
  !-----------------------MDF-----------------------!
  !zz  = LOG10(zmetarr/0.0190)
  !mdf = EXP(-(zz-zmet0)**2/2/zmets)
  !mdf = 1.0
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
  hess = getmodel(pos,im)

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
