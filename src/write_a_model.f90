PROGRAM WRITE_A_MODEL

  ! Routine to create one model Hess diagram convolved
  ! with an MDF, a SFH, and a P(Mpix)

  !syntax: write_a_model.exe Mpix SFH_type tau/log_age log(EBV) tag
  !e.g., write_a_model.exe 2.0 0 10.0 -> 10 Gyr SSP
  !e.g., write_a_model.exe 2.0 1 2.0  -> tau=2 Gyr SFH
  !   tau>10 Gyr -> constant SFH
  !   defaults: Mpix=2.0, constant SFH

  USE nrtype; USE pixcmd_utils; USE pixcmd_vars
  USE nr, ONLY : locate,ran1,gasdev
  IMPLICIT NONE
  
  INTEGER  :: i,j,sfhflag=1,ii
  REAL(SP) :: mpix=2.0,lebv=-5.0,zmet,tau=10.,dt,twgt=0.
  CHARACTER(10) :: char_mpix='2.0',char_flag='1',char_tau='10.0'
  CHARACTER(10) :: char_lebv='-5.0',char_zmet='0.0'
  CHARACTER(60) :: outfile='',tag=''
  REAL(SP), DIMENSION(npar)      :: pos
  REAL(SP), DIMENSION(nx,ny)     :: hess
  REAL(SP), DIMENSION(npix,npix) :: im
  REAL(SP), DIMENSION(nage)      :: sfh
  REAL(SP), DIMENSION(nage)      :: swgt

  !------------------------------------------------------------!

  !generate truly random models if fix_seed=0
  fix_seed = 0

  IF (IARGC().GT.0) THEN
     IF (IARGC().LT.4) THEN
        WRITE(*,*) 'ERROR: syntax: write_a_model.exe Mpix '//&
             'SFH_type tau/log_age log(EBV) [Z/H] [tag]'
        STOP
     ENDIF
     CALL GETARG(1,char_mpix)
     READ(char_mpix,'(F4.1)') mpix
     CALL GETARG(2,char_flag)
     READ(char_flag,'(I1)') sfhflag
     CALL GETARG(3,char_tau)
     READ(char_tau,'(F4.1)') tau
     CALL GETARG(4,char_lebv)
     READ(char_lebv,'(F4.1)') lebv
     CALL GETARG(5,char_zmet)
     READ(char_zmet,'(F4.1)') zmet
  ENDIF

  IF (IARGC().EQ.6) THEN
     tag(1:1)='_'
     CALL GETARG(6,tag(2:))
  ENDIF

  IF (tau.LT.6.0.AND.sfhflag.EQ.0) THEN
     WRITE(*,*) 'ERROR: age<6.0, dont forget to input *log* age for SFH=0!'
     STOP
  ENDIF

  IF (tau.LT.10.) THEN
     char_tau(2:4)=char_tau(1:3)
     char_tau(1:1)='0'
  ENDIF

  IF (lebv.GE.0.0) THEN
     char_lebv(2:4)=char_lebv(1:3)
     char_lebv(1:1)='+'
  ENDIF

  IF (zmet.GE.0.0) THEN
     char_zmet(2:4)=char_zmet(1:3)
     char_zmet(1:1)='+'
  ENDIF

  !output file name
  IF (sfhflag.EQ.1) THEN
     outfile = 'model_M'//TRIM(char_mpix(1:3))//'_SFH'//&
          TRIM(char_flag(1:1))//'_tau'//TRIM(char_tau(1:4))//&
          '_lebv'//TRIM(char_lebv(1:4))//'_Z'//&
          TRIM(char_zmet(1:4))//TRIM(tag)
  ELSE
     outfile = 'model_M'//TRIM(char_mpix(1:3))//'_SFH'//&
          TRIM(char_flag(1:1))//'_t'//TRIM(char_tau(1:4))//&
          '_lebv'//TRIM(char_lebv(1:4))//'_Z'//&
          TRIM(char_zmet(1:4))//TRIM(tag)
  ENDIF

  WRITE(*,'("  Output File =  ",A60)') outfile
  WRITE(*,'("  log(Mpix)   = ",F4.1)') mpix
  WRITE(*,'("  log(EBV)    = ",F4.1)') lebv
  WRITE(*,'("  [Z/H]       = ",F4.1)') zmet
  WRITE(*,'("  SFH flag    =  ",I1)') sfhflag
  IF (sfhflag.EQ.1) THEN
     WRITE(*,'("  tau         = ",F4.1)') tau
  ELSE
     WRITE(*,'("  log(age)    = ",F4.1)') tau
  ENDIF

  !initialize the random number generator
  CALL INIT_RANDOM_SEED()
  DO i=1,niso_max
     CALL RAN1(ranarr(:,i))
  ENDDO
  DO i=1,npix
     CALL GASDEV(gdev(:,i))
  ENDDO
  DO i=1,npix
     DO j=1,npix
        CALL RAN1(ebv_ran(i,j,:))
     ENDDO
  ENDDO
  ebv_ran = 10**ebv_ran

  !setup the model grid
  CALL SETUP_MODELS()
 
  !-----------------------SFH-----------------------!

  IF (sfhflag.EQ.1) THEN

     IF (tau.LT.10.) THEN
        !tau model
        sfh = EXP(-(10**10.201-10**agesarr)/1E9/tau) / tau/1E9
     ELSE
        !constant SFH
        sfh = 1/10**agesarr2(nage+1)
     ENDIF

     DO j=1,nage
        dt = (10**agesarr2(j+1)-10**agesarr2(j))
        swgt(j) = sfh(j)*dt
        twgt = twgt+swgt(j)
     ENDDO

     !ensure the weights sum to unity.  This is 
     !automatically the case for a constant SFH but not for
     !a tau model b/c of numerical errors in the integrals
     swgt = swgt/twgt

  ELSE

     !SSP at age=tau
     ii = locate(agesarr,tau)
     swgt     = 10**prlo_sfh
     swgt(ii) = 1.0

  ENDIF

  !transfer the parameters to the parameter array 
  pos(1) = lebv
  pos(2) = 1E-3
  pos(1+nxpar:nxpar+nage) = LOG10(swgt) + mpix
  pos(1+nxpar+nage)  = zmet

  IF (1.EQ.1) THEN
     pos(1) = -2.34+1
     pos(2) = -0.3
     pos(3) = -3.05
     pos(4) = -1.91
     pos(5) = -1.48
     pos(6) = -0.87
     pos(7) = -0.33
     pos(8) =  0.86
     pos(9) =  1.48
     pos(10) = -0.05
  ENDIF

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
