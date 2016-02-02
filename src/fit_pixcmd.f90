PROGRAM FIT_PIXCMD

  !To Do: 1) include E(B-V) and Z as free parameters
  !syntax: mpirun -np XX fit_pixcmd.exe input_data Mpix0 lebv0 [Z/H] seed tag

  USE pixcmd_utils; USE pixcmd_vars; USE nrtype
  USE nr, ONLY : ran1,gasdev; USE mpi
  USE ran_state, ONLY : ran_seed,ran_init

  IMPLICIT NONE

  !key emcee parameters
  INTEGER, PARAMETER :: nwalkers=256,nburn=1000,nmcmc=100
  !flag for testing clock time
  INTEGER, PARAMETER :: test_time=0
  !fit for tau-Mpix
  INTEGER, PARAMETER :: dotaufit=1
  !fix the SFH=const
  INTEGER, PARAMETER :: doinitsfh=0

  !down-sample the data by this factor
  REAL(SP), PARAMETER :: subsample=1.0
  !initialization width, and re-initializtion width
  REAL(SP), PARAMETER :: wdth0=1E-2,wdth1=1E-3

  INTEGER  :: i,j,k,ml,stat,iter=30,totacc=0,npos,ii
  REAL(SP) :: dt,cmin,cstd,minchi2=huge_number,ndat,tmp
  REAL(SP) :: time1,time2,twgt=0.0,zmet0
  CHARACTER(10) :: time,is,tmpstr
  CHARACTER(2)  :: tis
  CHARACTER(50) :: infile,tag=''
  REAL(SP), DIMENSION(2) :: dumt,dumt2
  REAL(SP), DIMENSION(nx,ny) :: bmodel=0.,imodel=0.
  REAL(SP), DIMENSION(nage)  :: sfh,wgt
  REAL(SP), DIMENSION(npix,npix) :: im

  !emcee variables
  REAL(SP), DIMENSION(npar) :: bpos=-99.,dumpos=-99.
  REAL(SP), DIMENSION(npar,nwalkers) :: pos_emcee_in,pos_emcee_out
  REAL(SP), DIMENSION(nwalkers)      :: lp_emcee_in,lp_emcee_out,lp_mpi,tchi2
  INTEGER,  DIMENSION(nwalkers)      :: accept_emcee
  REAL(SP), DIMENSION(npar,nwalkers) :: mpiposarr=0.0

  !MPI variables
  INTEGER :: ierr,taskid,ntasks,received_tag,status(MPI_STATUS_SIZE)
  INTEGER :: KILL=99,BEGIN=0
  LOGICAL :: wait=.TRUE.
  INTEGER, PARAMETER :: masterid=0

  !------------------------------------------------------------!

  !Initialize MPI, and get the total number of processes and
  !your process number
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ntasks,ierr)

  IF (IARGC().LT.1) THEN
     !infile='m31_bulge'
     infile='model_M2.0_SFH1_tau10.0'
     mpix0=2.0
  ELSE
     CALL GETARG(1,infile)
  ENDIF

  IF (IARGC().GE.2) THEN
     CALL GETARG(2,tmpstr)
     READ(tmpstr,'(F4.1)') mpix0
  ENDIF

  IF (IARGC().GE.3) THEN
     CALL GETARG(3,tmpstr)
     READ(tmpstr,'(F4.1)') lebv0
  ENDIF

  IF (IARGC().GE.4) THEN
     CALL GETARG(4,tmpstr)
     READ(tmpstr,'(F4.1)') zmet0
  ENDIF

  IF (IARGC().GE.5) THEN
     CALL GETARG(5,tmpstr)
     READ(tmpstr,'(I5)') fix_seed
  ENDIF

  IF (IARGC().EQ.6) THEN
     tag(1:1)='_'
     CALL GETARG(6,tag(2:))
  ENDIF


  IF (ntasks.EQ.1) THEN
     WRITE(*,*) 'ERROR: you are not using mpirun!'
     STOP
  ENDIF

  IF (taskid.EQ.masterid) THEN
     !write some important variables to screen
     WRITE(*,*)
     WRITE(*,'(" ************************************")')
     WRITE(*,'("   Mpix_init  = ",1x,F4.1)') mpix0
     WRITE(*,'("   lebv_init  = ",1x,F4.1)') lebv0
     WRITE(*,'("   [Z/H]_init = ",1x,F4.1)') zmet0
     WRITE(*,'("   Npix       = ",I5)') npix
     WRITE(*,'("   dotaufit   = ",I5)') dotaufit
     WRITE(*,'("   doinitsfh  = ",I5)') doinitsfh
     WRITE(*,'("   Nwalkers   = ",I5)') nwalkers
     WRITE(*,'("   Nburn      = ",I5)') nburn
     WRITE(*,'("   Nchain     = ",I5)') nmcmc
     WRITE(*,'("   Ntasks     = ",I5)') ntasks
     WRITE(*,'("   Rseed      = ",I5)') fix_seed
     WRITE(*,'("   filename   = ",A)') TRIM(infile)//TRIM(tag)
     WRITE(*,'(" ************************************")')
  ENDIF

  !transfer the priors to the prior array
  prlo(1)                  = prlo_lebv
  prlo(2)                  = prlo_lebvw
  prlo(1+nxpar:nxpar+nage) = prlo_sfh+mpix0
  prlo(1+nxpar+nage)       = prlo_zmet

  prhi(1)                  = prhi_lebv
  prhi(2)                  = prhi_lebvw
  prhi(1+nxpar:nxpar+nage) = prhi_sfh+mpix0
  prhi(1+nxpar+nage)       = prhi_zmet

  !initialize the random number generator
  !set each task to sleep for a different length of time
  !so that each task has its own unique random number seed
  IF (fix_seed.EQ.0) THEN
     CALL SLEEP(taskid)
  ENDIF
  CALL INIT_RANDOM_SEED()
  DO i=1,niso_max
     CALL RAN1(ranarr(:,i))
  ENDDO
  DO i=1,npix
     CALL GASDEV(gdev(:,i))
     CALL GASDEV(gdev2(:,i))
  ENDDO

  !now that the ranarr is identically initialized,
  !re-set the seed on each taskid for the emcee steps
  IF (fix_seed.NE.0) THEN
     fix_seed=0
     CALL SLEEP(taskid)
     CALL INIT_RANDOM_SEED()
  ENDIF

  !setup the model grid, PSF, etc.
  CALL SETUP_MODELS()

  !read in the Hess diagram for the data
  OPEN(1,IOSTAT=stat,FILE=TRIM(PIXCMD_HOME)//'/data/'//&
       TRIM(infile)//'.hess',FORM='UNFORMATTED',STATUS='OLD',&
       ACCESS='direct',recl=nx*ny*4,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'ERROR: input file not found:'
     WRITE(*,*) TRIM(infile)//'.hess'
     STOP
  ENDIF
  READ(1,rec=1) hess_data
  CLOSE(1)

  !down-sample the data
  IF (subsample.LT.1.0) THEN
     hess_data = hess_data * subsample
  ENDIF

  ndat = SUM(hess_data)

  !Poisson error at each CMD pixel
  hess_err = SQRT(hess_data)

  !test for bad values and inflate errors as needed
  DO i=1,nx
     DO j=1,ny
        IF (hess_data(i,j).LE.tiny_number) hess_err(i,j)=1.0
        IF (hess_data(i,j).LT.0.0) THEN
           WRITE(*,*) 'ERROR: input data has negative Hess entry'
           STOP
        ENDIF
        IF (hess_data(i,j).GT.tiny_number.AND.hess_data(i,j).LE.2.0) &
             hess_err(i,j)=hess_err(i,j)*10
     ENDDO
  ENDDO

  !normalize the data to unity
  hess_err  = hess_err  / ndat
  hess_data = hess_data / ndat

  ! The worker's only job is to calculate the value of a function
  ! after receiving a parameter vector.
  IF (taskid.NE.masterid) THEN
     
     ! Start event loop
     DO WHILE (wait)

        ! Look for data from the master. This call can accept up
        ! to nwalkers paramater positions, but it expects
        ! that the actual number of positions is smaller and is
        ! given by the MPI_TAG.  This call does not return until
        ! a set of parameter vectors is received
        CALL MPI_RECV(npos, 1, MPI_INTEGER, &
             masterid, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        received_tag = status(MPI_TAG)
        IF ((received_tag.EQ.KILL).OR.(npos.EQ.0)) EXIT
        CALL MPI_RECV(mpiposarr(1,1), npos*npar, MPI_REAL, &
             masterid, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
   
        CALL CPU_TIME(time1)
        !Calculate the probability for these parameter positions
        DO k=1,npos
           !note: func returns chi^2
           lp_mpi(k) = -0.5*func(mpiposarr(:,k))
        ENDDO
        CALL CPU_TIME(time2)

        IF (test_time.EQ.1) THEN
           WRITE(*,'(" Task ID ",I3": Elapsed Time: ",F6.2," s", '//&
                '", N=",I2,", chi^2=",99ES11.3)') &
                taskid,time2-time1,npos,-2.0*lp_mpi(1:npos)
           CALL FLUSH()
        ENDIF

        !Send it back to the master
        CALL MPI_SEND(lp_mpi(1), npos, MPI_REAL, &
             masterid, BEGIN, MPI_COMM_WORLD, ierr)

     ENDDO

  ENDIF

  !this is the master process
  IF (taskid.EQ.masterid) THEN

     CALL DATE_AND_TIME(TIME=time)
     WRITE(*,*) 'Start Time '//time(1:2)//':'//time(3:4)//':'//time(5:6)

     !----------------------Intial Best-Fit--------------------------!

     IF (dotaufit.EQ.1) THEN

        !fit in a grid of tau and Mpix
        WRITE(*,*) 'Running tau-mpix fitter'
        CALL FIT_TAU(bpos,zmet0)

     ELSE IF (doinitsfh.EQ.1) THEN

        !initialize with a constant SFH
        sfh = 1/10**agesarr2(nage+1)
        DO j=1,nage
           dt = (10**agesarr2(j+1)-10**agesarr2(j))
           wgt(j) = sfh(j)*dt
           twgt = twgt+wgt(j)
        ENDDO
        wgt = wgt/twgt
        !transfer the parameters to the parameter array
        bpos(1) = lebv0
        bpos(2) = -2.0
        bpos(1+nxpar:nxpar+nage) = LOG10(wgt)+mpix0
        bpos(1+nxpar+nage)  = zmet0

     ELSE

        bpos = -99.

     ENDIF
     
     !!test smoothness of chi^2 surface
     !immed  = 1E9
     !imodel = getmodel(bpos,im)
     !immed  = SUM(10**(-2./5*im))/npix**2
     !write(*,*) immed
     !ii = 1
     !tmp = bpos(ii)
     !write(*,*) tmp
     !DO i=1,20
     !   bpos(ii) = tmp + 0.01*i-0.1
     !   write(*,*) bpos(ii),func(bpos)
     !   write(tis,'(I2)') i+10
     ! !  bmodel = getmodel(bpos)
     ! !  OPEN(11,FILE=TRIM(PIXCMD_HOME)//'/results2/test'//TRIM(tis)&
     ! !       //'.hess',FORM='UNFORMATTED',STATUS='REPLACE',&
     ! !       access='DIRECT',recl=nx*ny*4)
     ! !  WRITE(11,rec=1) bmodel
     ! !  CLOSE(11)
     !ENDDO
     !STOP


     !-------------------------------------------------------------------!
     !---------------------------Run emcee-------------------------------!
     !-------------------------------------------------------------------!

     !setup the starting positions
     WRITE(*,*) 'initial parameters:'

     IF (dotaufit.EQ.1.OR.doinitsfh.EQ.1) THEN
        !initialize parameters near the best-fit intializer
        DO j=1,nwalkers
           DO i=1,npar
              pos_emcee_in(i,j) = bpos(i) + wdth0*(2.*myran()-1.0)
           ENDDO
           WRITE(*,'(30(F5.2,1x))') pos_emcee_in(:,j)
        ENDDO
     ELSE
        !initialize randomly across parameter space
        DO j=1,nwalkers
           DO i=1,npar
              pos_emcee_in(i,j) = myran()*(prhi(i)-prlo(i)-6*wdth0) + &
                   (prlo(i)+3*wdth0)
           ENDDO
           WRITE(*,'(30(F6.2,1x))') pos_emcee_in(:,j)
        ENDDO
     ENDIF

     !Compute the initial log-probability for each walker
     CALL FUNCTION_PARALLEL_MAP(npar,nwalkers,ntasks-1,&
          pos_emcee_in,lp_emcee_in)

     WRITE(*,*) 'chi^2 for initialized walkers:'
     WRITE(*,'(10(ES10.3,1x))') -2.0*lp_emcee_in

     !tchi2 = -2.0*lp_emcee_in !LOG10(-2.0*lp_emcee_in)
     !cstd  = SQRT(SUM( (tchi2-SUM(tchi2)/nwalkers)**2 )/(nwalkers-1))
     !write(*,*) cstd

     IF (-2.0*MAXVAL(lp_emcee_in).GE.huge_number) THEN
        WRITE(*,*) 'FIT_PIXCMD ERROR: initial parameters are out of bounds'
        STOP
     ENDIF


     !---------------------initial burn-in---------------------!

     WRITE(*,'(A)') ' initial burn-in...  '
     DO i=1,nburn/5
        IF (test_time.EQ.1) THEN
           WRITE(*,'("Iteration ",I4)') i
           CALL FLUSH()
        ENDIF
        CALL EMCEE_ADVANCE_MPI(npar,nwalkers,2.0,pos_emcee_in,&
             lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee,ntasks-1)
        pos_emcee_in = pos_emcee_out
        lp_emcee_in  = lp_emcee_out
     ENDDO

     WRITE(*,*) 'parameters after first-pass:'
     DO j=1,nwalkers
        WRITE(*,'(30(F5.2,1x))') pos_emcee_in(:,j)
     ENDDO

     WRITE(*,*) 'chi^2 after first-pass:'
     WRITE(*,'(10(ES10.3,1x))') -2.0*lp_emcee_in

     !take the min(chi^2) and re-initialize a ball around
     !the minimum, if nburn>100
     IF (nburn.GT.100) THEN

        ml    = MINLOC(tchi2,1)    
        bpos  = pos_emcee_in(:,ml) 

        DO j=1,nwalkers
           DO i=1,npar
              pos_emcee_in(i,j) = bpos(i)+wdth1*(2.*myran()-1.0)
              IF (pos_emcee_in(i,j).LT.prlo(i)) &
                   pos_emcee_in(i,j)=prlo(i)+2*wdth1
              IF (pos_emcee_in(i,j).GT.prhi(i)) &
                   pos_emcee_in(i,j)=prhi(i)-2*wdth1
           ENDDO
        ENDDO

        WRITE(*,*) 're-initialized parameters:'
        DO j=1,nwalkers
           WRITE(*,'(30(F5.2,1x))') pos_emcee_in(:,j)
        ENDDO

        !Compute the initial log-probability for each walker
        CALL FUNCTION_PARALLEL_MAP(npar,nwalkers,ntasks-1,&
             pos_emcee_in,lp_emcee_in)

        WRITE(*,*) 'chi^2 for re-initialized walkers:'
        WRITE(*,'(10(ES10.3,1x))') -2.0*lp_emcee_in

     ENDIF

     !-------------------second-pass burn-in-------------------!

     WRITE(*,'(A)',advance='no') ' second burn-in: '
     DO i=nburn/5,nburn
        IF (test_time.EQ.1) WRITE(*,'("Iteration ",I4)') i
        CALL EMCEE_ADVANCE_MPI(npar,nwalkers,2.0,pos_emcee_in,&
             lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee,ntasks-1)
        pos_emcee_in = pos_emcee_out
        lp_emcee_in  = lp_emcee_out
        IF (i.EQ.nburn/4.*1) THEN
           WRITE (*,'(A)',advance='no') '...25%'
           CALL FLUSH()
        ENDIF
        IF (i.EQ.nburn/4.*2) THEN
           WRITE (*,'(A)',advance='no') '...50%'
           CALL FLUSH()
        ENDIF
        IF (i.EQ.nburn/4.*3) THEN
           WRITE (*,'(A)',advance='no') '...75%'
           CALL FLUSH()
        ENDIF
     ENDDO
     WRITE (*,'(A)') '...100%'
     CALL FLUSH()

     WRITE(*,*) 'chi^2 after second-pass:'
     WRITE(*,'(10(ES10.3,1x))') -2.0*lp_emcee_in

     !-------------------production chain--------------------!

     OPEN(12,FILE=TRIM(PIXCMD_HOME)//'/results2/'//&
          TRIM(infile)//TRIM(tag)//'.mcmc',STATUS='REPLACE')

     !write a file header
     WRITE(12,'("#   Mpix_init  = ",1x,F4.1)') mpix0
     WRITE(12,'("#   lebv_init  = ",1x,F4.1)') lebv0
     WRITE(12,'("#   Npix       = ",I5)') npix
     WRITE(12,'("#   dotaufit   = ",I5)') dotaufit
     WRITE(12,'("#   doinitsfh  = ",I5)') doinitsfh
     WRITE(12,'("#   Nwalkers   = ",I5)') nwalkers
     WRITE(12,'("#   Nburn      = ",I5)') nburn
     WRITE(12,'("#   Nchain     = ",I5)') nmcmc
     WRITE(12,'("#   Ntasks     = ",I5)') ntasks
     WRITE(12,'(2I3)') nage, nzi
     WRITE(12,'(20(F5.2,1x))') agesarr


     WRITE(*,'(A)',advance='no') ' production run: '       
     DO i=1,nmcmc
        IF (test_time.EQ.1) WRITE(*,'("Iteration ",I4)') i        
        CALL EMCEE_ADVANCE_MPI(npar,nwalkers,2.0,pos_emcee_in,&
             lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee,ntasks-1)
        pos_emcee_in = pos_emcee_out
        lp_emcee_in  = lp_emcee_out
        totacc = totacc + SUM(accept_emcee)
        DO j=1,nwalkers
           !keep the model with the lowest chi2
           IF (-2.0*lp_emcee_in(j).LT.minchi2) THEN
              bpos = pos_emcee_in(:, j)
              minchi2 = -2.0*lp_emcee_in(j)
           ENDIF
           !write the chain elements to file
           IF (-2.0*lp_emcee_in(j).EQ.huge_number) THEN
              WRITE(12,'(F10.6,1x,999(F7.4,1x))') &
                   LOG10(-2.0*lp_emcee_in(j)),dumpos
           ELSE
              WRITE(12,'(F10.6,1x,999(F7.4,1x))') &
                   LOG10(-2.0*lp_emcee_in(j)),pos_emcee_in(:, j)
           ENDIF
        ENDDO
        IF (i.EQ.nmcmc/4.*1) THEN
           WRITE (*,'(A)',advance='no') '...25%'
           CALL FLUSH()
        ENDIF
        IF (i.EQ.nmcmc/4.*2) THEN
           WRITE (*,'(A)',advance='no') '...50%'
           CALL FLUSH()
        ENDIF
        IF (i.EQ.nmcmc/4.*3) THEN
           WRITE (*,'(A)',advance='no') '...75%'
           CALL FLUSH()
        ENDIF
     ENDDO
     WRITE (*,'(A)') '...100%'
     CALL FLUSH()

     CLOSE(12)

     WRITE(*,'("  Facc: ",F6.3)') REAL(totacc)/REAL(nmcmc*nwalkers)
     
     !write the best model to a binary file
     bmodel = getmodel(bpos)
     OPEN(11,FILE=TRIM(PIXCMD_HOME)//'/results2/'//TRIM(infile)&
          //TRIM(tag)//'.hess',FORM='UNFORMATTED',STATUS='REPLACE',&
          access='DIRECT',recl=nx*ny*4)
     WRITE(11,rec=1) bmodel
     CLOSE(11)
     
     CALL DATE_AND_TIME(TIME=time)
     WRITE(*,*) 'End Time '//time(1:2)//':'//time(3:4)//':'//time(5:6)

     !break the workers out of their event loops so they can close
     CALL FREE_WORKERS(ntasks-1)
     
  ENDIF

  CALL MPI_FINALIZE(ierr)
 

END PROGRAM FIT_PIXCMD
