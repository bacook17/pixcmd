PROGRAM FIT_PIXCMD

  !To Do: 1) include E(B-V) and Z as free parameters
  !syntax: mpirun -np XX fit_pixcmd.exe input_data Mpix_init tag

  USE pixcmd_utils; USE pixcmd_vars; USE nrtype
  USE nr, ONLY : powell,ran1; USE mpi
  USE ran_state, ONLY : ran_seed,ran_init

  IMPLICIT NONE

  !key emcee parameters
  INTEGER, PARAMETER :: nwalkers=64,nburn=10000,nmcmc=100

  !starting guess for the Mpix parameter
  REAL(SP) :: mpix=2.0

  !flag for testing clock time
  INTEGER, PARAMETER :: test_time=0
  !Powell minimization
  INTEGER, PARAMETER :: dopowell=0
  !fit for tau-Mpix
  INTEGER, PARAMETER :: dotaufit=0
  !fix the SFH=const
  INTEGER, PARAMETER :: doinitsfh=0

  INTEGER  :: i,j,k,ml,ndat,stat,iter=30,totacc=0,npos
  REAL(SP) :: fret,bret=huge_number,dt,cmin,cmean,cstd,minchi2=huge_number
  CHARACTER(10) :: time,is,tmpstr
  REAL(SP) :: time1,time2
  REAL(SP), DIMENSION(2) :: dumt,dumt2
  CHARACTER(50) :: infile,tag=''
  REAL(SP), DIMENSION(nx,ny) :: bmodel=0.,imodel=0.
  REAL(SP), DIMENSION(nage) :: sfh,wgt

  !emcee variables
  REAL(SP), DIMENSION(npar,nwalkers) :: pos_emcee_in,pos_emcee_out
  REAL(SP), DIMENSION(nwalkers)      :: lp_emcee_in,lp_emcee_out,lp_mpi,tchi2
  INTEGER,  DIMENSION(nwalkers)      :: accept_emcee
  REAL(SP), DIMENSION(npar,nwalkers) :: mpiposarr=0.0

  !Powell variables
  REAL(SP), PARAMETER :: ftol=10.
  REAL(SP), DIMENSION(npar,npar) :: xi=0.0
  REAL(SP), DIMENSION(npar)      :: pos=0.0,bpos=0.,dum9=-9.0

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
  ELSE
     CALL GETARG(1,infile)
  ENDIF

  IF (IARGC().GE.2) THEN
     CALL GETARG(2,tmpstr)
     READ(tmpstr,'(F4.1)') mpix
  ENDIF

  IF (IARGC().EQ.3) THEN
     tag(1:1)='_'
     CALL GETARG(3,tag(2:))
  ENDIF

  IF (ntasks.EQ.1) THEN
     WRITE(*,*) 'ERROR: you are not using mpirun!'
     STOP
  ENDIF

  IF (taskid.EQ.masterid) THEN
     !write some important variables to screen
     WRITE(*,*)
     WRITE(*,'(" ************************************")')
     WRITE(*,'("   Mpix_init  = ",1x,F4.1)') mpix
     WRITE(*,'("   dopowell   = ",I5)') dopowell
     WRITE(*,'("   dotaufit   = ",I5)') dotaufit
     WRITE(*,'("   doinitsfh  = ",I5)') doinitsfh
     WRITE(*,'("   Nwalkers   = ",I5)') nwalkers
     WRITE(*,'("   Nburn      = ",I5)') nburn
     WRITE(*,'("   Nchain     = ",I5)') nmcmc
     WRITE(*,'("   Ntasks     = ",I5)') ntasks
     WRITE(*,'("   filename   = ",A)') TRIM(infile)//TRIM(tag)
     WRITE(*,'(" ************************************")')
  ENDIF


  !initialize the random number generator
  !set each task to sleep for a different length of time
  !so that each task has its own unique random number seed
  CALL SLEEP(taskid)
  CALL INIT_RANDOM_SEED()
  CALL RAN1(ranarr)

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
  ndat = INT(SUM(hess_data))

  !Poisson error at each CMD pixel
  hess_err = SQRT(hess_data)
  DO i=1,nx
     DO j=1,ny
        IF (hess_data(i,j).LE.tiny_number) hess_err(i,j)=1.0
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
           WRITE(*,'(" Task ID ",I2": Elapsed Time: ",F6.2," s", '//&
                '", N=",I2,", chi^2=",20ES11.3)') &
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


     !----------------------Initialization--------------------------!

     IF (dopowell.EQ.1) THEN

        !Powell minimization
        WRITE(*,*) 'Running Powell minimization'
        
        DO j=1,10
           !setup params
           pos(1) = (myran()-0.5)+mpix
           DO i=1+nxpar,npar
              pos(i) = myran()*(prhi-prlo-3*wdth0) + (prlo+1.5*wdth0)
           ENDDO
           pos(1+nxpar:npar) = pos(1+nxpar:npar) - &
                LOG10(SUM(10**pos(1+nxpar:npar)))
           xi=0.0
           xi(1,1) = 0.2
           DO i=1+nxpar,npar
              xi(i,i) = 1.0
           ENDDO
           fret = huge_number
           CALL POWELL(pos,xi,ftol,iter,fret)
           WRITE(*,'(50F10.5)') LOG10(fret),pos
           IF (fret.LT.bret) THEN
              bret = fret
              bpos = pos
           ENDIF
        ENDDO
        WRITE(*,'(5F10.5)') LOG10(bret),log10(bret/(nx*ny-npar))

     ELSE IF (dotaufit.EQ.1) THEN

        !fit in a grid of tau and Mpix
        WRITE(*,*) 'Running tau-mpix fitter'
        CALL FIT_TAU(bpos,mpix)

     ELSE IF (doinitsfh.EQ.1) THEN

        !initialize with a constant SFH
        sfh = 1/10**agesarr(nage)
        DO j=1,nage
           IF (j.EQ.1) THEN
              dt = (10**agesarr(j)-10**(agesarr(j)-dage))
           ELSE
              dt = (10**agesarr(j)-10**agesarr(j-1))
           ENDIF
           wgt(j) = sfh(j)*dt
        ENDDO
        !transfer the parameters to the parameter array
        bpos(1)      = mpix
        bpos(1+nxpar:npar) = LOG10(wgt)

     ELSE

        bpos = -99.

     ENDIF
     
     !-------------------------------------------------------------------!
     !---------------------------Run emcee-------------------------------!
     !-------------------------------------------------------------------!

     !setup the starting positions
     WRITE(*,*) 'initial parameters:'

     IF (dotaufit.EQ.1.OR.dopowell.EQ.1.OR.doinitsfh.EQ.1) THEN
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
           pos_emcee_in(1,j) = (myran()-0.5)+mpix
           DO i=1+nxpar,npar
              pos_emcee_in(i,j) = myran()*(prhi-prlo-3*wdth0) + &
                   (prlo+1.5*wdth0)
           ENDDO
           WRITE(*,'(30(F5.2,1x))') pos_emcee_in(:,j)
        ENDDO
     ENDIF

     !Compute the initial log-probability for each walker
     CALL FUNCTION_PARALLEL_MAP(npar,nwalkers,ntasks-1,&
          pos_emcee_in,lp_emcee_in)

     WRITE(*,*) 'chi^2 for initialized walkers:'
     WRITE(*,'(10(ES10.3,1x))') -2.0*lp_emcee_in

     !tchi2 = LOG10(-2.0*lp_emcee_in)
     !cstd  = SQRT(SUM( (tchi2-SUM(tchi2)/nwalkers)**2 )/(nwalkers-1))
     !write(*,*) cstd

     IF (-2.0*MAXVAL(lp_emcee_in).EQ.huge_number) THEN
        WRITE(*,*) 'FIT_PIXCMD ERROR: initial parameters are out of bounds'
        STOP
     ENDIF


     !---------------------initial burn-in---------------------!

     WRITE(*,'(A)',advance='no') ' initial burn-in:  '
     DO i=1,nburn/10
        IF (test_time.EQ.1) THEN
           WRITE(*,'("Iteration ",I4)') i
           CALL FLUSH()
        ENDIF
        CALL EMCEE_ADVANCE_MPI(npar,nwalkers,2.0,pos_emcee_in,&
             lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee,ntasks-1)
        pos_emcee_in = pos_emcee_out
        lp_emcee_in  = lp_emcee_out
        IF (i.EQ.nburn/4.*1) WRITE (*,'(A)',advance='no') ' ...25%'
        IF (i.EQ.nburn/4.*2) WRITE (*,'(A)',advance='no') '...50%'
        IF (i.EQ.nburn/4.*3) WRITE (*,'(A)',advance='no') '...75%'
     ENDDO
     WRITE (*,'(A)') '...100%'
     CALL FLUSH()

     WRITE(*,*) 'parameters after first-pass:'
     DO j=1,nwalkers
        WRITE(*,'(30(F5.2,1x))') pos_emcee_in(:,j)
     ENDDO

     WRITE(*,*) 'chi^2 after first-pass:'
     WRITE(*,'(10(ES10.3,1x))') -2.0*lp_emcee_in

     tchi2 = LOG10(-2.0*lp_emcee_in)
     cmin  = MINVAL(tchi2) !min chi^2
     cmean = SUM(tchi2)/nwalkers !mean chi^2
     cstd  = SQRT( SUM( (tchi2-cmean)**2 )/(nwalkers-1) ) !stddev chi^2
     ml    = MINLOC(tchi2,1) !location of min chi^2
     bpos  = pos_emcee_in(:,ml) !params at min chi^2
     k=0
     DO j=1,nwalkers
        !if chi^2 > min+stddev, then re-initialize
        IF (tchi2(j).GT.(cmin+cstd)) THEN
           DO i=1,npar
              pos_emcee_in(i,j) = bpos(i)+wdth0/5.*(2.*myran()-1.0)
              IF (i.EQ.1) CYCLE
              IF (pos_emcee_in(i,j).LT.prlo) &
                   pos_emcee_in(i,j)=prlo+wdth0/5.
              IF (pos_emcee_in(i,j).GT.prhi) &
                   pos_emcee_in(i,j)=prhi-wdth0/5.
           ENDDO
           k=k+1
        ENDIF
     ENDDO

     WRITE(*,'("Reinitialized ",I2," out of ",I3," walkers")') k,nwalkers

     !Compute the initial log-probability for each walker
     CALL FUNCTION_PARALLEL_MAP(npar,nwalkers,ntasks-1,&
          pos_emcee_in,lp_emcee_in)

     WRITE(*,*) 'chi^2 for re-initialized walkers:'
     WRITE(*,'(10(ES10.3,1x))') -2.0*lp_emcee_in

     !-------------------second-pass burn-in-------------------!

     WRITE(*,'(A)',advance='no') ' second burn-in: '
     DO i=1,nburn
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
                   LOG10(-2.0*lp_emcee_in(j)),dum9
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
