PROGRAM FIT_PIXCMD

  !To Do: 1) include E(B-V) and Mpix as free parameters

  USE pixcmd_utils; USE pixcmd_vars; USE nrtype
  USE nr, ONLY : powell; USE mpi
  USE ran_state, ONLY : ran_seed,ran_init

  IMPLICIT NONE

  !emcee variables
  INTEGER, PARAMETER :: nwalkers=512,nburn=500,nmcmc=50
  REAL(SP), DIMENSION(npar,nwalkers) :: pos_emcee_in,pos_emcee_out
  REAL(SP), DIMENSION(nwalkers)      :: lp_emcee_in,lp_emcee_out,lp_mpi
  INTEGER,  DIMENSION(nwalkers)      :: accept_emcee
  REAL(SP), DIMENSION(npar,nwalkers) :: mpiposarr=0.0

  INTEGER  :: i,j,k,ndat,stat,i1,i2,iter=30,totacc=0,npos
  REAL(SP) :: fret,bret=huge_number
  CHARACTER(10) :: time
  CHARACTER(50) :: infile,tag=''
  REAL(SP), DIMENSION(nx,ny) :: bmodel=0.

  !Powell parameters
  INTEGER, PARAMETER  :: dopowell=1
  REAL(SP), PARAMETER :: ftol=0.1
  REAL(SP), DIMENSION(npar,npar) :: xi=0.0
  REAL(SP), DIMENSION(npar)      :: pos=0.0,bpos=0.,dum9=-9.0

  !variables for MPI
  INTEGER :: ierr,taskid,ntasks,received_tag,status(MPI_STATUS_SIZE)
  INTEGER :: KILL=99,BEGIN=0
  LOGICAL :: wait=.TRUE.
  INTEGER, PARAMETER :: masterid=0
  INTEGER, PARAMETER :: test_time=0

  !------------------------------------------------------------!

  ! Initialize MPI, and get the total number of processes and
  ! your process number
  CALL MPI_INIT( ierr )
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, taskid, ierr )
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, ntasks, ierr )

  IF (IARGC().LT.1) THEN
     infile='m31_bulge'
  ELSE
     CALL GETARG(1,infile)
  ENDIF

  IF (IARGC().GT.1) THEN
     tag(1:1)='_'
     CALL GETARG(2,tag(2:))
  ENDIF

  IF (taskid.EQ.masterid) THEN
     !write some important variables to screen
     WRITE(*,*)
     WRITE(*,'(" ************************************")')
     WRITE(*,'("  dopowell   = ",I5)') dopowell
     WRITE(*,'("  Nwalkers   = ",I5)') nwalkers
     WRITE(*,'("  Nburn      = ",I5)') nburn
     WRITE(*,'("  Nchain     = ",I5)') nmcmc
     WRITE(*,'("  Ntasks     = ",I5)') ntasks
     WRITE(*,'("  filename   = ",A)') TRIM(infile)//TRIM(tag)
     WRITE(*,'(" ************************************")')
  ENDIF


  !initialize the random number generator
  CALL INIT_RANDOM_SEED()

  !setup the model grid, PSF, etc.
  CALL SETUP_MODELS(1)

  !read in the Hess diagram for the data
  OPEN(1,FILE=TRIM(PIXCMD_HOME)//'/data/'//TRIM(infile)//'.hess',&
       FORM='UNFORMATTED',STATUS='OLD',access='direct',&
       recl=nx*ny*4,ACTION='READ')
  READ(1,rec=1) hess_data
  CLOSE(1)
  ndat = SUM(hess_data)

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
        ! to ``nwalkers`` paramater positions, but it expects
        ! that the actual number of positions is smaller and is
        ! given by the MPI_TAG.  This call does not return until
        ! a set of parameter vectors is received
        CALL MPI_RECV(npos, 1, MPI_INTEGER, &
             masterid, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        received_tag = status(MPI_TAG)
        IF ((received_tag.EQ.KILL).OR.(npos.EQ.0)) EXIT
        CALL MPI_RECV(mpiposarr(1,1), npos*npar, MPI_REAL, &
             masterid, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
   
        IF (taskid.EQ.1.AND.test_time.EQ.1) THEN
           CALL DATE_AND_TIME(TIME=time)
           WRITE(*,*) '1 Time '//time(1:2)//':'//time(3:4)//':'&
                //time(5:9),npos,taskid
        ENDIF

        !Calculate the probability for these parameter positions
        DO k=1,npos
           lp_mpi(k) = -0.5*func(mpiposarr(:,k))
        ENDDO

         IF (taskid.EQ.1.AND.test_time.EQ.1) THEN
           CALL DATE_AND_TIME(TIME=time)
           WRITE(*,*) '2 Time '//time(1:2)//':'//time(3:4)//':'&
                //time(5:9),npos,taskid
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

     IF (1.EQ.0) THEN
 
        !single age-Z minimization
        WRITE(*,*) 'Running single age-Z minimization'

        bpos=-8.0
        k=1
        DO j=1,nage
           DO i=1,nz
              bpos(k)=0.0
              fret = func(bpos)
              bpos(k)=-8.0
              k=k+1
              write(*,*) i,j,LOG10(fret)
           ENDDO
        ENDDO

     ELSE IF (dopowell.EQ.0) THEN

        !Powell minimization
        WRITE(*,*) 'Running Powell minimization'
        
        DO j=1,10
           !setup params
           DO i=1,npar
              pos(i) = LOG10(myran()/npar)
           ENDDO
           xi=0.0
           DO i=1,npar
              xi(i,i) = 1E-2
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

     ELSE
        
        !random initialization

        DO i=1,npar
           bpos(i) = myran()*(prhi-prlo-3*wdth0) + (prlo+1.5*wdth0)
        ENDDO
        
     ENDIF
     
     !-------------------------Run emcee---------------------------------!
     
     !initialize the lp
     DO j=1,nwalkers
        DO i=1,npar
           pos_emcee_in(i,j) = bpos(i) + wdth0*(2.*myran()-1.0)
        ENDDO
     ENDDO
     !Compute the initial log-probability for each walker
     CALL FUNCTION_PARALLEL_MAP(npar,nwalkers,ntasks-1,&
          pos_emcee_in,lp_emcee_in)

     !initial burn-in
     WRITE(*,'(A)',advance='no') ' first burn-in:  '
     DO i=1,nburn
        CALL EMCEE_ADVANCE_MPI(npar,nwalkers,2.0,pos_emcee_in,&
             lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee,ntasks-1)
        pos_emcee_in = pos_emcee_out
        lp_emcee_in  = lp_emcee_out
        IF (i.EQ.nburn/4.*1) THEN
           WRITE (*,'(A)',advance='no') ' ...25%'
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

     WRITE(*,*) 'log(chi^2) after first-pass:'
     WRITE(*,'(10(ES10.2,1x))') -2.0*lp_emcee_in

     !prune the walkers and re-initialize
     i    = MAXLOC(lp_emcee_in,1)
     bpos = pos_emcee_in(:,i)
     DO j=1,nwalkers
        DO i=1,npar
           pos_emcee_in(i,j) = bpos(i)+wdth0/5.*(2.*myran()-1.0)
           IF (pos_emcee_in(i,j).LT.(prlo+2*wdth0)) &
                pos_emcee_in(i,j)=prlo+2*wdth0
           IF (pos_emcee_in(i,j).GT.(prhi-wdth0)) &
                pos_emcee_in(i,j)=prhi-wdth0
        ENDDO
     ENDDO
     !Compute the initial log-probability for each walker
     CALL FUNCTION_PARALLEL_MAP(npar,nwalkers,ntasks-1,&
          pos_emcee_in,lp_emcee_in)

     WRITE(*,*) 'log(chi^2) for re-initialized walkers:'
     WRITE(*,'(10(ES10.2,1x))') -2.0*lp_emcee_in

     !second-pass burn-in
     WRITE(*,'(A)',advance='no') ' second burn-in: '
     DO i=1,nburn
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

     OPEN(12,FILE=TRIM(PIXCMD_HOME)//'/results2/'//&
          TRIM(infile)//TRIM(tag)//'.mcmc',STATUS='REPLACE')

     !production chain
     WRITE(*,'(A)',advance='no') ' production run: '       
     DO i=1,nmcmc
        CALL EMCEE_ADVANCE_MPI(npar,nwalkers,2.0,pos_emcee_in,&
             lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee,ntasks-1)
        pos_emcee_in = pos_emcee_out
        lp_emcee_in  = lp_emcee_out
        totacc = totacc + SUM(accept_emcee)
        !write the chain elements to file
        DO j=1,nwalkers
           IF (-2.0*lp_emcee_in(j).GT.1E31) THEN
              WRITE(12,'(ES12.5,1x,999(F9.4,1x))') &
                   -2.0*lp_emcee_in(j),dum9
           ELSE
              WRITE(12,'(ES12.5,1x,999(F9.4,1x))') &
                   -2.0*lp_emcee_in(j),pos_emcee_in(:, j)
           ENDIF
        ENDDO
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

     CLOSE(12)

     WRITE(*,'("  Facc: ",F5.2)') REAL(totacc)/REAL(nmcmc*nwalkers)
     
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
