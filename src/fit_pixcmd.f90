PROGRAM FIT_PIXCMD

  !To Do: 1) include E(B-V) as free parameter
  !       2) MPI parallelize

  USE pixcmd_utils; USE pixcmd_vars; USE nrtype
  USE nr, ONLY : powell; USE mpi
  USE ran_state, ONLY : ran_seed,ran_init

  IMPLICIT NONE

  !emcee variables
  INTEGER, PARAMETER :: nwalkers=256,nburn=10,nmcmc=1000
  REAL(SP), DIMENSION(npar,nwalkers) :: pos_emcee_in,pos_emcee_out
  REAL(SP), DIMENSION(nwalkers)      :: lp_emcee_in,lp_emcee_out
  INTEGER,  DIMENSION(nwalkers)      :: accept_emcee,lp_mpi
  REAL(SP), DIMENSION(npar,nwalkers) :: mpiposarr=0.0

  INTEGER :: i,j,k,ndat,stat,i1,i2,iter=30,totacc=0,npos
  REAL(SP) :: mpix,fret,bret=huge_number,wdth=0.1
  CHARACTER(10) :: time
  CHARACTER(50) :: infile
  REAL(SP), DIMENSION(nx,ny) :: bmodel=0.

  !Powell parameters
  INTEGER, PARAMETER :: dopowell=0
  REAL(SP), PARAMETER :: ftol=0.1
  REAL(SP), DIMENSION(npar,npar) :: xi=0.0
  REAL(SP), DIMENSION(npar)      :: pos=0.0,bpos=0.

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

  !initialize the random number generator
  CALL INIT_RANDOM_SEED()

  !setup the model grid
  CALL SETUP_MODELS()

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
        IF (hess_data(i,j).LE.tiny_number) hess_err(i,j)=1.0!huge_number
     ENDDO
  ENDDO

  !normalize the data to unity
  hess_err  = hess_err /ndat
  hess_data = hess_data/ndat


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
        CALL MPI_RECV(mpiposarr(1,1), npos*npar, MPI_DOUBLE_PRECISION, &
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
        CALL MPI_SEND(lp_mpi(1), npos, MPI_DOUBLE_PRECISION, &
             masterid, BEGIN, MPI_COMM_WORLD, ierr)

     ENDDO

  ENDIF

  !this is the master process
  IF (taskid.EQ.masterid) THEN

     CALL DATE_AND_TIME(TIME=time)
     WRITE(*,*) 'Start Time '//time(1:2)//':'//time(3:4)//':'//time(5:6)

     !----------------------single value fitter--------------------------!

     IF (1.EQ.0) THEN
        
        bpos=-8.0
        k=1
        DO j=1,nage
           DO i=1,nz
              bpos(k)=0.0
              fret = func(bpos)
              bpos(k)=-8.0
              k=k+1
              write(*,*) i,j,fret,LOG10(fret)
           ENDDO
        ENDDO
        
     ENDIF

     !---------------------Run powell minimization-----------------------!

     IF (dopowell.EQ.1) THEN
        
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
           WRITE(*,'(F10.5)') log10(fret/(nx*ny-npar))
           IF (fret.LT.bret) THEN
              bret = fret
              bpos = pos
           ENDIF
        ENDDO
        WRITE(*,'(5F10.5)') Log10(bret),log10(bret/(nx*ny-npar))
        
     ELSE
        
        !DO i=1,npar
        !   bpos(i) = LOG10(myran()/npar)
        !ENDDO
        bpos=-4.0
        k=1
        DO i=1,nage
           DO j=1,nz
              IF (i.GE.(nage-2).AND.j.GE.3) bpos(k) = LOG10(1/5.)
              k=k+1
           ENDDO
        ENDDO
        
     ENDIF
     
     !-------------------------Run emcee---------------------------------!
     
     !initialize the walkers
     DO j=1,nwalkers
        DO i=1,npar
           pos_emcee_in(i,j) = bpos(i) + wdth*(2.*myran()-1.0)
        ENDDO
        !Compute the initial log-probability for each walker
        lp_emcee_in(j) = -0.5*func(pos_emcee_in(:, j))
     ENDDO
  
     !initial burn-in the chain
     WRITE(*,*) 'first burn-in...'
     DO i=1,nburn/2-1
        CALL EMCEE_ADVANCE_MPI(npar,nwalkers,2.0,pos_emcee_in,&
             lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee,ntasks-1)
        pos_emcee_in = pos_emcee_out
        lp_emcee_in  = lp_emcee_out
     ENDDO

     !prune the walkers and re-initialize
     i = MAXLOC(lp_emcee_in,1)
     bpos = pos_emcee_in(:,i)
     DO j=1,nwalkers
        DO i=1,npar
           pos_emcee_in(i,j) = bpos(i) + wdth*(2.*myran()-1.0)
        ENDDO
        !Compute the initial log-probability for each walker
        lp_emcee_in(j) = -0.5*func(pos_emcee_in(:, j))
     ENDDO

     !second-pass burn-in the chain
     WRITE(*,*) 'second burn-in...'
     DO i=nburn/2,nburn
        CALL EMCEE_ADVANCE_MPI(npar,nwalkers,2.0,pos_emcee_in,&
             lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee,ntasks-1)
        pos_emcee_in = pos_emcee_out
        lp_emcee_in  = lp_emcee_out
     ENDDO
     

     OPEN(12,FILE=TRIM(PIXCMD_HOME)//'/results2/fit.mcmc',STATUS='REPLACE')

     !production chain
     WRITE(*,*) 'production run...'       
     DO i=1,nmcmc
        CALL EMCEE_ADVANCE_MPI(npar,nwalkers,2.0,pos_emcee_in,&
             lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee,ntasks-1)
        pos_emcee_in = pos_emcee_out
        lp_emcee_in  = lp_emcee_out
        totacc = totacc + SUM(accept_emcee)
        
        !write the chain elements to file
        DO j=1,nwalkers
           WRITE(12,'(ES12.5,1x,999(F9.4,1x))') &
                -2.0*lp_emcee_in(j),pos_emcee_in(:, j)
        ENDDO
        
     ENDDO

     WRITE(*,'("  Facc: ",F5.2)') REAL(totacc)/REAL(nmcmc*nwalkers)
     
     !write the best model to a binary file
     bmodel = getmodel(bpos)
     OPEN(11,FILE=TRIM(PIXCMD_HOME)//'/results2/fit.hess',&
          FORM='UNFORMATTED',STATUS='REPLACE',access='DIRECT',&
          recl=nx*ny*4)
     WRITE(11,rec=1) bmodel
     CLOSE(11)

     CLOSE(12)
     
     CALL DATE_AND_TIME(TIME=time)
     WRITE(*,*) 'End Time '//time(1:2)//':'//time(3:4)//':'//time(5:6)

     !break the workers out of their event loops so they can close
     CALL FREE_WORKERS(ntasks-1)
     
  ENDIF

  CALL MPI_FINALIZE(ierr)
 


END PROGRAM FIT_PIXCMD
