PROGRAM FIT_PIXCMD

  !To Do: 1) include E(B-V) as free parameter
  !       2) MPI parallelize

  !Note: 1E6 pix is not quite enough compared to 2E7.

  USE pixcmd_utils; USE pixcmd_vars; USE nrtype
  USE nr, ONLY : powell
  USE ran_state, ONLY : ran_seed,ran_init

  IMPLICIT NONE

  !emcee variables
  INTEGER, PARAMETER :: nwalkers=100,nburn=100,nmcmc=100
  REAL(SP), DIMENSION(npar,nwalkers) :: pos_emcee_in,pos_emcee_out
  REAL(SP), DIMENSION(nwalkers)      :: lp_emcee_in,lp_emcee_out
  INTEGER,  DIMENSION(nwalkers)      :: accept_emcee

  INTEGER :: i,j,ndat,stat,i1,i2,iter=30,totacc=0
  REAL(SP) :: mpix,fret,bret=huge_number,wdth=0.1
  CHARACTER(10) :: time
  CHARACTER(50) :: infile
  REAL(SP), DIMENSION(nx,ny) :: bmodel=0.

  !Powell parameters
  INTEGER, PARAMETER :: dopowell=0
  REAL(SP), PARAMETER :: ftol=0.01
  REAL(SP), DIMENSION(npar,npar) :: xi=0.0
  REAL(SP), DIMENSION(npar)      :: pos=0.0,bpos=0.

  !------------------------------------------------------------!

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
        IF (hess_data(i,j).LE.tiny_number) hess_err(i,j)=huge_number
     ENDDO
  ENDDO

  !normalize the data to unity
  hess_err  = hess_err /ndat
  hess_data = hess_data/ndat
  
  CALL DATE_AND_TIME(TIME=time)
  WRITE(*,*) 'Start Time '//time(1:2)//':'//time(3:4)//':'//time(5:6)
 
  !---------------------Run powell minimization-----------------------!

  IF (dopowell.EQ.1) THEN
  
     DO j=1,10
        !setup params
        pos(1) = myran()+1.5
        DO i=2,npar
           pos(i) = LOG10(myran()/npar)
        ENDDO
        xi=0.0
        DO i=1,npar
           xi(i,i) = 1E-2
        ENDDO
        fret = huge_number
        CALL POWELL(pos,xi,ftol,iter,fret)
        IF (fret.LT.bret) THEN
           bret = fret
           bpos = pos
        ENDIF
     ENDDO
     WRITE(*,'(200F10.5)') log10(bret/(nx*ny-npar))!,bpos

  ELSE

     DO i=1,npar
        bpos(i) = LOG10(myran()/npar)
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
     CALL EMCEE_ADVANCE(npar,nwalkers,2.0,pos_emcee_in,&
          lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee)
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
     CALL EMCEE_ADVANCE(npar,nwalkers,2.0,pos_emcee_in,&
          lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee)
     pos_emcee_in = pos_emcee_out
     lp_emcee_in  = lp_emcee_out
  ENDDO


  OPEN(12,FILE=TRIM(PIXCMD_HOME)//'/results2/fit.mcmc',STATUS='REPLACE')

  !production chain
  WRITE(*,*) 'production run...'       
  DO i=1,nmcmc
     CALL EMCEE_ADVANCE(npar,nwalkers,2.0,pos_emcee_in,&
          lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee)
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


END PROGRAM FIT_PIXCMD
