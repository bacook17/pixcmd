PROGRAM FIT_PIXCMD

  !To Do: 1) fit the 2D space of age-Z
  !       2) include E(B-V), Mpix as free parameters

  !Note: 1E6 pix is not quite enough compared to 2E7.

  USE pixcmd_utils; USE pixcmd_vars; USE nrtype
  USE nr, ONLY : locate,powell
  USE ran_state, ONLY : ran_seed,ran_init

  IMPLICIT NONE

  INTEGER, PARAMETER :: nwalkers=100,nburn=300,nmcmc=200
  INTEGER, PARAMETER :: dopowell=0
  INTEGER :: i,j,ndat,stat,i1,i2,iter=30,totacc=0
  REAL(SP) :: mpix,fret,bret=huge_number,wdth=0.1
  REAL(SP), DIMENSION(nz) :: zmet
  CHARACTER(10) :: time
  CHARACTER(6)  :: zstr
  CHARACTER(4)  :: mstr
  
  REAL(SP), DIMENSION(nx) :: xarr=0.
  REAL(SP), DIMENSION(ny) :: yarr=0.
  REAL(SP), DIMENSION(ndat_max) :: xdat=0.,ydat=0.
  REAL(SP), DIMENSION(nx,ny) :: bmodel=0.

  !Powell iteration tolerance
  REAL(SP), PARAMETER :: ftol=0.01
  REAL(SP), DIMENSION(npar,npar) :: xi=0.0
  REAL(SP), DIMENSION(npar)      :: pos=0.0,bpos=0.

  !emcee variables
  REAL(SP), DIMENSION(npar,nwalkers) :: pos_emcee_in,pos_emcee_out
  REAL(SP), DIMENSION(nwalkers)      :: lp_emcee_in,lp_emcee_out
  INTEGER,  DIMENSION(nwalkers)  :: accept_emcee
  !------------------------------------------------------------!

  !initialize the random number generator
  CALL INIT_RANDOM_SEED()

  mpix = 1E2
  WRITE(mstr,'(F4.2)') LOG10(mpix)
  zmet = (/0.0010,0.0025,0.0040,0.0080,0.0190,0.0290/)

  !read in the model Hess diagrams
  DO i=1,nz
     WRITE(zstr,'(F6.4)') zmet(i)
     OPEN(11,FILE='../hess/hess_M'//mstr//'_Z'//zstr//'.dat',&
          FORM='UNFORMATTED',STATUS='OLD',access='direct',&
          recl=nage*nx*ny*4,ACTION='READ')
     READ(11,rec=1) model(i,:,:,:)
     CLOSE(11)
     !normalize each model to unity and compute Poisson err
     DO j=1,nage
        model(i,j,:,:) = model(i,j,:,:)/npix**2
     ENDDO
  ENDDO

  !set up the Hess arrays
  DO i=1,nx
     xarr(i) = xmin+(i-1)*dx
  ENDDO
  DO i=1,ny
     yarr(i) = ymin+(i-1)*dy
  ENDDO

  !set up model ages array
  DO i=1,nage
     model_ages(i) = age0+(i-1)*dage
  ENDDO

  !read in the data
  OPEN(12,FILE='../data/m31_bulge.dat',STATUS='OLD',iostat=stat,&
       ACTION='READ')
  DO i=1,ndat_max
     READ(12,*,IOSTAT=stat) xdat(i),ydat(i)
     IF (stat.NE.0) GOTO 20
  ENDDO
  WRITE(*,*) 'FIT_PIXCMD ERROR: did not finish reading in data file'
  STOP
20 CONTINUE
  ndat = i-1
  CLOSE(12)

  !create a Hess diagram for the data
  hess_data=0.0
  DO i=1,ndat
     IF (xdat(i).LT.xarr(1).OR.xdat(i).GT.xarr(nx).OR.&
          ydat(i).LT.yarr(1).OR.ydat(i).GT.yarr(ny)) CYCLE
     i1 = locate(xarr,xdat(i))
     i2 = locate(yarr,ydat(i))
     hess_data(i1,i2) = hess_data(i1,i2)+1.
  ENDDO

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
        DO i=1,npar
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
  WRITE(*,*) 'burning in...'
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
  WRITE(*,*) 'burning in...'
  DO i=nburn/2,nburn
     CALL EMCEE_ADVANCE(npar,nwalkers,2.0,pos_emcee_in,&
          lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee)
     pos_emcee_in = pos_emcee_out
     lp_emcee_in  = lp_emcee_out
  ENDDO


  OPEN(12,FILE='../results2/fit.mcmc',STATUS='REPLACE')

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
  bmodel = get_model(bpos)
  OPEN(11,FILE='../results2/fit.hess',&
       FORM='UNFORMATTED',STATUS='REPLACE',access='DIRECT',&
       recl=nx*ny*4)
  WRITE(11,rec=1) bmodel
  CLOSE(11)

  CLOSE(12)

  CALL DATE_AND_TIME(TIME=time)
  WRITE(*,*) 'End Time '//time(1:2)//':'//time(3:4)//':'//time(5:6)


END PROGRAM FIT_PIXCMD
