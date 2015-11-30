PROGRAM CREATE_A_MODEL

  ! Routine to create one model Hess diagram convolved
  ! with an MDF, a SFH, and a P(Mpix)

  USE nrtype; USE pixcmd_utils; USE pixcmd_vars
  USE nr, ONLY : locate
  IMPLICIT NONE
  
  INTEGER  :: ilo,i,j,k
  REAL(SP) :: di,mpix,wgt,zmet0,zmets,tau,tot
  CHARACTER(50) :: infile
  REAL(SP), DIMENSION(nx,ny) :: onemodel,m1,m2
  REAL(SP), DIMENSION(nage)  :: sfh
  REAL(SP), DIMENSION(nz)    :: mdf,zz

  !------------------------------------------------------------!

  infile = 'model_tau2.0_mpix3.5_mdf0.0_0.05'
  mpix   = 3.5
  zmet0  = 0.0
  zmets  = 0.05
  tau    = 0.5

  !setup the model grid
  CALL SETUP_MODELS()

  !(nm,nz,nage,nx,ny)
  !mpixarr, zmetarr, agesarr

  !MDF
  zz  = LOG10(zmetarr/0.0190)
  mdf = EXP(-(zz-zmet0)**2/2/zmets)

  !SFH
  sfh = EXP(-(10**10.201-10**agesarr)/1E9/tau)
  sfh=0.
  sfh(20)=1.

  tot=0.
  DO i=1,nage
     DO j=1,nz
        tot = tot + mdf(j)*sfh(i)
     ENDDO
  ENDDO

  m1 = 0.0
  m2 = 0.0

  ilo = MAX(MIN(locate(mpixarr,mpix),nm-1),1)
  di  = (mpix-mpixarr(ilo))/(mpixarr(ilo+1)-mpixarr(ilo))

  DO j=1,nage
     DO i=1,nz
        wgt = mdf(i)*sfh(j)/tot
        m1  = m1+wgt*model(ilo,i,j,:,:)
        m2  = m2+wgt*model(ilo+1,i,j,:,:)
     ENDDO
  ENDDO
  
  onemodel = m1*(1-di) + m2*di

  onemodel=0.0
  DO i=1,nx
     DO j=1,ny
        onemodel(i,j) = TSUM(zz,mdf*model(ilo,:,20,i,j))
     ENDDO
  ENDDO


  !save the Hess diagram to file
  OPEN(1,FILE=TRIM(PIXCMD_HOME)//'data/'//TRIM(infile)//'.hess',&
       FORM='UNFORMATTED',STATUS='REPLACE',access='direct',&
       recl=nx*ny*4)
  WRITE(1,rec=1) onemodel
  CLOSE(1)



END PROGRAM CREATE_A_MODEL
