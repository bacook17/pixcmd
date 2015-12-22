PROGRAM WRITE_A_MODEL

  ! Routine to create one model Hess diagram convolved
  ! with an MDF, a SFH, and a P(Mpix)

  USE nrtype; USE pixcmd_utils; USE pixcmd_vars
  USE nr, ONLY : locate
  IMPLICIT NONE
  
  INTEGER  :: i,j,k,ii,jj,wgti,ilo=1,imax,jmax
  REAL(SP) :: mpix,zmet0,zmets,tau,tot=1.,dt,twgt=0.,wgtmax=0.
  CHARACTER(50) :: infile
  REAL(SP), DIMENSION(nx,ny) :: onemodel
  REAL(SP), DIMENSION(npix,npix,nfil) :: f1,cf1,of1
  REAL(SP), DIMENSION(nage)  :: sfh
  REAL(SP), DIMENSION(nz,nage)  :: wgt
  REAL(SP), DIMENSION(nz)    :: mdf,zz
  INTEGER, DIMENSION(npix,2)  :: sarr
  CHARACTER(10) :: time

  !------------------------------------------------------------!

  infile = 'model_M2.0_cSFH'
  ii = 1
  jj = 5

  !initialize the random number generator
  CALL INIT_RANDOM_SEED()

  !setup the model grid
  CALL SETUP_MODELS(1)

  !set up the random indices
  DO i=1,npix
     sarr(i,:) = i
  ENDDO
  CALL SHUFFLE(sarr(:,1))
  CALL SHUFFLE(sarr(:,2))

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

  f1 = 0.0
  DO j=1,nage
     DO i=1,nz
        IF (j.EQ.1) THEN
           dt = (10**agesarr(j)-10**(agesarr(j)-dage))
        ELSE
           dt = (10**agesarr(j)-10**agesarr(j-1))
        ENDIF
        wgt(i,j) = mdf(i)*sfh(j)*dt
        IF (wgt(i,j).GT.wgtmax) THEN
           wgtmax = wgt(i,j)
           imax = i
           jmax = j
        ENDIF
        twgt = twgt+wgt(i,j)
     ENDDO
  ENDDO

  wgt(imax,jmax)=0.0
  f1 = model(:,:,:,imax,jmax)

  write(*,*) wgt

stop
  DO j=1,nage
     DO i=1,nz
        wgti = INT(wgt(i,j)*npix**2)
        DO k=1,wgti
           ii = INT(myran()*npix)
           jj = INT(myran()*npix)
           f1(ii,jj,:) = model(ii,jj,:,i,j)
        ENDDO
     ENDDO
  ENDDO


  !convolve with PSF
  cf1 = convolve(f1)
  
  !add obs errors
  of1 = add_obs_err(cf1)

  !compute Hess diagram
  onemodel = hist_2d(of1(:,:,1)-of1(:,:,2),of1(:,:,2),&
       xhess,yhess,nx,ny,npix)

  !save the PSF-convolved, obs err-included Hess diagram
  OPEN(11,FILE=TRIM(PIXCMD_HOME)//'/data/'//TRIM(infile)//'.im',&
       FORM='UNFORMATTED',STATUS='REPLACE',access='direct',&
          recl=npix*npix*nfil*4)
  WRITE(11,rec=1) of1
  CLOSE(11)

  !save the Hess diagram to file
  OPEN(1,FILE=TRIM(PIXCMD_HOME)//'data/'//TRIM(infile)//'.hess',&
       FORM='UNFORMATTED',STATUS='REPLACE',access='direct',&
       recl=nx*ny*4)
  WRITE(1,rec=1) onemodel
  CLOSE(1)


END PROGRAM WRITE_A_MODEL
