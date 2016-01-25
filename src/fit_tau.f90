SUBROUTINE FIT_TAU(pos,zmet0)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : func

  IMPLICIT NONE

  INTEGER, PARAMETER :: ifprint=0

  REAL(SP), DIMENSION(npar), INTENT(inout) :: pos
  REAL(SP), INTENT(in) :: zmet0
  INTEGER :: i,j
  INTEGER, PARAMETER :: ntau=7,nmpix=10
  REAL(SP), DIMENSION(npar) :: tpos
  REAL(SP), DIMENSION(nage) :: sfh,wgt
  REAL(SP) :: bestchi2,chi2,tau,twgt,dt,btau,tmpix,bmpix,zmet
  REAL(SP), DIMENSION(ntau) :: tauarr

 !------------------------------------------------------------!

  bestchi2 = huge_number

  tauarr = (/1.0,2.0,3.0,5.0,7.0,10.0,20.0/)

  IF (ifprint.EQ.1) WRITE(*,*) '   chi2    tau  Mpix'

  DO i=1,ntau

     tau = tauarr(i)

     sfh = EXP(-(10**10.201-10**agesarr)/1E9/tau) / tau/1E9

     twgt=0.
     DO j=1,nage
        IF (j.EQ.1) THEN
           dt = (10**agesarr(j)-10**(agesarr(j)-dage))
        ELSE
           dt = (10**agesarr(j)-10**agesarr(j-1))
        ENDIF
        wgt(j) = sfh(j)*dt
        twgt = twgt+wgt(j)
     ENDDO

     !ensure the weights sum to unity.  This is 
     !automatically the case for a constant SFH but not for
     !a tau model b/c of numerical errors in the integrals
     wgt = wgt/twgt
     DO j=1,nage
        wgt(j) = MAX(wgt(j),10**(prlo_sfh+1))
     ENDDO

     DO j=1,nmpix

        tpos(1) = lebv0
        tmpix = mpix0+(REAL(j)-1.)/nmpix-0.5
        tpos(1+nxpar:nxpar+nage) = LOG10(wgt)+tmpix
        tpos(1+nxpar+nage:npar)  = zmet0

        chi2 = func(tpos)
        IF (chi2.LT.bestchi2) THEN
           bestchi2 = chi2
           pos      = tpos
           btau     = tau
           bmpix    = tmpix
        ENDIF

        IF (ifprint.EQ.1) &
             WRITE(*,'(ES10.3,50F6.2)') chi2,tau,tmpix

     ENDDO

  ENDDO

  
  WRITE(*,'("fit_tau: tau=",F4.1,", Mpix=",F4.1)') btau,bmpix


END SUBROUTINE FIT_TAU
