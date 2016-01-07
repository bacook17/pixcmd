SUBROUTINE FIT_TAU(pos,mpix)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : func

  IMPLICIT NONE

  REAL(SP), INTENT(in) :: mpix
  REAL(SP), DIMENSION(npar), INTENT(inout) :: pos
  INTEGER :: i,j
  INTEGER, PARAMETER :: ntau=15,nmpix=10
  REAL(SP), DIMENSION(npar) :: tpos
  REAL(SP), DIMENSION(nage) :: sfh,wgt
  REAL(SP) :: bestchi2,chi2,tau,twgt,dt,btau

 !------------------------------------------------------------!

  bestchi2 = huge_number

  DO i=1,ntau

     tau = REAL(i)

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
        wgt(j) = MAX(wgt(j),10**prlo)
     ENDDO

     tpos(1+nxpar:npar) = LOG10(wgt)

     DO j=1,nmpix

        tpos(1) = mpix+(real(j)-1.)/nmpix-0.5

        chi2 = func(tpos)
        IF (chi2.LT.bestchi2) THEN
           bestchi2 = chi2
           pos      = tpos
           btau     = tau
        ENDIF

        !WRITE(*,'(ES10.3,50F6.2)') chi2,tau,tpos(1)

     ENDDO

  ENDDO

  WRITE(*,'("fit_tau: tau=",F4.1,", Mpix=",F4.1)') btau,pos(1)


END SUBROUTINE FIT_TAU
