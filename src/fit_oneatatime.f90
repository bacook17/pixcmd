SUBROUTINE FIT_ONEATATIME(pos)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : func

  IMPLICIT NONE

  REAL(SP), DIMENSION(npar), INTENT(inout) :: pos
  INTEGER :: i,j,im,jbest
  INTEGER, PARAMETER :: ns=18
  REAL(SP), DIMENSION(npar) :: tpos
  REAL(SP) :: dpp=0.5,chimax, bestchi
  REAL(SP), DIMENSION(ns) :: pp, chi2

 !------------------------------------------------------------!

  tpos   = prlo
  chimax = func(tpos)

  DO i=1,npar

     tpos    = prlo
     bestchi = chimax
     jbest   = 1

     DO j=1,ns

        tpos(i) = prlo + j*dpp
        chi2(j) = func(tpos)
        IF (chi2(j).LT.bestchi) THEN
           bestchi = chi2(j)
           jbest = j
        ENDIF

     ENDDO

     pos(i) = prlo+jbest*dpp

     write(*,'(F5.2,1x,F5.2,ES12.3)') agesarr(i),pos(i),bestchi

  ENDDO
  


END SUBROUTINE FIT_ONEATATIME
