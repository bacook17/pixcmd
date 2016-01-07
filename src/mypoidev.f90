FUNCTION MYPOIDEV(xm)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : myran
  IMPLICIT NONE
  
  INTEGER :: i,j,k=1
  REAL(SP), INTENT(IN) :: xm
  REAL(SP) :: em,t,g
  REAL(SP), DIMENSION(npix,npix) :: mypoidev

  !------------------------------------------------------------!

  g = EXP(-xm)

  DO i=1,npix
     DO j=1,npix
        em = -1.0
        t  = 1.0
        DO
           em = em+1.0
           IF (true_poisson.EQ.1) THEN
              t  = t*myran()
           ELSE
              t = t*ranarr(k)
           ENDIF
           k = MOD((k+1),nran)
           IF (t.LE.g) EXIT
        ENDDO
        mypoidev(i,j)=em
     ENDDO
  ENDDO

END FUNCTION MYPOIDEV
