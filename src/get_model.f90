FUNCTION GET_MODEL(inpos)

  USE pixcmd_vars; USE nrtype
  USE nr, ONLY : locate
  IMPLICIT NONE

  REAL(SP), DIMENSION(npar) :: inpos
  REAL(SP), DIMENSION(nx,ny) :: get_model,m1,m2
  INTEGER :: ilo,i,j,k
  REAL(SP) :: di

  !------------------------------------------------------------!

  m1 = 0.0
  m2 = 0.0

  ilo = MAX(MIN(locate(mpixarr,inpos(1)),nm-1),1)
  di  = (inpos(1)-mpixarr(ilo))/(mpixarr(ilo+1)-mpixarr(ilo))

  k=1
  DO j=1,nage
     DO i=1,nz
        m1 = m1+10**inpos(k)*model(ilo,i,j,:,:)
        m2 = m2+10**inpos(k)*model(ilo+1,i,j,:,:)
        k=k+1
     ENDDO
  ENDDO
  
  get_model = m1*(1-di) + m2*di


END FUNCTION GET_MODEL
