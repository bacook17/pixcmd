FUNCTION GET_MODEL(inpos)

  USE pixcmd_vars; USE nrtype
  USE nr, ONLY : locate
  IMPLICIT NONE

  REAL(SP), DIMENSION(npar) :: inpos
  REAL(SP), DIMENSION(nx,ny) :: get_model,m1,m2
  INTEGER :: ilo
  REAL(SP) :: dt, ff

  !------------------------------------------------------------!

  ilo = MAX(MIN(locate(model_ages,inpos(6)),nage-1),1)
  dt  = (inpos(6)-model_ages(ilo))/&
       (model_ages(ilo+1)-model_ages(ilo))

  !the weight of the 6th metallicity point is defined
  !so that the total is one.
  ff = 1-SUM(10**inpos(1:5))
 
  m1 = 10**inpos(1)*model(1,ilo,:,:) + &
       10**inpos(2)*model(2,ilo,:,:) + &
       10**inpos(3)*model(3,ilo,:,:) + &
       10**inpos(4)*model(4,ilo,:,:) + &
       10**inpos(5)*model(5,ilo,:,:) + &
       ff*model(6,ilo,:,:) 

  m2 = 10**inpos(1)*model(1,ilo+1,:,:) + &
       10**inpos(2)*model(2,ilo+1,:,:) + &
       10**inpos(3)*model(3,ilo+1,:,:) + &
       10**inpos(4)*model(4,ilo+1,:,:) + &
       10**inpos(5)*model(5,ilo+1,:,:) + &
       ff*model(6,ilo+1,:,:) 

  get_model = m1*(1-dt) + m2*dt


END FUNCTION GET_MODEL
