FUNCTION GET_MODEL(inpos)

  USE pixcmd_vars; USE nrtype
  USE pixcmd_utils, ONLY : convolve,hist_2d,add_obs_err
  USE nr, ONLY : locate
  IMPLICIT NONE

  REAL(SP), DIMENSION(npar) :: inpos
  REAL(SP), DIMENSION(nx,ny) :: get_model,m1,m2
  REAL(SP), DIMENSION(npix,npix,nfil) :: f1,f2
  INTEGER :: ilo,i,j,k
  REAL(SP) :: di

  !------------------------------------------------------------!

  f1 = 0.0
  f2 = 0.0

  ilo = MAX(MIN(locate(mpixarr,inpos(1)),nm-1),1)
  di  = (inpos(1)-mpixarr(ilo))/(mpixarr(ilo+1)-mpixarr(ilo))

  k=1
  DO j=1,nage
     DO i=1,nz
        f1 = f1+10**inpos(k)*model(ilo,i,j,:,:,:)
        f2 = f2+10**inpos(k)*model(ilo+1,i,j,:,:,:)
        k=k+1
     ENDDO
  ENDDO
  
  !convolve with PSF
  f1 = -2.5*LOG10(convolve(f1))
  f2 = -2.5*LOG10(convolve(f2))

  !add obs errors
  f1 = add_obs_err(f1)
  f2 = add_obs_err(f2)

  !compute Hess diagram
  m1 = hist_2d(f1(:,:,1)-f1(:,:,2),f1(:,:,2),&
       xhess,yhess,nx,ny,npix)
  m2 = hist_2d(f2(:,:,1)-f2(:,:,2),f2(:,:,2),&
       xhess,yhess,nx,ny,npix)


  get_model = m1*(1-di) + m2*di

  !normalize each model to unity
  !DO j=1,nage
  !   model(m,i,j,:,:) = model(m,i,j,:,:)/npix**2
  !ENDDO



END FUNCTION GET_MODEL
