MODULE PIXCMD_UTILS

  INTERFACE
     FUNCTION ADD_OBS_ERR(flux,dm,exptime,zpt)
       USE pixcmd_vars; USE nrtype
       REAL(SP), DIMENSION(nfil,npix,npix), INTENT(in) :: flux
       REAL(SP), INTENT(in) :: dm
       REAL(SP), DIMENSION(nfil), INTENT(in) :: exptime,zpt
       REAL(SP), DIMENSION(nfil,npix,npix) :: add_obs_err
     END FUNCTION ADD_OBS_ERR
  END INTERFACE

  INTERFACE
     FUNCTION CONVOLVE(arr,psf,na,nb,np)
       USE nrtype
       INTEGER, INTENT(in) :: na,np
       REAL(SP), DIMENSION(nb,na,na), INTENT(in) :: arr
       REAL(SP), DIMENSION(np,np), INTENT(in) :: psf  
       REAL(SP), DIMENSION(nb,np,np) :: convolve
     END FUNCTION CONVOLVE
  END INTERFACE

  INTERFACE
     SUBROUTINE EMCEE_ADVANCE(ndim,nwalkers,a,pin,lpin,&
          pout,lpout,accept)
       USE pixcmd_vars; USE nrtype
       INTEGER, INTENT(in) :: ndim, nwalkers
       REAL(SP), INTENT(in) :: a
       REAL(SP), INTENT(in), DIMENSION(ndim,nwalkers) :: pin
       REAL(SP), INTENT(in), DIMENSION(nwalkers) :: lpin
       REAL(SP), INTENT(out), DIMENSION(ndim,nwalkers) :: pout
       REAL(SP), INTENT(out), DIMENSION(nwalkers) :: lpout
       INTEGER, INTENT(out), DIMENSION(nwalkers) :: accept
     END SUBROUTINE EMCEE_ADVANCE
  END INTERFACE

  INTERFACE
     FUNCTION FUNC(inpos)
       USE nrtype
       REAL(SP), DIMENSION(:) :: inpos
       REAL(SP) :: func
     END FUNCTION FUNC
  END INTERFACE

  INTERFACE
     FUNCTION GET_MODEL(inpos)
       USE pixcmd_vars; USE nrtype
       REAL(SP), DIMENSION(npar) :: inpos
       REAL(SP), DIMENSION(nx,ny) :: get_model
     END FUNCTION GET_MODEL
  END INTERFACE

  INTERFACE
     FUNCTION HIST_2D(xx,yy,xarr,yarr,nx,ny,npix)
       USE nrtype
       INTEGER, INTENT(in) :: nx,ny,npix
       REAL(SP), DIMENSION(npix,npix), INTENT(in) :: xx,yy
       REAL(SP), DIMENSION(nx), INTENT(in) :: xarr
       REAL(SP), DIMENSION(ny), INTENT(in) :: yarr
       REAL(SP), DIMENSION(nx,ny) :: hist_2d
     END FUNCTION HIST_2D
  END INTERFACE

  INTERFACE
     FUNCTION MYRAN()
       USE nr, ONLY : ran1; USE nrtype
       REAL(SP) :: myran
     END FUNCTION MYRAN
  END INTERFACE
  
END MODULE PIXCMD_UTILS

