MODULE PIXCMD_UTILS

  INTERFACE
     FUNCTION ADD_OBS_ERR(flux)
       USE pixcmd_vars; USE nrtype
       REAL(SP), DIMENSION(npix,npix,nfil), INTENT(in) :: flux
       REAL(SP), DIMENSION(npix,npix,nfil) :: add_obs_err
     END FUNCTION ADD_OBS_ERR
  END INTERFACE

  INTERFACE
     FUNCTION CONVOLVE(arr)
       USE pixcmd_vars; USE nrtype
       REAL(SP), DIMENSION(npix,npix,nfil), INTENT(inout) :: arr
       REAL(SP), DIMENSION(npix,npix,nfil) :: convolve
     END FUNCTION CONVOLVE
  END INTERFACE

  INTERFACE
     FUNCTION DRAWN(nn)
       USE pixcmd_vars; USE nrtype
       REAL(SP), INTENT(in) :: nn
       REAL(SP), DIMENSION(npix,npix) :: drawn
     END FUNCTION DRAWN
  END INTERFACE

  INTERFACE
     SUBROUTINE EMCEE_ADVANCE_MPI(ndim,nwalkers,a,pin,lpin,&
          pout,lpout,accept,nworkers)
       USE nrtype
       INTEGER, INTENT(in) :: ndim, nwalkers,nworkers
       REAL(SP), INTENT(in) :: a
       REAL(SP), INTENT(in), DIMENSION(ndim,nwalkers) :: pin
       REAL(SP), INTENT(in), DIMENSION(nwalkers) :: lpin
       REAL(SP), INTENT(out), DIMENSION(ndim,nwalkers) :: pout
       REAL(SP), INTENT(out), DIMENSION(nwalkers) :: lpout
       INTEGER, INTENT(out), DIMENSION(nwalkers) :: accept
     END SUBROUTINE EMCEE_ADVANCE_MPI
  END INTERFACE

  INTERFACE
     SUBROUTINE FIT_ONEATATIME(pos)
       USE pixcmd_vars; USE nrtype
       REAL(SP), DIMENSION(npar), INTENT(inout) :: pos
     END SUBROUTINE FIT_ONEATATIME
  END INTERFACE

  INTERFACE
     SUBROUTINE FREE_WORKERS(nworkers)
       INTEGER, INTENT(in) :: nworkers
     END SUBROUTINE FREE_WORKERS
  END INTERFACE

  INTERFACE
     FUNCTION FUNC(inpos)
       USE nrtype
       REAL(SP), DIMENSION(:) :: inpos
       REAL(SP) :: func
     END FUNCTION FUNC
  END INTERFACE

  INTERFACE
     SUBROUTINE FUNCTION_PARALLEL_MAP(ndim,nk,nworkers,pos,lnpout)
       USE nrtype
       INTEGER, INTENT(in) :: ndim, nk, nworkers
       REAL(SP), INTENT(in), DIMENSION(ndim,nk) :: pos
       REAL(SP), INTENT(out), DIMENSION(nk) :: lnpout
    END SUBROUTINE FUNCTION_PARALLEL_MAP
  END INTERFACE

  INTERFACE
     FUNCTION GETMODEL(inpos,im)
       USE pixcmd_vars; USE nrtype
       REAL(SP), DIMENSION(npar), INTENT(in) :: inpos
       REAL(SP), DIMENSION(npix,npix), OPTIONAL :: im
       REAL(SP), DIMENSION(nx,ny) :: getmodel
     END FUNCTION GETMODEL
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
  
  INTERFACE
     FUNCTION MYPOIDEV(xm)
       USE pixcmd_vars; USE nrtype
       REAL(SP), INTENT(IN) :: xm
       REAL(SP), DIMENSION(npix,npix) :: mypoidev
     END FUNCTION MYPOIDEV
  END INTERFACE

  INTERFACE
     SUBROUTINE SETUP_MODELS()
     END SUBROUTINE SETUP_MODELS
  END INTERFACE

END MODULE PIXCMD_UTILS

