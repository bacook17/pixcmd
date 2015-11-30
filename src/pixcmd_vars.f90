MODULE PIXCMD_VARS

  ! module to set up common arrays and variables

  USE nrtype
  IMPLICIT NONE
  SAVE

  !-------------------common parameters-------------------!

  !variables for CMD Hess diagram
  INTEGER, PARAMETER :: nx=121,ny=221,npix=1000
  REAL(SP) :: xmin=-1.5,ymin=-6.0,dx=0.05,dy=0.05
  REAL(SP), DIMENSION(nx) :: xhess=0.
  REAL(SP), DIMENSION(ny) :: yhess=0.

  !number of age and metallicity points in the model
  INTEGER, PARAMETER :: nage=22,nz=6,nm=30
  REAL(SP) :: dage=0.2,age0=6.0,mpix0=1.0,dmpix=0.1
  REAL(SP), DIMENSION(nage) :: agesarr=0.
  REAL(SP), DIMENSION(nm)   :: mpixarr=0.
  REAL(SP), DIMENSION(nz)   :: zmetarr=0.

  !number of free parameters
  INTEGER, PARAMETER :: npar=nage*nz+1

  !max size of array for data and isochrones
  INTEGER, PARAMETER :: ndat_max=3000000,niso_max=5000

  !number of filters, dimension of PSF array
  INTEGER, PARAMETER :: nfil=2,npsf=59

  !define small and large numbers
  REAL(SP), PARAMETER :: huge_number = 1E33
  REAL(SP), PARAMETER :: tiny_number = 1E-33

  CHARACTER(250) :: PIXCMD_HOME=''

  !---------------------common arrays---------------------!

  !array for model grids
  REAL(SP), DIMENSION(nm,nz,nage,nx,ny) :: model=0.
  !array for the data
  REAL(SP), DIMENSION(nx,ny)    :: hess_data=0., hess_err=0.


END MODULE PIXCMD_VARS
