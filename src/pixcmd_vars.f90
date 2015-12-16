MODULE PIXCMD_VARS

  ! module to set up common arrays and variables

  USE nrtype
  IMPLICIT NONE
  SAVE

  !-------------------common parameters-------------------!

  !flag for convolution (FFT=1, brute force=0)
  INTEGER, PARAMETER :: fft=1
 
  !variables for CMD Hess diagram
  INTEGER, PARAMETER :: nx=121,ny=221,npix=1024,nfil=2
  REAL(SP) :: xmin=-1.5,ymin=-6.0,dx=0.05,dy=0.05
  REAL(SP), DIMENSION(nx) :: xhess=0.
  REAL(SP), DIMENSION(ny) :: yhess=0.

  !data-specific parameters
  REAL(SP) :: dm=24.47  !M31 distance modulus
  !exposure times B,I
  REAL(SP), DIMENSION(nfil) :: exptime=(/1720.+1900.,1520.+1715./)
  !zero-points, B,I
  REAL(SP), DIMENSION(nfil) :: zpt=(/26.0593,25.9433/)
 
  !upper/lower priors
  REAL(SP), PARAMETER :: prlo=-9.0,prhi=0.75,wdth0=0.5

  !number of age and metallicity points in the model
  INTEGER, PARAMETER :: nage=15,nz=5,nm=1
  !parametres defining the age and mpix arrays
  REAL(SP) :: dage=0.2,age0=7.4,mpix0=2.0,dmpix=0.2
  REAL(SP), DIMENSION(nage) :: agesarr=0.
  REAL(SP), DIMENSION(nm)   :: mpixarr=0.
  REAL(SP), DIMENSION(nz)   :: zmetarr=0.

  !number of free parameters
  INTEGER, PARAMETER :: npar=nage*nz

  !max size of array for data and isochrones
  INTEGER, PARAMETER :: ndat_max=3000000,niso_max=5000

  !PSF array
  INTEGER, PARAMETER :: psf_step=4
  INTEGER, PARAMETER :: npsf=59,npsf2=npix/psf_step
  REAL(SP), DIMENSION(npsf,npsf) :: psfi=0.
  REAL(SP), DIMENSION(npsf2,npsf2,psf_step,psf_step) :: psf=0.

  !define small and large numbers
  REAL(SP), PARAMETER :: huge_number = 1E33
  REAL(SP), PARAMETER :: tiny_number = 1E-33

  CHARACTER(250) :: PIXCMD_HOME=''

  !---------------------common arrays---------------------!

  !array for model grids
  REAL(SP), DIMENSION(npix,npix,nfil,nz,nage) :: model=0.
  !array for the data
  REAL(SP), DIMENSION(nx,ny) :: hess_data=0.,hess_err=0.


END MODULE PIXCMD_VARS
