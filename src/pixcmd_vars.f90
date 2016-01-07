MODULE PIXCMD_VARS

  ! module to set up common arrays and variables

  USE nrtype
  IMPLICIT NONE
  SAVE

  !-------------------common parameters-------------------!

  CHARACTER(10) :: iso_tag='_x5FEWER'

  !flag for convolution (FFT=1, brute force=0)
  INTEGER, PARAMETER :: fft=1
 
  !variables for model image and CMD Hess diagram
  INTEGER, PARAMETER  :: nx=121,ny=301,npix=512,nfil=2
  REAL(SP), PARAMETER :: xmin=-1.5,ymin=-10.0,dx=0.05,dy=0.05
  REAL(SP), DIMENSION(nx) :: xhess=0.
  REAL(SP), DIMENSION(ny) :: yhess=0.

  !data-specific parameters
  REAL(SP) :: dm=24.47  !M31 distance modulus
  !exposure times: F475W, F814W
  REAL(SP), DIMENSION(nfil), PARAMETER :: exptime=(/3620.,3235./)
  !zero-points: F475W, F814W
  REAL(SP), DIMENSION(nfil), PARAMETER :: zpt=(/26.0593,25.9433/)
  !Reddening values from Schalfly & Finkbeiner 2011 for F475W, F814W
  REAL(SP), DIMENSION(nfil), PARAMETER :: red_per_ebv = (/3.268,1.526/)
 
  !upper/lower priors
  REAL(SP), PARAMETER :: prlo=-7.0,prhi=0.3,wdth0=0.5
  REAL(SP), PARAMETER :: prlo_m=0.5,prhi_m=7.0
  REAL(SP), PARAMETER :: prlo_lebv=-7.0,prhi_lebv=0.5

  !stellar mass below which the IMF is assumed to be fully populated
  REAL(SP), PARAMETER :: minmass=0.8

  !number of age and metallicity points in the model
  INTEGER, PARAMETER :: nage=22,nz=1,nzskip=5
  !parameters defining the age array
  REAL(SP) :: dage=0.2,age0=6.0
  REAL(SP), DIMENSION(nage) :: agesarr=0.
  REAL(SP), DIMENSION(nz)   :: zmetarr=0.

  !number of free parameters 
  INTEGER, PARAMETER :: nxpar = 1 !mpix
  INTEGER, PARAMETER :: npar=nage*nz+nxpar

  !max size of array for data and isochrones
  INTEGER, PARAMETER :: niso_max=5000

  !background flux level (less than 1 dM per pixel)
  REAL(SP), PARAMETER :: bkgnd = 1E-10

  !PSF array
  INTEGER, PARAMETER :: psf_step=4
  INTEGER, PARAMETER :: npsf=59,npsf2=npix/psf_step
  REAL(SP), DIMENSION(npsf,npsf) :: psfi=0.
  REAL(SP), DIMENSION(npsf2,npsf2,psf_step,psf_step) :: psf=0.

  !define small and large numbers
  REAL(SP), PARAMETER :: huge_number = 1E33
  REAL(SP), PARAMETER :: tiny_number = 1E-33

  CHARACTER(250) :: PIXCMD_HOME=''

  !define the random number array
  INTEGER, PARAMETER :: nran=npix*npix*256
  REAL(SP), DIMENSION(nran) :: ranarr


  !---------------------common arrays---------------------!

  !array for the data
  REAL(SP), DIMENSION(nx,ny) :: hess_data=0.,hess_err=0.

  !type structure for the isochrones
  TYPE TISO
     !band order is BVIJH
     REAL(SP), DIMENSION(nfil) :: bands=0.0
     REAL(SP) :: imf=-99.,mass=0.0,age=0.0
  END TYPE TISO

  !array holding the actual isochrones
  TYPE(TISO), DIMENSION(niso_max) :: iso

  !actual number of isochrone points
  INTEGER :: niso
  !location of age breaks in the isochrone file
  INTEGER, DIMENSION(nage+1) :: ageind

END MODULE PIXCMD_VARS
