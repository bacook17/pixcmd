MODULE PIXCMD_VARS

  ! module to set up common arrays and variables

  USE nrtype
  IMPLICIT NONE
  SAVE

  ! If you want to make sure that the results are insensitive
  ! to the various approximations that I have made, do the following:
  ! 1) set true_poisson=1 
  ! 2) set minmass=0.0
  ! 3) set minnum=huge_number
  ! 4) set iso_tag=''
  ! 5) set npix=2048 (or larger!)
  ! 6) set nran=1 (otherwise you'll run out of memory)
  !
  ! Be prepared for the code to run MUCH MUCH slower!

  !-------------------common parameters-------------------!

  !number of pixels for the model image.  computation time
  !scales as npix^2
  INTEGER, PARAMETER :: npix=256

  !isochrone flag.  default to using 5x fewer points
  CHARACTER(10) :: iso_tag='_x5FEWER'

  !flag for convolution (FFT=1, brute force=0)
  INTEGER, PARAMETER :: fft=1
  !if set, draw random numbers for each Poisson draw
  !if not set, use the master random number array
  INTEGER, PARAMETER :: true_poisson=0
  !if set, include observational errors
  INTEGER, PARAMETER :: incl_obs_err=1
  !if >0, then set the random seed to fix_seed
  INTEGER :: fix_seed=3300

  !variables for model image and CMD Hess diagram
  INTEGER, PARAMETER  :: nx=121,ny=351,nfil=2
  REAL(SP), PARAMETER :: xmin=-1.5,ymin=-12.0,dx=0.05,dy=0.05

  !data-specific parameters
  REAL(SP), PARAMETER :: dm=24.47  !M31 distance modulus
  !exposure times: F475W, F814W
  REAL(SP), DIMENSION(nfil), PARAMETER :: exptime=(/3620.,3235./)
  !zero-points: F475W, F814W
  REAL(SP), DIMENSION(nfil), PARAMETER :: zpt=(/26.0593,25.9433/)
  !Reddening values from Schalfly & Finkbeiner 2011 for F475W, F814W
  !should double check that this is the correct implementation
  REAL(SP), DIMENSION(nfil), PARAMETER :: red_per_ebv = (/3.268,1.526/)
 
  !upper/lower priors
  REAL(SP), PARAMETER :: prlo_sfh=-10.0,prhi_sfh=0.5,wdth0=1E-3
  REAL(SP), PARAMETER :: prlo_lebv=-5.0,prhi_lebv=0.0
  REAL(SP), PARAMETER :: prlo_zmet=-5.0,prhi_zmet=0.0
  
  !stellar mass below which the IMF is assumed to be fully populated
  REAL(SP), PARAMETER :: minmass=0.8
  !number below which we sample an isochrone point via Poisson/Gaussian,
  !above which we assume the pixel-to-pixel variance is 0.0
  REAL(SP), PARAMETER :: minnum=500.
  !N below which we use Poisson, above which we approx with a Gaussian
  REAL(SP), PARAMETER :: maxpoidev=40.
  
  !number of age and metallicity points in the model
  INTEGER, PARAMETER :: nage=7,nz=5,nzskip=3,nzi=4
  !parameters defining the age array
  REAL(SP) :: dage=0.5,age0=6.0
  REAL(SP), DIMENSION(nage)   :: agesarr=0.
  REAL(SP), DIMENSION(nage+1) :: agesarr2=0.
  REAL(SP), DIMENSION(nz)     :: zmetarr=0.

  !number of free parameters 
  INTEGER, PARAMETER :: nxpar = 1 !lebv
  INTEGER, PARAMETER :: npar=nage+nxpar+nz-1

  !max size of array for data and isochrones
  INTEGER, PARAMETER :: niso_max=4096

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
  INTEGER, PARAMETER :: nran=256*256
  INTEGER :: kran=1
  REAL(SP), DIMENSION(nran,niso_max) :: ranarr
  REAL(SP), DIMENSION(npix,npix) :: gdev,gdev2

  !---------------------common arrays---------------------!

  !array for the x and y axes of the Hess diagram
  REAL(SP), DIMENSION(nx) :: xhess=0.
  REAL(SP), DIMENSION(ny) :: yhess=0.

  !array holding the final priors
  REAL(SP), DIMENSION(npar) :: prlo=-99.,prhi=99.

  !array for the data
  REAL(SP), DIMENSION(nx,ny) :: hess_data=0.,hess_err=0.

  !first guess for Mpix, log(EBV)
  REAL(SP) :: mpix0=2.0,lebv0=-2.0

  !type structure for the isochrones
  TYPE TISO
     !band order is BI
     REAL(SP), DIMENSION(nfil) :: bands=0.0
     REAL(SP) :: imf=-99.,mass=0.0,age=0.0
     INTEGER  :: aind=1
  END TYPE TISO

  !array holding the actual isochrones
  TYPE(TISO), DIMENSION(nz,niso_max) :: iso

  !actual number of isochrone points
  INTEGER, DIMENSION(nz) :: niso=1
  !number of isochrones per age bin
  INTEGER, DIMENSION(nage) :: nisoage=1

END MODULE PIXCMD_VARS
