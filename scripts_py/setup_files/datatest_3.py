# setup_files/datatest_bulge.py
# Ben Cook (bcook@cfa.harvard.edu)

###############################################
# SETUP FILE for Data Test on Disk with tau model

import pcmdpy.instrument as ins
import pcmdpy.isochrones as iso
import pcmdpy.galaxy as gal
import pcmdpy.driver as driver
import pcmdpy.utils as utils
import pcmdpy.gpu_utils as gpu_utils
from astropy.io import ascii
import multiprocessing

import time

import numpy as np
import sys

###############################################
## IMPLEMENTATION SETTINGS

## Whether to use GPU acceleration
use_gpu = True

## Whether to output progress steps
verbose = True

## The number of parallel processes to run.
## Using more threads than available CPUs (or GPUs, if gpu=True) will not improve performance
N_threads = 1

## Setup the multiprocessing pool, for parallel evaluation
pool = None
if N_threads > 1:
    if use_gpu:
        pool = multiprocessing.Pool(processes=N_threads, initializer=gpu_utils.initialize_gpu)
        time.sleep(10)
    else:
        pool = multiprocessing.Pool(processes=N_threads)

if use_gpu:
    gpu_utils.initialize_gpu(n=0)

## Check to see if GPU is available. If not, exit
if use_gpu:
    if not gpu_utils._GPU_AVAIL:
        print('GPU NOT AVAILABLE, SEE ERROR LOGS. QUITTING')
        sys.exit(2)

## Whether to require CUDAC (fasetest) implementation
use_cudac = True

## Check to see if CUDAC is available. If not, exit
if use_cudac:
    if not gpu_utils._CUDAC_AVAIL:
        print('CUDAC NOT AVAILABLE, SEE ERROR LOGS. QUITTING')
        sys.exit(2)

## Whether to use a fixed random-number seed (decreases stochasticity of likelihood calls)
fixed_seed = True

##Add the binned hess values and the mean magnitude and color terms
like_mode = 2

## Cut out stars rarer than some limit (as fraction of total mass)
lum_cut = np.inf

##### TIMING NOTE:
## The evaluation time of the fitting process will scale as:
## N_walkers * (N_burn + N_sample) / N_threads

## Whether to use dynesty (otherwise, use Nestle)
use_dynesty = True
dynamic = False

## The number of dynesty live points
N_points = 50

## The number of burn-in iterations, per walker
N_burn = 0

## The number of max calls for dynesty
N_max = 500000

## The error tolerance for Nestle stopping criterion
dlogz = 0.5

###############################################
## MODELLING SETTINGS

## The size (N_scale x N_scale) of the simulated image
N_scale = 1024

## The filters (photometry bands) to model
## There should be at least 2 filters.
###### Using more than 2 filters is currently not implemented
dmod = 24.47
d_mpc = 10.**((dmod - 25.)/5.) #about 0.78
filters = np.array([ins.Filter.HST_F475W(d_mpc), ins.Filter.HST_F814W(d_mpc)])

## Initialize the isochrone models for the current set of filters
iso_model = iso.Isochrone_Model(filters,
                                conversions=np.array([-0.0978757217, 0.4236631949]))

## The galaxy class to use to model the data
model_class = gal.Tau_Model # Tau model

#### Initialize the emcee chains
# p0 = None #will initialize randomly over the prior space

## Initialize with a ball around a particular starting position
## for SSP mock model
#params_start = np.array([-0.2, -2., 2., 9.6])

## Constrain the prior volume to a small region around the correct answer?
small_prior = True

## for FULL mock model
## constant SFH, summing to Npix = 1e2
#Npix = 1e2
#age_edges = np.array([6., 7., 8., 8.5, 9.0, 9.5, 10., 10.2])
#bin_widths = 10.**age_edges[1:] - 10.**age_edges[:-1]
#logsfhs = np.log10(Npix * bin_widths / np.sum(bin_widths)) 
#params_start = np.append(np.array([-0.2, -2]), logsfhs)

#assert(len(params_start) == model_class._num_params)

## Initialize the ball with a particular width
#std = 0.1 * np.ones_like(params_start)
#p0 = sample_ball(params_start, std, size=N_walkers)

###############################################
## DATA SETTINGS

## Load observed data
#data_is_mock = False
#data_pcmd = ???????

## Mock tests
data_is_mock = False

## scale of mock image (N_mock x N_mock)
#N_mock = 256

## model of the mock galaxy

## SSP model
#model_mock = gal.Galaxy_SSP
#params_mock = np.array([-0.2, -2., 2., 9.6])

## Tau model with Npix = 1e4, tau=5 Gyr
#model_mock = gal.Constant_SFR
#Npix = 1e2
#age_edges = np.array([6., 7., 8., 8.5, 9.0, 9.5, 10., 10.2])
#bin_widths = 10.**age_edges[1:] - 10.**age_edges[:-1]
#logsfhs = np.log10(Npix * bin_widths / np.sum(bin_widths)) 
#params_mock = np.append(np.array([-0.2, -2]), logsfhs)

params_center = np.array([-0.5, -1., 1.5, 4.]) 
params_width = np.array([1., 1.,  1., 3.])

def prior_trans(normed_params):
    #+/- 0.5 around correct answer
    return params_center + params_width*(-1. + 2.*normed_params)

def lnprior_func(params):
    z, log_dust, log_Npix, tau = params
    if (z < -2.) or (z > 0.5):
        return -np.inf
    if (log_dust < -3.) or (log_dust > 0.5):
        return -np.inf
    if (log_Npix < -1.) or (log_Npix > 8):
        return -np.inf
    if (tau < .1) or (tau > 20.):
        return -np.inf
    return 0.

## Create the mock data
#driv = driver.Driver(iso_model, gpu=use_gpu) #temporary driver to model the data
#mags, _ = driv.simulate(galaxy_mock, N_mock, fixed_seed=fixed_seed)

## The mock data
data_pcmd = np.loadtxt('/n/home01/bcook/pixcmd/data/m31_b06-263.dat', unpack=True)

##############################################
## SAVE FILE SETTINGS

## Directory to save results to
results_dir = '/n/home01/bcook/pixcmd/scripts_py/results/'
## NAME OF THIS PARTICULAR RUN
name = "datatest_3"
## the file to save the data
output_file = results_dir + name + '.csv'
