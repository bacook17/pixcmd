# setup_files/TAU_mock.py
# Ben Cook (bcook@cfa.harvard.edu)

###############################################
# SETUP FILE for TAU Mock Test
# Mock data is a tau=1Gyr SFH

import pcmdpy.instrument as ins
import pcmdpy.isochrones as iso
import pcmdpy.galaxy as gal
import pcmdpy.driver as driver
import pcmdpy.utils as utils
import pcmdpy.gpu_utils as gpu_utils

from emcee.utils import sample_ball
import multiprocessing

import time

import numpy as np
import sys

###############################################
## IMPLEMENTATION SETTINGS

## Whether to use GPU acceleration
use_gpu = True

## The number of parallel processes to run.
## Using more threads than available CPUs (or GPUs, if gpu=True) will not improve performance
########## IMPORTANT NOTE:
##### Not currently implemented for N_threads > 1 and use_gpu = True, will fail
##### Hopefully this will be addressed soon
N_threads = 6

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

## Whether to add an additional likelihood term, a 2D gaussian fit of the data
add_total = True

## Cut out stars rarer than some limit (as fraction of total mass)
rare_cut = 0.

##### TIMING NOTE:
## The evaluation time of the fitting process will scale as:
## N_walkers * (N_burn + N_sample) / N_threads

## The number of emcee walkers
N_walkers = 256

## The number of burn-in iterations, per walker
N_burn = 50

## The number of sampling iterations, per walker
N_sample = 200

###############################################
## MODELLING SETTINGS

## The size (N_scale x N_scale) of the simulated image
N_scale = 1024

## The filters (photometry bands) to model
## There should be at least 2 filters.
###### Using more than 2 filters is currently not implemented
filters = np.array([ins.Filter.HST_F475W(1.0), ins.Filter.HST_F814W(1.0)])

## Initialize the isochrone models for the current set of filters
iso_model = iso.Isochrone_Model(filters)

## The galaxy class to use to model the data
#model_class = gal.Galaxy_SSP # simple stellar population (SSP)
model_class = gal.Galaxy_Model # 7-bin non-parametric SFH (FULL)

#### Initialize the emcee chains
# p0 = None #will initialize randomly over the prior space

## Initialize with a ball around a particular starting position
## for SSP mock model
#params_start = np.array([-0.2, -2., 2., 9.6])

## for TAU mock model
## tau=1Gyr, summing to Npix = 1e2
Npix = 1e2
tau = 1e9
age_edges = 10.**np.array([6., 7., 8., 8.5, 9.0, 9.5, 10., 10.2])
exp_factors = np.exp((age_edges - 10.**10.2) / tau)
logsfhs = np.log10(Npix * (exp_factors[1:] - exp_factors[:-1])) 
params_start = np.append(np.array([-0.2, -2]), logsfhs)

assert(len(params_start) == model_class._num_params)

## Initialize the ball with a particular width
std = 0.1 * np.ones_like(params_start)
p0 = sample_ball(params_start, std, size=N_walkers)

###############################################
## DATA SETTINGS

## Load observed data
#data_is_mock = False
#data_pcmd = ???????

## Mock tests
data_is_mock = True

## scale of mock image (N_mock x N_mock)
N_mock = 128

## model of the mock galaxy

## SSP model
#model_mock = gal.Galaxy_SSP
#params_mock = np.array([-0.2, -2., 2., 9.6])

## TAU model with Npix = 1e2, tau = 1 Gyr
model_mock = gal.Galaxy_Model
Npix = 1e2
tau = 1e9
age_edges = 10.**np.array([6., 7., 8., 8.5, 9.0, 9.5, 10., 10.2])
exp_factors = np.exp((age_edges - 10.**10.2) / tau)
exp_factors = np.exp(- age_edges / tau)
params_mock = np.append(np.array([-0.2, -2]), logsfhs)

galaxy_mock = model_mock(params_mock)

## Create the mock data
driv = driver.Driver(iso_model, gpu=use_gpu) #temporary driver to model the data
mags, _ = driv.simulate(galaxy_mock, N_mock, fixed_seed=fixed_seed)

## The mock data
data_pcmd = utils.make_pcmd(mags)


##############################################
## SAVE FILE SETTINGS

## Directory to save results to
results_dir = '/n/home01/bcook/pixcmd/scripts_py/results/'
## NAME OF THIS PARTICULAR RUN
name = "TAU_mock"
## the file to save the data
chain_file = results_dir + name + '.csv'
