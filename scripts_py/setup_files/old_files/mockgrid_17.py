# mockgrid_base_const.py configuration file
# Ben Cook (bcook@cfa.harvard.edu)

###############################################
# CONFIG FILE for mock test
# MOCK Galaxy: has constant SFH, single FeH, Npix=1e4
# MODEL Galaxy: has NonParam SFH, single MDF

from pcmdpy import instrument, galaxy, gpu_utils, priors
from pcmdpy import agemodels, dustmodels, metalmodels
# only required for creating a mock
from pcmdpy import driver, utils
from pcmdpy.isochrones import Isochrone_Model
import multiprocessing

import time

import numpy as np
import sys

###############################################
# IMPLEMENTATION SETTINGS

params = {}
sampler_params = {}
run_params = {}

# Whether to use GPU acceleration
params['use_gpu'] = True

# Whether to output progress steps
params['verbose'] = True

# The number of parallel processes to run. Using more threads than available
# CPUs (or GPUs, if gpu=True) will not improve performance
N_threads = 1

# Setup the multiprocessing pool, for parallel evaluation
pool = None
if N_threads > 1:
    if params['use_gpu']:
        pool = multiprocessing.Pool(processes=N_threads,
                                    initializer=gpu_utils.initialize_gpu)
        time.sleep(10)
    else:
        pool = multiprocessing.Pool(processes=N_threads)
sampler_params['pool'] = pool

# Initialize the GPU with pycuda
if params['use_gpu']:
    gpu_utils.initialize_gpu(n=0)

# Check to see if GPU is available. If not, exit
if params['use_gpu']:
    if not gpu_utils._GPU_AVAIL:
        print('GPU NOT AVAILABLE, SEE ERROR LOGS. QUITTING')
        sys.exit(2)

# Whether to require CUDAC (fasetest) implementation
params['use_cudac'] = True

# Check to see if CUDAC is available. If not, exit
if params['use_cudac']:
    if not gpu_utils._CUDAC_AVAIL:
        print('CUDAC NOT AVAILABLE, SEE ERROR LOGS. QUITTING')
        sys.exit(2)

###############################################
# DYNESTY SAMPLER SETTINGS
# These parameters are passed to initialization of
# Dynesty Sampler object

# Whether to use dynamic nested sampling
params['dynamic'] = False
DYNAMIC = params['dynamic']

# The number of dynesty live points
_nlive = 300
if DYNAMIC:
    run_params['nlive_init'] = _nlive
else:
    sampler_params['nlive'] = _nlive

# How to bound the prior
sampler_params['bound'] = 'multi'

# How to sample within the prior bounds
sampler_params['method'] = 'unif'

# Number of parallel processes
sampler_params['nprocs'] = N_threads

# Only update the bounding distribution after this many calls
sampler_params['update_interval'] = 1

# Compute multiple realizations of bounding objects
sampler_params['bootstrap'] = 0

# Enlarge volume of bounding ellipsoids
sampler_params['enlarge'] = 1.1

# When should sampler update bounding from unit-cube
sampler_params['first_update'] = {'min_eff': 30.}

###############################################
# DYNESTY RUN_NESTED SETTINGS

# The number of max calls for dynesty
run_params['maxcall'] = 120000

# The error tolerance for dynesty stopping criterion
_dlogz = 0.5
if DYNAMIC:
    run_params['dlogz_init'] = _dlogz
else:
    run_params['dlogz'] = _dlogz

if DYNAMIC:
    # How many batches?
    run_params['maxbatch'] = 0
    # How many live points per batch?
    run_params['nlive_batch'] = 0
    # weight function parameters
    run_params['wt_kwargs'] = {'pfrac': 1.0}
    # How many max calls per iteration?
    run_params['maxcall_per_it'] = 1000

###############################################
# PCMD MODELLING SETTINGS

# The size (N_im x N_im) of the simulated image
params['N_im'] = 1024

# The filters (photometry bands) to model
# There should be at least 2 filters.

######
# Using more than 2 filters is currently not implemented
dmod = 24.47    # distance modulus to M31
d_mpc = 10.**((dmod - 25.)/5.)   # about 0.78
params['filters'] = instrument.m31_filters(dist=d_mpc)

# Initialize the isochrone models for the current set of filters
params['iso_model'] = Isochrone_Model(params['filters'])

# The galaxy class to use to model the data
metals = metalmodels.SingleFeH  # Single Metallicity
dust = dustmodels.SingleDust  # single dust screen
age = agemodels.NonParam  # tau SFH
params['gal_class'] = galaxy.CustomGalaxy(metals, dust, age)

# Add the binned hess values and the mean magnitude and color terms
params['like_mode'] = 2

# Cut out stars rarer than some limit (as fraction of total mass)
params['lum_cut'] = np.inf

# Whether to use a fixed random-number seed
# (decreases stochasticity of likelihood calls)
params['fixed_seed'] = True

###############################################
# PRIOR SETTINGS

z_bound = [-1.5, 0.5]
dust_med_bound = [-2.0, 0.5]

lower_sfh = np.log10(agemodels.ConstantSFR(np.array([2.5]), iso_step=-1).SFH)
upper_sfh = np.log10(agemodels.ConstantSFR(np.array([5.0]), iso_step=-1).SFH)
SFH_bounds_arr = np.array([lower_sfh, upper_sfh]).T
SFH_bounds = list(list(bound) for bound in SFH_bounds_arr)

prior_bounds = {}
prior_bounds['feh_bounds'] = [z_bound]
prior_bounds['dust_bounds'] = [dust_med_bound]
prior_bounds['age_bounds'] = SFH_bounds

params['prior'] = params['gal_class'].get_flat_prior(**prior_bounds)

###############################################
# DATA / MOCK SETTINGS

params['data_is_mock'] = True

# scale of mock image (N_mock x N_mock)
N_mock = 256

# model of the mock galaxy
metals = metalmodels.SingleFeH  # single metallicity
dust = dustmodels.SingleDust  # single dust screen
age = agemodels.ConstantSFR  # Constant SFR
model_mock = galaxy.CustomGalaxy(metals, dust, age)

# Tau model with [Fe/H]=-0.2, log E(B-V) = -.5
# Npix = 1e4 constant SFR
gal_params = np.array([-0.2, -0.5, 4.])
galaxy_mock = model_mock.get_model(gal_params)

# Create the mock data
# temporary driver to make mock
driv = driver.Driver(params['iso_model'], gpu=params['use_gpu'])
# The mock data
params['data_pcmd'], _ = driv.simulate(galaxy_mock, N_mock,
                                       fixed_seed=params['fixed_seed'])

del driv
