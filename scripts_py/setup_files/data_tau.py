# data_tau.py configuration file
# Ben Cook (bcook@cfa.harvard.edu)

###############################################
# CONFIG FILE for Data analysis with Tau model

from pcmdpy import instrument, galaxy, gpu_utils, priors
# only required for creating a mock
# from pcmdpy import driver, utils
from pcmdpy.isochrones import Isochrone_Model
import multiprocessing

import time

import numpy as np
import sys

###############################################
# IMPLEMENTATION SETTINGS

params = {}

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
params['pool'] = pool

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
# DYNESTY SETTINGS

# Whether to use dynamic nested sampling
params['dynamic'] = False

# The number of dynesty live points
params['nlive'] = 50

# The number of max calls for dynesty
params['maxcall'] = 200000

# The error tolerance for dynesty stopping criterion
params['dlogz'] = 0.5

# How many max calls per iteration?
params['maxcall_per_it'] = 1000

# How many batches??
params['maxbatch'] = 0

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
params['gal_class'] = galaxy.TauModel

# Add the binned hess values and the mean magnitude and color terms
params['like_mode'] = 2

# Cut out stars rarer than some limit (as fraction of total mass)
params['lum_cut'] = np.inf

# Whether to use a fixed random-number seed
# (decreases stochasticity of likelihood calls)
params['fixed_seed'] = True

###############################################
# PRIOR SETTINGS

prior_args = {}
prior_args['z_bound'] = [-1.5, 0.5]
prior_args['dust_bound'] = [-2.5, -0.5]
prior_args['npix_bound'] = [1., 4.]
prior_args['tau_bound'] = [1., 9.]

params['prior'] = priors.TauFlatPrior(**prior_args)

###############################################
# DATA / MOCK SETTINGS

params['data_is_mock'] = False
