# M87 Model 1a (NonParam, Free, Single)
# Ben Cook (bcook@cfa.harvard.edu)

###############################################
# CONFIG FILE for m87 Model 1a
# MODEL Galaxy:
#    Single FeH
#    Single Dust
#    NonParam SFH
#    Distance Free

import pcmdpy_gpu as ppy
import multiprocessing

import time

import numpy as np

###############################################
# IMPLEMENTATION SETTINGS

params = {}  # arguments passed to pcmdpy_integrate
sampler_params = {}  # arguments passed to dynesty sampler initialization
run_params = {}  # arguments passed to sampler's run_nested()

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
                                    initializer=ppy.gpu_utils.initialize_gpu)
        time.sleep(10)
    else:
        pool = multiprocessing.Pool(processes=N_threads)
sampler_params['pool'] = pool

# Initialize the GPU with pycuda
# if params['use_gpu']:
#     ppy.gpu_utils.initialize_gpu(n=0)

# Check to see if GPU is available and properly initialized. If not, exit
# if params['use_gpu']:
#     assert ppy.gpu_utils._GPU_AVAIL, ('GPU NOT AVAILABLE, SEE ERROR LOGS. ',
#                                       'QUITTING')
#     assert ppy.gpu_utils._CUDAC_AVAIL, ('CUDAC COMPILATION FAILED, SEE ERROR ',
#                                         'LOGS. QUITTING')

###############################################
# DYNESTY SAMPLER SETTINGS
# These parameters are passed to initialization of
# Dynesty Sampler object

# Whether to use dynamic nested sampling
params['dynamic'] = DYNAMIC = False

# The number of dynesty live points
_nlive = 500
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
sampler_params['update_interval'] = 100

# Compute multiple realizations of bounding objects
sampler_params['bootstrap'] = 0

# Enlarge volume of bounding ellipsoids
sampler_params['enlarge'] = 1.1

# When should sampler update bounding from unit-cube
sampler_params['first_update'] = {'min_eff': 30.}

###############################################
# DYNESTY RUN_NESTED SETTINGS

# The number of max calls for dynesty
run_params['maxcall'] = 500000

# The error tolerance for dynesty stopping criterion
_dlogz = 0.5
if DYNAMIC:
    run_params['dlogz_init'] = _dlogz
else:
    run_params['dlogz'] = _dlogz
    sampler_params['add_live'] = True

if DYNAMIC:
    # How many batches?
    run_params['maxbatch'] = 10
    # How many live points per batch?
    run_params['nlive_batch'] = 100
    # weight function parameters
    run_params['wt_kwargs'] = {'pfrac': 1.0}
    # How many max calls per iteration?
    run_params['maxcall_per_iter'] = 100000
    # Don't keep boundaries
    run_params['save_bounds'] = False

###############################################
# PCMD MODELLING SETTINGS

# The size (Nim x Nim) of the simulated image
params['Nim'] = 512

# The filters (photometry bands) to model. There should be at least 2 filters.
# Default choice: F814W and F475W
params['filters'] = ppy.instrument.default_m87_filters()

# Alternative choice: F814W, F555W, and F435W
# params['filters'] = ppy.instrument.default_m51_filters()

# To manually set options:
# filters = []
# filters.append(ppy.instrument.ACS_WFC_F814W(exposure=8160., psf=....))
# filters.append(ppy.instrument.ACS_WFC_F475W(exposure=3120., psf=....))

# Initialize the isochrone models for the current set of filters
params['iso_model'] = ppy.isochrones.Isochrone_Model(params['filters'])

# Set a custom Galaxy Model with four parts

# Metallicity model
metalmodel = ppy.metalmodels.SingleFeH()  # Single Metallicity
# metalmodel = ppy.metalmodels.NormMDF()  # Gaussian MDF
# metalmodel = ppy.metalmodels.FixedWidthNormMDF(0.2)  # fixed width MDF

# Dust model
dustmodel = ppy.dustmodels.SingleDust()  # single dust screen
# dustmodel = ppy.dustmodels.LogNormDust()  # lognormal screen
# dustmodel = ppy.dustmodels.FixedWidthLogNormDust(0.3)  # fixed width lognorm

# Age model
sfhmodel = ppy.sfhmodels.NonParam()  # Fully non-parametric model
# sfhmodel = ppy.sfhmodels.ConstantSFR()  # constant Star Formation Rate
# sfhmodel = ppy.sfhmodels.TauModel()  # exponential SFR decline
# sfhmodel = ppy.sfhmodels.RisingTau()  # Linear x exponential decline
# sfhmodel = ppy.sfhmodels.SSPModel()  # single age SSP

# Distance model
distancemodel = ppy.distancemodels.FixedDistance(31.0)  # fixed dmod=31
# distancemodel = ppy.distancemodels.VariableDistance()  # dmod floats
params['gal_model'] = ppy.galaxy.CustomGalaxy(metalmodel, dustmodel, sfhmodel,
                                              distancemodel)

# Add the binned hess values and the mean magnitude and color terms
params['like_mode'] = 2

# The hess bins to compute the likelihood in
# The magnitude upper/lower bounds are very important to consider
# relative to distance
magbins = np.arange(10, 45, 0.05)
colorbins = np.arange(-1.5, 5.6, 0.05)  # fairly insensitive to distance
params['bins'] = [magbins, colorbins]

# Factor to downsample the isochrones
params['downsample'] = 5

# which magnitude system
params['mag_system'] = 'vega'

# Cut out stars brighter than some limit (of mean luminosity)
params['lum_cut'] = np.inf

# Whether to use a fixed random-number seed
# (decreases stochasticity of likelihood calls)
params['fixed_seed'] = True

# Average counts of "sky noise" to add in each band
params['sky_noise'] = [292.8, 377.1]

params['shot_noise'] = True

params['dither'] = False

###############################################
# PRIOR SETTINGS

# The bounds on the flat prior for each parameter
z_bound = [-1.0, 0.25]  # metallicity
dust_med_bound = [-2.0, -1.0]  # log dust median
# Only set the distance bounds if allowed to float
dmod_bound = None
# dmod_bound = [[30., 32.]]

# Compute the 7-param SFH bound using tau models to bound
# Npix_bound = [2., 5.]
# tau_bound = [0.1, 5.]
logNpix = 4.4
dN = 1.0
tau = 2.0

model = ppy.sfhmodels.TauModel(iso_step=-1)
model.set_params([logNpix - dN, tau])
lower_sfh = np.log10(model.SFH)

model.set_params([logNpix + dN, tau])
upper_sfh = np.log10(model.SFH)

SFH_bounds_arr = np.array([lower_sfh, upper_sfh]).T
SFH_bounds = list(list(bound) for bound in SFH_bounds_arr)

# Create a Prior object with given bounds
prior_bounds = {}
prior_bounds['feh_bounds'] = [z_bound]
prior_bounds['dust_bounds'] = [dust_med_bound]
prior_bounds['age_bounds'] = SFH_bounds
prior_bounds['dmod_bounds'] = dmod_bound

params['prior'] = params['gal_model'].get_flat_prior(**prior_bounds)

###############################################
# DATA / MOCK SETTINGS

# Is the data created manually, or should it be read from a file?
params['data_is_mock'] = False
