# Mock configuration file for Paper 1
# Ben Cook (bcook@cfa.harvard.edu)

###############################################
# CONFIG FILE for mock run #49
# MOCK Galaxy:
#    Metallicity Model: Fixed-Width (0.2) MDF
#            [Fe/H] = -0.25
#    Dust Model:        Fixed-Width (0.1) LogNormal (Fdust = 0.5)
#        log E(B-V) = -0.5
#    SFH Model: Tau
#              Npix = 6.0
#              tau  = 3.0
#    Distance
#              dmod = 26.0
#
# MODEL Galaxy: Matches input model
#       WITH Non-Param SFH
# Priors:
#           [Fe/H] : [-0.5, 0.25]
#       log E(B-V) : [-1.0, 0.0]
#            SFH_i : tau_model +/- 1 dex 
#         distance : [25.0, 30.0]

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
params['dynamic'] = DYNAMIC = True

# The number of dynesty live points
_nlive = 1000
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
run_params['maxcall'] = 200000
params['continue_run'] = False

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
    run_params['maxcall_per_iter'] = 10000
    # Don't keep boundaries
    run_params['save_bounds'] = False

###############################################
# PCMD MODELLING SETTINGS

# The size (Nim x Nim) of the simulated image
params['Nim'] = 512

# The filters (photometry bands) to model. There should be at least 2 filters.
# Default choice: F814W and F475W
params['filters'] = ppy.instrument.default_m31_filters()

# Initialize the isochrone models for the current set of filters
params['iso_model'] = ppy.isochrones.Isochrone_Model(params['filters'])

# Set a custom Galaxy Model with four parts

# Metallicity model
# metalmodel = ppy.metalmodels.SingleFeH()  # Single Metallicity
# metalmodel = ppy.metalmodels.NormMDF()  # Gaussian MDF
metalmodel = ppy.metalmodels.FixedWidthNormMDF(0.2)  # fixed width MDF

# Dust model
# dustmodel = ppy.dustmodels.SingleDust()  # single dust screen
# dustmodel = ppy.dustmodels.LogNormDust()  # lognormal screen
dustmodel = ppy.dustmodels.FixedWidthLogNormDust(0.1)  # fixed width lognorm

# Age model
sfhmodel = ppy.sfhmodels.NonParam()  # Fully non-parametric model
# sfhmodel = ppy.sfhmodels.ConstantSFR()  # constant Star Formation Rate
# sfhmodel = ppy.sfhmodels.TauModel()  # exponential SFR decline
# sfhmodel = ppy.sfhmodels.RisingTau()  # Linear x exponential decline
# sfhmodel = ppy.sfhmodels.SSPModel()  # single age SSP

# Distance model
# distancemodel = ppy.distancemodels.FixedDistance(26.0)  # fixed dmod=26.0 (1.6 Mpc)
distancemodel = ppy.distancemodels.VariableDistance()  # dmod floats

params['gal_model'] = ppy.galaxy.CustomGalaxy(
    metalmodel,
    dustmodel,
    sfhmodel,
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
params['sky_noise'] = None

params['shot_noise'] = True

###############################################
# PRIOR SETTINGS

# The bounds on the flat prior for each parameter
z_bound = [-0.5, 0.25]  # metallicity

dust_med_bound = [-1.0, 0.0]  # log dust

# Only set the distance bounds if allowed to float
# dmod_bound = None
dmod_bound = [[25.0, 30.0]]

# Compute the 7-param SFH bound using tau models to bound
Npix_low, tau = 5.0, 3.0
model = ppy.sfhmodels.TauModel(iso_step=-1)
model.set_params([Npix_low, tau])
lower_sfh = np.log10(model.SFH)

Npix_high = 7.0
model.set_params([Npix_high, tau])
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
params['data_is_mock'] = True

# scale of mock image (N_mock x N_mock)
N_mock = 256

# model of the mock galaxy
feh = -0.25
log_ebv = -0.5
log_npix = 6.0
tau = 3.0
dmod = 26.0

# Mock data is generated with same model as is fit (except Tau Model)
metalmodel = metalmodel
dustmodel = dustmodel
sfhmodel = ppy.sfhmodels.TauModel()
distancemodel = ppy.distancemodels.VariableDistance()  # dmod floats
model_mock = ppy.galaxy.CustomGalaxy(
    metalmodel,
    dustmodel,
    sfhmodel,
    distancemodel)

gal_params = np.array([feh, log_ebv, log_npix, tau, dmod])
model_mock.set_params(gal_params)

# Create the mock data
# temporary driver to make mock
driv = ppy.driver.Driver(params['iso_model'], gpu=True)
# The mock data
params['data_pcmd'], _ = driv.simulate(model_mock, N_mock,
                                       fixed_seed=params['fixed_seed'],
                                       shot_noise=params['shot_noise'],
                                       sky_noise=params['sky_noise'],
                                       downsample=params['downsample'],
                                       mag_system=params['mag_system'])

del driv
