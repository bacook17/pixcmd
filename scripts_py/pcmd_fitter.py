# pcmd_fitter.py
# Ben Cook (bcook@cfa.harvard.edu)

import numpy as np
from pcmdpy import fit_model

import pandas as pd
import sys

import imp

if __name__ == "__main__":

    #Setup from an external file
    setup_file = sys.argv[1]
    setup_mod = sys.argv[1].strip('.py').rpartition('/')[-1]
    print('Loading Setup File: %s'%setup_file)
    setup = imp.load_source(setup_mod, setup_file)

    N_scale = setup.N_scale
    N_walkers = setup.N_walkers
    N_burn = setup.N_burn
    N_sample = setup.N_sample
    N_threads = setup.N_threads
    pool = setup.pool 
    gpu = setup.use_gpu

    fixed_seed = setup.fixed_seed

    filters = setup.filters
    iso_model = setup.iso_model

    data_is_mock = setup.data_is_mock
    data_pcmd = setup.data_pcmd
    if data_is_mock:
        galaxy_mock = setup.galaxy_mock
        N_mock = setup.N_mock

    model_class = setup.model_class
    param_names = model_class._param_names
    N_params = len(param_names)
    p0 = setup.p0
    
    chain_file = setup.chain_file

    print('Running emcee')
    sampler = fit_model.sample_post(data_pcmd, filters, N_scale, N_walkers, N_burn, N_sample, fixed_seed=fixed_seed,
                                    gal_class=model_class, gpu=gpu, pool=pool, p0=p0)

    print('emcee complete, saving results')
    #Save results of the chain
    chain_df = pd.DataFrame()
    for d in range(N_params):
        chain_df[param_names[d]] = sampler.flatchain[:,d]
    chain_df['lnprob'] = sampler.flatlnprobability
    chain_df['walker'] = np.repeat(np.arange(N_walkers), N_sample)
    chain_df['accept_frac'] = np.repeat(sampler.acceptance_fraction, N_sample)

    chain_df.to_csv(chain_file, index=False, float_format='%.4f', compression='gzip')
