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

    #arguments
    args = {}
    args['pcmd'] = setup.data_pcmd
    args['filters'] = setup.filters
    args['im_scale'] = setup.N_scale
    args['N_walkers'] = setup.N_walkers
    args['N_burn'] = setup.N_burn
    args['N_sample'] = setup.N_sample

    #optional key-word arguments (defaults are set by fit_model.sample_post)
    args['pool'] = setup.pool 
    args['gpu'] = setup.use_gpu
    args['fixed_seed'] = setup.fixed_seed
    args['gal_class'] = setup.model_class
    args['p0'] = setup.p0
    try:
        args['like_mode'] = setup.like_mode
    except:
        pass
    try:
        args['rare_cut'] = setup.rare_cut
    except:
        pass
    #if setup.data_is_mock:
    #    args['galaxy_mock'] = setup.galaxy_mock
    #    args['N_mock'] = setup.N_mock

    print('Running emcee')
    sampler = fit_model.sample_post(**args)

    print('emcee complete, saving results')
    #to use in saving results
    param_names = setup.model_class._param_names
    N_params = len(param_names)
    chain_file = setup.chain_file

    #Save results of the chain
    chain_df = pd.DataFrame()
    for d in range(N_params):
        chain_df[param_names[d]] = sampler.flatchain[:,d]
    chain_df['lnprob'] = sampler.flatlnprobability
    chain_df['walker'] = np.repeat(np.arange(N_walkers), N_sample)
    chain_df['accept_frac'] = np.repeat(sampler.acceptance_fraction, N_sample)

    chain_df.to_csv(chain_file, index=False, float_format='%.4e', compression='gzip')
