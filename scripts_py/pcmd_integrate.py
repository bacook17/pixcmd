# pcmd_integrate.py
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
    args['N_points'] = setup.N_points

    #optional key-word arguments (defaults are set by fit_model.nested_integrate)
    args['max_call'] = setup.N_max
    args['gpu'] = setup.use_gpu
    args['fixed_seed'] = setup.fixed_seed
    try:
        args['like_mode'] = setup.like_mode
    except:
        pass
    try:
        args['rare_cut'] = setup.rare_cut
    except:
        pass
    args['gal_class'] = setup.model_class
    args['verbose'] = setup.verbose
    
    print('Running Nestle')
    sampler = fit_model.nested_integrate(**args)

    print('Nestle complete, saving results')
    #Used for saving output
    param_names = model_class._param_names
    N_params = len(param_names)
    chain_file = setup.chain_file
    
    #Save results of the chain
    chain_df = pd.DataFrame()
    for d in range(N_params):
        chain_df[param_names[d]] = sampler.samples[:,d]
    chain_df['lnlike'] = sampler.logl
    chain_df['weights'] = sampler.weights
    chain_df['logvol'] = sampler.logvol

    chain_df['niter'] = sampler.niter
    chain_df['log_evidence'] = sampler.logz
    chain_df['error_log_evidence'] = sampler.logzerr
    chain_df['information'] = sampler.h

    chain_df.to_csv(chain_file, index=False, float_format='%.4f', compression='gzip')
