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
    try:
        args['N_batch'] = setup.N_batch
    except:
        pass
    #optional key-word arguments (defaults are set by fit_model.nested_integrate)
    args['max_call'] = setup.N_max
    args['gpu'] = setup.use_gpu
    args['fixed_seed'] = setup.fixed_seed
    args['like_mode'] = setup.like_mode
    args['small_prior'] = setup.small_prior
    try:
        args['lum_cut'] = setup.lum_cut
    except:
        pass
    try:
        args['dlogz'] = setup.dlogz
    except:
        pass
    try:
        args['use_dynesty'] = setup.use_dynesty
    except:
        args['use_dynesty'] = False

    try:
        args['dynamic'] = setup.dynamic
    except:
        args['dynamic'] = False
    
    args['gal_class'] = setup.model_class
    args['verbose'] = setup.verbose
    
    print('Running Nested Sampling')
    results = fit_model.nested_integrate(**args)

    print('Nested Sampling Complete, saving results')
    #Used for saving output
    param_names = setup.model_class._param_names
    N_params = len(param_names)
    chain_file = setup.chain_file
    
    results_dict = {'nlive':'nlive', 'niter':'niter', 'ncall':'ncall', 'eff':'eff',
                     'logl':'lnlike', 'logwt':'logwt', 'logz':'log_evidence',
                     'logzerr':'e_log_evidence', 'h':'information', 'weights':'weights'}

    #Save results of the chain
    chain_df = pd.DataFrame()
    for d in range(N_params):
        chain_df[param_names[d]] = results.samples[:,d]
    #######Fix nomenclature
    for their_name, our_name in results_dict.iteritems():
        try:
            chain_df[our_name] = getattr(results, their_name)
        except:
            print('%s not found among result keys'%(their_name))
    
    chain_df.to_csv(chain_file, index=False, float_format='%.3e', compression='gzip')
