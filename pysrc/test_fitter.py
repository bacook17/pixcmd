# test_fitter.py
# Ben Cook (bcook@cfa.harvard.edu)

import numpy as np
import instrument as ins
import isochrones as iso
import galaxy as gal
import driver
import fit_model
import utils
import pandas as pd
import os
import sys, getopt

if __name__ == "__main__":
    
    N_scale = 1024
    N_walkers = 10
    N_burn = 20
    N_sample = 500
    gpu=True

    #Take in optional arguments from command line
    print('Loading command line arguments')
    usage_message = 'usage: test_fitter.py [--N_scale <N_scale>] [--N_walkers <N_walkers>] [--N_burn <N_burn>] '\
                    +'[--N_sample <N_sample>]'
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'aohn:', ['N_scale=', 'N_walkers=', 'N_burn=', 'N_sample=', 'no_gpu'])
    except getopt.GetoptError:
        print(usage_message)
        sys.exit(2)

    for opt, arg in opts:
        if opt == '--N_scale':
            print('Setting N_scale to %d'%int(arg))
            N_scale = int(arg)
        elif opt == '--N_walkers':
            print('Setting N_walkers to %d'%int(arg))
            N_walkers = int(arg)
        elif opt == '--N_burn':
            print('Setting N_burn to %d'%int(arg))
            N_burn = int(arg)
        elif opt == '--N_sample':
            print('Setting N_sample to %d'%int(arg))
            N_sample = int(arg)
        elif opt == '--no_gpu':
            print('Disabling GPU')
            gpu=False
        else:
            print(opt, arg)
            sys.exit(2)
    
    print('---Setting up model galaxy') 
    filters = np.array([ins.Filter.HST_F475W(1.0), ins.Filter.HST_F814W(1.0)])
    iso_model = iso.Isochrone_Model(filters)
    driv = driver.Driver(iso_model, gpu=gpu)
    
    true_params = [-0.2, -2, 2, 9.6]
    model_galaxy = gal.Galaxy_SSP(true_params)
    print('---Simulating model galaxy')
    _, mags, _, _ = driv.simulate(model_galaxy, N_scale)
    pcmd_model = utils.make_pcmd(mags)

    params = ['logz', 'logdust', 'logNpix', 'logage']

    print('---Running emcee')
    sampler = fit_model.sample_post(pcmd_model, filters, N_scale, N_walkers, N_burn, N_sample,
                                    gal_class=gal.Galaxy_SSP, gpu=gpu)

    print('---Emcee done, saving results')
    chain_df = pd.DataFrame()
    for w in np.arange(N_walkers):
        for d in np.arange(4):
            chain_df[params[d] + '_%d'%w] = sampler.chain[w,:,d]
        chain_df['lnprob_%d'%w] = sampler.lnprobability[w]

    if 'MacBook' in os.uname()[1]:
        pcmd_dir = '/Users/bcook/pCMDs/pixcmd/'
    else:
        pcmd_dir = '/n/home01/bcook/pixcmd/'
        
    chain_file = pcmd_dir + 'pysrc/results/test_chain.csv'
    chain_df.to_csv(chain_file, index=False)
            
    accept_df = pd.DataFrame()
    accept_df['acceptance'] = sampler.acceptance_fraction
    accept_file = pcmd_dir + 'pysrc/results/test_accept.csv'
    accept_df.to_csv(accept_file, index=False)
    
