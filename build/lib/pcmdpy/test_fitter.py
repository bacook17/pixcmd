# test_fitter.py
# Ben Cook (bcook@cfa.harvard.edu)

import numpy as np
import instrument as ins
import isochrones as iso
import galaxy as gal
import driver
import fit_model
import utils
import gpu_utils
import pandas as pd
import os
import sys, getopt
import multiprocessing
import emcee

if __name__ == "__main__":
    
    N_scale = 1024
    N_walkers = 10
    N_burn = 20
    N_sample = 500
    N_threads = 1
    pool = None
    gpu=True
    force_gpu=False
    offset=False
    ssp=False
    fixed_seed=False
    ball=False
    append = ''

    #Take in optional arguments from command line
    print('Loading command line arguments')
    usage_message = 'usage: test_fitter.py [--N_scale=<N_scale>] [--N_walkers=<N_walkers>] [--N_burn=<N_burn>] '\
                    +'[--N_sample=<N_sample>] [--no_gpu] [--require_gpu] [--require_cudac] [--SSP] [--append=<append>] [--N_threads=<N_threads>]'
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'aohn:', ['N_scale=', 'N_walkers=', 'N_burn=', 'N_sample=', 'fixed_seed', 'no_gpu', 'require_gpu', 'require_cudac', 'SSP',
                                                           'append=', 'N_threads=', 'sample_ball', 'offset_ball'])
    except getopt.GetoptError:
        print(usage_message)
        sys.exit(2)

    for opt, arg in opts:
        if opt == '--N_scale':
            print('Setting N_scale to %d'%int(arg))
            N_scale = int(arg)
        elif opt == '--N_threads':
            print('Using %d threads'%int(arg))
            N_threads = int(arg)
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
        elif opt == '--require_gpu':
            gpu=True
            print('Requiring GPU')
            if not gpu_utils._GPU_AVAIL:
                print('GPU NOT AVAILABLE. QUITTING')
                sys.exit(2)
        elif opt == '--require_cudac':
            gpu=True
            print('Requiring CUDAC')
            if not gpu_utils._CUDAC_AVAIL:
                print('CUDAC NOT AVAILABLE. QUITTING')
                sys.exit(2)
        elif opt == '--SSP':
            print('Using SSP')
            ssp = True
        elif opt == '--fixed_seed':
            print('Using fixed seed')
            fixed_seed = True
        elif opt == '--append':
            print('Appending %s to filenames'%arg)
            append = arg
        elif opt == '--sample_ball':
            print('Setting initial state as ball')
            ball = True
        elif opt == '--offset_ball':
            print('Initial chain ball will be slightly offset')
            ball = True
            offset = True
        else:
            print(opt, arg)
            sys.exit(2)
    
    print('---Setting up model galaxy') 
    filters = np.array([ins.Filter.HST_F475W(1.0), ins.Filter.HST_F814W(1.0)])
    iso_model = iso.Isochrone_Model(filters)
    driv = driver.Driver(iso_model, gpu=gpu)
    
    SSP_params = np.array([-0.2, -2., 2., 9.6])
    log_SFH_1e2 = np.log10(1e2 / 7.) 
    full_params = np.array([-0.2, -2, log_SFH_1e2, log_SFH_1e2,log_SFH_1e2,log_SFH_1e2,log_SFH_1e2,log_SFH_1e2,log_SFH_1e2,]) 

    if ssp:
        model_galaxy = gal.Galaxy_SSP(SSP_params)
    else:
        model_galaxy = gal.Galaxy_Model(full_params)
    print('---Simulating model galaxy')
    N_data = 128
    _, mags, _, _ = driv.simulate(model_galaxy, N_data, fixed_seed=fixed_seed)
    pcmd_model = utils.make_pcmd(mags)

    if ssp:
        params = ['logz', 'logdust', 'logNpix', 'logage']
    else:
        params = ['logz', 'logdust', 'logSFH0', 'logSFH1', 'logSFH2', 'logSFH3', 'logSFH4', 'logSFH5', 'logSFH6']

    if N_threads > 1:
        print('Setting up multiple GPUs')
        pool = multiprocessing.Pool(processes=N_threads, initializer=gpu_utils.initialize_gpu)

    if ball:
        if ssp:
            std = 0.1 * np.ones_like(SSP_params)
            if offset:
                SSP_params += 0.25 * np.random.random(size=len(SSP_params))
            p0 = emcee.utils.sample_ball(SSP_params, std, size=N_walkers)
        else:
            std = 0.1 * np.ones_like(full_params)
            if offset:
                full_params += 0.25 * np.random.random(size=len(full_params))
            p0 = emcee.utils.sample_ball(full_params, std, size=N_walkers)
    else:
        p0 = None

    print('---Running emcee')
    if ssp:
        sampler = fit_model.sample_post(pcmd_model, filters, N_scale, N_walkers, N_burn, N_sample, fixed_seed=fixed_seed,
                                        gal_class=gal.Galaxy_SSP, gpu=gpu, pool=pool, p0=p0)
    else:
        sampler = fit_model.sample_post(pcmd_model, filters, N_scale, N_walkers, N_burn, N_sample, fixed_seed=fixed_seed,
                                        gal_class=gal.Galaxy_Model, gpu=gpu, pool=pool, p0=p0)

    print('---Emcee done, saving results')
    chain_df = pd.DataFrame()
    for d in np.arange(len(params)):
        chain_df[params[d]] = sampler.flatchain[:,d]
    chain_df['lnprob'] = sampler.flatlnprobability
    chain_df['N_walkers'] = N_walkers

    if 'MacBook' in os.uname()[1]:
        pcmd_dir = '/Users/bcook/pCMDs/pixcmd/'
    else:
        pcmd_dir = '/n/home01/bcook/pixcmd/'
        
    chain_file = pcmd_dir + 'pysrc/results/chain%s.csv'%append
    chain_df.to_csv(chain_file, index=False, float_format='%.4f')
    
    accept_df = pd.DataFrame()
    accept_df['acceptance'] = sampler.acceptance_fraction
    accept_file = pcmd_dir + 'pysrc/results/accept%s.csv'%append
    accept_df.to_csv(accept_file, index=False, float_format='%.4f')
    
