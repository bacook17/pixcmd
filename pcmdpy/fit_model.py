# fit_model.py
# Ben Cook (bcook@cfa.harvard.edu)

import numpy as np
import instrument as ins
import isochrones as iso
import galaxy as gal
import sys
import driver
import utils
import emcee
import nestle
try:
    import dynesty
except:
    _DYNESTY_INSTALLED = False
else:
    _DYNESTY_INSTALLED = True
from datetime import datetime

def lnprior_ssp(gal_params):
    z, log_dust, log_Npix, age = gal_params
    #Flat priors
    # log (z/z_solar) between -2 and 0.5
    if (z < -2.) or (z>0.5):
        return -np.inf
    # E(B-V) between 1e-3 and 3
    if (log_dust < -3) or (log_dust > 0.5):
        return -np.inf
    #Npix between 0.1 and 1e6
    if (log_Npix < -1) or (log_Npix > 6):
        return -np.inf
    #age between 6 (1 Myr) and 10.3 (50 Gyr)
    if (age < 6.) or (age > 10.3):
        return -np.inf
    return 0.

def lnprior_transform_ssp(normed_params):
    results = np.zeros(len(normed_params))
    #Flat priors
    # log (z/z_solar) between -2 and 0.5
    results[0] = -2 + 2.5*normed_params[0]
    # E(B-V) between 1e-3 and 3
    results[1] = -3 + 3.5*normed_params[1]
    #log Npix between -1 and 6
    results[2] = -1 + 7*normed_params[2]
    #age between 6 (1 Myr) and 10.3 (50 Gyr)
    results[3] = 6 + 4.3*normed_params[3]
    return results

def lnprior_transform_ssp_small(normed_params):
    results = np.zeros(len(normed_params))
    #Flat priors
    # log (z/z_solar) between -0.5 and 0.0
    results[0] = -.5 + .5*normed_params[0]
    # E(B-V) between -2.5 and -1.5
    results[1] = -2.5 + normed_params[1]
    #log Npix between 1.5 and 2.5
    results[2] = 1.5 + normed_params[2]
    #age between 9.5 (3 Gyr) and 10.0 (10 Gyr)
    results[3] = 9.5 + 0.5*normed_params[3]
    return results

def lnprior(gal_params):
    z, log_dust = gal_params[:2]
    log_SFH = gal_params[2:]
    log_Npix = np.log10(np.sum(10.**log_SFH))
    #Flat priors
    # log Npix between -1 and 6
    if (log_Npix < -1.) or (log_Npix > 6.):
        return -np.inf
    # log (z/z_solar) between -2 and 0.5
    if (z < -2.) or (z > 0.5):
        return -np.inf
    # E(B-V) between 1e-3 and 3
    if (log_dust < -3) or (log_dust > 0.5):
        return -np.inf
    # log M_i / Mtot between 1e-6 and 1
    for log_m in log_SFH:
        if (log_m < -10 + log_Npix) or (log_m > 0 + log_Npix):
            return -np.inf
    return 0.

def lnprior_transform(normed_params):
    #HUGE HACK: returns more dimensions than used in model!
    results = np.zeros(len(normed_params))
    #Flat priors
    # log (z/z_solar) between -2 and 0.5
    results[0] = -2. + 2.5*normed_params[0]
    # E(B-V) between 1e-3 and 3
    results[1] = -3 + 3.5*normed_params[1]
    # log M_i between -6 and 0
    for i in range(2, len(normed_params)-1):
        results[i] = -6 + 6*normed_params[i]
    log_total = np.log10(np.sum(10.**results[2:-1]))
    #HACKHACKHACK
    #log Npix between -1 and 6
    log_Npix = -1 + 7*normed_params[2]
    results[-1] = log_Npix
    #Normalize the mass bins to sum to log_Npix
    results[2:-1] += log_Npix - log_total
    return results

def lnprior_transform_small(normed_params):
    results = np.zeros(len(normed_params))
    #Flat priors
    # log (z/z_solar) between -0.5 and 0.0
    results[0] = -0.5 + 0.5*normed_params[0]
    # E(B-V) between 3e-3 and 3e-2
    results[1] = -2.5 + normed_params[1]
    # log M_i between +/- 0.5 of truth
    appx_truth = np.array([-1.25, -0.25,  0.135,  0.635,
                            1.135, 1.635,  1.57])
    for i in range(2, len(normed_params)):
        results[i] = appx_truth[i-2] - 0.5 + normed_params[i]
    return results

def lnlike(gal_params, driv, im_scale, gal_class=gal.Galaxy_Model, **kwargs):
    if (gal_class is gal.Galaxy_SSP):
        pri = lnprior_ssp(gal_params)
    else:
        pri = lnprior(gal_params)
    if np.isinf(pri):
        return -np.inf
    gal_model = gal_class(gal_params)
    mags, _ = driv.simulate(gal_model, im_scale, **kwargs)
    pcmd = utils.make_pcmd(mags)
    like = driv.loglike(pcmd, **kwargs)

    return like

def lnprob(gal_params, driv, im_scale, gal_class=gal.Galaxy_Model, **kwargs):
    if (gal_class is gal.Galaxy_SSP):
        pri = lnprior_ssp(gal_params)
    else:
        pri = lnprior(gal_params)
    if np.isinf(pri):
        return -np.inf
    like = lnlike(gal_params, driv, im_scale, gal_class=gal_class, **kwargs)
    return pri + like



def nested_integrate(pcmd, filters, im_scale, N_points, method='multi', max_call=100000, gal_class=gal.Galaxy_Model, gpu=True,
                     bins=None, verbose=False, small_prior=False, dlogz=None, use_dynesty=False, **kwargs):
    if (not _DYNESTY_INSTALLED) and use_dynesty:
        raise ImportError('Dynesty not installed correctly')
    print('-initializing models')
    n_filters = len(filters)
    assert(pcmd.shape[0] == n_filters)
    n_dim = gal_class._num_params
    
    iso_model = iso.Isochrone_Model(filters)
    driv = driver.Driver(iso_model, gpu=gpu)
    if bins is None:
        assert(n_filters == 2)
        xbins = np.arange(-1.5, 4.6, 0.05)
        ybins = np.arange(-12, 15.6, 0.05)
        bins = np.array([xbins,ybins])
    driv.initialize_data(pcmd,bins)

    if gal_class is gal.Galaxy_Model:
        if small_prior:
            this_pri_transform = lnprior_transform_small
        else:
            ndim += 1 #HACKHACKHACK
            this_pri_transform = lnprior_transform
    else:
        if small_prior:
            this_pri_transform = lnprior_transform_ssp_small
        else:
            this_pri_transform = lnprior_transform_ssp

    def this_lnlike(gal_params):
        #HACKHACKHACK to remove trailing zero
        if (not small_prior) and (gal_class is gal.Galaxy_Model):
            gal_params = gal_params[:-1]
        return lnlike(gal_params, driv, im_scale, gal_class=gal_class, **kwargs)

    callback = None
    if verbose:
        def my_progress(callback_info):
            it = callback_info['it']
            logz = callback_info['logz']
            n_calls = driv.num_calls
            print('----------------')
            if (np.abs(logz) <= 1e5):
                print('Iteration Number: %d, Likelihood Calls: %d, logz: %.2f'%(it, n_calls, logz))
            else:
                print('Iteration Number: %d, Likelihood Calls: %d, logz: %.4e'%(it, n_calls, logz))
            print('Current time: %s'%(str(datetime.now())))
            sys.stdout.flush()
        callback = my_progress

    #Initialize the nestle sampler with a different random state than global
    #This is important because the driver resets the global seed
    rstate = np.random.RandomState(1234)

    if use_dynesty:
        sampler = dynesty.NestedSampler(this_lnlike, this_pri_transform, ndim=n_dim, bound=method, sample='unif', nlive=N_points,
                                        update_interval=1, rstate=rstate)
        print('-Running dynesty sampler')
        for it, results in enumerate(sampler.sample(dlogz=dlogz,maxcall=max_call)):
            (worst, ustar, vstar, loglstar, logvol, logwt, logz, logzerr, h, nc) = results
            #compute delta_logz
            logz_remain = np.max(sampler.live_logl) + logvol
            delta_logz = np.logaddexp(logz, logz_remain) - logz
            message = 'iteration: %d | ncalls: %d | logz: %6.3f +/- %6.3f | dlogz: %6.3f'%(it, nc, logz, logzerr, delta_logz)
            message += '\n Current time: ' + '%s'%(str(datetime.now()))
            message += '\n --------------------------'
            print(message)
        results = sampler.results
    else:
        print('-Running nestle sampler')
        results = nestle.sample(this_lnlike, this_pri_transform, n_dim, method=method, npoints=N_points, maxcall=max_call, callback=callback,
                                update_interval=1, rstate=rstate, dlogz=dlogz)

    if driv.num_calls >= (max_call - 1):
        print('Terminated after surpassing max likelihood calls')
    else:
        print('Reached desired convergence')

    return results

def sample_post(pcmd, filters, im_scale, N_walkers, N_burn, N_sample, 
                p0=None, gal_class=gal.Galaxy_Model, gpu=True, bins=None, threads=1, fixed_seed=True,
                rare_cut=0., like_mode=0, 
                **kwargs):

    print('-initializing models')
    N_filters = len(filters)
    assert(pcmd.shape[0] == N_filters)
    N_dim = gal_class._num_params
    
    iso_model = iso.Isochrone_Model(filters)
    driv = driver.Driver(iso_model, gpu=gpu)
    if bins is None:
        assert(N_filters == 2)
        xbins = np.arange(-1.5, 4.6, 0.05)
        ybins = np.arange(-12, 15.6, 0.05)
        bins = np.array([xbins,ybins])
    driv.initialize_data(pcmd,bins)

    print('-Setting up emcee sampler')
    
    sampler = emcee.EnsembleSampler(N_walkers, N_dim, lnprob, args=[driv, im_scale], kwargs={'gal_class':gal_class, 'fixed_seed':fixed_seed,'rare_cut':rare_cut, "like_mode":like_mode},
                                    threads=threads, **kwargs)

    if p0 is None:
        if (gal_class is gal.Galaxy_Model):
            np.random.seed(0)
            z0 = np.random.uniform(-2., 0.5, N_walkers)
            np.random.seed(0)
            dust0 = np.random.uniform(-6, 0, N_walkers)
            np.random.seed(0)
            npix0 = 10.**np.random.uniform(-1, 6, N_walkers)
            sfh0 = 10.**np.random.uniform(-10, 0, (7, N_walkers))
            sfh0 *= npix0 / np.sum(sfh0, axis=0)
            sfh0 = np.log10(sfh0)
            p0 = np.array([z0, dust0])
            p0 = np.concatenate([p0, sfh0]).T
        else:
            np.random.seed(0)
            z0 = np.random.uniform(-2., 0.5, N_walkers)
            np.random.seed(0)
            dust0 = np.random.uniform(-6, 0, N_walkers)
            np.random.seed(0)
            npix0 = np.random.uniform(-1, 6, N_walkers)
            np.random.seed(0)
            age0 = np.random.uniform(6, 10.3, N_walkers)
            p0 = np.array([z0, dust0, npix0, age0]).T
    assert(p0.shape == (N_walkers, N_dim))

    if N_burn > 0:
        print('-emcee burn-in')
        pos,prob,state = sampler.run_mcmc(p0, N_burn)

        print('-emcee sampling')
        sampler.reset()
        sampler.run_mcmc(pos, N_sample)

    else:
        sampler.run_mcmc(p0, N_sample)

    return sampler

