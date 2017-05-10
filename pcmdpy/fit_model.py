# fit_model.py
# Ben Cook (bcook@cfa.harvard.edu)

import numpy as np
import instrument as ins
import isochrones as iso
import galaxy as gal
import driver
import utils
import emcee
import nestle

def lnprior_ssp(gal_params):
    z, log_dust, log_Npix, age = gal_params
    #Flat priors
    # log (z/z_solar) between -2 and 0.5
    if (z < -2.) or (z>0.5):
        return -np.inf
    # E(B-V) between 1e-2 and 3
    if (log_dust < -2) or (log_dust > 0.5):
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
    results[0] = -2. + 2.5*normed_params[0]
    # E(B-V) between 1e-2 and 3
    results[1] = -2 + 2.5*normed_params[1]
    #log Npix between -1 and 6
    results[2] = -1 + 7*normed_params[2]
    #age between 6 (1 Myr) and 10.3 (50 Gyr)
    results[3] = 6 + 4.3*normed_params[3]
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
    # E(B-V) between 1e-2 and 3
    if (log_dust < -2) or (log_dust > 0.5):
        return -np.inf
    # log M_i / Mtot between 1e-6 and 1
    for log_m in log_SFH:
        if (log_m < -10 + log_Npix) or (log_m > 0 + log_Npix):
            return -np.inf
    return 0.

def lnprior_transform(normed_params):
    results = np.zeros(len(normed_params))
    #Flat priors
    # log (z/z_solar) between -2 and 0.5
    results[0] = -2. + 2.5*normed_params[0]
    # E(B-V) between 1e-2 and 3
    results[1] = -2 + 2.5*normed_params[1]
    #log Npix between -1 and 6
    log_Npix = -1 + 7*normed_params[2]
    total = 0.
    # log M_i between -6 and 0
    for i in range(2, len(normed_params)-1):
        log_mi = -6 + 6*normed_params[i+1]
        results[i] = log_Npix + log_mi 
        total += 10.**results[i]
    results[-1] = np.log10(10.**log_Npix - total)
    return results

def lnlike(gal_params, driv, im_scale, gal_class=gal.Galaxy_Model, **kwargs):
    if (gal_class is gal.Galaxy_SSP):
        pri = lnprior_ssp(gal_params)
    else:
        pri = lnprior(gal_params)
    if np.isinf(pri):
        return -np.inf
    gal_model = gal_class(gal_params)
    _, mags, _, _ = driv.simulate(gal_model, im_scale, **kwargs)
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

def nested_integrate(pcmd, filters, im_scale, n_points=200, method='multi', max_iters=100000, gal_class=gal.Galaxy_Model, bins=None, **kwargs):
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
        this_pri_transform = lnprior_transform
    else:
        this_pri_transform = lnprior_transform_ssp

    def this_lnlke(gal_params):
        return lnlike(gal_params, driv, im_scale, gal_class=gal_class, **kwargs)

    print('-Running nestle sampler')
    sampler = nestle.sample(this_lnlike, this_pri_transform, ndim, method=method, npoints=n_points, maxiter=max_iters)
    return sampler

def sample_post(pcmd, filters, im_scale, n_walkers, n_burn, n_sample, 
                p0=None, gal_class=gal.Galaxy_Model, gpu=True, bins=None, threads=1, fixed_seed=True, add_total=True,
                rare_cut=0.,
                **kwargs):

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

    print('-Setting up emcee sampler')
    
    sampler = emcee.EnsembleSampler(n_walkers, n_dim, lnprob, args=[driv, im_scale], kwargs={'gal_class':gal_class, 'fixed_seed':fixed_seed,'add_total':add_total, 'rare_cut':rare_cut},
                                    threads=threads, **kwargs)

    if p0 is None:
        if (gal_class is gal.Galaxy_Model):
            np.random.seed(0)
            z0 = np.random.uniform(-2., 0.5, n_walkers)
            np.random.seed(0)
            dust0 = np.random.uniform(-6, 0, n_walkers)
            np.random.seed(0)
            npix0 = 10.**np.random.uniform(-1, 6, n_walkers)
            sfh0 = 10.**np.random.uniform(-10, 0, (7, n_walkers))
            sfh0 *= npix0 / np.sum(sfh0, axis=0)
            sfh0 = np.log10(sfh0)
            p0 = np.array([z0, dust0])
            p0 = np.concatenate([p0, sfh0]).T
        else:
            np.random.seed(0)
            z0 = np.random.uniform(-2., 0.5, n_walkers)
            np.random.seed(0)
            dust0 = np.random.uniform(-6, 0, n_walkers)
            np.random.seed(0)
            npix0 = np.random.uniform(-1, 6, n_walkers)
            np.random.seed(0)
            age0 = np.random.uniform(6, 10.3, n_walkers)
            p0 = np.array([z0, dust0, npix0, age0]).T
    assert(p0.shape == (n_walkers, n_dim))

    if n_burn > 0:
        print('-emcee burn-in')
        pos,prob,state = sampler.run_mcmc(p0, n_burn)

        print('-emcee sampling')
        sampler.reset()
        sampler.run_mcmc(pos, n_sample)

    else:
        sampler.run_mcmc(p0, n_sample)

    return sampler

