import pcmdpy as ppy

import numpy as np, pandas as pd
import sys
import time

ppy.gpu_utils.initialize_gpu(n=0)

filters = np.array([ppy.instrument.Filter.HST_F475W(1.0), ppy.instrument.Filter.HST_F814W(1.0)])
iso_model = ppy.isochrones.Isochrone_Model(filters)
iso_model_x5 = ppy.isochrones.Isochrone_Model(filters, MIST_path = '/n/home01/bcook/pixcmd/pcmdpy/isoc_MIST_v0.29_x5FEWER/')

driv_full = ppy.driver.Driver(iso_model, gpu=True)
driv_x5 = ppy.driver.Driver(iso_model_x5, gpu=True)

gal1 = ppy.galaxy.Constant_SFR(np.array([-0.5, -.3, 4.]))
gal2 = ppy.galaxy.Constant_SFR(np.array([-0.5, -.3, 2.]))
gal3 = ppy.galaxy.Constant_SFR(np.array([-0.5, -.3, 0.]))

N_im = [512, 1024, 2048]
N_sim = 50

mags, _ = driv_full.simulate(gal1, 512, fixed_seed=True)
p = ppy.utils.make_pcmd(mags)
driv_full.initialize_data(p)

mags, _ = driv_x5.simulate(gal1, 512, fixed_seed=True)
p = ppy.utils.make_pcmd(mags)
driv_x5.initialize_data(p)

print('Starting timing tests\n\n')

print('5x Undersampled')
driv = driv_x5
for n in N_im:
    print('N_im: %d'%n)
    start = time.time()
    for i in range(N_sim):
        mags, _ = driv.simulate(gal1, n, fixed_seed=False)
        p = ppy.utils.make_pcmd(mags)
        _ = driv.loglike(p)
    ave_t = (time.time() - start) / float(N_sim)
    print('---Gal1 (Npix = 1e4): %.2f sec'%(ave_t)) 
    start = time.time()
    for i in range(N_sim):
        mags, _ = driv.simulate(gal2, n, fixed_seed=False)
        p = ppy.utils.make_pcmd(mags)
        _ = driv.loglike(p)
    ave_t = (time.time() - start) / float(N_sim)
    print('---Gal2 (Npix = 1e2): %.2f sec'%(ave_t)) 
    start = time.time()
    for i in range(N_sim):
        mags, _ = driv.simulate(gal3, n, fixed_seed=False)
        p = ppy.utils.make_pcmd(mags)
        _ = driv.loglike(p)
    ave_t = (time.time() - start) / float(N_sim)
    print('---Gal3 (Npix = 1): %.2f sec'%(ave_t)) 

print('Full MIST')
driv = driv_full
for n in N_im:
    print('N_im: %d'%n)
    start = time.time()
    for i in range(N_sim):
        mags, _ = driv.simulate(gal1, n, fixed_seed=False)
        p = ppy.utils.make_pcmd(mags)
        _ = driv.loglike(p)
    ave_t = (time.time() - start) / float(N_sim)
    print('---Gal1 (Npix = 1e4): %.2f sec'%(ave_t)) 
    start = time.time()
    for i in range(N_sim):
        mags, _ = driv.simulate(gal2, n, fixed_seed=False)
        p = ppy.utils.make_pcmd(mags)
        _ = driv.loglike(p)
    ave_t = (time.time() - start) / float(N_sim)
    print('---Gal2 (Npix = 1e2): %.2f sec'%(ave_t)) 
    start = time.time()
    for i in range(N_sim):
        mags, _ = driv.simulate(gal3, n, fixed_seed=False)
        p = ppy.utils.make_pcmd(mags)
        _ = driv.loglike(p)
    ave_t = (time.time() - start) / float(N_sim)
    print('---Gal3 (Npix = 1): %.2f sec'%(ave_t)) 
                                
