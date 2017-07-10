import pcmdpy.instrument as ins
import pcmdpy.isochrones as iso
import pcmdpy.galaxy as gal
import pcmdpy.driver as driver
import pcmdpy.utils as utils
import pcmdpy.gpu_utils as gpu_utils

import numpy as np
import pandas as pd
import sys

#gpu_utils.initialize_gpu(n=0)

filters = np.array([ins.Filter.HST_F475W(1.0), ins.Filter.HST_F814W(1.0)])

iso_model = iso.Isochrone_Model(filters)
iso_model_x5 = iso.Isochrone_Model(filters, MIST_path='/n/home01/bcook/pixcmd/pcmdpy/isoc_MIST_v0.29_x5FEWER/')

gal_model = gal.Galaxy_Model(np.array([-0.2, -2., -1., -1., -1., -1., -1., -1., -1.]))

driv = driver.Driver(iso_model, gpu=False)
driv_x5 = driver.Driver(iso_model_x5, gpu=False)

N_scale = []
rare_cut = []
full_MIST = []
like_mode = []
log_like = []

N_data = 128
N_loops = 2

xbins = np.arange(-1.5, 4.6, 0.05)
ybins = np.arange(-12, 15.6, 0.05)
bins = np.array([xbins,ybins])

for r_cut in [0.]:
    for full in [True, False]:
        if full:
            print('Running rare_cut=%.1e, full MIST models'%(r_cut))
            d = driv
        else:
            print('Running rare_cut=%.1e, 5xFEWER MIST models'%(r_cut))
            d = driv_x5
        #Simulate mock data
        mags, _ = d.simulate(gal_model, N_data, rare_cut=r_cut, fixed_seed=True)
        data_pcmd = utils.make_pcmd(mags)
        d.initialize_data(data_pcmd, bins)

        for n in [128]:
            print('---N_scale = %d'%(n))
            for i in range(N_loops):
                print('       %d of %d'%(i, N_loops))
                mags, _ = d.simulate(gal_model, n, fixed_seed=False, rare_cut=r_cut)
                pcmd = utils.make_pcmd(mags)
                for l_mode in [0, 1, 2, 3]: 
                    N_scale.append(n)
                    rare_cut.append(r_cut)
                    full_MIST.append(full)
                    like_mode.append(l_mode)
                    log_like.append(d.loglike(pcmd, like_mode=l_mode))
        
df = pd.DataFrame({'N_scale':N_scale, 'rare_cut':rare_cut, 'full_MIST':full_MIST,
                   'log_like':log_like, 'like_mode':like_mode})
df.to_csv('/n/home01/bcook/pixcmd/scripts_py/results/loglike_tests_cpu.csv', index=False, 
          float_format='%.4e')
