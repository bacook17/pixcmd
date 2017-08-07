import pcmdpy as ppy

import numpy as np, pandas as pd
import sys

ppy.gpu_utils.initialize_gpu(n=0)

filters = np.array([ppy.instrument.Filter.HST_F475W(1.0), ppy.instrument.Filter.HST_F814W(1.0)])
iso_model = ppy.isochrones.Isochrone_Model(filters)

driv = ppy.driver.Driver(iso_model, gpu=True)

data_gals = []
data_gals.append(ppy.galaxy.Constant_SFR(np.array([-0.5, -.3, 4.])))
data_gals.append(ppy.galaxy.Constant_SFR(np.array([-1., -.1, 2.])))
data_gals.append(ppy.galaxy.Tau_Model(np.array([-.5, -.3, 4., 4.])))
data_gals.append(ppy.galaxy.Tau_Model(np.array([-.5, -.3, 2., 4.])))
data_gals.append(ppy.galaxy.Tau_Model(np.array([-1., -.1, 4., 1.])))
data_gals.append(ppy.galaxy.Tau_Model(np.array([-1., -.1, 2., 1.])))

lum_cuts = [np.inf, 1e4, 1e3]
N_im = [256, 512, 1024, 2048]
like_mode = [0,1,2]

vary_fracs = np.linspace(-.5, .5, 11)
N_sim = 1
#N_sim = 10
N_data = 512

cols = np.array(['gal_num', 'N_im', 'param_num', 'vary_frac', 'like_mode', 'lum_cut', 'log_like', 'data_params', 'data_meta', 'model_params', 'model_meta'])
results = pd.DataFrame(columns=cols)

for g, gal in enumerate(data_gals):
    print('Simulating Dataset %d'%g)
    row = {}
    row['gal_num'] = g
    row['data_params'] = gal._params
    row['data_meta'] = gal._meta_params

    for l_cut in lum_cuts:
        print('--Luminosity Cut: %.1e'%l_cut)
        row['lum_cut'] = l_cut
        mags, _ = driv.simulate(gal, N_data, lum_cut=l_cut, fixed_seed=True)
        data_pcmd = ppy.utils.make_pcmd(mags)
        driv.initialize_data(data_pcmd)

        for n in N_im:
            print('----N_im: %d'%n)
            row['N_im'] = n
            #vary the parameters
            for p in np.arange(-1, 9):
                print('------Varying parameter %d'%p)
                row['param_num'] = p
                for vf in vary_fracs:
                    if p == -1:
                        row['vary_frac'] = 0.
                        model_gal = gal
                    else:
                        row['vary_frac'] = vf
                        params = np.copy(gal._params)
                        params[p] += vf
                        model_gal = ppy.galaxy.Galaxy_Model(params)
                    row['model_params'] = model_gal._params
                    row['model_meta'] = model_gal._meta_params
                    for i in range(N_sim):
                        mags, _ = driv.simulate(model_gal, n, lum_cut=l_cut, fixed_seed=False)
                        model_pcmd = ppy.utils.make_pcmd(mags)
                        for l_mode in like_mode:
                            row['like_mode'] = l_mode
                            row['log_like'] = driv.loglike(model_pcmd, like_mode=l_mode)
                            results = results.append(row, ignore_index=True)

    results.to_csv('/n/home01/bcook/pixcmd/scripts_py/results/varations.csv', index=False, float_format='%.4e')
