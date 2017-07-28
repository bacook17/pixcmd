import pcmdpy as ppy

import numpy as np
import pandas as pd
import sys

ppy.gpu_utils.initialize_gpu(n=0)

filters = np.array([ppy.instrument.Filter.HST_F475W(1.0), ppy.instrument.Filter.HST_F814W(1.0)])
iso_model = ppy.isochrones.Isochrone_Model(filters)
iso_model_x5 = ppy.isochrones.Isochrone_Model(filters, MIST_path = '/n/home01/bcook/pixcmd/pcmdpy/isoc_MIST_v0.29_x5FEWER/')

full_gal = ppy.galaxy.Galaxy_Model
ssp_gal = ppy.galaxy.Galaxy_SSP

driv = ppy.driver.Driver(iso_model, gpu=True)
driv_x5 = ppy.driver.Driver(iso_model_x5, gpu=True)


#Create a list of random (deterministic) galaxy models to use as "data" sources
data_gals = []
r = np.random.RandomState(seed=0)
#for i in range(2):
for i in range(10):
    u = r.rand(10)
    v = ppy.fit_model.lnprior_transform(u)[:-1] #HACKHACKHACK
    gal = full_gal(v)
    gal.u = u #keep track of random param
    data_gals.append(gal)

    u = r.rand(4)
    v = ppy.fit_model.lnprior_transform_ssp(u)
    gal = ssp_gal(v)
    gal.u = u 
    data_gals.append(gal)

#lum_cuts = [np.inf, 1e2]
lum_cuts = [np.inf, 1e4, 3e3, 1e3, 3e2, 1e2, 3e1, 1e1, 3e0]
#N_scale = [128, 256, 512]
N_scale = [128, 256, 512, 1024, 2048]
like_mode = [0,1,2,3]

#N_sim = 1
N_sim = 10
N_data = 512
close_width = 0.02
similar_width = 0.1

cols = np.array(['gal_type', 'compare_type', 'N_scale', 'MIST_full', 'like_mode', 'lum_cut', 'log_like', 'data_params', 'model_params'])

results = pd.DataFrame(columns=cols)

for i, gal in enumerate(data_gals):
    print('Simulating dataset %d'%i)
    row = {}
    if isinstance(gal, full_gal):
        row['gal_type'] = 'Full'
        galaxy_model = full_gal
        transform = ppy.fit_model.lnprior_transform
    else:
        row['gal_type'] = 'SSP'
        galaxy_model = ssp_gal
        transform = ppy.fit_model.lnprior_transform_ssp
    row['data_params'] = gal._params
    
    for l_cut in lum_cuts:
        print('--Luminosity cut: %.1e'%l_cut)
        row['lum_cut'] = l_cut
        for MIST_full in [True, False]:
            print('----MIST_full: ' + str(MIST_full))
            row['MIST_full'] = MIST_full
            if MIST_full:
                d = driv
            else:
                d = driv_x5
            #simulate mock data
            mags, _ = d.simulate(gal, N_data, lum_cut=l_cut, fixed_seed=True)
            data_pcmd = ppy.utils.make_pcmd(mags)
            d.initialize_data(data_pcmd)

            for n in N_scale:
                #reset index to redraw same random galaxies
                r = np.random.RandomState(seed=0)
                print('------N_scale: %d'%n)
                row['N_scale'] = n
                for j in range(N_sim):
                    print('--------Iteration: %d'%j)
                    for compare in ['same', 'close', 'similar', 'rand']:
                        print('----------Compare Type: %s'%compare)
                        row['compare_type'] = compare
                        #compare same model to itself
                        if compare is 'same':
                            new_gal = gal

                        #compare same model to very similar ones
                        elif compare is 'close':
                            u = np.copy(gal.u)
                            u += (r.rand(len(u))-0.5) * close_width
                            #reflect around [0,1]
                            u[u<0.] = -u[u<0.]
                            u[u>1.] = 2 - u[u>1.]
                            v = transform(u)
                            if isinstance(gal, full_gal):
                                v = v[:-1]
                            new_gal = galaxy_model(v)

                        #compare same model to pretty similar ones
                        elif compare is 'similar':
                            u = np.copy(gal.u)
                            u += (r.rand(len(u))-0.5) * similar_width
                            #reflect around [0,1]
                            u[u<0.] = -u[u<0.]
                            u[u>1.] = 2 - u[u>1.]
                            v = transform(u)
                            if isinstance(gal, full_gal):
                                v = v[:-1]
                            new_gal = galaxy_model(v)

                        #compare to random model
                        elif compare is 'rand':
                            u = r.rand(len(gal.u))
                            v = transform(u)
                            if isinstance(gal, full_gal):
                                v = v[:-1]
                            new_gal = galaxy_model(v)

                        #simulate the model, and store the comparison
                        row['model_params'] = new_gal._params
                        mags, _ = d.simulate(new_gal, n, lum_cut=l_cut, fixed_seed=False)
                        model_pcmd = ppy.utils.make_pcmd(mags)
                        for l_mode in like_mode:
                            row['like_mode'] = l_mode
                            row['log_like'] = d.loglike(model_pcmd, like_mode=l_mode)
                            print('               appending like_mode: %d'%l_mode)
                            results = results.append(row, ignore_index=True)

results.to_csv('/n/home01/bcook/pixcmd/scripts_py/results/full_loglike_test.csv', index=False, float_format='%.4e')
