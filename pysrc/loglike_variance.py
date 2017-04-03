import numpy as np, pandas as pd
import isochrones as iso, instrument as ins
import galaxy as gal
import driver, utils, gpu_utils

print('initializing')
filters = [ins.Filter.HST_F475W(1.), ins.Filter.HST_F814W(1.)]
iso_model = iso.Isochrone_Model(filters)
driver_gpu = driver.Driver(iso_model, gpu=True)

xbins, ybins = np.arange(-1.5, 4.6, 0.05), np.arange(-12, 15.6, 0.05)
bins = np.array([xbins, ybins])

loglikes = {}
N_samples = 1000
for log_Npix in [0., 2., 4.]:
    model = gal.Galaxy_SSP(np.array([0.,-2.,log_Npix,10.]))
    _, mags, _, _ = driver_gpu.simulate(model, 128, fixed_seed=True)
    pcmd = utils.make_pcmd(mags)
    driver_gpu.initialize_data(pcmd, bins)
    for N_scale in [128, 256, 512, 1024, 2048]:
        print((N_scale, log_Npix))
        temp = []
        for i in range(N_samples):
            _, mags, _, _ = driver_gpu.simulate(model, N_scale, fixed_seed=True)
            pcmd = utils.make_pcmd(mags)
            temp.append(driver_gpu.loglike(pcmd))
        label = "N_scale: %d, log_Npix: %d"%(N_scale, log_Npix)
        loglikes[label] = np.array(temp)

df = pd.DataFrame(loglikes)
df.to_csv('results/delta_chi2_fixed.csv', index=False)
