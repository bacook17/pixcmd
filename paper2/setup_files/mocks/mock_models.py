import pcmdpy as ppy
import numpy as np

models = {}
run_names = {}
results = {}
results_dir = '~/pCMDs/pixcmd/paper2/results/'

models['mock_1'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.25),
        ppy.dustmodels.SingleDust(-1.0),
        ppy.sfhmodels.TauModel(np.array([3.0, 2.0])),
        ppy.distancemodels.VariableDistance(29.0)
    )
run_names['mock_1'] = 'Tau Model, Npix=1e3'

models['mock_2'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.25),
        ppy.dustmodels.SingleDust(-1.0),
        ppy.sfhmodels.TauModel(np.array([5.0, 2.0])),
        ppy.distancemodels.VariableDistance(29.0)
    )
run_names['mock_2'] = 'Tau Model, Npix=1e5'

tau_SFHbins = ppy.sfhmodels.TauModel(np.array([3.0, 2.0]), iso_step=-1).logSFH
models['mock_3'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.25),
        ppy.dustmodels.SingleDust(-1.0),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_3'] = r'NonParam Model, Npix=1e3'

tau_SFHbins = ppy.sfhmodels.TauModel(np.array([5.0, 2.0]), iso_step=-1).logSFH
models['mock_4'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.25),
        ppy.dustmodels.SingleDust(-1.0),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_4'] = r'NonParam Model, Npix=1e5'

models['mock_5'] = models['mock_1'].copy()
run_names['mock_5'] = 'Tau Model, Npix=1e3, LL=5'

models['mock_6'] = models['mock_2'].copy()
run_names['mock_6'] = 'Tau Model, Npix=1e5, LL=5'

models['mock_7'] = models['mock_3'].copy()
run_names['mock_7'] = r'NonParam Model, Npix=1e3, LL=5'

models['mock_8'] = models['mock_4'].copy()
run_names['mock_8'] = r'NonParam Model, Npix=1e5, LL=5'


for key in models.keys():
    res_file = results_dir + 'paper2_' + key + '.csv'
    live_file = res_file.replace('.csv', '_live.csv')
    try:
        results[key] = ppy.results.ResultsPlotter(
            res_file, live_file=live_file, run_name=run_names[key],
            gal_model=models[key], model_is_truth=('mock in key'))
    except (FileNotFoundError, ValueError, AttributeError) as e:
        print(key, e)
    except Exception as e:
        print(key, e)
        raise e
        
assert np.all([k in models for k in run_names.keys()]), str([k for k in run_names.keys() if k not in models])
assert np.all([k in run_names for k in models.keys()]), str([k for k in models.keys() if k not in run_names])
