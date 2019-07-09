import pcmdpy as ppy
import numpy as np

models = {}
run_names = {}
results = {}
nlive = {}
results_dir = '~/pCMDs/pixcmd/paper2/results/'

models['mock_1'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.25),
        ppy.dustmodels.SingleDust(-1.0),
        ppy.sfhmodels.TauModel(np.array([3.0, 2.0])),
        ppy.distancemodels.VariableDistance(29.0)
    )
run_names['mock_1'] = 'Tau Model, Npix=1e3'
nlive['mock_1'] = 300

models['mock_2'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.25),
        ppy.dustmodels.SingleDust(-1.0),
        ppy.sfhmodels.TauModel(np.array([5.0, 2.0])),
        ppy.distancemodels.VariableDistance(29.0)
    )
run_names['mock_2'] = 'Tau Model, Npix=1e5'
nlive['mock_2'] = 300

tau_SFHbins = ppy.sfhmodels.TauModel(np.array([3.0, 2.0]), iso_step=-1).logSFH
models['mock_3'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.25),
        ppy.dustmodels.SingleDust(-1.0),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(29.0)
    )
run_names['mock_3'] = r'NonParam Model, Npix=1e3'
nlive['mock_3'] = 300

tau_SFHbins = ppy.sfhmodels.TauModel(np.array([5.0, 2.0]), iso_step=-1).logSFH
models['mock_4'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.25),
        ppy.dustmodels.SingleDust(-1.0),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(29.0)
    )
run_names['mock_4'] = r'NonParam Model, Npix=1e5'
nlive['mock_4'] = 300

models['mock_5'] = models['mock_1'].copy()
run_names['mock_5'] = 'Tau Model, Npix=1e3, LL=5'
nlive['mock_5'] = nlive['mock_1']

models['mock_6'] = models['mock_2'].copy()
run_names['mock_6'] = 'Tau Model, Npix=1e5, LL=5'
nlive['mock_6'] = nlive['mock_2']

models['mock_7'] = models['mock_3'].copy()
run_names['mock_7'] = r'NonParam Model, Npix=1e3, LL=5'
nlive['mock_7'] = nlive['mock_3']

models['mock_8'] = models['mock_4'].copy()
run_names['mock_8'] = r'NonParam Model, Npix=1e5, LL=5'
nlive['mock_8'] = nlive['mock_4']

models['mock_9'] = models['mock_1'].copy()
run_names['mock_9'] = 'Tau Model, Npix=1e3, LL=2'
nlive['mock_9'] = nlive['mock_1']

models['mock_10'] = models['mock_2'].copy()
run_names['mock_10'] = 'Tau Model, Npix=1e5, LL=2'
nlive['mock_10'] = nlive['mock_2']

models['mock_11'] = models['mock_3'].copy()
run_names['mock_11'] = r'NonParam Model, Npix=1e3, LL=2'
nlive['mock_11'] = nlive['mock_3']

models['mock_12'] = models['mock_4'].copy()
run_names['mock_12'] = r'NonParam Model, Npix=1e5, LL=2'
nlive['mock_12'] = nlive['mock_4']

for k in ['mock_13', 'mock_14', 'mock_15', 'mock_16', 'mock_17', 'mock_18', 'mock_19', 'mock_20']:
    models[k] = models['mock_1'].copy()
    run_names[k] = models['mock_1']
    nlive[k] = nlive['mock_1']

nlive['mock_20'] = 100

for key in models.keys():
    res_file = results_dir + 'paper2_' + key + '.csv'
    live_file = res_file.replace('.csv', '_live.csv')
    try:
        results[key] = ppy.results.ResultsPlotter(
            res_file, live_file=live_file, run_name=run_names[key],
            gal_model=models[key], model_is_truth=('mock in key'), nlive=nlive[key],
            nlive_batch=100)
    except (FileNotFoundError, ValueError, AttributeError) as e:
        print(key, e)
    except Exception as e:
        print(key, e)
        raise e
        
assert np.all([k in models for k in run_names.keys()]), str([k for k in run_names.keys() if k not in models])
assert np.all([k in run_names for k in models.keys()]), str([k for k in models.keys() if k not in run_names])
