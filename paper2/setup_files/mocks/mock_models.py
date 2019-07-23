import pcmdpy as ppy
import numpy as np

models = {}
run_names = {}
results = {}
nlive = {}
results_dir = '~/pCMDs/pixcmd/paper2/results/'

models['mock_1'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.3),
        ppy.dustmodels.SingleDust(-1.4),
        ppy.sfhmodels.TauModel(np.array([2.0, 2.0])),
        ppy.distancemodels.VariableDistance(25.0)
    )
run_names['mock_1'] = 'Tau Model, 1 Mpc'

models['mock_2'] = models['mock_1'].copy()
models['mock_2'].set_params(np.array([-0.3, -1.4, 3.0, 2.0, 27.5]))
run_names['mock_2'] = 'Tau Model, 3 Mpc'

models['mock_3'] = models['mock_1'].copy()
models['mock_3'].set_params(np.array([-0.3, -1.4, 3.0, 2.0, 30.0]))
run_names['mock_3'] = 'Tau Model, 10 Mpc'

models['mock_4'] = models['mock_1'].copy()
models['mock_4'].set_params(np.array([-0.3, -1.4, 5.0, 2.0, 32.5]))
run_names['mock_4'] = 'Tau Model, 30 Mpc'

models['mock_5'] = models['mock_1'].copy()
models['mock_5'].set_params(np.array([-0.3, -1.4, 6.0, 2.0, 35.0]))
run_names['mock_5'] = 'Tau Model, 100 Mpc'

tau_SFHbins = ppy.sfhmodels.TauModel(np.array([3.0, 2.0]), iso_step=-1).logSFH
models['mock_6'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.3),
        ppy.dustmodels.SingleDust(-1.4),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(30.0)
    )
run_names['mock_6'] = 'NonParam Model, Npix=1e3'

tau_SFHbins = ppy.sfhmodels.TauModel(np.array([5.0, 2.0]), iso_step=-1).logSFH
models['mock_7'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.3),
        ppy.dustmodels.SingleDust(-1.4),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(30.0)
    )
run_names['mock_7'] = 'NonParam Model, Npix=1e3'

for i in range(1, 8):
    models[f'mock_2{i:d}'] = models[f'mock_{i:d}'].copy()
    run_names[f'mock_2{i:d}'] = run_names[f'mock_{i:d}'] + ' , Hess Binned'

for key in models.keys():
    nlive[key] = 500
    res_file = results_dir + 'paper2_newmock_' + key.strip('mock_') + '.csv'
    live_file = res_file.replace('.csv', '_live.csv')
    try:
        results[key] = ppy.results.ResultsPlotter(
            res_file, live_file=live_file, run_name=run_names[key],
            gal_model=models[key], model_is_truth=('mock in key'), nlive=nlive[key],
            max_logl=np.inf,
            nlive_batch=100)
    except (FileNotFoundError, ValueError, AttributeError) as e:
        print(key, e)
    except Exception as e:
        print(key, e)
        raise e
        
assert np.all([k in models for k in run_names.keys()]), str([k for k in run_names.keys() if k not in models])
assert np.all([k in run_names for k in models.keys()]), str([k for k in models.keys() if k not in run_names])
