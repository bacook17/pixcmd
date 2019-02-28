import pcmdpy_gpu as ppy
import numpy as np
import pandas as pd
from importlib import util
from setup_files.mocks_paper1.mock_models import models, run_names, results


def get_max_logl(key, n_samples=100, verbose=True):
    config_mod = 'setup_files/mocks_paper1'
    config_file = f'setup_files/mocks_paper1/{key}.py'
    spec = util.spec_from_file_location(config_mod, config_file)
    config = util.module_from_spec(spec)
    spec.loader.exec_module(config)
    params = config.params
    
    assert params['data_is_mock'], "Have not implemented data runs yet"
    
    driv.filters = params['filters']
    driv.iso_model = params['iso_model']
    driv.initialize_data(params['data_pcmd'], bins=params['bins'])
    best_params = results[key].best_params
    Nim = params['Nim']
    gal_model = params['gal_model']
    prior = params.get('prior', gal_model.get_flat_prior())
    
    logl_kwargs = {'like_mode': params['like_mode'],
                   'downsample': params['downsample'],
                   'mag_system': params['mag_system'],
                   'sky_noise': params['sky_noise'],
                   'shot_noise': params['shot_noise']}
    
    logls = [ppy.fit_model.lnlike(best_params, driv, Nim,
                                  prior.lnprior, gal_model,
                                  fixed_seed=False,
                                  **logl_kwargs) for _ in range(n_samples)]
    return np.median(logls)


def add_max_logl(key, max_logl):
    filename = f'results/paper1_{key}.csv'
    with open(filename, 'r') as f:
        text = f.read()
        if 'max_logl' in text:
            text = text.partition('\n')[-1]
    with open(filename, 'w') as f:
        f.write(f'# max_logl : {max_logl:.5e}\n')
        f.write(text)


if __name__ == '__main__':
    f_mock = ppy.instrument.default_m31_filters()
    iso_model = ppy.isochrones.Isochrone_Model(f_mock, mag_system='vega')
    driv = ppy.driver.Driver(iso_model, gpu=True)

    max_logls = {}
    for i, k in enumerate(['mock_11']):
        print(f'{k}, {i} of {len(results)}')
        max_logl = get_max_logl(k, n_samples=100, verbose=False)
        max_logls[k] = max_logl
        print(k, max_logl)
        add_max_logl(k, max_logl)

