__all__ = ['models', 'results', 'pcmds', 'data', 'radii_am', 'radii_kpc',
           'dmods', 'regions']

try:
    import pcmdpy_gpu as ppy
except:
    import pcmdpy as ppy
import numpy as np
import pandas as pd
from os.path import expanduser

models = {}
results = {}
pcmds = {}
data = {}
results_dir = expanduser('~/pCMDs/pixcmd/paper2/results/')
data_dir = expanduser('~/pCMDs/pixcmd/data/')
df_radii = pd.read_csv(results_dir + 'regions_radii.csv')
radii_am = {}
radii_kpc = {}
regions = {}

dmods = {
    'NGC3377': 30.18,
    'NGC4993': 33.05,
    'M87': 31.02,
    'M49': 31.0,
    'DF2': 31.505,
}


base_models = {}
base_models[1] = ppy.galaxy.CustomGalaxy(
    ppy.metalmodels.SingleFeH(),
    ppy.dustmodels.SingleDust(),
    ppy.sfhmodels.NonParam(),
    ppy.distancemodels.VariableDistance()
)
base_models[2] = ppy.galaxy.CustomGalaxy(
    ppy.metalmodels.SingleFeH(),
    ppy.dustmodels.SingleDust(),
    ppy.sfhmodels.NonParam(),
    ppy.distancemodels.FixedDistance()
)
base_models[3] = ppy.galaxy.CustomGalaxy(
    ppy.metalmodels.FixedWidthNormMDF(0.3),
    ppy.dustmodels.SingleDust(),
    ppy.sfhmodels.NonParam(),
    ppy.distancemodels.VariableDistance()
)
base_models[4] = ppy.galaxy.CustomGalaxy(
    ppy.metalmodels.SingleFeH(),
    ppy.dustmodels.SingleDust(),
    ppy.sfhmodels.TauModel(),
    ppy.distancemodels.VariableDistance()
)
base_models[5] = ppy.galaxy.CustomGalaxy(
    ppy.metalmodels.SingleFeH(),
    ppy.dustmodels.SingleDust(),
    ppy.sfhmodels.SSPModel(),
    ppy.distancemodels.VariableDistance()
)


def add_set(galaxy, mnum, region, key,
            colors='z_gz'):
    data_file = data_dir + f'{galaxy.lower()}/pcmds/{galaxy}_{colors}_{region}.pcmd'
    res_file = results_dir + f'{galaxy}_r{region}_m{mnum}.csv'
    live_file = res_file.replace('.csv', '_live.csv')
    pcmd_file = res_file.replace('.csv', '.pcmd')
    regions[key] = region
    models[key] = base_models[mnum].copy()
    results[key] = ppy.results.ResultsPlotter(
        res_file, live_file=live_file, dmod_true=dmods[galaxy],
        gal_model=models[key], model_is_truth=False)
    data[key] = np.loadtxt(data_file, unpack=True)
    radii_am[key] = df_radii[galaxy][region] * 0.05 / 60.
    radii_kpc[key] = radii_am[key] * (np.pi/(180.*60.)) * 1e3 * ppy.distancemodels.dmod_to_mpc(dmods[galaxy])
    try:
        pcmds[key] = np.loadtxt(pcmd_file, unpack=True)
    except:
        pass

# M87
print('M87')
for m in range(1, 6):
    add_set('M87', m, 204, f'M87_a1_m{m}', colors='I_VI')
    add_set('M87', m, 128, f'M87_b1_m{m}', colors='I_VI')
    add_set('M87', m, 44, f'M87_c1_m{m}', colors='I_VI')

print('M49')
for m in range(1, 6):
    add_set('M49', m, 204, f'M49_a1_m{m}')
    add_set('M49', m, 124, f'M49_b1_m{m}')
    add_set('M49', m, 40, f'M49_c1_m{m}')
    
print('NGC3377')
for m in range(1, 6):
    add_set('NGC3377', m, 173, f'NGC3377_a1_m{m}')
    add_set('NGC3377', m, 97, f'NGC3377_b1_m{m}')
    add_set('NGC3377', m, 41, f'NGC3377_c1_m{m}')
