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
df_radii = pd.read_csv(results_dir + 'regions_radii.csv', index_col=0)
radii_am = {}
radii_kpc = {}
regions = {}

dmods = {
    'NGC3377': 30.18,
    'NGC4993': 33.05,
    'M87': 31.02,
    'M87v2': 31.02,
    'M49': 31.0,
    'DF2': 31.505,
    'M31': 24.42,
    'M51': 29.67
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
base_models[6] = base_models[1].copy()
base_models[7] = base_models[5].copy()


def add_set(galaxy, mnum, region, key,
            colors='z_gz', model=None):
    g_orig = galaxy.replace('v2', '')
    data_file = data_dir + f'{galaxy.lower()}/pcmds/{g_orig}_{colors}_{region}.pcmd'
    res_file = results_dir + f'{galaxy}_r{region}_m{mnum}.csv'
    live_file = res_file.replace('.csv', '_live.csv')
    pcmd_file = res_file.replace('.csv', '.pcmd')
    regions[key] = region
    models[key] = model or base_models[mnum].copy()
    results[key] = ppy.results.ResultsPlotter(
        res_file, live_file=live_file, dmod_true=dmods[galaxy],
        gal_model=models[key], model_is_truth=False)
    if key.replace(f'_m{mnum}', '_m1') in data:
        data[key] = data[key.replace(f'_m{mnum}', '_m1')]
    elif key not in data:
        data[key] = np.loadtxt(data_file, unpack=True)
    radii_am[key] = df_radii[g_orig][region] * 0.05 / 60.
    radii_kpc[key] = radii_am[key] * (np.pi/(180.*60.)) * 1e3 * ppy.distancemodels.dmod_to_mpc(dmods[galaxy])
    try:
        pcmds[key] = np.loadtxt(pcmd_file, unpack=True)
    except:
        pass

    
def add_set_v2(galaxy, mnum, region, key, data_name, model=None):
    data_file = data_dir + f'{galaxy.lower()}/pcmds/{data_name}.pcmd'
    res_file = results_dir + f'{galaxy}_{region}_m{mnum}.csv'
    live_file = res_file.replace('.csv', '_live.csv')
    pcmd_file = res_file.replace('.csv', '.pcmd')
    regions[key] = region
    models[key] = model or base_models[mnum].copy()
    results[key] = ppy.results.ResultsPlotter(
        res_file, live_file=live_file, dmod_true=dmods[galaxy],
        gal_model=models[key], model_is_truth=False)
    if key.replace(f'_m{mnum}', '_m1') in data:
        data[key] = data[key.replace(f'_m{mnum}', '_m1')]
    else:
        data[key] = np.loadtxt(data_file, unpack=True)
    try:
        pcmds[key] = np.loadtxt(pcmd_file, unpack=True)
    except:
        pass


def load_M87():
    # M87
    print('M87')
    for m in range(1, 6):
        add_set('M87', m, 204, f'M87_a1_m{m}', colors='I_VI')
        add_set('M87', m, 128, f'M87_b1_m{m}', colors='I_VI')
        add_set('M87', m, 44, f'M87_c1_m{m}', colors='I_VI')
    for i, r in enumerate([204, 201, 202, 203]):
        add_set('M87', 6, r, f'M87_a1_m6', colors='I_VI')
    for i, r in enumerate([128, 125, 126, 127]):
        add_set('M87', 6, r, f'M87_b1_m6', colors='I_VI')
    for i, r in enumerate([44, 41, 42, 43]):
        add_set('M87', 6, r, f'M87_c1_m6', colors='I_VI')

        
def load_M87v2():
    print('M87v2')
    for m in range(1, 6):
        add_set('M87v2', m, 204, f'M87v2_a1_m{m}', colors='I_gI')
        add_set('M87v2', m, 128, f'M87v2_b1_m{m}', colors='I_gI')
        add_set('M87v2', m, 44, f'M87v2_c1_m{m}', colors='I_gI')
    for i, r in enumerate([204, 201, 202, 203]):
        add_set('M87v2', 6, r, f'M87v2_a1_m6', colors='I_gI')
    for i, r in enumerate([128, 125, 126, 127]):
        add_set('M87v2', 6, r, f'M87v2_b1_m6', colors='I_gI')
    for i, r in enumerate([44, 41, 42, 43]):
        add_set('M87v2', 6, r, f'M87v2_c1_m6', colors='I_gI')

        
def load_M49():
    print('M49')
    for m in range(1, 6):
        add_set('M49', m, 204, f'M49_a1_m{m}')
        add_set('M49', m, 124, f'M49_b1_m{m}')
        add_set('M49', m, 40, f'M49_c1_m{m}')
    for i, r in enumerate([204, 201, 202, 203]):
        add_set('M49', 6, r, f'M49_a1_m6')
    for i, r in enumerate([124, 121, 122, 123]):
        add_set('M49', 6, r, f'M49_b1_m6')
    for i, r in enumerate([40, 37, 38, 39]):
        add_set('M49', 6, r, f'M49_c1_m6')

        
def load_NGC3377():
    print('NGC3377')
    for m in range(1, 6):
        add_set('NGC3377', m, 173, f'NGC3377_a1_m{m}')
        add_set('NGC3377', m, 97, f'NGC3377_b1_m{m}')
        add_set('NGC3377', m, 41, f'NGC3377_c1_m{m}')
    for i, r in enumerate([173, 174, 175, 176]):
        add_set('NGC3377', 6, r, f'NGC3377_a1_m6')
    for i, r in enumerate([97, 98, 99, 100]):
        add_set('NGC3377', 6, r, f'NGC3377_b1_m6')
    for i, r in enumerate([41, 42, 43, 44]):
        add_set('NGC3377', 6, r, f'NGC3377_c1_m6')
    for m in range(7, 8):
        add_set('NGC3377', m, 173, f'NGC3377_a1_m{m}')
        add_set('NGC3377', m, 97, f'NGC3377_b1_m{m}')
        add_set('NGC3377', m, 41, f'NGC3377_c1_m{m}')

        
def load_NGC4993():
    print('NGC4993')
    for m in range(1, 6):
        add_set('NGC4993', m, 203, f'NGC4993_a1_m{m}')
        add_set('NGC4993', m, 143, f'NGC4993_b1_m{m}')
        add_set('NGC4993', m, 83, f'NGC4993_c1_m{m}')
    for m in [6]:
        add_set('NGC4993', m, 203, f'NGC4993_a1_m{m}', model=base_models[5])
        add_set('NGC4993', m, 143, f'NGC4993_b1_m{m}', model=base_models[5])
        add_set('NGC4993', m, 83, f'NGC4993_c1_m{m}', model=base_models[5])

        
def load_M31():
    print('M31')
    m31_regions = {
        'a': 'summer_r1',
        'b': 'winter_r9',
        'c': 'winter_r1',
        'd': 'winter_r2',
        'e': 'winter_r10'
    }
    m31_rpix = {
        'a': 6250,
        'b': 3550,
        'c': 2050,
        'd': 1350,
        'e': 700
    }
    for m in range(1, 7):
        for i, r in enumerate(['a', 'b', 'c', 'd', 'e']):
            k = f'M31_{r}_m{m}'
            add_set_v2('M31', m, r, f'M31_{r}_m{m}', 'm31_'+m31_regions[r])
            radii_am[k] = m31_rpix[r] * 0.05 / 60.
            radii_kpc[k] = radii_am[k] * (np.pi/(180.*60.)) * 1e3 * ppy.distancemodels.dmod_to_mpc(dmods['M31'])

            
def load_M51():
    print('M51')
    m51_rpix = {
        'a': 1745,
        'b': 2776,
        'c': 2471,
        'd': 2709,
        'e': 678
    }
    for m in range(1, 4):
        model = base_models[m].copy()
        model.dust_model = ppy.dustmodels.FixedWidthLogNormDust(0.1)
        for i, r in enumerate(['a', 'b', 'c', 'd', 'e']):
            k = f'M51_{r}_m{m}'
            add_set_v2('M51', m, r, f'M51_{r}_m{m}', f'M51_I_BI_{i+1}',
                       model=model)
            radii_am[k] = m51_rpix[r] * 0.05 / 60.
            radii_kpc[k] = radii_am[k] * (np.pi/(180.*60.)) * 1e3 * ppy.distancemodels.dmod_to_mpc(dmods['M51'])


load_M87()
load_M87v2()
load_M49()
load_NGC3377()
load_NGC4993()
load_M31()
load_M51()
