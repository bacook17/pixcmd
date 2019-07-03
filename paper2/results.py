__all__ = ['models', 'results', 'pcmds', 'data', 'radii_am', 'radii_kpc',
           'dmods', 'regions']

try:
    import pcmdpy_gpu as ppy
except:
    import pcmdpy as ppy
import numpy as np
import pandas as pd
from os.path import expanduser, isfile

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
    'M87': 30.9,
    'M87v2': 30.9,
    'M49': 31.12,
    'DF2': 31.505,
    'M31': 24.44,
    'M31d': 24.44,
    'M51': 29.67
}
m51_rpix = {
    'a': 1745,
    'b': 2776,
    'c': 2471,
    'd': 2709,
    'e': 678
}
m31_rpix = {
    'e': 6250,
    'd': 3550,
    'c': 2050,
    'b': 1350,
    'a': 700
}
m31d_rpix = {
    'a': 21397,
    'b': 23423,
    'c': 25217
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
base_models[7] = base_models[1].copy()
base_models[8] = base_models[5].copy()
base_models[9] = base_models[1].copy()
base_models[10] = base_models[2].copy()
base_models[11] = ppy.galaxy.CustomGalaxy(
    ppy.metalmodels.ClosedBoxMDF(),
    ppy.dustmodels.SingleDust(),
    ppy.sfhmodels.NonParam(),
    ppy.distancemodels.VariableDistance()
)
base_models[12] = base_models[4].copy()
base_models[13] = base_models[9].copy()
base_models[14] = base_models[9].copy()


def add_set(galaxy, mnum, region, key,
            colors='z_gz', model=None,
            dataname=None):
    g_orig = galaxy.replace('v2', '').replace('M31d', 'M31')
    g_dir = galaxy.replace('d','').lower()
    if dataname is not None:
        data_file = data_dir + f'{g_dir}/pcmds/{dataname}.pcmd'
    else:
        data_file = data_dir + f'{g_dir}/pcmds/{g_orig}_{colors}_{region}.pcmd'
    res_file = results_dir + f'{galaxy}_r{region}_m{mnum}.csv'
    live_file = res_file.replace('.csv', '_live.csv')
    pcmd_file = res_file.replace('.csv', '.pcmd')
    if not isfile(res_file):
        print(f'Skipping {key}')
        return
    regions[key] = region
    models[key] = model or base_models[mnum].copy()
    try:
        results[key] = ppy.results.ResultsPlotter(
            res_file, live_file=live_file, dmod_true=dmods[galaxy],
            gal_model=models[key], model_is_truth=False)
    except Exception as e:
        print('Error loading ', key)
        print(e)
        return
    ks = [d for d in data.keys() if key.replace(f'_m{mnum}', '') in d]
    if len(ks) > 0:
        data[key] = data[ks[0]]
    elif key not in data:
        data[key] = np.loadtxt(data_file, unpack=True)
    if g_orig in df_radii:
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
    try:
        results[key] = ppy.results.ResultsPlotter(
            res_file, live_file=live_file, dmod_true=dmods[galaxy],
            gal_model=models[key], model_is_truth=False)
    except Exception as e:
        print('Error loading ', key)
        print(e)
        return
    ks = [d for d in data.keys() if key.replace(f'_m{mnum}', '') in d]
    if len(ks) > 0:
        data[key] = data[ks[0]]
    else:
        data[key] = np.loadtxt(data_file, unpack=True)
    try:
        pcmds[key] = np.loadtxt(pcmd_file, unpack=True)
    except:
        pass


def load_model(m, all_quads=False):
    print('M87')
    for let, rs in zip(['a','b','c'], [[204,201,202,203], [128,125,126,127], [44,41,42,43]]):
        for i, r in enumerate(rs):
            if (i > 0) and not all_quads:
                continue
            add_set('M87', m, r, f'M87_{let}{i+1}_m{m}', colors='I_VI')

    print('M87v2')
    for let, rs in zip(['a','b','c'], [[204,201,202,203], [128,125,126,127], [44,41,42,43]]):
        for i, r in enumerate(rs):
            if (i > 0) and not all_quads:
                continue
            add_set('M87v2', m, r, f'M87v2_{let}{i+1}_m{m}', colors='I_gI')

    print('M49')
    for let, rs in zip(['a','b','c'], [[204,201,202,203], [124,121,122,123], [40,37,38,39]]):
        for i, r in enumerate(rs):
            if (i > 0) and not all_quads:
                continue
            add_set('M49', m, r, f'M49_{let}{i+1}_m{m}')

    print('NGC3377')
    for let, rs in zip(['a','b','c'], [[173,174,175,176], [97,98,99,100], [41,42,43,44]]):
        for i, r in enumerate(rs):
            if (i > 0) and not all_quads:
                continue
            if r == 175:
                continue
            add_set('NGC3377', m, r, f'NGC3377_{let}{i+1}_m{m}')

    print('NGC4993')
    for let, rs in zip(['a','b','c'], [[203,204,201,202], [143,144,141,142], [83,84,81,82]]):
        for i, r in enumerate(rs):
            if (i > 0) and not all_quads:
                continue
            if r == 201:
                continue
            model = (base_models[5].copy() if m==6 else None)
            add_set('NGC4993', m, r, f'NGC4993_{let}{i+1}_m{m}', model=model)

    print('M31 Bulge')
    for i, let in enumerate(['e','d','c','b','a']):
        k = f'M31_{let}_m{m}'
        add_set('M31', m, i+1, k,
                dataname=f'm31_bulge_r{i+1}')
        radii_am[k] = m31_rpix[let] * 0.05/60.
        radii_kpc[k] = radii_am[k] * (np.pi/(180.*60.)) * 1e3 * ppy.distancemodels.dmod_to_mpc(dmods['M31'])

    # if m in [7]:
    #     print('M31 Disk')
    #     for i, let in enumerate(['a','b','c']):
    #         k = f'M31d_{let}_m{m}'
    #         add_set('M31d', m, i+1, k,
    #                 dataname=f'm31_disk_r{i+1}')
    #         radii_am[k] = m31d_rpix[let] * 0.05/60.
    #         radii_kpc[k] = radii_am[k] * (np.pi/(180.*60.)) * 1e3 * ppy.distancemodels.dmod_to_mpc(dmods['M31'])

    # print('M51')
    # for i, let in enumerate(['a','b','c','d','e']):
    #     k = f'M51_{let}_m{m}'
    #     model = base_models[m].copy()
    #     model.dust_model = ppy.dustmodels.FixedWidthLogNormDust(0.1)
    #     add_set('M51', m, i+1, k, colors='I_BI', model=model)
    #     radii_am[k] = m51_rpix[let] * 0.05/60.
    #     radii_kpc[k] = radii_am[k] * (np.pi/(180.*60.)) * 1e3 * ppy.distancemodels.dmod_to_mpc(dmods['M51'])
    

#load_M87()
#load_M87v2()
#load_M49()
#load_NGC3377()
#load_NGC4993()
#load_M31()
#load_M51()
