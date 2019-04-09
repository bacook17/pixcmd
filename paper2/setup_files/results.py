__all__ = ['models', 'results', 'pcmds', 'data']

try:
    import pcmdpy_gpu as ppy
except:
    import pcmdpy as ppy
import numpy as np
from os.path import expanduser

models = {}
run_names = {}
results = {}
pcmds = {}
data = {}
results_dir = expanduser('~/pCMDs/pixcmd/paper2/results/')
data_dir = expanduser('~/pCMDs/pixcmd/data/')

model_nonparam = ppy.galaxy.CustomGalaxy(
    ppy.metalmodels.SingleFeH(),
    ppy.dustmodels.SingleDust(),
    ppy.sfhmodels.NonParam(),
    ppy.distancemodels.VariableDistance()
)
model_fixeddist = ppy.galaxy.CustomGalaxy(
    ppy.metalmodels.SingleFeH(),
    ppy.dustmodels.SingleDust(),
    ppy.sfhmodels.NonParam(),
    ppy.distancemodels.FixedDistance()
)
model_tau = ppy.galaxy.CustomGalaxy(
    ppy.metalmodels.SingleFeH(),
    ppy.dustmodels.SingleDust(),
    ppy.sfhmodels.TauModel(),
    ppy.distancemodels.VariableDistance()
)
model_ssp = ppy.galaxy.CustomGalaxy(
    ppy.metalmodels.SingleFeH(),
    ppy.dustmodels.SingleDust(),
    ppy.sfhmodels.SSPModel(),
    ppy.distancemodels.VariableDistance()
)
model_ssp_mdf = ppy.galaxy.CustomGalaxy(
    ppy.metalmodels.FixedWidthNormMDF(0.2),
    ppy.dustmodels.SingleDust(),
    ppy.sfhmodels.SSPModel(),
    ppy.distancemodels.VariableDistance()
)
model_ssp_fixed = ppy.galaxy.CustomGalaxy(
    ppy.metalmodels.SingleFeH(),
    ppy.dustmodels.SingleDust(),
    ppy.sfhmodels.SSPModel(),
    ppy.distancemodels.FixedDistance()
)
custom_sfh = ppy.sfhmodels.NonParam()
custom_sfh.update_sfh_edges(np.array([9.5, 9.75, 10.0, 10.2]))
custom_sfh.update_edges(np.arange(9.0, 10.3, 0.1))
model_df2_nonparam = ppy.galaxy.CustomGalaxy(
    ppy.metalmodels.SingleFeH(),
    ppy.dustmodels.SingleDust(),
    custom_sfh,
    ppy.distancemodels.VariableDistance()
)


def add_set(galaxy, mnum, region, key, model=model_nonparam,
            colors='z_gz', run_name=None):
    data_file = data_dir + f'{galaxy.lower()}/pcmds/{galaxy}_{colors}_{region}.pcmd'
    run_names[key] = run_name
    res_file = results_dir + f'{galaxy}_m{mnum}_r{region}.csv'
    live_file = res_file.replace('.csv', '_live.csv')
    pcmd_file = res_file.replace('.csv', '.pcmd')
    models[key] = model.copy()
    results[key] = ppy.results.ResultsPlotter(
        res_file, live_file=live_file, run_name=run_name,
        gal_model=models[key], model_is_truth=False)
    data[key] = np.loadtxt(data_file, unpack=True)
    try:
        pcmds[key] = np.loadtxt(pcmd_file, unpack=True)
    except:
        pass

    
# M87
print('M87')
add_set('M87', 3, 44, 'M87_m3', colors='I_VI')
add_set('M87', 4, 104, 'M87_m4', colors='I_VI')
add_set('M87', 4, 101, 'M87_m4_q1', colors='I_VI')
add_set('M87', 4, 102, 'M87_m4_q2', colors='I_VI')
add_set('M87', 4, 103, 'M87_m4_q3', colors='I_VI')
add_set('M87', 5, 204, 'M87_m5', colors='I_VI')
add_set('M87', 6, 264, 'M87_m6', colors='I_VI')
add_set('M87', 7, 104, 'M87_m7', model=model_fixeddist, colors='I_VI')
add_set('M87', 8, 104, 'M87_m8', model=model_tau, colors='I_VI')
add_set('M87', 9, 104, 'M87_m9', model=model_ssp, colors='I_VI')
add_set('M87', 10, 104, 'M87_m10', model=model_ssp, colors='I_VI')
add_set('M87', 11, 104, 'M87_m11', model=model_ssp, colors='I_VI')
add_set('M87', 12, 104, 'M87_m12', model=model_ssp_fixed, colors='I_VI')
add_set('M87', 13, 104, 'M87_m13', model=model_ssp, colors='I_VI')
add_set('M87', 14, 104, 'M87_m14', model=model_ssp, colors='I_VI')
add_set('M87', 15, 104, 'M87_m15', model=model_ssp_mdf, colors='I_VI')
add_set('M87', 16, 104, 'M87_m16', model=model_ssp, colors='I_VI')
add_set('M87', 17, 104, 'M87_m17', model=model_ssp, colors='I_VI')

# M49
print('M49')
add_set('M49', 3, 40, 'M49_m3')
add_set('M49', 4, 100, 'M49_m4')
add_set('M49', 4, 97, 'M49_m4_q1')
add_set('M49', 4, 98, 'M49_m4_q2')
add_set('M49', 4, 99, 'M49_m4_q3')
add_set('M49', 5, 204, 'M49_m5')
add_set('M49', 6, 256, 'M49_m6')
add_set('M49', 7, 100, 'M49_m7', model=model_fixeddist)
add_set('M49', 8, 100, 'M49_m8', model=model_tau)
add_set('M49', 9, 100, 'M49_m9', model=model_ssp)
add_set('M49', 10, 100, 'M49_m10', model=model_ssp)

# NGC 3377
print('NGC3377')
add_set('NGC3377', 3, 41, 'NGC3377_m3')
add_set('NGC3377', 4, 97, 'NGC3377_m4')
add_set('NGC3377', 4, 98, 'NGC3377_m4_q1')
add_set('NGC3377', 4, 99, 'NGC3377_m4_q2')
add_set('NGC3377', 4, 100, 'NGC3377_m4_q3')
add_set('NGC3377', 5, 173, 'NGC3377_m5')
add_set('NGC3377', 6, 241, 'NGC3377_m6')
add_set('NGC3377', 7, 97, 'NGC3377_m7', model=model_fixeddist)
add_set('NGC3377', 8, 97, 'NGC3377_m8', model=model_tau)
add_set('NGC3377', 9, 97, 'NGC3377_m9', model=model_ssp)
add_set('NGC3377', 10, 97, 'NGC3377_m10', model=model_ssp)

# NGC 4993
print('NGC4993')
add_set('NGC4993', 3, 35, 'NGC4993_m3')
add_set('NGC4993', 4, 83, 'NGC4993_m4')
add_set('NGC4993', 4, 81, 'NGC4993_m4_q1')
add_set('NGC4993', 4, 82, 'NGC4993_m4_q2')
add_set('NGC4993', 4, 84, 'NGC4993_m4_q3')
add_set('NGC4993', 5, 103, 'NGC4993_m5')
# add_set('NGC4993', 6, 241, 'NGC4993_m6')
add_set('NGC4993', 7, 83, 'NGC4993_m7', model=model_fixeddist)
add_set('NGC4993', 8, 83, 'NGC4993_m8', model=model_tau)
add_set('NGC4993', 9, 83, 'NGC4993_m9', model=model_ssp)

# DF2
print('DF2')
for i in range(1, 5):
    df2_res = results_dir + f'DF2_m{i}.csv'
    df2_live = df2_res.replace('.csv', '_live.csv')
    df2_data = data_dir + 'DF2/pcmds/DF2_I_VI_1.pcmd'
    if i in [2, 4]:
        model = model_df2_nonparam.copy()
    else:
        model = model_ssp.copy()
    results[f'DF2_m{i}'] = ppy.results.ResultsPlotter(
        df2_res, live_file=df2_live, run_name=f'DF2, model {i}',
        gal_model=model, model_is_truth=False)
    data[f'DF2_m{i}'] = np.loadtxt(df2_data, unpack=True)
    try:
        pcmds[f'DF2_m{i}'] = np.loadtxt(df2_res.replace('.csv', '.pcmd'), unpack=True)
    except:
        pass
