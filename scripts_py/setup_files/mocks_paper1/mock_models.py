import pcmdpy as ppy
import numpy as np

models = {}
run_names = {}
results = {}
results_dir = '~/pCMDs/pixcmd/scripts_py/results/'

models['mock_1'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.25),
        ppy.dustmodels.SingleDust(-0.5),
        ppy.sfhmodels.TauModel(np.array([2.0, 3.0])),
        ppy.distancemodels.FixedDistance(26.0)
    )
run_names['mock_1'] = 'Tau Model 1 (Distance Fixed)'

models['mock_2'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.25),
        ppy.dustmodels.SingleDust(-0.5),
        ppy.sfhmodels.TauModel(np.array([4.5, 1.2])),
        ppy.distancemodels.FixedDistance(26.0)
    )
run_names['mock_2'] = 'Tau Model 2 (Distance Fixed)'

models['mock_3'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.25),
        ppy.dustmodels.SingleDust(-0.5),
        ppy.sfhmodels.TauModel(np.array([2.0, 3.0])),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_3'] = 'Tau Model 1 (Distance Free)'

models['mock_4'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.25),
        ppy.dustmodels.SingleDust(-0.5),
        ppy.sfhmodels.TauModel(np.array([4.5, 1.2])),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_4'] = 'Tau Model 2 (Distance Free)'

models['mock_5'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([2.0, 3.0])),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_5'] = 'Tau-MDF Model (Distance Free)'

models['mock_6'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.NormMDF([-0.25, 0.2]),
        ppy.dustmodels.LogNormDust([-0.5, 0.1]),
        ppy.sfhmodels.TauModel(np.array([2.0, 3.0])),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_6'] = r'Tau-MDF Model (Distance Free, $\sigma$ Free)'

tau_SFHbins = ppy.sfhmodels.TauModel(np.array([2.0, 3.0]), iso_step=-1).logSFH

models['mock_7'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.FixedDistance(26.0)
    )
run_names['mock_7'] = r'NonParam Model (Distance Fixed)'

models['mock_8'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_8'] = r'NonParam Model (Distance Free)'

models['mock_9'] = models['mock_1'].copy()
run_names['mock_9'] = 'Tau Model 1 (Distance Fixed, LikeMode=3)'

models['mock_10'] = models['mock_5'].copy()
run_names['mock_10'] = 'Tau-MDF Model (Distance Free, Nim=256)'

models['mock_11'] = models['mock_5'].copy()
run_names['mock_11'] = 'Tau-MDF Model (Distance Free, Nim=1024)'

models['mock_12'] = models['mock_2'].copy()
run_names['mock_12'] = 'Tau Model 2 (Distance Fixed, Nim=1024)'

models['mock_13'] = models['mock_4'].copy()
run_names['mock_13'] = 'Tau Model 2 (Distance Free, Nim=1024)'

models['mock_14'] = models['mock_2'].copy()
run_names['mock_14'] = 'Tau Model 2 (Distance Fixed, more live points)'

models['mock_15'] = models['mock_4'].copy()
run_names['mock_15'] = 'Tau Model 2 (Distance Free, more live points)'

models['mock_16'] = models['mock_2'].copy()
run_names['mock_16'] = 'Tau Model 2 (Distance Fixed, Narrow Npix prior)'

models['mock_17'] = models['mock_4'].copy()
run_names['mock_17'] = 'Tau Model 2 (Distance Free, Narrow Npix prior)'


# initial conditions, same realization
models['mock_18'] = models['mock_5'].copy()
run_names['mock_18'] = 'Tau-MDF Model (Distance Fixed, Dynesty Seed 2)'
models['mock_19'] = models['mock_5'].copy()
run_names['mock_19'] = 'Tau-MDF Model (Distance Fixed, Dynesty Seed 3)'
models['mock_20'] = models['mock_5'].copy()
run_names['mock_20'] = 'Tau-MDF Model (Distance Fixed, Dynesty Seed 4)'
models['mock_21'] = models['mock_5'].copy()
run_names['mock_21'] = 'Tau-MDF Model (Distance Fixed, Dynesty Seed 5)'
models['mock_22'] = models['mock_5'].copy()
run_names['mock_22'] = 'Tau-MDF Model (Distance Fixed, Dynesty Seed 6)'
models['mock_23'] = models['mock_5'].copy()
run_names['mock_23'] = 'Tau-MDF Model (Distance Fixed, Dynesty Seed 7)'
models['mock_24'] = models['mock_5'].copy()
run_names['mock_24'] = 'Tau-MDF Model (Distance Fixed, Dynesty Seed 8)'

# Npix
models['mock_25'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([3.0, 3.0])),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_25'] = 'Tau-MDF Model (Distance Free, Npix=1e3)'

models['mock_26'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([4.0, 3.0])),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_26'] = 'Tau-MDF Model (Distance Free, Npix=1e4)'

models['mock_27'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([5.0, 3.0])),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_27'] = 'Tau-MDF Model (Distance Free, Npix=1e5)'

models['mock_28'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([6.0, 3.0])),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_28'] = 'Tau-MDF Model (Distance Free, Npix=1e6)'

# Hess Binning
models['mock_29'] = models['mock_5'].copy()
run_names['mock_29'] = 'Tau-MDF Model (Distance Free, Wide Hess Bins)'

models['mock_30'] = models['mock_5'].copy()
run_names['mock_30'] = 'Tau-MDF Model (Distance Free, Narrow Hess Bins)'

models['mock_31'] = models['mock_8'].copy()
run_names['mock_31'] = 'NonParam Model (Distance Free, Wide Hess Bins)'

models['mock_32'] = models['mock_8'].copy()
run_names['mock_32'] = 'NonParam Model (Distance Free, Narrow Hess Bins)'

# Fixed Distance Model 5
models['mock_33'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([2.0, 3.0])),
        ppy.distancemodels.FixedDistance(26.0)
    )
run_names['mock_33'] = 'Tau-MDF Model (Distance Fixed)'

# Different realizations
models['mock_34'] = models['mock_5'].copy()
run_names['mock_34'] = 'Tau-MDF Model (Distance Free, Realization 2)'

models['mock_35'] = models['mock_5'].copy()
run_names['mock_35'] = 'Tau-MDF Model (Distance Free, Realization 3)'

models['mock_36'] = models['mock_5'].copy()
run_names['mock_36'] = 'Tau-MDF Model (Distance Free, Realization 4)'

models['mock_37'] = models['mock_5'].copy()
run_names['mock_37'] = 'Tau-MDF Model (Distance Free, Realization 5)'

models['mock_38'] = models['mock_5'].copy()
run_names['mock_38'] = 'Tau-MDF Model (Distance Free, Realization 6)'

models['mock_39'] = models['mock_5'].copy()
run_names['mock_39'] = 'Tau-MDF Model (Distance Free, Realization 7)'

models['mock_40'] = models['mock_5'].copy()
run_names['mock_40'] = 'Tau-MDF Model (Distance Free, Realization 8)'

models['mock_41'] = models['mock_5'].copy()
run_names['mock_41'] = 'Tau-MDF Model (Distance Free, Nim = 2048, Npix=1e2)'

models['mock_42'] = models['mock_25'].copy()
run_names['mock_42'] = 'Tau-MDF Model (Distance Free, Nim = 648, Npix=1e3)'

models['mock_43'] = models['mock_26'].copy()
run_names['mock_43'] = 'Tau-MDF Model (Distance Free, Nim = 205, Npix=1e4)'

models['mock_44'] = models['mock_27'].copy()
run_names['mock_44'] = 'Tau-MDF Model (Distance Free, Nim = 65, Npix=1e5)'

models['mock_45'] = models['mock_28'].copy()
run_names['mock_45'] = 'Tau-MDF Model (Distance Free, Nim = 21, Npix=1e6)'
 
tau_SFHbins = ppy.sfhmodels.TauModel(np.array([3.0, 3.0]), iso_step=-1).logSFH
models['mock_46'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_46'] = r'NonParam Model (Distance Free, Npix=1e3)'

tau_SFHbins = ppy.sfhmodels.TauModel(np.array([4.0, 3.0]), iso_step=-1).logSFH
models['mock_47'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_47'] = r'NonParam Model (Distance Free, Npix=1e4)'

tau_SFHbins = ppy.sfhmodels.TauModel(np.array([5.0, 3.0]), iso_step=-1).logSFH
models['mock_48'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_48'] = r'NonParam Model (Distance Free, Npix=1e5)'

tau_SFHbins = ppy.sfhmodels.TauModel(np.array([6.0, 3.0]), iso_step=-1).logSFH
models['mock_49'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_49'] = r'NonParam Model (Distance Free, Npix=1e6)'

tau_SFHbins = ppy.sfhmodels.ConstantSFR(np.array([2.0]), iso_step=-1).logSFH
models['mock_50'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_50'] = r'NonParam Model (Distance Free, Npix=1e2, Constant SFR)'

tau_SFHbins = ppy.sfhmodels.ConstantSFR(np.array([3.0]), iso_step=-1).logSFH
models['mock_51'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_51'] = r'NonParam Model (Distance Free, Npix=1e3, Constant SFR)'

tau_SFHbins = ppy.sfhmodels.ConstantSFR(np.array([4.0]), iso_step=-1).logSFH
models['mock_52'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_52'] = r'NonParam Model (Distance Free, Npix=1e4, Constant SFR)'

tau_SFHbins = ppy.sfhmodels.ConstantSFR(np.array([5.0]), iso_step=-1).logSFH
models['mock_53'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_53'] = r'NonParam Model (Distance Free, Npix=1e5, Constant SFR)'

tau_SFHbins = ppy.sfhmodels.ConstantSFR(np.array([6.0]), iso_step=-1).logSFH
models['mock_54'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.NonParam(tau_SFHbins),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mock_54'] = r'NonParam Model (Distance Free, Npix=1e6, Constant SFR)'

models['mock_55'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([2.0, 3.0])),
        ppy.distancemodels.VariableDistance(25.0)
    )
run_names['mock_55'] = 'Tau-MDF Model (Distance Free, Nim = 2048, Npix=1e2, Dmod=25)'

models['mock_56'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([3.0, 3.0])),
        ppy.distancemodels.VariableDistance(27.5)
    )
run_names['mock_56'] = 'Tau-MDF Model (Distance Free, Nim = 658, Npix=1e3, Dmod=27.5)'

models['mock_57'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([4.0, 3.0])),
        ppy.distancemodels.VariableDistance(30.0)
    )
run_names['mock_57'] = 'Tau-MDF Model (Distance Free, Nim = 205, Npix=1e4, Dmod=30)'

models['mock_58'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([5.0, 3.0])),
        ppy.distancemodels.VariableDistance(32.5)
    )
run_names['mock_58'] = 'Tau-MDF Model (Distance Free, Nim = 66, Npix=1e5, Dmod=32.5)'

models['mock_59'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([6.0, 3.0])),
        ppy.distancemodels.VariableDistance(35.0)
    )
run_names['mock_59'] = 'Tau-MDF Model (Distance Free, Nim = 21, Npix=1e6, Dmod=35)'

models['mock_60'] = models['mock_5'].copy()
run_names['mock_60'] = 'Tau-MDF Model (Poisson Likelihood)'

#######

models['mismatch_1'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),  # model SingleFeH with FixedWidthNormMDF
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([2.0, 3.0])),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mismatch_1'] = 'Model: MDF, Truth: Single [Fe/H]'
 
models['mismatch_2'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.SingleFeH(-0.25),  # model FixedWidthNormMDF with SingleFeH
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([2.0, 3.0])),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mismatch_2'] = 'Model: Single [Fe/H], Truth: MDF'

models['mismatch_3'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),  # model SingleDust (Fdust=0.5) with FixedWidthLogNormDust
        ppy.sfhmodels.TauModel(np.array([2.0, 3.0])),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mismatch_3'] = 'Model: LogNorm Dust, Truth: Single Dust'

models['mismatch_4'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.SingleDust(-0.5, dust_frac=0.5),  # model FixedWidthLogNormDust with SingleDust (Fdust=0.5)
        ppy.sfhmodels.TauModel(np.array([2.0, 3.0])),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mismatch_4'] = 'Model: Single Dust, Truth: LogNorm Dust'

models['mismatch_5'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5, dust_frac=1.0),  # model Fdust=0.5 with Fdust=1.0
        ppy.sfhmodels.TauModel(np.array([2.0, 3.0])),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mismatch_5'] = 'Model: DustFrac = 1, Truth: DustFrac = 0.5'

models['mismatch_6'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2, -0.25),
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([2.0, 3.0])),
        ppy.distancemodels.FixedDistance(28.0)  # model has wrong distance
    )
run_names['mismatch_6'] = 'Model: Dmod=28, Truth: Dmod=26'

models['mismatch_7'] = models['mock_5'].copy()  # Exposure overestimated in model
run_names['mismatch_7'] = 'Model: Overestimates Exposure Time by 5x'

models['mismatch_8'] = models['mock_5'].copy()  # mock has PSF 10% narrower
run_names['mismatch_8'] = 'Model: Overestimates PSF width by 10%'

models['mismatch_9'] = models['mock_5'].copy()  # mock has PSF 10% narrower in F814W
run_names['mismatch_9'] = 'Model: Overestimates F814W PSF width by 10%'

models['mismatch_10'] = ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.3, -0.25),  # model sig=0.1 with sig=0.3
        ppy.dustmodels.FixedWidthLogNormDust(0.1, -0.5),
        ppy.sfhmodels.TauModel(np.array([2.0, 3.0])),
        ppy.distancemodels.VariableDistance(26.0)
    )
run_names['mismatch_10'] = 'Model: Overestimates MDF width (0.3 vs 0.1)'

models['mismatch_11'] = models['mock_5'].copy()  # Mock has no downsample
run_names['mismatch_11'] = 'Mock: no isochrone downsampling'

for key in models.keys():
    res_file = results_dir + 'paper1_' + key + '.csv'
    live_file = res_file.replace('.csv', '_live.csv')
    try:
        results[key] = ppy.results.ResultsPlotter(
            res_file, live_file=live_file, run_name=run_names[key],
            gal_model=models[key], model_is_truth=('mock in key'))
    except (FileNotFoundError, ValueError, AttributeError):
        pass
    except Exception as e:
        print(key, e)
        raise e
        
assert np.all([k in models for k in run_names.keys()]), str([k for k in run_names.keys() if k not in models])
assert np.all([k in run_names for k in models.keys()]), str([k for k in models.keys() if k not in run_names])
