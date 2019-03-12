import pcmdpy as ppy
import numpy as np
from astropy.coordinates import SkyCoord

models = {}
run_names = {}
results = {}
results_dir = '~/pCMDs/pixcmd/scripts_py/results/'

bulge_seasons = ['s1', 'w1', 'w2', 's8', 's9',
                 'w7', 'w8', 'w9', 'w10']
disk_seasons = ['w3', 'w4', 'w5', 's2', 's3',
                'w6', 's4', 's5', 's6', 's7']

m31_kpc = 7.66e2

coord_center = SkyCoord("0h42m44.341s", "+41d16m08.72s")
coords_disk = []
coords_disk.append(SkyCoord("0h44m14.282s", "+41d21m53.59s"))
coords_disk.append(SkyCoord("0h44m19.127s", "+41d24m02.26s"))
coords_disk.append(SkyCoord("0h44m32.039s", "+41d21m46.35s"))
coords_disk.append(SkyCoord("0h44m41.015s", "+41d20m41.62s"))
coords_disk.append(SkyCoord("0h44m46.214s", "+41d22m43.44s"))
coords_disk.append(SkyCoord("0h44m25.489s", "+41d19m32.61s"))
coords_disk.append(SkyCoord("0h44m53.254s", "+41d18m27.43s"))
coords_disk.append(SkyCoord("0h44m59.235s", "+41d20m40.15s"))
coords_disk.append(SkyCoord("0h44m56.598s", "+41d22m41.50s"))
coords_disk.append(SkyCoord("0h44m34.199s", "+41d20m31.13s"))

disk_arcmins = [coord_center.separation(c).arcmin for c in coords_disk]
disk_kpc = [coord_center.separation(c).radian*m31_kpc for c in coords_disk]

bulge_pix = [6225., 2050., 1350., 6225., 6225., 3550., 3550., 3550., 705.]
bulge_arcmins = [bp * 0.05 / 60. for bp in bulge_pix]
bulge_kpc = [am * (np.pi / (180.*60.)) * m31_kpc for am in bulge_arcmins]

base_models = {
    1: ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2),
        ppy.dustmodels.SingleDust(),
        ppy.sfhmodels.TauModel(),
        ppy.distancemodels.FixedDistance(24.42)),
    2: ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2),
        ppy.dustmodels.SingleDust(),
        ppy.sfhmodels.TauModel(),
        ppy.distancemodels.VariableDistance()),
    3: ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2),
        ppy.dustmodels.SingleDust(),
        ppy.sfhmodels.NonParam(),
        ppy.distancemodels.FixedDistance(24.42)),
    4: ppy.galaxy.CustomGalaxy(
        ppy.metalmodels.FixedWidthNormMDF(0.2),
        ppy.dustmodels.SingleDust(),
        ppy.sfhmodels.NonParam(),
        ppy.distancemodels.VariableDistance()),
}

for m in range(1, 5):
    for bulge_id, season_id in enumerate(bulge_seasons):
        key = f'b{bulge_id+1}-m{m}'
        models[key] = base_models[m].copy()
        run_names[key] = f'Bulge Region {bulge_id+1} (Model {m})'
        season = 'summer' if 's' in season_id else 'winter'
        r_id = int(season_id[1:])
        res_file = results_dir + f'm31_{season}_model{m}_r{r_id}.csv'
        live_file = res_file.replace('.csv', '_live.csv')
        results[key] = ppy.results.ResultsPlotter(
            res_file,
            live_file=live_file,
            run_name=run_names[key])
        results[key].r_kpc = bulge_kpc[bulge_id]
        results[key].r_arcmin = bulge_arcmins[bulge_id]
        results[key].season = season
    for disk_id, season_id in enumerate(disk_seasons):
        try:
            key = f'd{disk_id+1}-m{m}'
            models[key] = base_models[m].copy()
            run_names[key] = f'Disk Region {disk_id+1} (Model {m})'
            season = 'summer' if 's' in season_id else 'winter'
            r_id = int(season_id[1:])
            res_file = results_dir + f'm31_{season}_model{m}_r{r_id}.csv'
            live_file = res_file.replace('.csv', '_live.csv')
            results[key] = ppy.results.ResultsPlotter(
                res_file, live_file=live_file, run_name=run_names[key],
                gal_model=models[key], model_is_truth=False)
            results[key].r_kpc = disk_kpc[disk_id]
            results[key].r_arcmin = disk_arcmins[disk_id]
            results[key].season = season
        except FileNotFoundError:
            models.pop(key)
            run_names.pop(key)

assert np.all([k in models for k in run_names.keys()]), str([k for k in run_names.keys() if k not in models])
assert np.all([k in run_names for k in models.keys()]), str([k for k in models.keys() if k not in run_names])
assert np.all([k in results for k in models.keys()]), str([k for k in models.keys() if k not in results])
