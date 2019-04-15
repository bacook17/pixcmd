try:
    import pcmdpy_gpu as ppy
except:
    import pcmdpy as ppy

filters = {}
filters['M87'] = ppy.instrument.default_m87_filters()
filters['M87v2'] = ppy.instrument.m87_filters_v2()
filters['M49'] = ppy.instrument.default_m49_filters()
filters['NGC3377'] = ppy.instrument.default_ngc3377_filters()
filters['NGC4993'] = ppy.instrument.default_ngc4993_filters()
filters['DF2'] = ppy.instrument.default_df2_filters()
filters['M31'] = ppy.instrument.m31_winter_filters()
filters['M51'] = ppy.instrument.default_m51_filters()

sky = {
    'M87': [292.8, 377.1],
    'M87v2': [292.8, 43.4],
    'M49': [56.9, 50.1],
    'M31': [178.1, 103.5],
    'M51': [72.9, 40.5],
    'NGC3377': [126.9, 76.5],
    'NGC4993': [62.0, 175.3],
    'DF2': [163.7, 207.1]
}

print('Loading Isochrone Models')
iso_models = {
    'M87': ppy.isochrones.Isochrone_Model(filters['M87']),  # 814, 606
    'NGC3377': ppy.isochrones.Isochrone_Model(filters['NGC3377']),  # 850, 475
    'M51': ppy.isochrones.Isochrone_Model(filters['M51']),  # 814, 435
    'M31': ppy.isochrones.Isochrone_Model(filters['M31']),  # 814, 475
}
iso_models['M87v2'] = iso_models['M31']
iso_models['NGC4993'] = iso_models['M49'] = iso_models['NGC3377']
iso_models['DF2'] = iso_models['M87']
drivers = {}
for k, i in iso_models.items():
    drivers[k] = ppy.driver.Driver(i)
    drivers[k].filters = filters[k]

