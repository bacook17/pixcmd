try:
    import pcmdpy_gpu as ppy
except:
    import pcmdpy as ppy

filters = {}
filters['M87'] = ppy.instrument.default_m87_filters()
filters['M49'] = ppy.instrument.default_m49_filters()
filters['NGC3377'] = ppy.instrument.default_ngc3377_filters()
filters['NGC4993'] = ppy.instrument.default_ngc4993_filters()
filters['DF2'] = ppy.instrument.default_df2_filters()

print('Loading Isochrone Models')
iso_models = {}
drivers = {}
for k, f in filters.items():
    iso_models[k] = ppy.isochrones.Isochrone_Model(f)
    drivers[k] = ppy.driver.Driver(iso_models[k])
