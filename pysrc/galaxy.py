# galaxy.py
# Ben Cook (bcook@cfa.harvard.edu)

"""Define the Galaxy_Model class"""

import numpy as np

class Galaxy_Model:

    age_arr = np.array([6.5, 7.5, 8.25, 8.75, 9.25, 9.75, 10.1])
    _num_ages = len(age_arr)
    _num_params = _num_ages + 2
    
    def __init__(self, gal_params):
        """
        gal_params:
           0 -- log (z / z_solar) metallicity
           1 -- log E(B-V) dust extinction
           2... -- log SFH in age bin
        """
        self.ages = self.age_arr
        assert(len(gal_params) == self._num_params)

        self._params = gal_params
        self.z = gal_params[0]
        self.dust = 10.**gal_params[1]
        self.SFH = 10.**gal_params[2:]
        self.Npix = np.sum(self.SFH)
        
class Galaxy_SSP:
    _num_params = 4

    def __init__(self, gal_params):
        """
        gal_params:
           0 -- log (z / z_solar) metallicity
           1 -- log E(B-V) dust extinction
           2 -- log Npix
           3 -- log age (in yrs)
        """
        assert(len(gal_params) == self._num_params)
        self._params = gal_params
        self.z = gal_params[0]
        self.dust = 10.**gal_params[1]
        Npix = 10.**gal_params[2]
        self.SFH = np.array([Npix])
        self.ages = np.array([gal_params[3]])
        self._num_ages = 1
    
