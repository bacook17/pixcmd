# galaxy.py
# Ben Cook (bcook@cfa.harvard.edu)

"""Define the Galaxy_Model class"""

import numpy as np

class Galaxy_Model:

    def __init__(self, Npix, age_arr, age_weights, z, dust):
        self.Npix = Npix
        self.ages = age_arr
        self._num_ages = len(self.ages)
        assert(len(age_weights) == self._num_ages)
        weights = age_weights/np.sum(age_weights)
        self.SFH = self.Npix * weights
        self.z = z
        self.dust = dust
        
