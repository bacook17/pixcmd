# galaxy.py
# Ben Cook (bcook@cfa.harvard.edu)

"""Define the Galaxy_Model class"""

import numpy as np

class Galaxy_Model:

    def __init__(self, age_arr, log_weights, z, dust):
        self.ages = age_arr
        self._num_ages = len(self.ages)
        assert(len(age_weights) == self._num_ages)
        self.SFH = 10.**log_weights
        self.Npix = np.sum(self.SFH)
        self.z = z
        self.dust = dust
        
class Simple_Galaxy:

    def __init__(self,gal_params):
        assert(len(gal_params) == 4)
        self.Npix = 10.**gal_params[0]
        self.ages = [gal_params[1]]
        self._num_ages = len(self.ages)
        self.SFH = [self.Npix]
        self.z = gal_params[2]
        self.dust = gal_params[3]
