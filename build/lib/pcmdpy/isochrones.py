# isochrones.py
# Ben Cook (bcook@cfa.harvard.edu)

"""Define the Isocrhone_Model class"""

import numpy as np
import pandas as pd

# The pre-computed MIST model metallicities
_z_arr_default = np.array([-2.15, -1.13, -0.73, -0.52, -0.22, 0., 0.3, 0.5])

class Isochrone_Model:
    """Models Isochrones (IMF, and magnitudes in particular Filters) using
       linear interpolation of MIST models

    An Isocrhone_Model incorporates a collection of MIST models, and
       allows for interpolating the IMF and magnitudes (for given Filter
       objects) at any arbitrary metallicity and mass 

    Attributes: 
       MIST_df-- A pandas Dataframe containing all pre-computed MIST datapoints
       ages -- An array of ages (in log years) which are valid for the model
    Methods:
       get_magnitudes -- Pass a Galaxy_Model object, return IMF and magnitudes for each mass, age, metallicity bin
    Constructors:
       __init__ -- Pass a list of Filter objects, path to MIST model files, and array of metallicities.
    """
    # CHANGE PATH 
    def __init__(self, filters, path='/n/home01/bcook/pixcmd/pcmdpy/isoc_csv/', z_arr=_z_arr_default):
        """Creates a new Isochrone_Model, given a list of Filter objects
        
        Arguments:
           filters -- list of Filter objects
        Keyword Arguments:
           path -- directory containing MIST model files
           z_arr -- array of MIST metallicity values to use
        """

        #Import all MIST model files into Pandas dataframe
        self.MIST_df = pd.DataFrame()
        self._z_arr = z_arr
        self.num_filters = len(filters)
        self.filters = filters
        self.filter_names = [f.tex_name for f in self.filters]
        #MIST files are organized by metallicity
        for z in self._z_arr:
            MIST_doc = path + 'MIST_v29_Z'+ self._z_to_str(z) + '_x5FEWER.csv'
            try:
                new_df = pd.read_csv(MIST_doc)
                new_df['z'] = z
                self.MIST_df = self.MIST_df.append([new_df], ignore_index=True)
            except IOError:
                raise IOError('No MIST file found for z=%.2f'%z)

        self.ages = self.MIST_df.age.unique()
            
        #The MIST columns that will be interpolated (mass, logIMF, and all input filters)
        self._interp_cols = ['logIMF']
        for f in filters:
            c = f.MIST_column
            if c in self.MIST_df.columns:
                self._interp_cols.append(f.MIST_column)
            else:
                raise ValueError('Filter input does not have a valid MIST_column')


            
    @staticmethod
    def _interp_arrays(arr1, arr2, f):
        """Linearly interpolate between two (potentially unequal length) arrays 
        
        Arguments:
           arr1 -- first (lower) array (len N1 or N1xD)
           arr2 -- second (upper) array (len N2 or N2xD, N2 doesn't have to equal N1)
           f -- linear interpolation fraction (float between 0 and 1)
        Output: interpolated array (len max(N1,N2) or max(N1,N2)xD)
        """
        assert(arr1.ndim == arr2.ndim)
        
        l1,l2 = len(arr1), len(arr2)
        #If arrays are unequal length, extrapolate shorter using trend of longer
        if (l1 < l2):
            delta = arr2[l1:] - arr2[l1-1]
            added = arr1[-1] + delta
            arr1 = np.append(arr1, added, axis=0)
        elif (l1 > l2):
            delta = arr1[l2:] - arr1[l2-1]
            added = arr2[-1] + delta
            arr2 = np.append(arr2, added, axis=0)
        return (1-f)*arr1 + f*arr2

    @staticmethod
    def _z_to_str(z):
        """Converts a metallicity value to MIST string
        Example Usage: 
           _z_to_str(-0.5313) -> "m0.53"
           _z_to_str(1.326)   -> "p1.33"

        Arguments: 
           z -- metallicity (float)
        Output: string representing metallicity
        """
        result = ''
        if (z < 0):
            result += 'm'
        else:
            result += 'p'
        result += '%1.2f'%(np.abs(z))
        return result

    """
    def _interp_df(self, df, col, z):
        if z in self._z_arr:
            return df[df.z == z][col]
        else:
            i = self._z_arr.searchsorted(z)
            if (i == 0):
                i = 1 #will extrapolate low
            elif (i == len(self._z_arr)):
                i = -1 #will extrapolate high  
            zlow, zhigh = self._z_arr[i-1:i+1] #bounding metallicities
            frac_between = (z - zlow) / (zhigh - zlow)
            if (frac_between >= 2) or (frac_between <= -1):
                raise ValueError('Extrapolating metallicity more than one entire metallicity bin')
            dflow, dfhigh = df[df.z == zlow][col], df[df.z == zhigh][col]
            return self._interp_arrays(dflow.values, dfhigh.values, frac_between).T
    """
    
    def get_isochrone(self, age, z, norm_IMF=True, rare_cut=0., **kwargs):
        """Interpolate MIST isochrones for given age and metallicity
        
        Arguments:
           age --- 
           z ---
        Output:
           IMF -- array of IMF points (len N)
           mags -- 2D array of magnitudes (DxN, where D is number of filters the model was initialized with)
        """
        #Find closest age in MIST database
        if age not in self.ages:
            age = self.ages[np.abs(self.ages - age).argmin()]
        this_age = self.MIST_df[self.MIST_df.age == age]
        #Output MIST values for known metallicities
        if z in self._z_arr:
            inter = this_age[this_age.z == z][self._interp_cols].values
        #Interpolate/extrapolate for other metallicities
        else:
            i = self._z_arr.searchsorted(z)
            if (i == 0):
                i = 1 #will extrapolate low
            elif (i == len(self._z_arr)):
                i = -1 #will extrapolate high  
            zlow, zhigh = self._z_arr[i-1:i+1] #bounding metallicities
            frac_between = (z - zlow) / (zhigh - zlow)
            if (frac_between >= 2) or (frac_between <= -1):
                raise ValueError('Extrapolating metallicity more than one entire metallicity bin')
            dflow, dfhigh = this_age[this_age.z == zlow][self._interp_cols], this_age[this_age.z == zhigh][self._interp_cols]
            inter = self._interp_arrays(dflow.values, dfhigh.values, frac_between)
            
        IMF = 10.**inter[:,0]
        if norm_IMF:
            IMF /= np.sum(IMF)
        mags = inter[:,1:].T

        #remove stars that are extremely rare
        to_keep = (IMF >= rare_cut)
        
        return IMF[to_keep], mags[:,to_keep]

    def model_galaxy(self, galaxy, **kwargs):
        weights = np.empty((1,0),dtype=float)
        mags = np.empty((self.num_filters, 0), dtype=float)
        for sfh, age in zip(galaxy.SFH, galaxy.ages):
            imf,m = self.get_isochrone(age, galaxy.z, **kwargs)
            weights = np.append(weights, imf*sfh)
            mags = np.append(mags, m, axis=-1)
        return weights, mags

    def plot_isochrone(self, galaxy, ax=None, **kwargs):
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots()
        for age in galaxy.ages:
            _, mags = self.get_isochrone(age, galaxy.z)
            ax.plot(mags[0]-mags[1], mags[1], 'k-', label='age: %d'%(age),**kwargs)
        names = self.filter_names
        ax.set_ylabel(names[1],fontsize='x-large')
        ax.set_xlabel('%s - %s'%(names[0], names[1]), fontsize='x-large')
        return ax
