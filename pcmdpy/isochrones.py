# isochrones.py
# Ben Cook (bcook@cfa.harvard.edu)

"""Define the Isocrhone_Model class"""

import numpy as np
import pandas as pd
try:
    from pkg_resources import resource_filename
except ImportError:
    pass

#-----------------
# Useful Utilities
#-----------------

def salpeter_IMF(mass, cutoff=0.08, normed=True):
    imf = np.power(mass, -2.35)
    imf[mass < cutoff] = 0.
    if normed:
        imf /= np.sum(imf)
    return imf

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


# The pre-computed MIST model metallicities
_z_arr_default = np.array([-4.0, -3.5, -3.0, -2.5, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5])
#_z_arr_default = np.array([-2.15, -1.13, -0.73, -0.52, -0.22, 0., 0.3, 0.5])

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
    def __init__(self, filters, MIST_path=None, z_arr=_z_arr_default):
        """Creates a new Isochrone_Model, given a list of Filter objects
        
        Arguments:
           filters -- list of Filter objects
        Keyword Arguments:
           MIST_path -- directory containing MIST model files
           z_arr -- array of MIST metallicity values to use
        """

        #Locate MIST files
        if MIST_path is None:
            try:
                MIST_path = resource_filename('pcmdpy', 'isoc_MIST_v1.0/')
            except:
                MIST_path = '/n/home01/bcook/pixcmd/pcmdpy/isoc_MIST_v1.0/'
        
        #Import all MIST model files into Pandas dataframe
        self.MIST_df = pd.DataFrame()
        self._z_arr = z_arr
        self.num_filters = len(filters)
        self.filters = filters
        self.filter_names = [f.tex_name for f in self.filters]
        self.colnames = pd.read_table(MIST_path + 'columns.txt', delim_whitespace=True).columns
        #MIST files are organized by metallicity
        for z in self._z_arr:
            MIST_doc = MIST_path + 'MIST_v1.0_feh_'+ _z_to_str(z) + '_afe_p0.0_vvcrit0.0_HST_ACSWF.iso.cmd'
            try:
                new_df = pd.read_table(MIST_doc, names=self.colnames, comment='#', delim_whitespace=True)
                new_df['z'] = z
                self.MIST_df = self.MIST_df.append([new_df], ignore_index=True)
            except IOError:
                raise IOError('No MIST file found for z=%.2f, tried: %s'%(z, MIST_doc))

        self.MIST_df.rename(columns={'log10_isochrone_age_yr':'age'}, inplace=True)
        self.ages = self.MIST_df.age.unique()
        #The MIST columns that will be interpolated (mass, logIMF, and all input filters)
        self._interp_cols = ['initial_mass']
        for f in self.filters:
            c = f.MIST_column
            if c in self.MIST_df.columns:
                self._interp_cols.append(f.MIST_column)
            else:
                raise ValueError('Filter input does not have a valid MIST_column')

    
    def get_isochrone(self, age, z, imf_func=salpeter_IMF, rare_cut=0., **kwargs):
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
            inter = _interp_arrays(dflow.values, dfhigh.values, frac_between)
            
        IMF = imf_func(inter[:,0], **kwargs)
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


