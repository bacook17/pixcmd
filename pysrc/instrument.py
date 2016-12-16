# instrument.py
# Ben Cook (bcook@cfa.harvard.edu)

"""Define classes for Filters and other similar objects"""

import numpy as np
from scipy.signal import fftconvolve

class Filter:
    """Models observations in a single band
    
    A Filter specifies the conversions between magnitudes and counts, and PSF convolution
    
    Attributes:
       name -- descriptive name of the filter (string, eg: F475W)
       tex_name -- LaTeX formatted name of the filter, eg for use in plotting (string, eg: r"g$_{475}$")
       MIST_column -- column name in MIST tables (string, eg: vmag, fuv, or f275)
       props -- dictionary of other properties
    Methods:
       mag_to_counts -- convert absolute magnitudes to photon counts
       counts_to_mag -- convert photon counts to absolute magnitudes
       psf_convolve -- convolve (2D) array with the instrumental PSF
    Constructors:
       __init__ -- default, manual entry of all parameters
       HST_F475W -- the Hubble F475W filter (only free parameter is distance)
       HST_F814W -- the Hubble F814W filter (only free parameter is distance)
    """

    
    def __init__(self, exposure, zero_point, d_mpc, psf, name="", tex_name="", MIST_column="", **kwargs):
        """Create a new Filter, given input properties of observation

        Arguments:
           exposure -- exposure time of the observation, in seconds (int or float)
           zero_point -- apparent magnitude corresponding to 1 count / second (int or float)
                   this value is affected by telescope aperture, sensitivity, etc.
           d_mpc -- the assumed distance to the source in Mpc (int or float)
           psf -- the PSF kernel, should be normalized to one (2D square array of floats)
        Keyword Argments:
           name -- descriptive name of the filter (string)
           tex_name -- LaTeX formatted name of the filter, eg for use in plotting (string, eg: r"g$_{475}$")
           MIST_column -- column name in MIST tables (string)
           **kwargs -- all other keyword arguments will be saved as a dictionary
        """

        #validate and initialize internal attributes
        try:
            self._exposure = float(exposure)
            self._zero_point = float(zero_point)
            self._dmod = 25. + 5.*np.log10(d_mpc) #distance modulus
        except TypeError:
            print('First three arguments must each be either a float or integer')
            raise
        if np.isnan(self._dmod):
            raise ValueError('The third argument (d_mpc) must be greater than zero')
        if not isinstance(psf, np.ndarray):
            psf = np.array(psf)
        if (psf.ndim != 2) or (psf.shape[0] != psf.shape[1]) or (psf.dtype != float):
            raise TypeError('The fourth argument (psf) must be a 2x2 square array of floats')
        else:
            self._psf = psf

        #initialize public attributes
        self.name = name
        self.tex_name = tex_name
        self.MIST_column = MIST_column
        self.props = kwargs

    ##############################
    # Alternative constructors

    @classmethod
    def HST_F475W(cls, d_mpc):
        """Return a Filter with HST F475W default params
        
        Example usage:
           >>>my_filter = Filter.HST_F475W(1.0)

        Arguments:
           d_mpc -- the assumed distance to the source
        Output: Filter with the following attributes:
           exposure = 3620.
           zero_point = 26.0593
           d_mpc : set by argument
           psf : loaded from file
           name = "F475W"
           tex_name = r"g$_{475}$"
           MIST_column = "bmag"
        """

        assert(isinstance(d_mpc, int) or isinstance(d_mpc, float)) #d_mpc must be real number
        if (d_mpc < 0.):
            raise ValueError('Argument (d_mpc) must be greater than zero')
        
        exposure = 3620.
        zero_point = 26.0593
        psf_file = "../psf/f475w_00.psf"
        psf = 10.**np.loadtxt(psf_file)
        name= "F475W"
        tex_name = r"g$_{475}$"
        MIST_column = "bmag"

        return cls(exposure, zero_point, d_mpc, psf, name=name, tex_name=tex_name, MIST_column=MIST_column)

    @classmethod
    def HST_F814W(cls, d_mpc):
        """Return a Filter with HST F814W default params

        Example usage:
           >>>my_filter = Filter.HST_F814W(1.0)

        Arguments:
           d_mpc -- the assumed distance to the source
        Output: Filter with the following attributes:
           exposure = 3235.
           zero_point = 25.9433
           d_mpc : set by argument
           psf : loaded from file
           name = "F814W"
           tex_name = r"I$_{814}$"
           MIST_column = "imag"
        """

        assert(isinstance(d_mpc, int) or isinstance(d_mpc, float)) #d_mpc must be real number
        if (d_mpc < 0.):
            raise ValueError('Argument (d_mpc) must be greater than zero')
        
        exposure = 3235.
        zero_point = 25.9433
        psf_file = "../psf/f814w_00.psf"
        psf = 10.**np.loadtxt(psf_file)
        name= "F814W"
        tex_name = r"I$_{814}$"
        MIST_column = "imag"

        return cls(exposure, zero_point, d_mpc, psf, name=name, tex_name=tex_name, MIST_column=MIST_column)

    #########################
    # Filter methods
    
    def mag_to_counts(self, mags):
        """Convert absolute magnitudes to photon counts

        Arguments:
           mags -- absolute magnitudes (int or float or array or ndarray)
        Output:
           counts -- photon counts (same type as input)
        """

        return 10.**(-0.4 * (mags + self._dmod - self._zero_point)) * self._exposure

    def counts_to_mag(self, counts):
        """Convert photon counts to absolute magnitudes

        Arguments:
           counts -- photon counts (int or float or array or ndarray)
        Output:
           mags -- absolute magnitudes (same type as input)
        """

        return -2.5*np.log10(counts / self._exposure) + self._zero_point - self._dmod

    def psf_convolve(self, image, convolve_func=None, **kwargs):
        """Convolve image with instrumental PSF
        
        Arguments:
           image -- counts, or flux, in each pixel of image (2D array of integers or floats)
        Keyword Arguments:
           convolve_func -- function to convolve the image and PSF (default: scipy.signal.fftconvolve)
           **kwargs -- any additional keyword arguments will be passed to convolve_func
        Output:
           convolved_image -- image convolved with PSF (2D array of floats;
                                       guaranteed same shape as input if default convolve_func used)
        """

        if not isinstance(image, np.ndarray):
            image = np.array(image)
        if (image.ndim != 2):
            raise TypeError('The first argument (image) must be a 2D array of integers or floats')    

        if convolve_func is None:
            return fftconvolve(image, self._psf, mode='valid')
        else:
            return convolve_func(image, self._psf, **kwargs)

