            # instrument.py
# Ben Cook (bcook@cfa.harvard.edu)

"""Define classes for Filters and other similar objects"""

import numpy as np
import utils
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

    
    def __init__(self, exposure, zero_point, d_mpc, red_per_ebv, psf,  name="", tex_name="", MIST_column="", **kwargs):
        """Create a new Filter, given input properties of observation

        Arguments:
           exposure -- exposure time of the observation, in seconds (int or float)
           zero_point -- apparent magnitude corresponding to 1 count / second (int or float)
                   this value is affected by telescope aperture, sensitivity, etc.
           d_mpc -- the assumed distance to the source in Mpc (int or float)
           red_per_ebv -- the Reddening value [A_x / E(B-V)], such as from Schlafly & Finkbeiner 2011, Table 6 (float)
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
            self._red_per_ebv = float(red_per_ebv)
        except TypeError:
            print('First four arguments must each be either a float or integer')
            raise
        if np.isnan(self._dmod):
            raise ValueError('The third argument (d_mpc) must be greater than zero')
        if not isinstance(psf, np.ndarray):
            psf = np.array(psf)
        if (psf.shape[-2] != psf.shape[-1]) or (psf.dtype != float):
            raise TypeError('The fifth argument (psf) must be a square array (or 2D-array of square arrays) of floats')
        else:
            try:
                assert((psf.ndim == 2) or (psf.ndim == 4))
            except:
                raise TypeError('The fifth argument (psf) must be 2 or 4-dimensional (square array, or 2D-array of square arrays)')
            if (psf.ndim == 2):
                psf /= np.sum(psf)
            else:
                psf = np.array([[psf[i,j] / np.sum(psf[i,j]) for j in range(psf.shape[1])] for i in range(psf.shape[0])])
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
           red_per_ebv = 3.248 (Schlafly & Finkbeiner 2011, Table 6)
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
        red_per_ebv = 3.248
        psf_file = "../psf/f475w_%d%d.psf"
        psf = np.array([[10.**np.loadtxt(psf_file%(i,j)) for i in range(0,4)] for j in range(0,4)]) #4x4x73x73
        name= "F475W"
        tex_name = r"g$_{475}$"
        MIST_column = "bmag"

        return cls(exposure, zero_point, d_mpc, red_per_ebv, psf, name=name, tex_name=tex_name, MIST_column=MIST_column)

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
           red_per_ebv = 1.536 (Schlafly & Finkbeiner 2011, Table 6)
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
        red_per_ebv = 1.536
        psf_file = "../psf/f814w_%d%d.psf"
        psf = np.array([[10.**np.loadtxt(psf_file%(i,j)) for i in range(0,4)] for j in range(0,4)]) #4x4x73x73
        name= "F814W"
        tex_name = r"I$_{814}$"
        MIST_column = "imag"

        return cls(exposure, zero_point, d_mpc, red_per_ebv, psf, name=name, tex_name=tex_name, MIST_column=MIST_column)

    #########################
    # Filter methods
    
    def mag_to_counts(self, mags):
        """Convert absolute magnitudes to photon counts (no reddening assumed)

        Arguments:
           mags -- absolute magnitudes (int or float or array or ndarray)
        Output:
           counts -- photon counts (same type as input)
        """

        return 10.**(-0.4 * (mags + self._dmod - self._zero_point)) * self._exposure

    def counts_to_mag(self, counts, E_BV=0):
        """Convert photon counts to absolute magnitudes (assuming reddening)

        Arguments:
           counts -- photon counts (int or float or array or ndarray)
           E_BV -- E(B-V) attenuation factor (float)
        Output:
           mags -- absolute magnitudes (same type as input)
        """

        extinct = E_BV * self._red_per_ebv #magnitudes of extinction
        
        return -2.5*np.log10(counts / self._exposure) + self._zero_point - self._dmod + extinct

    def psf_convolve(self, image, multi_psf=True, convolve_func=None, **kwargs):
        """Convolve image with instrumental PSF
        
        Arguments:
           image -- counts, or flux, in each pixel of image (2D array of integers or floats)
        Keyword Arguments:
           multi_psf -- set to TRUE if 
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

        if self._psf.ndim != 4:
            multi_psf = False
        
        if convolve_func is None:
            N = image.shape[0]
            p = self._psf.shape[-1]
            if (p != self._psf.shape[-2]):
                message = 'each psf must be a square array'
                raise NotImplementedError(message)
            if multi_psf:
                assert(self._psf.ndim == 4)
                d_sub = self._psf.shape[0]
                assert(d_sub == self._psf.shape[1])
                #add border and subdivide
                sub_im_matrix = utils.subdivide_image(image, d_sub, w_border=p-1)
                convolved_matrix = np.array([[fftconvolve(sub_im_matrix[i,j], self._psf[i,j], mode='valid') for j in range(d_sub)] for i in range(d_sub)])
                im_convolved = np.concatenate(np.concatenate(convolved_matrix, axis=-2), axis=-1)
            else:
                #add border
                im_new = self._wrap_border(image, p-1)
                if self._psf.ndim == 2:
                    im_convolved = fftconvolve(im_new, self._psf, mode='valid')
                else:
                    im_convolved = fftconvolve(im_new, self._psf[0,0], mode='valid')
        else:
            im_convolved = convolve_func(image, self._psf, **kwargs)

        try:
            assert(im_convolved.shape == image.shape)
        except:
            print('Image shape has changed: ')
            print(image.shape)
            print('to ')
            print(im_convolved.shape)
            
        return im_convolved

