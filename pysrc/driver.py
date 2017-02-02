# driver.py
# Ben Cook (bcook@cfa.harvard.edu)

import numpy as np
from scipy.stats import poisson, norm
import instrument as ins
import isochrones as iso
import utils
import galaxy as gal
import warnings

#Try to import pycuda, initialize the GPU
GPU_AVAIL = False
try:
    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule
    import pycuda.curandom as curand
except ImportError as e:
    mess = e.__str__() #error message
    if 'No module named pycuda' in mess:
        warnings.warn('pycuda not installed.',ImportWarning)
    elif 'libcuda' in mess:
        warnings.warn('libcuda not found, likely because no GPU available.', RuntimeWarning)
    else:
        warnings.warn(mess, ImportWarning)
else:
    if (cuda.Device.count() > 0):
        GPU_AVAIL = True

if GPU_AVAIL:
    print('GPU acceleration enabled')
else:
    print('GPU acceleration not available')

class Driver:

    def __init__(self, iso_model, gpu=False):
        self.iso_model = iso_model
        self.filters = iso_model.filters
        self.n_filters = len(self.filters)

        if gpu:
            if GPU_AVAIL:
                self.gpu_on = True
            else:
                warnings.warn('GPU acceleration not available. Continuing without.', RuntimeWarning)
                self.gpu_on = False
        else:
            #No GPU acceleration
            self.gpu_on = False
        #No data has been initialized
        self._data_init = False
        
    def initialize_data(self, pcmd, bins, charlie=True):
        self.hess_bins = bins
        self.n_data = pcmd.shape[1]
        
        counts, hess, err = utils.bin_pcmd(pcmd, bins, charlie=charlie)
        self.counts_data = counts
        self.hess_data = hess
        self.err_data = err
        self._data_init = True

    def loglike(self, pcmd, use_gaussian=True, charlie=True):
        try:
            assert(self._data_init)
        except AssertionError:
            print('Cannot evaluate, as data has not been initialized (use driver.initialize_data)')
            return

        counts_model, hess_model, err_model = self._bin_pcmd(pcmd, self.hess_bins, charlie=charlie)
        n_model = pcmd.shape[1]

        if use_gaussian:
            #add error in quadrature
            combined_var = (self.err_data**2 + err_model**2) 
            hess_diff = (self.hess_data - hess_model)
            log_like = -np.sum(hess_diff**2 / combined_var)

        else:
            #Poisson Likelihood
            counts_model += 1e-3 #get NANs if model has zeros
            counts_model *= float(self.n_data) / n_model #normalize to same number of pixels as data
            log_like = np.sum(poisson.logpmf(self.counts_data, counts_model))

        return log_like

    def simulate(self, gal_model, im_scale, **kwargs):
        if self.gpu_on:
            return self._simulate_gpu(gal_model, im_scale, **kwargs)
        else:
            return self._simulate_cpu(gal_model, im_scale, **kwargs)

    def _simulate_cpu(self, gal_model, im_scale, psf='complex', fixed_seed=False, **kwargs):
        IMF, mags = self.iso_model.model_galaxy(gal_model)
        n_bins = len(IMF)
        psf_simple = (psf == 'simple')
        counts_per = np.array([f.mag_to_counts(m) for f,m in zip(self.filters, mags)])
        expected_num = IMF
        if fixed_seed:
            np.random.seed(0)
        realiz_num = np.random.poisson(lam=expected_num, size=(im_scale,im_scale, n_bins))
        
        raw_images = np.dot(realiz_num, counts_per.T).T
        raw_images += 1e-10
        convolved_images = np.array([f.psf_convolve(im, simple=psf_simple, **kwargs) for f,im in zip(self.filters,raw_images)])
        raw_mags = np.array([f.counts_to_mag(im, E_BV=gal_model.dust).flatten() for f,im in zip(self.filters, raw_images)])
        conv_mags = np.array([f.counts_to_mag(im, E_BV=gal_model.dust).flatten() for f,im in zip(self.filters, convolved_images)])
        return raw_mags, conv_mags, raw_images, convolved_images

    def _simulate_gpu(self, gal_model, im_scale, psf='complex', fixed_seed=False, **kwargs):
        assert(self.gpu_on)
        IMF, mags = self.iso_model.model_galaxy(gal_model)
        n_bins = len(IMF)
        psf_simple = (psf == 'simple')
        counts_per = np.array([f.mag_to_counts(m) for f,m in zip(self.filters, mags)])
        expected_num = IMF
        if fixed_seed:
            rng = curand.XORWOWRandomNumberGenerator(seed_getter=curand.seed_getter_uniform)
        else:
            rng = curand.XORWOWRandomNumberGenerator()
        
        raw_images = np.zeros((self.n_filters,im_scale*im_scale),dtype=float)
        
        for b in np.arange(n_bins):
            n_expected = expected_num[b]
            counts = counts_per[:,b]
            n_stars = rng.gen_poisson(im_scale*im_scale, np.uint32, n_expected).get()
            raw_images += np.array([c * n_stars for c in counts])
            
        raw_images += 1e-10
        raw_images = raw_images.reshape((self.n_filters, im_scale, im_scale))
        convolved_images = np.array([f.psf_convolve(im, simple=psf_simple, **kwargs) for f,im in zip(self.filters,raw_images)])
        raw_mags = np.array([f.counts_to_mag(im, E_BV=gal_model.dust).flatten() for f,im in zip(self.filters, raw_images)])
        conv_mags = np.array([f.counts_to_mag(im, E_BV=gal_model.dust).flatten() for f,im in zip(self.filters, convolved_images)])
        return raw_mags, conv_mags, raw_images, convolved_images
