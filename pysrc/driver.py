# driver.py
# Ben Cook (bcook@cfa.harvard.edu)

import numpy as np
import instrument as ins
import isochrones as iso
import galaxy as gal
import warnings

class Driver:

    def __init__(self, iso_model, gpu=False):
        self.iso_model = iso_model
        self.filters = iso_model.filters
        self.n_filters = len(self.filters)

        if gpu:
            #Try to import pycuda, initialize the GPU
            try:
                import pycuda.driver as cuda
                import pycuda.autoinit
                from pycuda.compiler import SourceModule
                import pycuda.curandom as curand
            except ImportError as e:
                mess = e.__str__() #error message
                if 'No module named pycuda' in mess:
                    warnings.warn('pycuda not installed.',ImportWarning)
                    self.gpu_on = False
                    print('Continuing without GPU acceleration')
                elif 'libcuda' in mess:
                    warnings.warn('libcuda not found, likely because no GPU available.', RuntimeWarning)
                    self.gpu_on = False
                    print('Continuing without GPU acceleration')
                else:
                    raise ImportError(mess) 
            else:
                assert(cuda.Device.count() > 0)
                self.gpu_on = True
        else:
            #No GPU acceleration
            self.gpu_on = False

    def simulate(self, gal_model, im_scale):
        if self.gpu_on:
            return self._simulate_gpu(gal_model, im_scale)
        else:
            return self._simulate_cpu(gal_model, im_scale)

    def _simulate_cpu(self, gal_model, im_scale):
        IMF, mags = self.iso_model.model_galaxy(gal_model)
        n_bins = len(IMF)
        counts_per = np.array([f.mag_to_counts(m) for f,m in zip(self.filters, mags)])
        expected_num = IMF
        realiz_num = np.random.poisson(lam=expected_num, size=(im_scale,im_scale, n_bins))
        
        raw_images = np.dot(realiz_num, counts_per.T).T
        raw_images += 1e-10
        convolved_images = np.array([f.psf_convolve(im) for f,im in zip(self.filters,raw_images)])
        raw_mags = np.array([f.counts_to_mag(im).flatten() for f,im in zip(self.filters, raw_images)])
        conv_mags = np.array([f.counts_to_mag(im).flatten() for f,im in zip(self.filters, convolved_images)])
        return raw_mags, conv_mags, raw_images, convolved_images

    def _simulate_gpu(self, gal_model, im_scale):
        assert(self.gpu_on)
        IMF, mags = self.iso_model.model_galaxy(gal_model)
        n_bins = len(IMF)
        counts_per = np.array([f.mag_to_counts(m) for f,m in zip(self.filters, mags)])
        expected_num = IMF
        rng = curand.XORWOWRandomNumberGenerator()
        
        raw_images = np.zeros((self.n_filters,im_scale*im_scale),dtype=float)
        
        for b in np.arange(n_bins):
            n_expected = expected_num[b]
            counts = counts_per[:,b]
            n_stars = rng.gen_poisson(im_scale*im_scale, np.uint32, n_expected).get()
            raw_images += np.array([c * n_stars for c in counts])
            
        raw_images += 1e-10
        convolved_images = np.array([f.psf_convolve(im) for f,im in zip(self.filters,raw_images)])
        raw_mags = np.array([f.counts_to_mag(im).flatten() for f,im in zip(self.filters, raw_images)])
        conv_mags = np.array([f.counts_to_mag(im).flatten() for f,im in zip(self.filters, convolved_images)])
        return raw_mags, conv_mags, raw_images, convolved_images
