# testing

import numpy as np
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import pycuda.curandom
import pycuda.gpuarray as gpuarray

code = """
   #include <curand_kernel.h>
   
   extern "C"
   {
   __global__ void my_random(curandState *global_state, const float *exp_nums, const float *fluxes, const int num_bins, const int num_bands, const int N, float *pixels) 
   {
      /* Initialize variables */
      int id_imx = blockIdx.x*blockDim.x + threadIdx.x;
      int id_imy = blockIdx.y*blockDim.y + threadIdx.y;
      int id_pix = id_imy + N*id_imx; 

      curandState local_state = global_state[threadIdx.x];
      float results[10] = {0};

      float flux, exp_num;
      int count;

      if (num_bands > 10) {
         /*Problem !!*/
      }


      if ((id_imx > N) || (id_imy > N)) {
      /* Out of image bound, don't do anything */
      }
      
      /* Draw stars from each mass bin, add up fluxes */
      for (int i = 0; i < num_bins, i++) {
         exp_num = exp_nums[i];
         count = curand_poisson(&local_state, exp_num);
         for (int f = 0; f < num_bands, f++) {
            flux = fluxes[i + (f*num_bins)];
            results[f] += count * flux;
         }
      }
      
      /* Save back state */
      global_state[threadIdx.x] = local_state;
      /* Save results into each pixel */
      for (int f = 0; f < num_bands, f++) {
          pixels[id_pix + (N*N)*f] = results[f];
      }
   }
   }
"""

code_simple = """
   #include <curand_kernel.h>

   extern "C"
   {
   __global__ void my_poisson(curandState *global_state, const float *exp_nums, const int N, const int num_lams, float *output) 
   {
      /* Initialize variables */
      int idx = blockIdx.x*blockDim.x + threadIdx.x;
      int idy = blockIdx.y*blockDim.y + threadIdx.y;

      int id_result = (num_lams * idy) + idx;

      curandState local_state = global_state[idx];

      if ((idx > N) || (idy > num_lams)) {
      /* Out of image bound, don't do anything */
      }

      output[id_result] = curand_poisson(&local_state, exp_nums[idy]);
      
      /* Save back state */
      global_state[threadIdx.x] = local_state;
   }

   }
"""

mod = SourceModule(code_simple, keep=False, no_extern_c=True)

my_poiss = mod.get_function('my_poisson')

def poiss_gpu(lam_arr, N):
    N_bins = len(lam_arr)
    generator = pycuda.curandom.XORWOWRandomNumberGenerator()
    result = np.zeros((N, N_bins), dtype=int)

    block_dim = (4,4,1)
    grid_dim = (N/4 + 1, N_bins/4 + 1)
    my_poiss(generator.state, cuda.In(lam_arr), N, len(lam_arr), cuda.Out(result), block=block_dim, grid=grid_dim)
    return result
