# testing

import numpy as np
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import pycuda.curandom
import pycuda.gpuarray as gpuarray

code_complex = """
   #include <curand_kernel.h>

   extern "C"
   {
   __global__ void poisson_sum(curandState *global_state, const float *exp_nums, const float *fluxes, const int num_bands, const int num_bins, const int N, float *pixels)
   {
      /* Initialize variables */
      int id_imx = blockIdx.x*blockDim.x + threadIdx.x;
      int id_imy = blockIdx.y*blockDim.y + threadIdx.y;
      int id_pix = (id_imx) + N*id_imy;
      int id_within_block = threadIdx.y + (blockDim.x * threadIdx.x);
      int block_id = blockIdx.x*gridDim.x + blockIdx.y;

      curandState local_state = global_state[id_within_block];
      float results[10] = {0.0};

      float flux;
      int count, skip;


      if ((id_imx < N) && (id_imy < N)) {
          /* Update local_state, to make sure values are very random */
          skip = 128 + block_id;
          skipahead(skip, &local_state);
          for (int i = 0; i < num_bins; i++){
             count = curand_poisson(&local_state, exp_nums[i]);
             for (int f = 0; f < num_bands; f++){
                flux = fluxes[i + (f*num_bins)];
                results[f] += count * flux;
             }
          }
          /* Save results for each band */
          for (int f = 0; f < num_bands; f++){
             pixels[id_pix + (N*N)*f] = results[f];
          }
      }

      /* Save back state */
      global_state[id_within_block] = local_state;
   }

   }
   }
"""

mod_2 = SourceModule(code_complex, keep=False, no_extern_c=True)

poiss_sum = mod_2.get_function('poisson_sum')

def poiss_sum_gpu(expected_nums, fluxes, N_scale):
    N_bins = len(expected_nums)
    N_bands = fluxes.shape[0]
    assert(N_bins == fluxes.shape[1])
    generator = pycuda.curandom.XORWOWRandomNumberGenerator(seed_getter=pycuda.curandom.seed_getter_unique)
    result = np.zeros((N_bands, N_scale, N_scale), dtype=np.float32)
    
    d_block = 32
    
    block_dim = (d_block, d_block,1)
    grid_dim = (N_scale/d_block + 1, N_scale/d_block + 1)
    poiss_sum(generator.state, cuda.In(expected_nums), cuda.In(fluxes), np.int32(N_bands), np.int32(N_bins), np.int32(N_scale), cuda.Out(result), block=block_dim, grid=grid_dim)
    return result

code_simple = """
   #include <curand_kernel.h>

   extern "C"
   {
   __global__ void my_poisson(curandState *global_state, const float *exp_nums, const int N, const int num_lams, int *output)
   {
      /* Initialize variables */
      int idx = blockIdx.x*blockDim.x + threadIdx.x;
      int idy = blockIdx.y*blockDim.y + threadIdx.y;

      int id_result = (num_lams * idx) + idy;

      curandState local_state = global_state[idx];

      if ((idx < N) && (idy < num_lams)) {
          output[id_result] = curand_poisson(&local_state, exp_nums[idy]);
      }

      /* Save back state */
      global_state[idx] = local_state;
   }

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
