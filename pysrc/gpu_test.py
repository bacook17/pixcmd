'''
Sample code showing running multiple processes including multiple GPUs
'''

import numpy as np
import multiprocessing, sys
import pycuda.driver as drv
import pycuda.curandom as curand
import pycuda.gpuarray as gpuarray
from pycuda.compiler import SourceModule
import os
import time

drv.init()

_code = """
   #include <curand_kernel.h>

   extern "C"
   {
   __global__ void poisson_sum(curandState *global_state, const float *exp_nums, const float *fluxes, const int num_bands, const int num_bins, const int N, float *pixels)
   {
      /* Initialize variables */
      int id_imx = blockIdx.x*blockDim.x + threadIdx.x;
      int id_imy = blockIdx.y*blockDim.y + threadIdx.y;
      int id_pix = (id_imx) + N*id_imy;
      int id_within_block = threadIdx.x + (blockDim.x * threadIdx.y);
      int block_id = blockIdx.y*gridDim.x + blockIdx.x;

      curandState local_state = global_state[id_within_block];
      float results[10] = {0.0};

      float flux;
      int count, skip;

      if ((id_imx < N) && (id_imy < N)) {
          /* Update local_state, to make sure values are very random */
          skip = 20 * block_id;
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
      /*global_state[id_within_block] = local_state;*/
   }
   }
"""

def _draw_image_cudac(N_bins, N_scale=64):

    """
    upper_lim = tolerance**-2.
    use_poisson = (expected_nums <= upper_lim)
    use_fixed = ~use_poisson

    #total flux from all "fully populated" bins above upper_lim
    fixed_fluxes = np.sum(expected_nums[use_fixed]*fluxes[:,use_fixed], axis=1)

    #remove fixed bins, and set proper byte size for cuda
    expected_nums = expected_nums[use_poisson].astype(np.float32)
    fluxes = fluxes[:,use_poisson].astype(np.int32)
    """

    expected_nums = np.random.uniform(0., 100., size=N_bins).astype(np.float32)

    fluxes = np.random.uniform(0., 100., size=N_bins).astype(np.float32)

    N_scale = np.int32(N_scale)

    seed_getter = curand.seed_getter_uniform
    generator = pycuda.curandom.XORWOWRandomNumberGenerator(seed_getter=seed_getter)
    result = np.zeros((N_scale, N_scale), dtype=np.float32)

    block_dim = (d_block, d_block,1)
    grid_dim = (N_scale/d_block + 1, N_scale/d_block + 1)
    _func(generator.state, cuda.In(expected_nums), cuda.In(fluxes), N_bands, N_bins, N_scale,
              cuda.Out(result), block=block_dim, grid=grid_dim)

    #Add on flux from fully-populated bins
    #result = np.array([result[i] + fixed_fluxes[i] for i in range(N_bands)]).astype(float)
    print(result)
    return result

def set_gpu_device():
    '''
    This function makes pycuda use GPU number n in the system.
    '''
    n = multiprocessing.current_process()._identity[0] - 1
    os.environ['CUDA_DEVICE'] = '%d'%n

    import pycuda.autoinit

    #drv.init()
    #global device
    #device = drv.Device(n)
    #global context
    #context = device.make_context()
    #context.push()
    #print(device)
    bid_current = pycuda.autoinit.device.pci_bus_id()
    bid_n = drv.Device(n).pci_bus_id()
    print('active PCI_BUS_ID matches nth: ' + str(str(bid_current) == str(bid_n)))
    global _mod
    _mod = SourceModule(_code, keep=False, no_extern_c=True)
    print('through here')
    global _func
    _func = _mod.get_function('possion_sum')

    #print('Setting context: %d'%n)
    #print(context)
    
def func(x):
    n = multiprocessing.current_process()._identity[0]-1
    #print('process %d, %d, %s'%(x, n, os.environ.get('CUDA_DEVICE')))

    gen = curand.XORWOWRandomNumberGenerator()
    a = np.empty(2, dtype=np.uint32)
    a_gpu = gpuarray.to_gpu(a.astype(np.uint32))
    gen.fill_poisson(a_gpu, x+1)
    a = a_gpu.get()

    return(a)

if __name__ == '__main__':
    # If we have GPUs we use them, otherwise we use CPUs
    numprocesses = drv.Device.count() # number of GPUs present in the system
    pool = multiprocessing.Pool(processes=numprocesses, initializer=set_gpu_device)
    args = range(50,100)
    sys.stdout.flush()
    sys.stderr.flush()
    try:
        #results = pool.map(_draw_image_cudac, args) # launches multiple processes
        results = pool.map(func, args)
    except drv.LogicError:
        print('Logic Error')
        raise
    print(results)
