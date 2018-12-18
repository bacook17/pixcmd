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

#import pycuda.autoinit
#drv.init()
#device = drv.Device(0)
#context = device.make_context()
#context.push()

_code = """
   #include <curand_kernel.h>

   extern "C"
   {
   __global__ void poisson_sum(curandState *global_state, const float *exp_nums, const int N, float *results)
   {
      /* Draw poisson random values and sum */

      /* Initialize variables */
      int idx = blockIdx.x*blockDim.x + threadIdx.x;
      int id_within_block = threadIdx.x;

      curandState local_state = global_state[id_within_block];
      int count;

      if (idx < N) {
          for (int i = 0; i < N; i++){
             count = curand_poisson(&local_state, exp_nums[i]);
             results[idx] += count * 1.0;
          }
      }

      /* Save back state */
      /*global_state[id_within_block] = local_state;*/
   }
   }
"""

#print('Unloading pycuda contexts')
#context.pop()
#context = None
#from pycuda.tools import clear_context_caches
#clear_context_caches()

print('Done')

def func_complex(x):

    N = np.int32(128)
    expected_nums = x * np.random.random(size=N).astype(np.float32)

    generator = curand.XORWOWRandomNumberGenerator()
    result = np.zeros(N, dtype=np.float32)

    d_block = 16
    block_dim = (d_block, 1, 1)
    grid_dim = (N / d_block + 1, 1)
    
    _poisson_sum(generator.state, drv.In(expected_nums), N,
          drv.Out(result), block=block_dim, grid=grid_dim)

    return result

def set_gpu_device():
    '''
    This function makes pycuda use GPU number n in the system.
    '''
    n = multiprocessing.current_process()._identity[0] - 1
    print('n: %d'%n)
    os.environ['CUDA_DEVICE'] = '%d'%n

    import pycuda.autoinit

    global _mod
    print('Starting SourceModule Code')
    _mod = SourceModule(_code, keep=False, no_extern_c=True)
    print('Getting function')
    global _poisson_sum
    _poisson_sum = _mod.get_function('poisson_sum')
    print('Past the SourceModule code')

    #if (n > 0):
        #reload(pycuda.autoinit)
        #Alternative version...also seems to work
        #drv.init()
        #global device
        #device = drv.Device(n)
        #global context
        #context = device.make_context()
        #context.push()
    print('n done: %d'%n)
    
def func(x):
    n = multiprocessing.current_process()._identity[0]-1
    print('process %d, arg %d, CUDA_DEVICE %s'%(n, x, os.environ.get('CUDA_DEVICE')))

    gen = curand.XORWOWRandomNumberGenerator()
    a = np.empty(2, dtype=np.uint32)
    a_gpu = gpuarray.to_gpu(a.astype(np.uint32))
    gen.fill_poisson(a_gpu, x+1)
    a = a_gpu.get()

    return(a)

print('Past defining functions')

if __name__ == '__main__':
    # If we have GPUs we use them, otherwise we use CPUs
    numprocesses = 2
    pool = multiprocessing.Pool(processes=numprocesses, initializer=set_gpu_device)
    print('Past creating pool and initalizing')
    args = range(32)
    sys.stdout.flush()
    sys.stderr.flush()
    try:
        results = pool.map(func_complex, args) # launches multiple processes
        #results = pool.map(func, args)
    except drv.LogicError:
        print('Logic Error')
        raise
    print(results)
