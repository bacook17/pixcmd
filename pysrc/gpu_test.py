'''
Sample code showing running multiple processes including multiple GPUs
'''

import numpy as np
import multiprocessing, sys
import pycuda.driver as drv
import pycuda.curandom as curand
from pycuda.compiler import SourceModule

import os
import time

drv.init()

def set_gpu_device(q):
    '''
    This function makes pycuda use GPU number n in the system.
    '''
    n = q.get()
    n2 = multiprocessing.current_process()._identity[0]
    os.environ['CUDA_DEVICE'] = '%d'%n

    import pycuda.autoinit

    #Alternative version...also seems to work
    """
    drv.init()
    global device
    device = drv.Device(n)
    global context
    context = device.make_context()
    context.push()
    """
    
def func(x):
    n = multiprocessing.current_process()._identity[0]
    print('process %d, %d'%(n, x))
    time.sleep(0.1)

    #a = np.array([np.random.randn() * n])
    #a = a.astype(np.float32)
    #a_gpu = drv.mem_alloc(a.nbytes)
    #drv.memcpy_htod(a_gpu, a)#
    #
    #mod = SourceModule("""
    #    __global__ void doublify(float *a)
    #   {
    #       int idx = threadIdx.x + threadIdx.y*4;
    #       a[idx] *= 2;
    #   }""")
    #func = mod.get_function("doublify")
    #func(a_gpu, block=(1,1,1))
    #a_doubled = np.empty_like(a)
    #drv.memcpy_dtoh(a_doubled, a_gpu)
    #print a_doubled
    #print a
    gen = curand.XORWOWRandomNumberGenerator()
    r = gen.gen_uniform(10, np.float32)
    print(r)
    return(r)

if __name__ == '__main__':
    # If we have GPUs we use them, otherwise we use CPUs
    numprocesses = drv.Device.count() # number of GPUs present in the system
    gpunums = multiprocessing.Queue()
    for n in range(numprocesses):
        gpunums.put(n)
    #print('done making queue')
    pool = multiprocessing.Pool(processes=numprocesses, initializer=set_gpu_device, initargs=(gpunums,))
    #print('done making pool')
    # args consists of pairs (process_n, x) where process_n is the number of the
    # process, which is used to initialiase the correct GPU if one is present,
    # and x is the argument to the function (see notes in module docstring
    # about this)
    args = range(10)
    sys.stdout.flush()
    sys.stderr.flush()
    try:
        results = pool.map_async(func, args) # launches multiple processes
        results.wait(10)
    except KeyboardInterrupt:
        pass
    except:
        print('other error')
    print(results.get())
