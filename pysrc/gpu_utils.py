import numpy as np
import warnings
import multiprocessing
_GPU_AVAIL = True
try:
    import pycuda
    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule
    import pycuda.curandom
except ImportError as e:
    mess = e.__str__() #error message
    if 'No module named pycuda' in mess:
        warnings.warn('pycuda not installed.',ImportWarning)
        print('pycuda not installed.')
    elif 'libcuda' in mess:
        warnings.warn('libcuda not found, likely because no GPU available.', RuntimeWarning)
        print('libcuda not found, likely because no GPU available.')
    else:
        warnings.warn(mess, ImportWarning)
        print(mess)
    _GPU_AVAIL = False
else:
    if (cuda.Device.count() > 0):
        _GPU_AVAIL = True

if _GPU_AVAIL:
    print('GPU acceleration enabled')
else:
    print('GPU acceleration not available, sorry')

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
      int id_within_block = threadIdx.y + (blockDim.x * threadIdx.x);
      int block_id = blockIdx.x*gridDim.x + blockIdx.y;

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
      global_state[id_within_block] = local_state;
   }
   }
"""

_CUDAC_AVAIL = False
if _GPU_AVAIL:
    try:
        _mod = SourceModule(_code, keep=False, no_extern_c=True)
        _func = _mod.get_function('poisson_sum')
    except:
        _CUDAC_AVAIL = False
    else:
        _CUDAC_AVAIL = True

def draw_image(expected_nums, fluxes, N_scale, gpu=_GPU_AVAIL, cudac=_CUDAC_AVAIL, **kwargs):
    if gpu:
        if cudac:
            func = _draw_image_cudac
        else:
            func = _draw_image_pycuda
    else:
        func = _draw_image_numpy
    return func(expected_nums, fluxes, N_scale, **kwargs)
        
def _draw_image_cudac(expected_nums, fluxes, N_scale, tolerance=1e-5, d_block=32):
    assert(_CUDAC_AVAIL)
    assert(_GPU_AVAIL)

    assert(len(expected_nums) == fluxes.shape[1])
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
    
    expected_nums = expected_nums.astype(np.float32)
    fluxes = fluxes.astype(np.float32)

    N_scale = np.int32(N_scale)

    N_bins = np.int32(len(expected_nums))
    N_bands = np.int32(fluxes.shape[0])
    
    generator = pycuda.curandom.XORWOWRandomNumberGenerator()
    result = np.zeros((N_bands, N_scale, N_scale), dtype=np.float32)
    
    block_dim = (d_block, d_block,1)
    grid_dim = (N_scale/d_block + 1, N_scale/d_block + 1)
    _func(generator.state, cuda.In(expected_nums), cuda.In(fluxes), N_bands, N_bins, N_scale,
              cuda.Out(result), block=block_dim, grid=grid_dim)

    #Add on flux from fully-populated bins
    #result = np.array([result[i] + fixed_fluxes[i] for i in range(N_bands)]).astype(float)
    return result

def _draw_image_pycuda(expected_nums, fluxes, N_scale, tolerance=-1., **kwargs):
    assert(_GPU_AVAIL)

    N_bins = len(expected_nums)
    N_bands = fluxes.shape[0]
    assert(N_bins == fluxes.shape[1])
    if (tolerance < 0.):
        upper_lim = np.inf
    else:
        upper_lim = tolerance**-2.
    
    generator = pycuda.curandom.XORWOWRandomNumberGenerator()
    result = np.zeros((N_bands, N_scale*N_scale), dtype=float)

    #Draw stars and cumulate flux for each mass bin
    for b in np.arange(N_bins):
        n_expected = expected_nums[b]
        counts = fluxes[:,b]
        if (n_expected <= upper_lim):
            n_stars = generator.gen_poisson(N_scale*N_scale, np.uint32, n_expected).get()
        #assume no poisson variance
        else:
            n_stars = n_expected
        result += np.array([c * n_stars for c in counts])

    return result.reshape([N_bands, N_scale, N_scale])

def _draw_image_numpy(expected_nums, fluxes, N_scale, fixed_seed=False, tolerance=-1., **kwargs):
    N_bins = len(expected_nums)
    assert(N_bins == fluxes.shape[1])
    if (tolerance < 0.):
        upper_lim = np.inf
    else:
        upper_lim = tolerance**-2.
    if fixed_seed:
        np.random.seed(0)

    realiz_num = np.zeros((N_scale, N_scale, N_bins))
    if not np.isinf(upper_lim):
        realiz_num = np.random.poisson(lam=expected_nums, size=(N_scale, N_scale, N_bins))
    else:
        use_poisson = (expected_nums <= upper_lim)
        use_fixed = ~use_poisson #Assume no poisson variance
        num_poisson = np.sum(use_poisson)
    
        realiz_num[:,:,use_fixed] = expected_nums[use_fixed]
        realiz_num[:,:,use_poisson] = np.random.poisson(lam=expected_nums[use_poisson], size=(N_scale, N_scale, num_poisson))
        
    return np.dot(realiz_num, fluxes.T).T
    
def set_gpu_device(n):
    """
    This function makes pycuda use GPU number n in the system.
    """
    cuda.init()
    assert(n < cuda.Device.count())
    try:
        pycuda.context.detach()
    except:
        pass
    pycuda.context = cuda.Device(n).make_context()

def initialize_process():
    cpu_id = multiprocessing.current_process()
    set_gpu_device(cpu_id)
