# pcmd_fitter.py
# Ben Cook (bcook@cfa.harvard.edu)

import numpy as np
import instrument as ins
import isochrones as iso
import galaxy as gal
import driver
import fit_model
import utils
import gpu_utils
import pandas as pd
import os
import sys, getopt
import multiprocessing
import emcee
import importlib

if __name__ == "__main__":

    #Setup from an external file
    setup_mod = sys.argv[1].strip('.py')
    setup = importlib.import_module(setup_mod)

    N_scale = setup.N_scale
    N_walkers = setup.N_walkers
    N_burn = setup.N_burn
    N_sample = setup.N_sample
    N_threads = setup.N_threads
    pool = None
    gpu = setup.gpu
    force_gpu = setup.force_gpu

    fixed_seed = setup.fixed_seed
