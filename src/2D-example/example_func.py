#!/usr/bin/env python
from multiprocessing import Pool
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from functools import partial
import time
import copy

from scipy.stats import multivariate_normal
from scipy import stats
# from scipy.optimize import root
from scipy.optimize import bisect

from sklearn.gaussian_process.kernels import RBF, \
    Matern

import matplotlib as mpl
mpl.rcParams['font.size'] = 16
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['text.usetex'] = True  # use latex for all text handling
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['savefig.format'] = 'pdf'  # gives best resolution plots
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['legend.fontsize'] = 16
# print mpl.rcParams.keys()
mpl.rcParams['text.latex.preamble'] = \
    r'\usepackage{siunitx}\usepackage{amsmath}\usepackage{amssymb}'

def rosenbrock_function(x):
    assert x.shape[0] == 2, "The first dimension of x should be 2 being equale to the number of parameters."
    # x = 4*x-2
    # vals = ((1.-x[0, :])**2+100*(x[1, :]-x[0, :]**2)**2)[:, np.newaxis]
    # vals = ((1.-x[0,:])**2+1*(x[1,:]-x[0,:]**2)**2)[:,np.newaxis]
    vals = np.zeros(shape=(2, x.shape[1]))
    vals[0] = 4 * x[0] - 2
    vals[1] = 4 * x[1] - 2 - (4 * x[0] - 2)**2
    return vals


def non_identifiable_example(x):
    assert x.shape[0], "The first dimension of x should be 2 being equale to the number of parameters."
    # x = 4 * x - 2
    # vals = ((1 - (x[0, :] + x[1, :])) ** 2 + 25 * (x[0, :] + x[1, :]) ** 2)[:, np.newaxis]
    vals = np.zeros(shape=(2, x.shape[1]))
    vals[0] = 4 * x[0] - 2 + 4 * x[1] - 2
    vals[1] = 1 /2 * (4 * x[0] - 2 + 4 * x[1] - 2)
    return vals

def wrap_function(ident):
    """
    Used to decide the function for use.
    Parameters:
    ===========
    ident: Boolean
    """
    if ident:
        return rosenbrock_function
    else:
        return non_identifiable_example

def call_functions(x, ident=True):
    """
    Parameters:
    ============
    x: np.ndarray, model inputs of the dimension (D, N) where D is the number of parameters and N is the sample size.
    ident: Boolean, use rosenbrock_function if True; non-identifiable_example if False.

    return:
    =========
    vals: model outputs
    """
    vals = wrap_function(ident)(x)
    
    return vals


if __name__ == '__main__':
    output_file = 'example_output.txt'
    print('Read Parameters')
    parameters = pd.read_csv('parameters.csv')
    parameters_vals = np.zeros(shape=(parameters.shape[0], 1))
    for i,j in parameters.iterrows():
        scaled_value = (j.upper - j.lower) * j.value/100 + j.lower 
        parameters_vals[i] = scaled_value

    # call the function
    vals = call_functions(parameters_vals)

    # write outputs to output.txt
    with open(output_file, 'w') as f:
        f.write('---- Obs 1----  \n')
        f.write(str(vals[0]) + '\n')
        
        f.write('---- Obs 2 ----  \n')
        f.write(str(vals[1]) + '\n') 