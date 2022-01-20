#!/usr/bin/env python
from multiprocessing import Pool
import numpy as np
import pandas as pd
import os

def identifiable_example(x):
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
    vals[1] = 0.9 * (4 * x[0] - 2 + 4 * x[1] - 2)
    return vals

def wrap_function(ident):
    """
    Used to decide the function for use.
    Parameters:
    ===========
    ident: Boolean
    """
    if ident:
        return identifiable_example
    else:
        return non_identifiable_example

def call_functions(x, ident=True):
    """
    Parameters:
    ============
    x: np.ndarray, model inputs of the dimension (D, N) where D is the number of parameters and N is the sample size.
    ident: Boolean, use rosenbrock_function if True; non-identifiable_example if False.

    return:
    ========
    vals: model outputs
    """
    vals = wrap_function(ident)(x)
    vals = np.squeeze(vals)
    
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
    vals = call_functions(parameters_vals, ident=False)

    # write outputs to output.txt
    with open(output_file, 'w') as f:
        f.write('---- Obs 1----  \n')
        f.write(str(vals[0]) + '\n')
        
        f.write('---- Obs 2 ----  \n')
        f.write(str(vals[1]) + '\n') 

