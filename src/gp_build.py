#!/usr/bin/env python
from multiprocessing import Pool
import numpy as np
import os
import matplotlib.pyplot as plt
from functools import partial
import copy
import pandas as pd
import pickle

from scipy import stats
# from scipy.optimize import root
from scipy.optimize import bisect

from sklearn.gaussian_process.kernels import RBF, \
    Matern

from pyapprox.density import tensor_product_pdf
from pyapprox.gaussian_process import CholeskySampler, AdaptiveGaussianProcess, generate_gp_candidate_samples
from pyapprox.low_discrepancy_sequences import transformed_halton_sequence
from pyapprox.utilities import compute_f_divergence, \
    get_tensor_product_quadrature_rule
from pyapprox.probability_measure_sampling import rejection_sampling
from pyapprox.visualization import get_meshgrid_function_data
from pyapprox import generate_independent_random_samples
from pyapprox.variables import IndependentMultivariateRandomVariable
from pyapprox.variable_transformations import AffineRandomVariableTransformation

import matplotlib as mpl
from matplotlib import rc
import spotpy as sp

from funcs.read_data import variables_prep, file_settings


mpl.rcParams['font.size'] = 16
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['text.usetex'] = False  # use latex for all text handling
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

def read_candidate_set(outpath, num_files):
    def rescale_parameters(outpath, pars):
        """
        Transform parameters into the original ranges
        """
        # read files containing parameter ranges
        par_range = pd.read_csv(f'{outpath}parameters.tpl', skiprows=1, index_col = 'parameter')
        for col in list(par_range.index):
            val_low, val_up = par_range.loc[col, 'lower'], par_range.loc[col, 'upper']
            pars[col.lower()] = pars[col.lower()]*(val_up - val_low) / 100 + val_low

        return pars

    for n in range(num_files[0], num_files[1]):
        par_temp = pd.read_csv(f'{outpath}126001A.{n}.par.csv', index_col = 'real_name', skiprows=[100])
        obs_temp = pd.read_csv(f'{outpath}126001A.{n}.obs.csv', index_col = 'real_name', skiprows=[100])
        try:
            pars = pd.concat([pars, par_temp])
            obs = pd.concat([obs, obs_temp])
        except NameError:
            pars = pd.read_csv(f'{outpath}126001A.{n}.par.csv', index_col = 'real_name', skiprows=[100])
            obs = pd.read_csv(f'{outpath}126001A.{n}.obs.csv', index_col = 'real_name', skiprows=[100])

    par_samples = rescale_parameters(outpath, pars)
    par_samples.drop(columns=['lcf'], inplace=True)
    true_meas = [52093.389, 99477.940, 44063.700, 57936.470, 
        53449.050, 21858.007, 38560.992, 51843.258, 14176.304]
    obs.reset_index(inplace=True)
    for ii in range(obs.shape[0]):
        obs.loc[ii, 'din_mse'] = (4.36e-05**2)*sp.objectivefunctions.mse(true_meas, obs.loc[:, 'din_2009':'din_2017'].values[ii])
    obs['din_pbias'] = 0.01*obs['din_pbias'] **2
    return pars, obs
    # END read_candidate_set()

outpath = '../output/work_run_0520/'
global year
global num_files
year = 'mse'
num_files = [0, 9]
par_samples, obs_samples = read_candidate_set(outpath, num_files=num_files)
par_vali, obs_vali = read_candidate_set(outpath, num_files=[9, 12])
values_vali = obs_vali[f'din_{year}'].values

if not os.path.exists(f'year_{year}/'):
    os.mkdir(f'year_{year}/')
par_samples[f'din_{year}'] = obs_samples[f'din_{year}'].values

def calculate_objectives(vars, pars_obs = par_samples):
    """
    Script used to run_source and return the output file.
    The function is called by AdaptiveLejaPCE.
    """

    pars_set = pars_obs.values[:, :-1].T
    assert pars_set.shape[0] == vars.shape[0], \
        "The number of parameters is the same with the first dimension of vars."

    idx = np.zeros(vars.shape[1])
    for ii in range(vars.shape[1]):
        idx[ii] = np.argwhere((pars_set.T == vars.T[ii]).all(axis=1)).ravel()

    y = pars_obs.values[:, -1:][idx.astype(int)]
    print(f'--------------Output values: {y}-------------')
    return y

def convergence_study(kernel, function, sampler,
                      num_vars, num_new_samples,
                      update_kernel_scale_num_samples,
                      noise_level=0, return_samples=False,
                      norm=np.linalg.norm, callback=None, gp_kernel=None, candidate_set=None):

    # dirty hack to include two GP kernel types (for IVAR)
    if hasattr(kernel, "__len__"):
        # in this case, kernel is an array and we assume to have received
        # two kernels
        sampler_kernel = kernel[1]
        kernel = kernel[0]
    else:
        sampler_kernel = kernel

    # Instantiate a Gaussian Process model
    if gp_kernel is None:
        gp_kernel = kernel
    

    gp = AdaptiveGaussianProcess(
        gp_kernel, n_restarts_optimizer=2, alpha=1e-12)
    gp.setup(function, sampler)
    if hasattr(sampler, "set_gaussian_process"):
        sampler.set_gaussian_process(gp)
    
    print('sampler kernel', kernel, 'gp kernel', gp_kernel)

    # Mesh the input space for evaluations of the real function,
    # the prediction and its MSE

    if candidate_set is not None:
        sampler.candidate_samples = candidate_set

    num_samples = np.cumsum(num_new_samples)
    num_steps = num_new_samples.shape[0]
    errors = np.empty(num_steps, dtype=float)
    nsamples = np.empty(num_steps, dtype=int)
    sample_step = 0
    optimizer_step = 0
    # breakpoint()
    while sample_step < num_steps:
        if hasattr(gp, 'kernel_'):
            # if using const * rbf + noise kernel
            # kernel.theta = gp.kernel_.k1.k2.theta
            # if using const * rbf
            # kernel.theta = gp.kernel_.k2.theta
            # if using rbf
            kernel.theta = gp.kernel_.theta

        # Fit to data using Maximum Likelihood Estimation of the parameters
        # if True:
        if ((optimizer_step >= update_kernel_scale_num_samples.shape[0]) or
            (sampler.ntraining_samples <
             update_kernel_scale_num_samples[optimizer_step])):
            gp.optimizer = None
        else:
            gp.optimizer = "fmin_l_bfgs_b"
            optimizer_step += 1

        flag = gp.refine(np.sum(num_new_samples[:sample_step+1]))

        # allow points to be added to gp more often than gp is evaluated for
        # validation
        if sampler.ntraining_samples >= num_samples[sample_step]:
            if callback is not None:
                callback(gp)

            print(gp.kernel_)
             
            if sample_step >=1:
                # Compute error
                gp_load = pickle.load(open(f'year_{year}/gp_{np.mod(sample_step - 1, 2)}_{year}.pkl', "rb"))
                validation_sub = par_vali.values.T
                pred_values = gp_load(validation_sub, return_cov=False).squeeze()
                values_sub = values_vali
                error_gp_comp = norm(pred_values-values_sub)/norm(values_sub)
                print('-----------error_gp_comp---------', error_gp_comp)

                print('N', sampler.ntraining_samples, 'Error', error_gp_comp)

            if sample_step >= 1:
                errors[sample_step -1] = error_gp_comp
                nsamples[sample_step - 1] = num_samples[sample_step -1]

            pickle.dump(gp, open(f'year_{year}/gp_{np.mod(sample_step, 2)}_{year}.pkl', "wb"))
            sample_step += 1            

        if flag > 0:
            errors, nsamples = errors[:sample_step], nsamples[:sample_step]
            print('Terminating study. Points are becoming ill conditioned')
            break

    if return_samples:
        return errors, nsamples, sampler.training_samples[:, 0:num_samples[sample_step - 1]]

    return errors, nsamples


def unnormalized_posterior(gp, prior_pdf, samples, temper_param=1):
    prior_vals = prior_pdf(samples).squeeze()
    gp_vals = gp.predict(samples.T).squeeze()
    # calculate the difference between gp values and the true model outputs
    true_meas = [52093.389, 99477.940, 44063.700, 57936.470, 
        53449.050, 21858.007, 38560.992, 51843.258, 14176.304]
    if isinstance(year, int):
        error_diff = np.square(gp_vals - true_meas[int(year - 2009)])
    else:

        error_diff = np.abs(gp_vals)
    unnormalized_posterior_vals = prior_vals*np.exp(-error_diff)**temper_param
    return unnormalized_posterior_vals


class BayesianInferenceCholeskySampler(CholeskySampler):
    def __init__(self, prior_pdf, num_vars,
                 num_candidate_samples, variables,
                 max_num_samples=None, generate_random_samples=None,
                 temper=True, true_nll=None):
        self.prior_pdf = prior_pdf
        if not temper:
            self.temper_param = 1
        else:
            self.temper_param = 0
        self.true_nll = true_nll
        self.gp = None

        super().__init__(num_vars, num_candidate_samples, variables,
                         None, generate_random_samples)

    def set_gaussian_process(self, gp):
        self.gp = gp

    # Qian: understand the purpose of function increment_temper_param()
    def increment_temper_param(self, num_training_samples):

        # samples = np.random.uniform(0, 1, (self.nvars, 1000))
        samples = generate_independent_random_samples(self.variables, 1000)
        density_vals_prev = self.weight_function(samples)

        def objective(beta):
            new_weight_function = partial(
                unnormalized_posterior, self.gp, self.prior_pdf,
                temper_param=beta)
            density_vals = new_weight_function(samples)

            # breakpoint()
            II = np.where(density_vals_prev > 1e-15)[0]
            JJ = np.where(density_vals_prev < 1e-15)[0]
            assert len(np.where(density_vals[JJ] > 1e-15)[0]) == 0
            ratio = np.zeros(samples.shape[1])
            ratio[II] = density_vals[II]/density_vals_prev[II]
            obj = ratio.std()/ratio.mean()
            return obj
        print('temper parameter', self.temper_param)
        x0 = self.temper_param+1e-4
        # result = root(lambda b: objective(b)-1, x0)
        # x_opt = result.x
        # breakpoint()
        
        x_opt = bisect(lambda b: objective(b)-1, x0, 1)
        # if not optimize temper_param
        # x_opt = self.temper_param + 1e-4
        self.temper_param = x_opt

    def __call__(self, num_samples):
        if self.gp is None:
            raise ValueError("must call self.set_gaussian_process()")
        
        if self.ntraining_samples > 0 and self.temper_param < 1:
            self.increment_temper_param(self.training_samples)
        assert self.temper_param <= 1
        if self.ntraining_samples == 0:
            weight_function = self.prior_pdf
        else:
            if self.true_nll is not None:
                def weight_function(x): return self.prior_pdf(x)*np.exp(
                    -self.true_nll(x)[:, 0])**self.temper_param
            else:
                weight_function = partial(
                    unnormalized_posterior, self.gp, self.prior_pdf,
                    temper_param=self.temper_param)

        self.set_weight_function(weight_function)

        samples, flag = super().__call__(num_samples)
        return samples, flag


def get_prior_samples(num_vars, variables, nsamples):
    rosenbrock_samples = generate_independent_random_samples(variables, nsamples)

    return rosenbrock_samples

def bayesian_inference_example():
    # read parameter distributions
    datapath = file_settings()[1]
    para_info = pd.read_csv(datapath + 'Parameters-PCE.csv')
    # define the variables for PCE
    param_file = file_settings()[-1]

    outpath = '../output/work_run_0520/'
    par_samples, obs_samples = read_candidate_set(outpath, num_files)
    par_samples[f'din_{year}'] = obs_samples[f'din_{year}'].values
    
    # Must set variables if not using uniform prior on [0,1]^D
    # variables = None
    ind_vars, variables = variables_prep(param_file, product_uniform='uniform', dummy=False)
    var_trans = AffineRandomVariableTransformation(variables, enforce_bounds=True)
    init_scale = 50# used to define length_scale for the kernel
    num_vars = variables.nvars
    num_candidate_samples = 20000
    num_new_samples = np.asarray([20]+[10]*10+[20]*9+[50]*8)


    prior_pdf = partial(tensor_product_pdf, 
        univariate_pdfs=[partial(stats.beta.pdf, a=1, b=1, scale=ind_vars[ii].args[1]) for ii in range(num_vars)])

# END calculate_objectives()

    def get_filename(method, fixed_scale):
        filename = 'bayes-example-%s-d-%d-n-%d.npz' % (
            method, num_vars, num_candidate_samples)
        if not fixed_scale:
            filename = filename[:-4]+'-opt.npz'
        return filename

    # defining kernel
    length_scale = [init_scale, init_scale, *(3*np.ones(num_vars -2, dtype=float))]
    kernel = RBF(length_scale, [(5e-2, 200), (5e-2, 100), (5e-2, 20), (5e-2, 15),
        (5e-2, 20), (5e-2, 15), (5e-2, 20), (5e-2, 15), (5e-2, 20), 
        (5e-2, 15), (5e-2, 20), (5e-2, 15), (5e-2, 20)])

    # this is the one Qian should use. The others are for comparision only
    adaptive_cholesky_sampler = BayesianInferenceCholeskySampler(
        prior_pdf, num_vars, num_candidate_samples, variables,
        max_num_samples=num_new_samples.sum(),
        generate_random_samples=None)
    adaptive_cholesky_sampler.set_kernel(copy.deepcopy(kernel))

    samplers = [adaptive_cholesky_sampler]
    methods = ['Learning-Weighted-Cholesky-b']
    labels = [r'$\mathrm{Adapted\;Weighted\;Cholesky}$']
    fixed_scales = [False]

    for sampler, method, fixed_scale in zip(samplers, methods, fixed_scales):
        filename = get_filename(method, fixed_scale)
        print(filename)
        if os.path.exists(filename):
            continue

        if fixed_scale:
            update_kernel_scale_num_samples = np.empty(0)
        else:
            update_kernel_scale_num_samples = np.cumsum(num_new_samples)

        cond_nums = []
        temper_params = []

        def callback(gp):
            cond_nums.append(np.linalg.cond(gp.L_.dot(gp.L_.T)))
            if hasattr(sampler, 'temper_param'):
                temper_params.append(sampler.temper_param)
                print(temper_params)

        errors, nsamples, samples = convergence_study(
            kernel, calculate_objectives, sampler, num_vars,
            num_new_samples, update_kernel_scale_num_samples, callback=callback,
            return_samples=True, candidate_set=par_samples.values[:, :-1].T)

        np.savez(f'year_{year}/{filename}', nsamples=nsamples, errors=errors,
                 cond_nums=np.asarray(cond_nums), samples=samples,
                 temper_params=np.asarray(temper_params))

    fig, axs = plt.subplots(1, 3, figsize=(3*8, 6), sharey=False)
    styles = ['-']
    # styles = ['k-','r-.','b--','g:']
    for method, label, ls, fixed_scale in zip(
            methods, labels, styles, fixed_scales):
        filename = get_filename(method, fixed_scale)
        data = np.load(f'year_{year}/{filename}')
        nsamples, errors = data['nsamples'][:-1], data['errors'][:-1]
        temper_params, cond_nums = data['temper_params'][1:-1], data['cond_nums'][:-1]
        axs[0].loglog(nsamples, errors, ls=ls, label=label)
        axs[1].loglog(nsamples, cond_nums, ls=ls, label=label)
        axs[2].semilogy(np.arange(1, nsamples.shape[0]),
                    temper_params, 'k-o')
        axs[2].set_xlabel(r'$\mathrm{Iteration}$ $j$')
        axs[2].set_ylabel(r'$\beta_j$')

    for ii in range(2):
        axs[ii].set_xlabel(r'$m$')
        axs[ii].set_xlim(10, 1000)
    axs[0].set_ylabel(r'$\tilde{\epsilon}_{\omega,2}$', rotation=90)
    ylim0 = axs[0].get_ylim()
    ylim1 = axs[1].get_ylim()
    ylim = [min(ylim0[0], ylim1[0]), max(ylim0[1], ylim1[1])]
    axs[0].set_ylim(ylim)
    axs[1].set_ylim(ylim)
    axs[1].set_ylabel(r'$\kappa$', rotation=90)

    figname = 'bayes_example_comparison_%d.pdf' % num_vars
    axs[0].legend()
    plt.savefig(f'year_{year}/{figname}') 

if __name__ == '__main__':
    try:
        import sklearn
    except:
        msg = 'Install sklearn using pip install sklearn'
        raise Exception(msg)

    bayesian_inference_example()
