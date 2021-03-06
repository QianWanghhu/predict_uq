# import packages
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pyapprox as pya
import json

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


# read files
f_names = [f'126001A.{i}.par.csv' for i in range(12)]
outpath = '../output/work_run_0520/'
fig_path = '../output/figs/'

def pair_plot():
    for fn in f_names:
        par_temp = pd.read_csv(outpath + fn, index_col = 'real_name')
        try:
            pars = pd.concat([pars, par_temp])
        except NameError:
            pars = pd.read_csv(outpath + fn, index_col = 'real_name')

    # transform parameters into the original ranges
    # read files containing parameter ranges
    par_range = pd.read_csv(f'{outpath}parameters.tpl', skiprows=1, index_col = 'parameter')
    for col in list(par_range.index):
        val_low, val_up = par_range.loc[col, 'lower'], par_range.loc[col, 'upper']
        pars[col.lower()] = pars[col.lower()]*(val_up - val_low) / 100 + val_low

    # select the parameter columns to plot
    joint_columns1 = [*pars.columns[:3], *pars.columns[4:5]]
    joint_columns2 = pars.columns[5:9]
    joint_columns3 = pars.columns[9:]
    # plot the realizations of parameter ensemble
    k = 0
    for cols in [joint_columns1, joint_columns2, joint_columns3]:
        fig = plt.figure(figsize=(10, 10))
        ax = sns.pairplot(pars[cols], kind='scatter', diag_kind='hist')
        plt.savefig(f'{fig_path}pars_{k}_pairplot_run_0520', dpi=300)
        k += 1

#---------calculate parameter sensitivity-----------------#
from funcs.read_data import file_settings, variables_prep
from funcs.utils import partial_rank
from gp_build import *
fpath = 'year_pbias/'
gp = pickle.load(open(f'{fpath}gp_1_pbias.pkl', "rb"))
x_training = gp.X_train_
y_training = gp.y_train_


# Resample in the ranges where the objective values are above 0
# x_select = x_training[np.where(y_training>0)[0], :]
# x_range = x_select.max(axis=0)
# univariable_temp = [stats.uniform(0, x_range[ii]) for ii in range(0, x_range.shape[0])]
# variable_temp = pya.IndependentMultivariateRandomVariable(univariable_temp)

# visualization the effects of factor fixing
# define the variables for PCE
param_file = file_settings()[-1]
ind_vars, variables = variables_prep(param_file, product_uniform='uniform', dummy=False)
var_trans = AffineRandomVariableTransformation(variables, enforce_bounds=True)
filename = f'{fpath}rankings.json'
param_names = pd.read_csv(param_file, usecols=[2]).values.flatten()
if not os.path.exists(filename):
    order = 2
    interaction_terms = pya.compute_hyperbolic_indices(len(ind_vars), order)
    interaction_terms = interaction_terms[:, np.where(
    interaction_terms.max(axis=0) == 1)[0]]
    sa = pya.analytic_sobol_indices_from_gaussian_process(gp, variables, 
        interaction_terms, ngp_realizations=100, ninterpolation_samples=500, 
            use_cholesky=True, ncandidate_samples=10000, nvalidation_samples=200)
    # np.savez(filename, total_effects=sa['total_effects']['mean'])
    ST = sa['total_effects']['values']
    ST_mean = sa['total_effects']['mean']
    ST_mean = pd.DataFrame(data = ST_mean, index = param_names)
    
    # save the results of parameter rankings
    ST_mean.to_csv(f'{fpath}ST.csv')
    index_sort = partial_rank(ST, ST.shape[1], conf_level=0.95)
    with open(filename, 'w') as fp:
        json.dump(index_sort, fp, indent=2)

else:
    with open(filename, 'r') as fp:
        index_sort_load = json.load(fp)
    index_sort = {}
    for k, v in index_sort_load.items():
        index_sort[int(k)] = index_sort_load[k]

