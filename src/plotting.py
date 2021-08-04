# import packages
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

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
