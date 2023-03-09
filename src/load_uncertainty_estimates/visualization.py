#!/usr/bin/python 

import pandas as pd
import numpy as np
import matplotlib as mpl
import seaborn as sns
import matplotlib.pyplot as plt
import os
import matplotlib.ticker as ticker
from scipy.stats import lognorm, norm
import spotpy
mpl.rc('xtick', labelsize=16) 
mpl.rc('ytick', labelsize=16)
mpl.rc('axes', labelsize=20)
mpl.rc('axes', titlesize=20) 

def calc_unc_min_max(fpath, fpath_fre):
    def plot_unc(df, height, unc_cov, ylabels):
        fig, axes = plt.subplots(nrows=1, ncols=2, sharey=False, sharex=True, figsize=(8*2, 6))
        axes[0].bar(x=df.index, height=height, yerr=unc_dict[unc_cov])
        axes[0].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
        axes[1].plot(df.index, unc_dict[f'{unc_cov}_ratio'][0], linestyle='', color='b', marker='o')
        axes[1].plot(df.index, unc_dict[f'{unc_cov}_ratio'][1], linestyle='', color='r', marker='o')
        axes[1].hlines(1.0, xmin=df.index[0]-0.5, xmax=df.index[-1], linestyle='--', color='grey')
        axes[0].set_title('(a)')
        axes[1].set_title('(b)')
        axes[0].set_ylim([0, 160])
        axes[1].set_ylim([0.4, 1.6])

        axes[1].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
        axes[0].set_ylabel(ylabels[0])
        axes[1].set_ylabel(ylabels[1])
        return fig

    df = pd.DataFrame(index = np.arange(2009, 2018))
    for fn in os.listdir(fpath+fpath_fre):
        """
        This script is used to visualize the outputs of load uncertainty estimation.
        """
        if '100' in fn:
            print(fn)
            df_temp = pd.read_csv(f'{fpath}{fpath_fre}{fn}', index_col='Unnamed: 0')
            df.loc[:, f'{fn[:-4]}_min'] = df_temp.min().values
            df.loc[:, f'{fn[:-4]}_max'] = df_temp.max().values

    df.loc[:, 'DIN_base'] = pd.read_csv(fpath+'DIN_base.csv', usecols=[3]).values
    df.loc[:, 'NOX_base'] = pd.read_csv(fpath+'NOX_base.csv', usecols=[3]).values
    cols = df.columns   
    tick_spacing = 1
    unc_dict = {}
    unc_settings = ['din_nocov', 'din_cov', 'nox_nocov', 'nox_cov'] # ['din_cov'] #

    if not (os.path.exists(f'{fpath}{fpath_fre}uncertainty_ranges.csv')):
        ##---------------------------------Prepare load uncertainty for plot----------------------------##
        unc_dict['din_cov'] = [(df.DIN_base - df[cols[0]]).values, (df[cols[5]] - df.DIN_base).values]
        unc_dict['din_cov_ratio'] = [df[cols[0]].values, df[cols[5]].values] / df.DIN_base.values
        unc_dict['din_nocov'] = [(df.DIN_base - df[cols[2]]).values, (df[cols[7]] - df.DIN_base).values]
        unc_dict['din_nocov_ratio'] = [df[cols[2]].values, df[cols[7]].values] / df.DIN_base.values

        unc_dict['nox_cov'] = [(df.NOX_base - df[cols[8]]).values, (df[cols[13]] - df.NOX_base).values]
        unc_dict['nox_cov_ratio'] = [df[cols[8]].values, df[cols[13]].values] / df.NOX_base.values
        unc_dict['nox_nocov'] = [(df.NOX_base - df[cols[10]]).values, (df[cols[15]] - df.NOX_base).values]
        unc_dict['nox_nocov_ratio'] = [df[cols[10]].values, df[cols[15]].values] / df.NOX_base.values

        unc_df = pd.DataFrame(index=np.arange(2009, 2018))
        for k, v in unc_dict.items():
            unc_df.loc[:, [k+'_lower', k+'_upper']] = np.array(v).T
        unc_df.index.name = 'Year'
        unc_df.to_csv(f'{fpath}{fpath_fre}uncertainty_ranges.csv')
    else:
        unc_df = pd.read_csv(f'{fpath}{fpath_fre}uncertainty_ranges.csv', index_col='Year')
        for key in unc_settings:
            breakpoint()
            unc_dict[key] = [unc_df[key+'_lower'], unc_df[key+'_upper']]
            unc_dict[key+'_ratio'] = [unc_df[key+'_ratio_lower'], unc_df[key+'_ratio_upper']]

    ##---------------------------------Plot load uncertainty----------------------------##
    sns.set_style('whitegrid')
    fig_names = ['DIN_load_nocov2.png', 'DIN_load_cov2.png', 'NOX_load_nocov2.png', 'NOX_load_cov2.png'] #['DIN_load_cov2.png'] # 
    for unc_setting, fig_name in zip(unc_settings, fig_names):
        fig = plot_unc(df, height=df.DIN_base, unc_cov= unc_setting, ylabels=['DIN Load (t/yr)', 'Load ratio (Bounds/Best)'])
        plt.savefig(fpath+fig_name, format='png', dpi=300)

fpath = 'output/'
fpath_fre = 'resample_freq100/'
calc_unc_min_max(fpath, fpath_fre)