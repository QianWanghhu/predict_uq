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

        axes[1].xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
        axes[0].set_ylabel(ylabels[0])
        axes[1].set_ylabel(ylabels[1])
        return fig

    df = pd.DataFrame(index = np.arange(2009, 2018))
    for fn in os.listdir(fpath+fpath_fre):
        """
        This script is used to visualize the outputs of load uncertainty estimation.
        """
        if '70' in fn:
            print(fn)
            df_temp = pd.read_csv(f'{fpath}{fpath_fre}{fn}', index_col='Unnamed: 0')
            df.loc[:, f'{fn[:-4]}_min'] = df_temp.min().values
            df.loc[:, f'{fn[:-4]}_max'] = df_temp.max().values

    df.loc[:, 'DIN_base'] = pd.read_csv(fpath+'DIN_base.csv', usecols=[3]).values
    df.loc[:, 'NOX_base'] = pd.read_csv(fpath+'NOX_base.csv', usecols=[3]).values
    cols = df.columns   
    tick_spacing = 1
    unc_dict = {}
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

    ##---------------------------------Plot load uncertainty----------------------------##
    sns.set_style('whitegrid')
    unc_settings = ['din_nocov', 'din_cov', 'nox_nocov', 'nox_cov']
    fig_names = ['DIN_load_nocov.png', 'DIN_load_cov.png', 'NOX_load_nocov.png', 'NOX_load_cov.png']
    for unc_setting, fig_name in unc_settings, fig_names:
        fig = plot_unc(df, height=df.DIN_base, unc_cov= unc_setting, ylabels=['DIN Load (t/yr)', 'Load ratio (Bounds/Best)'])
        plt.savefig(fpath+fig_name, format='png', dpi=300)

##--------------------------------Generate observation ensembles--------------------------##
def produce_meas_ensemble(fpath, fpath_freq, fname, n_sigma):
    def generate_measurement_realization(upper, lower, n_sample, n_sigma):
        """
        The function produces the measurement realizations using the uncertainty ranges.
        A log-normal distribution is assumed.   
        Parameters:
        ===========
        upper, lower: float, the upper and lower bounds of uncertainty as percentages.
        n_sample: int, the number of realizations to generate.
        meas_best: float, the best measurement estimate.

        Returns:
        =========== 
        meas_realizations: np.ndarray, of size (n_samples, 1)
        """
        
        log_sigma = (np.log(upper) - np.log(lower)) / n_sigma
        mu = 1
        log_mu = np.log(mu)
        lognorm_fitted = lognorm(s=log_sigma, loc=0, scale=mu)
        norm_dist = norm(loc = log_mu, scale = log_sigma)
        x = norm_dist.rvs(size = n_sample, random_state=612)
        meas_realizations = np.exp(x)
        return meas_realizations

    uncertain_df = pd.read_csv(fpath+fpath_freq+fname, index_col='Year')
    unc_bounds = uncertain_df.loc[:, 'din_cov_ratio_lower':'din_cov_ratio_upper'].values
    meas_base = pd.read_csv(fpath+'DIN_base.csv', usecols=[3]).values

    cols = [f'DIN_{year}' for year in range(2009, 2018)]
    unc_realizations = pd.DataFrame(index=np.arange(100), columns=cols)
    meas_realizations = pd.DataFrame(index=np.arange(100), columns=cols)
    for ii in np.arange(9):
        unc_realizations.loc[:, cols[ii]] = generate_measurement_realization(unc_bounds[ii, 1], \
            unc_bounds[ii, 0], 100, n_sigma=n_sigma)
        meas_realizations.loc[:, cols[ii]] = unc_realizations.loc[:, cols[ii]] * meas_base[ii]
    meas_realizations.loc[:, 'pbias'] = 100 * (meas_realizations.sum(axis=1) - np.sum(meas_base)) / np.sum(meas_base)
    meas_realizations = meas_realizations[['pbias', *cols]]

    fig, axes = plt.subplots(3, 3, figsize=(6*3, 5*3))
    k = 0
    for ii in range(3):
        for jj in range(3):
            meas_realizations[cols[k]].hist(ax=axes[ii, jj])
            axes[ii, jj].set_title(f'year {cols[k]}')
            k += 1
    plt.savefig(fpath+ fpath_freq + 'hist_meas_5.6_sigma.pdf', format='pdf', dpi=300)
    meas_stats = np.round(np.array([meas_realizations.min().values, meas_realizations.max().values]).T, 2)
    np.savetxt(fpath + fpath_fre +'meas_5.6_sigma.txt', meas_stats)
    
    # save the measurement realizations
    unc_realizations.to_csv(fpath+fpath_freq+'uncertainty_realizations.csv')
    meas_realizations.to_csv(fpath+fpath_freq+'measurement_realizations.csv')
##--------------------------------Generate observation ensembles--------------------------##
fpath = 'output/'
fpath_fre = 'resample_freq70/'
produce_meas_ensemble(fpath, fpath_fre, fname = 'uncertainty_ranges.csv', n_sigma=5.6)