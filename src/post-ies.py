#! usr/bin/env python

import numpy as np
from numpy.lib import quantile

def coverage_ratio(x_obs, x_pred, quantile_bool=False, quantiles = None):
    """
    This calculates how many observation measurements are contained by prediction 
    using the data for a year.
    Parameters:
    ============
    x_obs: np.ndarray, containing the measurement realizations.
    x_pred: np.ndarray, of predictions.

    return:
    ============
    cr: float, the coverage ratio.
    """
    if quantile_bool:
        assert quantiles != None, "when quantile_bool is True, must provide the quantiles."
        x_pred_min, x_pred_max = np.quantile(x_pred, quantiles)
    else:
        x_pred_min = x_pred.min()
        x_pred_max = x_pred.max()
    n_in = x_obs[(x_obs>=x_pred_min) & (x_obs <= x_pred_max)].shape[0]
    cr = n_in / x_pred.shape[0]
    return cr


def discri_diagram(x_obs, x_pred, thsds):
    """
    This calculates the discrimination diagram score.
    Using the data for a year.
    Parameters:
    ============
    x_obs: np.ndarray, containing the measurement realizations.
    x_pred: np.ndarray, of predictions.
    thsds: np.ndarray, thresholds for categorizing the output.

    return:
    ============
    dd_score: float, the discrimination diagram score.
    """
    x_obs_category = np.zeros_like(x_obs)
    x_pred_category = np.zeros_like(x_pred)
    dd_values = np.zeros_like(x_obs)
    x_obs_category[(x_obs>thsds[0]) & (x_obs<=thsds[1])] = 1
    x_pred_category[(x_pred>thsds[0]) & (x_pred<=thsds[1])] = 1
    x_obs_category[x_obs>thsds[1]] = 2
    x_pred_category[x_pred>thsds[1]] = 2
    dd_values[x_pred_category == x_obs_category] = 1
    dd_score = dd_values.sum()/dd_values.shape[0]
    return dd_score


def average_width(x_obs_all, x_pred_all, quantile_bool=True, quantiles = None):
    """
    This calculates the uncertainty width index. Using the data of all years.
    Parameters:
    ============
    x_obs: np.ndarray, containing the measurement realizations.
    x_pred: np.ndarray, of predictions.

    return:
    ============
    awi: float, the coverage ratio.
    """
    if quantile_bool:
        assert quantiles != None, "when quantile_bool is True, must provide the quantiles." 
        x_pred_min, x_pred_max = np.quantile(x_pred_all, quantiles, axis=0)
        x_obs_min, x_obs_max = np.quantile(x_obs_all, quantiles)
    else:
        x_pred_min = x_pred_all.min(axis=0)
        x_pred_max = x_pred_all.max(axis=0)
        x_obs_min = x_obs_all.min()
        x_obs_max = x_obs_all.max()
    # breakpoint()
    awi = 1 - np.mean(x_pred_max - x_pred_min) / (x_obs_max - x_obs_min)
    awi_annual = 1 - (x_pred_max - x_pred_min) / (x_obs_max - x_obs_min)
    return awi, awi_annual
    

def median_absolute_deviation(x_obs, x_pred):
    """
    Using the data for a year.
    """
    mad = np.round(np.median(np.abs(x_obs - x_pred)), 2)
    return mad

def interval_skill_score(x_obs, x_pred, quantile_bool=True, quantiles = None):
    if quantile_bool:
        assert quantiles != None, "when quantile_bool is True, must provide the quantiles."
        x_pred_min, x_pred_max = np.quantile(x_pred, quantiles)
        x_obs_min, x_obs_max = np.quantile(x_obs, quantiles)
    else:
        x_pred_min = x_pred.min()
        x_pred_max = x_pred.max()
        x_obs_min = x_obs.min()
        x_obs_max = x_obs.max()
        quantiles = [0.001, 0.999]
        
    si = np.zeros_like(x_obs)
    q_width = quantiles[1] - quantiles[0]
    for ii in range(x_obs.shape[0]):
        if (x_obs[ii] >= x_pred_min) & (x_obs[ii] <= x_pred_max):
            si[ii] = x_pred_max - x_pred_min
        elif x_obs[ii] < x_pred_min:
            si[ii] = (x_pred_max - x_pred_min) + 2 * (x_pred_min - x_obs[ii])
        else:
            si[ii] = (x_pred_max - x_pred_min) + 2 * (x_obs[ii] - x_pred_max)

    si_year = si.mean()
    return si_year

def average_interval_skill_score(x_obs_all, x_pred_all, quantile_bool=True, quantiles = None):
    si_all, si_all_clim = np.zeros(x_obs_all.shape[1]), np.zeros(x_obs_all.shape[1])
    
    for year in range(x_obs_all.shape[1]):
        si_all[year] = interval_skill_score(x_obs_all[:, year], x_pred_all[:, year], \
            quantile_bool=quantile_bool, quantiles = quantiles)
        
        si_all_clim[year] = interval_skill_score(x_obs_all[:, year], x_obs_all[:, :], \
            quantile_bool=quantile_bool, quantiles = quantiles)

    iss = 1 - np.mean(si_all) / np.mean(si_all_clim)
    return iss


def inter_qunatile_range(param_value, x_range, quantiles):
    """
    Using all the parameter values.
    """
    param_quantile = np.quantile(param_value, quantiles, axis=0) / x_range
    return param_quantile

import numpy as np
import pandas as pd
import math

# import the annual loads
file_date = '20220119'
fpath = f'../output/work_run_{file_date}/'
fn = '126001A.11.obs.csv'
fn_pars = '126001A.11.par.csv'
fn_meas = '126001A.base.obs.csv'
log_load = False

df = pd.read_csv(fpath + fn, index_col = 'real_name')
df_pars = pd.read_csv(fpath + fn_pars, index_col = 'real_name')
# select results of which the pbias is with 15%
df_meas = pd.read_csv(fpath + fn_meas, index_col = 'real_name')
if log_load:
    df_meas.loc[:, 'din_2009':] = 10**(df_meas.loc[:, 'din_2009':])
    df.loc[:, 'din_2009':] = 10**(df.loc[:, 'din_2009':])
df['average'] = df.loc[:, 'din_2009':'din_2017'].mean(axis=1).values
df_meas['average'] = df_meas.loc[:, 'din_2009':'din_2017'].mean(axis=1).values
df_meas = df_meas.loc[df.index, :]
# obs data
obs_annual = [52.093, 99.478, 44.064, 57.936, 53.449, 21.858, 38.561, 51.843, 14.176]
obs_annual.append(np.round(np.mean(obs_annual), 2))
obs_df = pd.DataFrame(data=np.log10(obs_annual), index = [*np.arange(2009, 2018), 'average'], columns=['Annual loads'])

quantiles = [0.05, 0.95]
awi, awi_annual = average_width(df_meas.values[:, 1:10], df.values[:, 1:10], \
    quantile_bool=True, quantiles = quantiles)
iss = average_interval_skill_score(df_meas.values[:, 1:10], df.values[:, 1:10], \
    quantile_bool=True, quantiles = quantiles)


cr_vals, dd_vals, mad_vals = np.zeros(shape=(9)), np.zeros(shape=(9)), np.zeros(shape=(9))
cols = df_meas.columns
for i in range(9):
    col = cols[i+1]
    cr_vals[i] = coverage_ratio(df_meas.loc[:, col].values, df.loc[:, col].values, \
        quantile_bool=True, quantiles = quantiles)
    dd_vals[i] = discri_diagram(df_meas.loc[:, col].values, df.loc[:, col].values, \
        thsds=np.quantile(df_meas.loc[:, col].values,  [0.25, 0.75]))
    mad_vals[i] = median_absolute_deviation(df_meas.loc[:, col].values, df.loc[:, col].values)

metric = pd.DataFrame(index=obs_df.index[:-1], \
    columns=['Coverage ratio', 'Discrimination diagram', 'Median absolute deviation', 'Width index'], \
        data=np.array([cr_vals, dd_vals, mad_vals, awi_annual]).T)
# metric.to_csv(fpath+'metric.csv')

##-----------------Check the residuals with weights------------------##
weight = np.array([0.01, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13])
weight2 = np.array([0.15, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03])
resd = (df_meas - df)
lsq = (resd.loc[:, 'din_pbias':'din_2017']*weight)**2
lsq.loc[:, 'din_2009':'din_2017'].sum(axis=1).mean()
lsq.sum(axis=1).mean()
lsq.sum(axis=1).std()
lsq.sum(axis=1).max()
lsq.sum(axis=1).min()


import spotpy
# df_temp  = df_meas.copy(deep=True)
# for i in range(9):
#     df_temp.loc[:, cols[i+1]] = 10**(df_meas.loc[:, cols[i+1]])

# for i in df_temp.index:
#     df_meas.loc[i, 'din_pbias'] = spotpy.objectivefunctions.pbias(10**(obs_df.iloc[0:9, :].values.flatten()), df_temp.loc[i, 'din_2009':'din_2017'].values)

# df_meas.to_csv('measurement_ensemble_log.csv')