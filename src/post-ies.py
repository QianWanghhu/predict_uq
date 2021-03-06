#! usr/bin/env python
from copy import deepcopy
from matplotlib.pyplot import axis
import numpy as np
from numpy.lib import quantile
import pandas as pd
import math
import spotpy
from spotpy import objectivefunctions

def coverage_ratio(x_obs, x_pred, quantile_bool=False, quantiles = None, raw_cr=False):
    """
    This calculates how many observation measurements are contained by prediction 
    using the data for a year.
    Parameters:
    ============
    x_obs: np.ndarray, containing the measurement realizations.
    x_pred: np.ndarray, of predictions.
    raw_cr: Bool, if true, calculate the coverage ratio, else produce overlap ratio.

    return:
    ============
    cr: float, the coverage ratio.
    """
    if quantile_bool:
        assert quantiles != None, "when quantile_bool is True, must provide the quantiles."
        x_pred_min, x_pred_max = np.quantile(x_pred, quantiles)
    else:
        x_std = x_pred.std()
        x_pred_min = x_pred.mean() - 1.96*x_std
        x_pred_max = x_pred.mean() + 1.96*x_std

    if raw_cr: # calculate coverage ratio
        if (x_obs>=x_pred_min) & (x_obs <= x_pred_max):
            cr = 1 / 9
        else:
            cr = 0
    else: # calculate overlap fraction
        n_in = x_obs[(x_obs>=x_pred_min) & (x_obs <= x_pred_max)].shape[0]
        cr = n_in / (9 * x_obs.shape[0])
    return cr


def average_width(x_obs, x_pred_all, quantile_bool=True, quantiles = None):
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
    obs_sigma = x_obs.std()
    if quantile_bool:
        assert quantiles != None, "when quantile_bool is True, must provide the quantiles." 
        x_pred_min, x_pred_max = np.quantile(x_pred_all, quantiles, axis=0)
        awi = np.mean(x_pred_max - x_pred_min) / obs_sigma
    else:
        x_std = x_pred_all.std(axis=0)
        awi = (3.92 * x_std.mean()) / obs_sigma
    
    return awi

def median_absolute_deviation(x_obs, x_pred_all):
    """
    This calculates the median absolute deviation of predictive uncertainty. Using the data of all years.
    """
    abs_deviation = np.zeros_like(x_pred_all)
    pred_median = np.median(x_pred_all, axis=0)
    for ii in range(x_pred_all.shape[1]):
        abs_deviation[:, ii] = np.abs(x_pred_all[:, ii] - pred_median[ii])
    
    obs_sigma = x_obs.std()
    mad = (np.median(abs_deviation, axis=0)).mean()
    mad = mad / obs_sigma
    return mad
   

def nse_cal(x_obs, x_pred):
    nse_result = np.zeros(x_obs.shape[0])
    for ii in range(nse_result.shape[0]):
        nse_result[ii] = objectivefunctions.nashsutcliffe(x_obs[ii, :], x_pred[ii, :])
    return nse_result


def pbias_cal(x_obs, x_pred):
    pbias_result = np.zeros(x_obs.shape[0])
    for ii in range(pbias_result.shape[0]):
        pbias_result[ii] = objectivefunctions.pbias(x_obs[ii, :], x_pred[ii, :])
    return pbias_result


def r2_cal(x_obs, x_pred):
    r2_result = np.zeros(x_obs.shape[0])
    for ii in range(r2_result.shape[0]):
        r2_result[ii] = objectivefunctions.rsquared(x_obs[ii, :], x_pred[ii, :])
    return r2_result


def measure_uncertainty(x_obs_log, std_log, nsample=10000):
    from SALib.sample import latin
    problem = {
        'num_vars': x_obs_log.shape[0],
        'names': [f'din_{year}' for year in np.arange(2009, 2018)],
        'bounds': np.array([[x_obs_log[ii], std_log] for ii in range(x_obs_log.shape[0])]),
        'dists': ['norm']*9
        }
    obs_ensemble = latin.sample(problem, nsample, seed=88)
    obs_ensemble = 10**obs_ensemble
    for ii in range(x_obs_log.shape[0]):
        obs_ensemble[:, ii] = obs_ensemble[:, ii] - 10 ** (x_obs_log[ii])

    return obs_ensemble

def pred_uncertainty(obs_ensemble, x_pred):
    total_uncertainty = np.zeros(shape=(100 * x_pred.shape[0], x_pred.shape[1]))
    for ii in range(x_pred.shape[0]):
        for jj in range(x_pred.shape[1]):
            total_uncertainty[ii*100 : (ii+1)*100, jj] = \
                obs_ensemble[ii*100 : (ii+1)*100, jj] + x_pred[ii, jj]
    
    return total_uncertainty


def inter_qunatile_range(param_value, quantiles):
    """
    Using all the parameter values.
    """
    param_quantiles = np.quantile(param_value, quantiles, axis=0)
    iqr = []
    for ii in range(param_value.shape[1]):
        iqr.append((param_quantiles[1, ii] - param_quantiles[0, ii]) / 100)

    return iqr


def par_range(param_value):
    """
    Using all the parameter values.
    """
    param_range = (param_value.max(axis=0) - param_value.min(axis=0)) / 100 
    return param_range

def direction_change_pars(par_df_pre, par_df_fix, pars_fixed):
    """
    Calculate the direction of parameters variation when fixing a parameter.
    Parameters:
    ===========
    par_df_pre: pd.DataFrame, of parameter estimates for the previous iteration.
    par_df_fix: pd.DataFrame, of parameter estimates for the current iteration.
    Return:
    ===========
    direction_df, pd.DataFrame, values indcate the direction that a parameter moves towards.
    """
    index_com = [ii for ii in par_df_pre.index if ii in par_df_fix.index]
    delta_pars = (par_df_fix.loc[index_com, :] - par_df_pre.loc[index_com, :]).values
    delta_direction = np.where(delta_pars>0, 1, 0)
    direction_df = pd.DataFrame(index = index_com, columns = par_df_fix.columns, data=delta_direction)
    direction_df[pars_fixed] = 0
    direction_df['total'] = direction_df.sum(axis=1)
    direction_df[pars_fixed] = (par_df_fix.loc[index_com, pars_fixed] - par_df_pre.loc[index_com, pars_fixed]).values
    return direction_df
    

# import the annual loads
# file_dates = ['20220118_full', '20220123', '20220125', '20220126', '20220128', '20220129', '20220130', '20220201', '20220203', '20220204']
# fn_index = [26, 12, 24, 8, 12, 12, 8, 11, 9, 7]
# pars_fix = ['odwc', 'gfdwc', 'drp', 'cdwc', ['godwc', 'fdwc'], 
#     ['gfemc', 'dwc'], 'oemc', 'cemc', ['drf', 'femc']]

file_dates = ['20220118_full', '20220123', '20220310', '20220311', '20220312', '20220313', '20220314', '20220315', '20220316']#, '20220126', '20220128', '20220129', '20220130', '20220201', '20220203', '20220204']
fn_index = [26, 12, 7, 7, 12, 7, 11, 11, 7]#, 8, 12, 12, 8, 11, 9, 7]
pars_fix = ['odwc', 'drp', 'godwc', 'gfdwc', 'dwc', ['fdwc', 'oemc'], ['gfemc', 'cdwc'], ['femc', 'drf', 'cemc']]
log_load = True
par_columns = pd.read_csv('parameter_ensemble.csv', index_col = 'real_name').columns
iqr_metrics_each_run = pd.DataFrame(index=np.arange(len(file_dates)), columns=par_columns)
quantile25_metrics_each_run = pd.DataFrame(index=np.arange(len(file_dates)), columns=par_columns)
quantile75_metrics_each_run = pd.DataFrame(index=np.arange(len(file_dates)), columns=par_columns)
range_metrics_each_run = pd.DataFrame(index=np.arange(len(file_dates)), columns=par_columns)


for run_id in range(len(file_dates)):
    print(f"---------------Read data for run_id: {run_id}---------------------")
    file_date = file_dates[run_id]
    fn_ind = fn_index[run_id]
    fpath = f'../output/test/work_run_{file_date}/'
    fn = f'126001A.{fn_ind}.obs.csv'
    fn_pars = f'126001A.{fn_ind}.par.csv'
    fn_meas = '126001A.base.obs.csv'
    df = pd.read_csv(fpath + fn, index_col = 'real_name')
    df_pars = pd.read_csv(fpath + fn_pars, index_col = 'real_name')
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
    obs_df = pd.DataFrame(data=np.log10(obs_annual), 
        index = [*np.arange(2009, 2018), 'average'], columns=['Annual loads'])

    ##-------------------Calculate metrics using pars uncertainty: Average width----------------##
    quantiles = [0.025, 0.975]
    awi = average_width(np.array(obs_annual)[:-1], df.values[:, 1:10], \
        quantile_bool=True, quantiles = quantiles)

    ##-------------------Calculate metrics using pars uncertainty: MAD and 95%CIs from std----------------##
    mad = median_absolute_deviation(np.array(obs_annual)[:-1], df.values[:, 1:10])
    std_ci = average_width(np.array(obs_annual)[:-1], df.values[:, 1:10], \
        quantile_bool=False)

    ##-------------------Calculate metrics using pars uncertainty: Coverage Ratio----------------##
    print("---------------------Calculate metrics using pars uncertainty: Coverage Ratio---------------------")
    cr_vals = np.zeros(shape=(9, 4))
    cols = df_meas.columns
    for i in range(9):
        col = cols[i+1]
        cr_vals[i, 0] = coverage_ratio(np.array(obs_annual)[i], df.loc[:, col].values, \
            quantile_bool=True, quantiles = quantiles, raw_cr=True)
        cr_vals[i, 1] = coverage_ratio(df_meas.loc[:, col].values, df.loc[:, col].values, \
            quantile_bool=True, quantiles = quantiles, raw_cr=False)
        cr_vals[i, 2] = coverage_ratio(np.array(obs_annual)[i], df.loc[:, col].values, \
            quantile_bool=False, raw_cr=True)
        cr_vals[i, 3] = coverage_ratio(df_meas.loc[:, col].values, df.loc[:, col].values, \
            quantile_bool=False, raw_cr=False)

    metric = pd.DataFrame(index=np.arange(1), \
        columns=['CR_perc', 'OLF_perc', 'CR_std', 'OLF_std', 'AWI_perc', 'MAD', 'AWI_std'], \
            data=np.array([*cr_vals.sum(axis=0), awi, mad, std_ci]).reshape(1, 7))
    metric.to_csv(fpath+'metric_param_unc.csv')

    ##-------------------Calculate metrics using pars uncertainty: Objective functions----------------##
    print("---------------------Calculate metrics using pars uncertainty: Objective functions---------------------")
    objfunc_results = np.zeros(shape=(df.shape[0], 3))
    obj_funcs = [nse_cal, pbias_cal, r2_cal]
    for ii in range(len(obj_funcs)):
        objfunc_results[:, ii] = obj_funcs[ii](df_meas.values, df.values)
    objfunc_results_df = pd.DataFrame(objfunc_results, index=np.arange(df.shape[0]), columns = ['NSE', 'PBIAS', 'R2'])
    objfunc_results_df.to_csv(fpath+'objective_functions.csv')

    ##-------------------Calculate total uncertainty: Predictive uncertainty----------------##
    print("---------------------Calculate total uncertainty---------------------")
    obs_ensemble = measure_uncertainty(obs_df.values[:-1].flatten(), np.round(1/11.13, 2), nsample=10000)
    total_uncertainty = pred_uncertainty(obs_ensemble, df.loc[:, cols[1:-1]].values)
    total_df = pd.DataFrame(data=total_uncertainty, index=np.arange(total_uncertainty.shape[0]), columns = cols[1:-1])
    total_df.index.name = 'real_name'
    total_df.to_csv(fpath+'total_uncertainty.csv')

    ##-------------------Calculate metrics using total uncertainty: Coverage Ratio----------------##
    awi_total = average_width(np.array(obs_annual)[:-1], total_uncertainty, \
        quantile_bool=True, quantiles = quantiles)
    
    ##-------------------Calculate metrics using pars uncertainty: MAD and 95%CIs from std----------------##
    mad_total = median_absolute_deviation(np.array(obs_annual)[:-1], total_uncertainty)
    std_ci_total = average_width(np.array(obs_annual)[:-1], total_uncertainty, \
        quantile_bool=False)

    ##-------------------Calculate metrics using total uncertainty: Coverage Ratio----------------##
    print("---------------------Calculate metrics using total uncertainty: Coverage Ratio---------------------")
    cr_vals_total = np.zeros(shape=(9, 4))
    cols = total_df.columns
    for i in range(9):
        col = cols[i]
        cr_vals_total[i, 0] = coverage_ratio(np.array(obs_annual)[i], total_df.loc[:, col].values, \
            quantile_bool=True, quantiles = quantiles, raw_cr=True)
        cr_vals_total[i, 1] = coverage_ratio(df_meas.loc[:, col].values, total_df.loc[:, col].values, \
            quantile_bool=True, quantiles = quantiles, raw_cr=False)
        cr_vals_total[i, 2] = coverage_ratio(np.array(obs_annual)[i], total_df.loc[:, col].values, \
            quantile_bool=False, raw_cr=True)
        cr_vals_total[i, 3] = coverage_ratio(df_meas.loc[:, col].values, total_df.loc[:, col].values, \
            quantile_bool=False, raw_cr=False)

    metric_total = pd.DataFrame(index=np.arange(1), \
        columns=['CR_perc', 'OLF_perc', 'CR_std', 'OLF_std', 'AWI_perc', 'MAD', 'AWI_std'], \
            data=np.array([*cr_vals_total.sum(axis=0), awi_total, mad_total, std_ci_total]).reshape(1, 7))
    metric_total.to_csv(fpath+'metric_total_unc.csv')

    ##-------------------Calculate metrics using total uncertainty: Coverage Ratio----------------##
    # Calculate the parameter changes
    print(f"---------------------Calculate the parameter changes---------------------")
    par_ranges = pd.read_csv('parameters.tpl', skiprows=1, index_col='parameter')
    df_pars_convert = df_pars.copy(deep=True)
    for par in par_ranges.index:
        df_pars_convert[par.lower()] = (par_ranges.loc[par, 'upper'] - par_ranges.loc[par, 'lower']) * df_pars[par.lower()] / 100 + par_ranges.loc[par, 'lower']
    
    iqr_metrics_each_run.loc[run_id, :] = inter_qunatile_range(df_pars.values, [0.25, 0.75]) 
    range_metrics_each_run.loc[run_id, :] = par_range(df_pars.values)
    quantile25_metrics_each_run.loc[run_id, :] = df_pars_convert.quantile(0.25).values
    quantile75_metrics_each_run.loc[run_id, :] = df_pars_convert.quantile(0.75).values

    ##---------------Calculate the direction of parameters changing--------------------##
    if run_id > 0:
        print(f"---------------------Calculate the direction of parameters changing---------------------")
        direction_df = direction_change_pars(df_pars_convert_pre, df_pars_convert, pars_fix[run_id-1])
        direction_df.to_csv(fpath+'parameter_direction.csv')

    df_pars_convert_pre = df_pars_convert
fp_fig = '../output/test/figs/'
iqr_metrics_each_run.to_csv(f'{fp_fig}parameter_iqr.csv')
range_metrics_each_run.to_csv(f'{fp_fig}parameter_range.csv')
quantile25_metrics_each_run.to_csv(f'{fp_fig}parameter_quantile25.csv')
quantile75_metrics_each_run.to_csv(f'{fp_fig}parameter_quantile75.csv')