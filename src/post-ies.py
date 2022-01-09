#! usr/bin/env python

import numpy as np

def coverage_ratio(x_obs, x_pred, quantile_bool=True, quantiles = None):
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
        x_pred_min, x_pred_max = np.quantile(x_pred, 0.025, 0.975)
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
        x_pred_min, x_pred_max = np.quantile(x_pred_all, quantiles)
        x_obs_min, x_obs_max = np.quantile(x_obs_all, quantiles)
    else:
        x_pred_min = x_pred_all.min()
        x_pred_max = x_pred_all.max()
        x_obs_min = x_obs_all.min()
        x_obs_max = x_obs_all.max()
    
    awi = 1 - np.mean(x_pred_max - x_pred_min) / np.mean(x_obs_max - x_obs_min)
    return awi
    

def mean_absolute_deviation(x_obs, x_pred):
    """
    Using the data for a year.
    """
    mad = np.mean(np.abs(x_obs - x_pred))
    return mad

def interval_skill_score(x_obs, x_pred, x_best, quantile_bool=True, quantiles = None):
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
    for ii in range(x_obs.shape[0]):
        if (x_obs[ii] >= x_pred_min) & (x_obs[ii] < x_pred_max):
            si[ii] = x_pred_max - x_pred_min
        elif x_obs[ii] < x_pred_min:
            si[ii] = x_pred_max - x_pred_min + 2 / (1 - quantiles[0]) * (x_pred_min - x_obs[ii])
        else:
            si[ii] = x_pred_max - x_pred_min + 2 / (1 - quantiles[0]) * (x_obs[ii] - x_pred_max)

    si_year = si.mean()
    si_obs = x_obs_max - x_obs_min
    return si_year, si_obs

def average_interval_skill_score(x_obs_all, x_pred_all, x_best_all, quantile_bool=True, quantiles = None):
    si_all, si_obs_all = np.zeros_like(x_best_all), np.zeros_like(x_best_all)
    
    for year in x_obs_all.shape[1]:
        si_all[year], si_obs_all[year] = interval_skill_score(x_obs_all[:, year], x_pred_all[:, year], \
            x_best_all[:, year], quantile_bool=quantile_bool, quantiles = quantiles)

    iss = 1 - np.mean(si_all) / np.mean(si_obs_all)
    return iss


def inter_qunatile_range(param_value, x_range, quantiles):
    """
    Using all the parameter values.
    """
    param_quantile = np.quantile(param_value, quantiles, axis=0) / x_range
    return param_quantile

