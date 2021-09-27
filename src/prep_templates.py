"""
The script is used to run commands for OED including the run for parameter,observation ensemble,
adding noise and running IES.
"""
# import packages
import os
import numpy as np
import pandas as pd
import spotpy as sp
from funcs.prep_template import pre_ies_obs_template, format_obs_ensemble, add_noise_obs, generate_obsnme

# format ies inputs
datapath = '../../data/'
obgnme = 'DIN'

fname = f'{datapath}observation.csv'
fname_obs_pst = 'obs_pst_noise10pct.csv'


# if not os.path.exists(fname_obs_pst):
#     print('Preparing the observation file for ies')
#     load = pd.read_csv(fname, index_col='Unnamed: 0')
#     obs_pst_name = f'{datapath}{fname_obs_pst}'
#     load = add_noise_obs(load, noise_level=0.1)
# # End if-else
def observation_ensemble(datapath, ies_obs_name):
    """Prepare the observation ensemble file for ies"""
    if not os.path.exists(ies_obs_name):
        obs_ensemble = pd.read_csv(f'{datapath}din_daily_0918.csv', index_col='Unnamed: 0')
        obs_ensemble.index = pd.to_datetime(obs_ensemble.index)
        obs_ensemble_month = obs_ensemble.resample('M').sum()
        obs_ensemble_annual = pd.DataFrame(index=np.arange(obs_ensemble.index[0].year, obs_ensemble.index[-1].year),
            columns = obs_ensemble.columns)
        for i in range(obs_ensemble_annual.shape[0]):
            obs_ensemble_annual.iloc[i] = obs_ensemble_month.iloc[i*12:(i+1)*12, :].sum(axis=0)
     
        obs_ensemble_ies = format_obs_ensemble(obs_ensemble_annual, obgnme)

        return obs_ensemble_ies, obs_ensemble
    else:
        print('The file exists.')
    # End if-else
    # End observation_ensemble()

def process_monitor(datapath):
    """
    convert the timeseries (monitoring and modeling) to a required temporal scale.
    read the monitoring data 
    (Note that the daily loads is of unit t and need to be converted to kg).
    return the observation section in the control file for PEST.
    """
    fname = datapath + 'low_interp_flow.csv'
    date_range = pd.to_datetime(['2009/07/01', '2018/06/30'])
    df = pd.read_csv(fname, index_col='Date')
    df.index = pd.to_datetime(df.index)
    df = df.loc[date_range[0]:date_range[1], :].filter(items=[df.columns[0]]).apply(lambda x: 1000 * x)
    # resample to obtain monthly loads
    df_month = df.resample('M').sum()
    df_annual = pd.DataFrame(index = ['pbias', *np.arange(date_range[0].year, date_range[1].year)], columns=df.columns)
    df_annual.iloc[0, :] = 0.0
    for i in range(df_annual.shape[0]-1):
        df_annual.iloc[i+1, :] = df_month.iloc[i*12: (i+1)*12, :].sum()
    # End for

    obs_index = list(df_annual.index)
    obsnme = generate_obsnme(obs_index, obgnme)
    obsval = df_annual.values
    fname_obs_pst = 'monitor_annual.csv'

    weight = np.array([1 / (df_annual.shape[0] - 1)]).repeat(df_annual.shape[0])
    weight[0] = 1.5
    pre_ies_obs_template(obsval, obsnme, obgnme, weight).to_csv(fname_obs_pst, sep=' ')
    return obsnme, df
# End process_monitor()

def write_instruction(fname, obsnme):
    with open(fname, 'w') as f:
        f.write('pif # \n')
        f.write(f'l2    [{obsnme[0]}]1:25 \n')
        f.write(f'l3    [{obsnme[1]}]1:25 \n')
        for i in obsnme[2:]:
            f.write(f'l1    [{i}]1:25 \n')
            
        f.close()

def calculate_obs_objetive(monitor, obs_ensemble_ies, obs_ensemble):
    """
    Calculate objectives with monitoring data and observation_ensemble.
    """
    objective_values = np.array([])
    # import pdb; pdb.set_trace()
    monitor_comp = monitor[monitor.columns[0]]
    for col in obs_ensemble.columns:
        objective_values = np.append(objective_values, sp.objectivefunctions.pbias(obs_ensemble[col], monitor_comp))
    
    obs_ensemble_ies.insert(0, 'DIN_pbias', objective_values)
    return obs_ensemble_ies
# write the observation section
obsnme, df = process_monitor(datapath)

# write the corresponding out instruction file for PEST
write_instruction('output_0918.ins', obsnme)

# write observation ensemble for ies
ies_obs_name = 'observation_ensemble_0918.csv'
obs_ensemble_ies, obs_ensemble = observation_ensemble(datapath, ies_obs_name)

obs_ensemble_ies = calculate_obs_objetive(df, obs_ensemble_ies, obs_ensemble)
obs_ensemble_ies.to_csv(ies_obs_name)