"""
RUN SOURCE to generate observation_ensemble.csv
"""
import pandas as pd
import numpy as np
import veneer
from veneer.pest_runtime import *
from veneer.manage import start,kill_all_now
import time
import os
import spotpy as sp

from funcs.modeling_funcs import vs_settings, obtain_initials,\
    change_param_values, generate_observation_ensemble, generate_parameter_ensemble,\
        generate_observation_default, modeling_settings

# check if there is a file containing synthetic observations
def run_default_obs(vs, criteria, start_date, end_date, 
    obs_name, retrieve_time, datapath):
    if not os.path.exists(f'{datapath}{obs_name}.csv'):
        load = generate_observation_default(vs_list[0], criteria, start_date, 
            end_date, retrieve_time)
        load.to_csv(f'{datapath}{obs_name}.csv')
    else:
        print(f'{datapath}.csv exists.')

# generate the observation ensemble
def run_obs_ensemble(vs, criteria, start_date, end_date, parameters, 
    obs_ensemble_name, retrieve_time, datapath):
    if not os.path.exists(f'{datapath}{obs_ensemble_name}.csv'):
        load = generate_observation_ensemble(vs_list, criteria, 
            start_date, end_date, parameters, retrieve_time)
        load.to_csv(f'{datapath}{obs_ensemble_name}.csv')
    else:
        print(f'{obs_ensemble_name}.csv exists.')

first_port=15000
num_copies = 8
# define catchment_project and veneer_exe
project_path = 'pest_source/'
# project_path = 'pest_source/'
catchment_project= project_path + '/MW_BASE_RC10.rsproj'

# Setup Veneer
# define paths to veneer command and the catchment project
veneer_path = project_path + 'vcmd45/FlowMatters.Source.VeneerCmd.exe'
#Now, go ahead and start source
processes, ports = start(catchment_project,
                        n_instances=num_copies,
                        ports=first_port,
                        debug=True,
                        veneer_exe=veneer_path,
                        remote=False,
                        overwrite_plugins=True)

NODEs, things_to_record, criteria, start_date, end_date = modeling_settings()
vs_list = vs_settings(ports, things_to_record)

# generate parameter emsenble
datapath = '../../data/'
# nsample = 512
# param_ensemble = 'parameter_ensemble_test.csv'
# generate_parameter_ensemble(nsample, param_ensemble, datapath, seed=88)

# obtain the initial values of parameter 
initial_values = obtain_initials(vs_list[0])

# run to generate observation with default parameter values in the model
print('------------------Generate observation with default parameter values-----------------')
obs_name = 'observation_0918'
retrieve_time = [pd.Timestamp('2009-07-01'), pd.Timestamp('2018-06-30')]
# end_date = '30/06/2018'
run_default_obs(vs_list[0], criteria, start_date, end_date, 
    obs_name, retrieve_time, datapath)

# run to generate observation ensemble with parameter ensemble
print('------------------Generate observation ensemble-----------------')
obs_ensemble_name = 'din_daily_0918'   
parameters = pd.read_csv('parameter_ensemble.csv', index_col='real_name')#[0:8] #OED/src/DIN/
run_obs_ensemble(vs_list, criteria, start_date, end_date, parameters, 
    obs_ensemble_name, retrieve_time, datapath)

# set parameter to the initial values
for vs in vs_list:
    vs = change_param_values(vs, initial_values, fromList=True)

kill_all_now(processes)


