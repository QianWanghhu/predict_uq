#!/usr/bin/env python
import subprocess
import shutil
import os
import numpy as np

parent_dir = os.getcwd()  
job_name = 'work'
num_copies = 6
pest_path= parent_dir + '\\pest_source' 
print('pest path ',pest_path) 

python_path = 'C:\\UserData\\Qian\\anaconda'
os.environ['PATH'] = os.environ['PATH'] + ';' + pest_path
os.environ['PATH'] = os.environ['PATH'] + ';' + python_path
print(os.environ['PATH'])  


pest_master_version = 'ipestpp-ies.exe'
pest_slave_version = 'ipestpp-ies.exe'  

# Set up job directories
instruction_file = 'example_output.ins'
template_file = 'example_parameters.tpl'
pst_file = 'example.pst'

supp_files = ['example_func.py', 'example_ensemble.csv', 'parameter_ensemble.csv']
            
parent_dir = os.getcwd()
os.makedirs(job_name,exist_ok=True)
shutil.copyfile(instruction_file, job_name + '/' + instruction_file)
shutil.copyfile(template_file, job_name + '/' + template_file)
shutil.copyfile(pst_file, job_name + '/' + pst_file)
for file in supp_files:
    shutil.copyfile(file, job_name + '/' + file)

os.chdir(job_name)
for port in np.arange(num_copies) :
    dir_name = 'Slave_' + str(port)
    os.makedirs(dir_name,exist_ok=True)
    shutil.copyfile(instruction_file, dir_name + '/' + instruction_file)
    shutil.copyfile(template_file, dir_name + '/' + template_file)
    shutil.copyfile(pst_file, dir_name + '/' + pst_file)
    for file in supp_files:
        shutil.copyfile(file, dir_name + '/' + file)

os.chdir(parent_dir)


# Start PEST  
os.chdir(job_name)
slave_processes = {}
pest_master_string = pest_master_version + '  ' +pst_file+ ' /h :4001 '
pest_slave_string = pest_slave_version + '  ' +pst_file+ ' /h localhost:4001'
current_dir = os.getcwd()
for port in np.arange(num_copies):
    dir_name = 'Slave_' + str(port)
    os.chdir(dir_name)
    slave_processes[port] = subprocess.Popen(["start", "cmd", "/k",pest_slave_string], shell=True  )
    # slave_processes[port] = Popen(pest_slave_string,shell=True, stdout=PIPE, stderr=PIPE)
    os.chdir(current_dir)
#subprocess.Popen(["start", "cmd", "/k",pest_master_string], shell=True  )
os.system(pest_master_string )
os.chdir(parent_dir)
# subprocess.run(pest_master_string)