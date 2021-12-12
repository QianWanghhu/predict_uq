#!/usr/bin/env python
from veneer.manage import start
import subprocess
import shutil
import os

parent_dir = os.getcwd()  
job_name = 'work'
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

supp_files = ['example_func.py']
            
# parent_dir = os.getcwd()
os.makedirs(job_name,exist_ok=True)
shutil.copyfile(instruction_file, job_name + '/' + instruction_file)
shutil.copyfile(template_file, job_name + '/' + template_file)
shutil.copyfile(pst_file, job_name + '/' + pst_file)
for file in supp_files:
    shutil.copyfile(file, job_name + '/' + file)

os.chdir(parent_dir)
# Start PEST  
pest_master_string = pest_master_version + '  ' +pst_file+ ' /h :4001 '
current_dir = os.getcwd()
os.system(pest_master_string )
# subprocess.run(pest_master_string)