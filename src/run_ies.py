from veneer.manage import start
import subprocess
import shutil
import os

parent_dir = os.getcwd()  
job_name = 'work'
pst_file = '126001A.pst'
catchment_project= parent_dir + '\\pest_source\\MW_BASE_RC10.rsproj'
pest_path= parent_dir + '\\pest_source' 
print('pest path ',pest_path) 

python_path = 'C:\\UserData\\Qian\\anaconda'
os.environ['PATH'] = os.environ['PATH'] + ';' + pest_path
os.environ['PATH'] = os.environ['PATH'] + ';' + python_path
print(os.environ['PATH'])  

# Setup Veneer
# define paths to veneer command and the catchment project
veneer_path = 'pest_source\\vcmd45\\FlowMatters.Source.VeneerCmd.exe'


# Number of instances to open
num_copies=4     # Important - set this to be a number ~ the number of CPU cores in your system!
first_port=15000

#Now, go ahead and start source
processes,ports = start(catchment_project,
                        n_instances=num_copies,
                        ports=first_port,
                        debug=True,
                        veneer_exe=veneer_path,
                        remote=False,
                        overwrite_plugins=True)

pest_master_version = 'ipestpp-ies.exe'
pest_slave_version = 'ipestpp-ies.exe'  

# Set up job directories
instruction_file = 'output.ins'
template_file = 'parameters.tpl'

supp_files = ['126001A.pst','126001A.csv', 'observation_ensemble_0918.csv' , 'parameter_ensemble.csv',
			  'run_source.py', 'Plugins.xml']
            
# parent_dir = os.getcwd()
os.makedirs(job_name,exist_ok=True)
shutil.copyfile(instruction_file, job_name + '/' + instruction_file)
shutil.copyfile(template_file, job_name + '/' + template_file)
shutil.copyfile(pst_file, job_name + '/' + pst_file)
for file in supp_files:
    shutil.copyfile(file, job_name + '/' + file)

model_func_file = 'modeling_funcs.py'
shutil.copyfile('funcs/' + model_func_file, job_name + '/' + model_func_file)

os.chdir(job_name)
for port in ports :
    dir_name = 'Slave_' + str(port)
    os.makedirs(dir_name,exist_ok=True)
    shutil.copyfile(instruction_file, dir_name + '/' + instruction_file)
    shutil.copyfile(template_file, dir_name + '/' + template_file)
    shutil.copyfile(pst_file, dir_name + '/' + pst_file)
    shutil.copyfile(model_func_file, dir_name + '/' + model_func_file)
    for file in supp_files:
        shutil.copyfile(file, dir_name + '/' + file)
    veneer_connection = dir_name + '/veneer_connection.txt'
    with open(veneer_connection, 'w') as f:
        f.write(str(port))

os.chdir(parent_dir)


# Start PEST  
os.chdir(job_name)
slave_processes = {}
pest_master_string = pest_master_version + '  ' +pst_file+ ' /h :4001 '
pest_slave_string = pest_slave_version + '  ' +pst_file+ ' /h localhost:4001'
current_dir = os.getcwd()
for port in ports:
    dir_name = 'Slave_' + str(port)
    os.chdir(dir_name)
    slave_processes[port] = subprocess.Popen(["start", "cmd", "/k",pest_slave_string], shell=True  )
    # slave_processes[port] = Popen(pest_slave_string,shell=True, stdout=PIPE, stderr=PIPE)
    os.chdir(current_dir)
#subprocess.Popen(["start", "cmd", "/k",pest_master_string], shell=True  )
os.system(pest_master_string )
os.chdir(parent_dir)
# subprocess.run(pest_master_string)