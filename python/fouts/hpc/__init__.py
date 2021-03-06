"""
  Contains functions to set up the environment on the a HPC server to run scripts in the 'simulate' 
  package
  
  @author Tom Close 
"""
import sys
if float(sys.version[0:3]) < 2.7:
    raise Exception("This script requires python version 2.7 or greater, you should add python 2.7 to your path (i.e. PATH=$PATH:/apps/python/272/bin)")
import os
import time
import shutil
import subprocess
from copy import copy, deepcopy
from collections import defaultdict

def get_project_dir():
    """
    Returns the root directory of the project
    """
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')) # Root directory of the project code

def create_seed(seed):
    if not seed:
        seed = long(time.time() * 256)
    else:
        seed = int(seed)
    return seed

def combo_params(args, ranging_params, combo=False):
    """ 
    Permutes the params (stored in the args variable) by the other variables provided or if combo is false replicate
    singleton dimensions to match the multiple parameters.
    
    @param args [namespace]: a namespace containing the members referenced in the 'ranging_params' list, usually the args passed by the ArgumentParser
    @param ranging_params [list(str)]: A list of strings that match selected fields in the 'args' namespace
    @param combo [bool]: A flag to decide whether the parameters with multiple values are combod or whether they are checked to be either the same length or singletons
    """
    args = deepcopy(args) # Deep copy so as not to affect the inputted args
    ranged_params = []
    if combo:
        ranged_params.append(getattr(args, ranging_params[0]))
        for par_name in ranging_params[1:]:
            param = getattr(args, par_name)
            curr_len = len(ranged_params[0]) # All parameters should have the same length
            for val in ranged_params:
                val *= len(param)
            new_param = []
            for p in param:
                new_param += [p] * curr_len
            ranged_params.append(new_param)
    else:
        num_param_sets = max([len(getattr(args, p)) for p in ranging_params])
        for par_name in ranging_params:
            param = getattr(args, par_name)
            if len(param) == 1:
                param *= num_param_sets
            elif len(param) != num_param_sets:
                raise Exception('Number of permuations for ''{param}'' ({num}) does not match number of params \
({total_num})'.format(param=par_name, num=len(param), total_num=num_param_sets))
            ranged_params.append(param)
    return ranged_params

def create_work_dir(script_name, output_dir_parent=None, required_dirs=[]):
    """
    Generates unique paths for the work and output directories, creating the work directory in the 
    process.
    
    @param script_name: The name of the script, used to name the directories appropriately
    @param output_dir_parent: The name of the parent directory in which the output directory will be created (defaults \
to $HOME/Output).
    @param required_dirs: The sub-directories that need to be copied into the work directory    
    """
    if not output_dir_parent:
        output_dir_parent = os.path.realpath(os.path.join(os.environ['HOME'], 'output'))
    work_dir_parent = os.path.realpath(os.path.join(os.environ['HOME'], 'work'))
    if not os.path.exists(work_dir_parent):
        raise Exception("Symbolic link to your work directory is missing from your home directory \
(i.e. $HOME/work). A symbolic link should be created that points to an appropriate \
directory in your units sub-directory of '/work' (i.e. ln -s /work/<unit-name>/<user-name> $HOME/work)")
    if not work_dir_parent.startswith('/work'):
        raise Exception("$HOME/work must be a symbolic link to a sub-directory of the high-performance \
filesystem mounted at '/work' (typically /work/<unit-name>/<user-name>).")
    # Automatically generate paths
    time_str = time.strftime('%Y-%m-%d-%A_%H-%M-%S', time.localtime()) # Unique time for distinguishing runs    
    work_dir = os.path.join(work_dir_parent, script_name + "." + time_str + ".1") # Working directory path
    #Ensure that working directory is unique
    created_work_dir = False
    count = 1
    while not created_work_dir:
        try:
            created_work_dir = not os.makedirs(work_dir)
        except OSError as e:
            count += 1
            if count > 1000:
                print "Somethinnessg has gone wrong, can't create directory '%s' after 1000 attempts" % work_dir
                raise e
            work_dir = '.'.join(work_dir.split('.')[:-1] + [str(count)]) # Replace old count at the end of work directory with new count
    output_dir = os.path.join(output_dir_parent, os.path.split(work_dir)[1])
    init_work_dir(work_dir, required_dirs, time_str)
    return work_dir, output_dir

def init_work_dir(work_dir, required_dirs, time_str):
    """
    Copies directories from the project directory to the work directory
    
    @param work_dir: The destination work directory
    @param required_dirs: The required sub-directories to be copied to the work directory
    """
    # Copy snapshot of selected subdirectories to working directory
    for directory in required_dirs:
        shutil.copytree(os.path.join(get_project_dir(), directory), os.path.join(work_dir, directory), symlinks=True)
    # Make output directory for the generated files
    os.mkdir(os.path.join(work_dir, 'output'))
    # Save the git revision in the output folder for reference
    subprocess.call('cd %s; git rev-parse HEAD > %s' % (get_project_dir(), os.path.join(work_dir, 'output', 'git_revision')), shell=True)
    # Write time string to file for future reference
    f = open(os.path.join(work_dir, 'output', 'time_stamp'), 'w')
    f.write(time_str + '\n')
    f.close()



def create_env(work_dir):
    """
    Creates a dictionary containing the appropriate environment variables
    
    @param work_dir: The work directory to set the envinroment variables for
    """
    env = os.environ.copy()
    env['PATH'] = (
        env['PATH'] + os.pathsep +
        os.path.join(env['HOME'], 'git', 'FouTS', 'bin') + os.pathsep +
        os.path.join(env['HOME'], 'git', 'MRtrix', 'bin'))
    env['LD_LIBRARY_PATH'] = (
        '/opt/mpi/gnu/openmpi-1.6.3/lib' + os.pathsep +
        os.path.join(env['HOME'], 'git', 'FouTS', 'lib') + os.pathsep +
        os.path.join(env['HOME'], 'git', 'MRtrix', 'lib'))
    return env


def submit_job(script_name, cmds, np, work_dir, output_dir, que_name='longP', env=None, copy_to_output=[],
               max_memory='4G', virtual_memory='4G', time_limit=None, dry_run=False):
    """
    Create a jobscript in the work directory and then submit it to the HPC que
    
    @param script_name: The name of the script (used to give a meaningful name to the job)
    @param cmds: The commands to run on the cluster
    @param np: The number of processors to request for the job
    @param work_dir: The working directory to run the script from
    @param output_dir: The output directory to copy the results to
    @param env: The required environment variables (defaults to those generated by 'create_env(work_dir)')    
    @param copy_to_output: Directories to copy into the output directory
    @param strip_build_from_copy: Removes all files and directories to be copied that have the name 'build'
    """
    if not env:
        env = create_env(work_dir)
    else:
        env = copy(env)
    if time_limit:
        if type(time_limit) != str or len(time_limit.split(':')) != 3:
            raise Exception("Poorly formatted time limit string '{}' passed to submit job"
                            .format(time_limit))
        time_limit_option= "\n# Set the maximum run time\n#$ -l h_rt {}\n".format(time_limit)
    else:
        time_limit_option=''
    copy_cmd = ''
    for to_copy in copy_to_output:
        origin = work_dir + os.path.sep + to_copy
        destination = output_dir + os.path.sep + to_copy
        base_dir = os.path.dirname(destination[:-1] if destination.endswith('/') else destination)
        copy_cmd += """
mkdir -p {base_dir}
cp -r {origin} {destination}
""".format(base_dir=base_dir, origin=origin, destination=destination)
    #Create jobscript
    jobscript_path = os.path.join(work_dir, script_name + '.job')
    f = open(jobscript_path, 'w')
    f.write("""#!/usr/bin/env sh
    
# Parse the job script using this shell
#$ -S /bin/bash

# Send stdout and stderr to the same file.
#$ -j y

# Standard output and standard error files
#$ -o {work_dir}/output_stream
#$ -e {work_dir}/output_stream

# Name of the queue
#$ -q {que_name}

# use OpenMPI parallel environment with {np} processes
#$ -pe openmpi {np}

# Export the following env variables:
#$ -v HOME
#$ -v PATH
#$ -v LD_LIBRARY_PATH

# Set the memory limits for the script
#$ -l h_vmem={max_mem}
#$ -l virtual_free={virt_mem}

{time_limit}

###################################################
### Copy the model to all machines we are using ###
###################################################

# Set up the correct paths 
export PATH={path}:$PATH
export PYTHONPATH={pythonpath}
export LD_LIBRARY_PATH={ld_library_path}

echo "============== Starting mpirun ===============" 

cd {work_dir}
{cmds}

echo "============== Mpirun has ended =============="

echo "Copying files to output directory '{output_dir}'"
cp -r {work_dir}/output {output_dir}
cp {jobscript_path} {output_dir}/job
{copy_cmd}
echo "============== Done ===============" 
cp {work_dir}/output_stream {output_dir}/output
""".format(work_dir=work_dir, path=env['PATH'], pythonpath=env['PYTHONPATH'],
      ld_library_path=env['LD_LIBRARY_PATH'], np=np,
      que_name=que_name, cmds=cmds, output_dir=output_dir, copy_cmd=copy_cmd,
      jobscript_path=jobscript_path, max_mem=max_memory, virt_mem=virtual_memory, 
      time_limit=time_limit_option))
    f.close()
    # Submit job
    print "Submitting job '%s' to que" % jobscript_path
    if not dry_run:
        subprocess.check_call('qsub %s' % jobscript_path, shell=True)
    print "Your job '%s' has been submitted" % jobscript_path
    print "The output stream can be viewed by:"
    print "less " + os.path.join(work_dir, 'output_stream')
    print "Once completed the output files (including the output stream and job script) of this job will be copied to:"
    print output_dir
