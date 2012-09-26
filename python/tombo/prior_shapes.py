#!/usr/bin/env python
"""
 This script prepares the environment for the simulate/fabios_network.py script on the Sun Grid 
 Engine batch script to run in by making a snapshot of the code base at the start of the run, 
 generating a jobscript and sending it to the job que

 Author: Tom Close (tclose@oist.jp)
 Created: 6/8/2012
"""

#Name of the script for the output directory and submitted mpi job
SCRIPT_NAME = 'prior_shapes'

import tombo
import argparse
import os.path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output_dir', default=None, type=str, help='The parent directory in which the output directory will be created (defaults to $HOME/Output)')
parser.add_argument('--num_runs', default=1, type=int, help='The number of runs to submit to the que')
parser.add_argument('--np', type=int, default=1, help='The the number of processes to use for the simulation (default: %(default)s)')
parser.add_argument('--que_name', type=str, default='short', help='The the que to submit the job to(default: %(default)s)')
args = parser.parse_args()

for i in xrange(args.num_runs):
    required_dirs = [os.path.join('fibre', 'tract'), os.path.join('image', 'noisy')]

    # Create work directory and get path for output directory
    work_dir, output_dir = tombo.create_work_dir(SCRIPT_NAME, args.output_dir, required_dirs=required_dirs)

    # Set up command to run the script
    cmd_line = "time metropolis {work_dir}/data/image/noisy/prior_background.mif \
{work_dir}/data/fibre/tract/single/x.tct".format(work_dir=work_dir)

    copy_to_output = required_dirs
    # Submit job to que
    tombo.submit_job(SCRIPT_NAME, cmd_line, args.np, work_dir, output_dir, copy_to_output=copy_to_output,
                                                                                que_name=args.que_name, dry_run=True)
