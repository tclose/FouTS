#!/usr/bin/env python
"""
 This script submits scripts to generate random prior shapes from noisy images and submit them to the job que on a sun
 grid engine

 Author: Tom Close (tclose@oist.jp)
 Created: 6/8/2012
"""

#Name of the script for the output directory and submitted mpi job
SCRIPT_NAME = 'prior_shapes'

import hpc
import argparse
import os.path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--img_dims', default=10, type=int, help='The size of the noisy image to fit against')
parser.add_argument('--prior_acs', nargs=2, default=(1000, 0.2), type=float, help='The standard deviation and mean of the ACS (Apparent Connection Strength) prior.')
parser.add_argument('--img_snr', default=5, type=float, help='The snr to used in the noisy image')
parser.add_argument('--like_snr', default=20, type=float, help='The assumed snr to used in the likelihood function in the metropolis sampling')
parser.add_argument('--output_dir', default=None, type=str, help='The parent directory in which the output directory will be created (defaults to $HOME/Output)')
parser.add_argument('--num_runs', default=1, type=int, help='The number of runs to submit to the que')
parser.add_argument('--dry_run', action='store_true', help='Only perform a dry run (create jobscript then quit)')
parser.add_argument('--np', type=int, default=1, help='The the number of processes to use for the simulation (default: %(default)s)')
parser.add_argument('--que_name', type=str, default='short', help='The the que to submit the job to(default: %(default)s)')
args = parser.parse_args()

for i in xrange(args.num_runs):
    # Create work directory and get path for output directory
    work_dir, output_dir = hpc.create_work_dir(SCRIPT_NAME, args.output_dir)
    # Set up command to run the script
    cmd_line = """
generate_image dummy.tct {work_dir}/output/noise.mif --empty --img_dims "{dim} {dim} {dim}" --noise_ref_signal 1 --noise_snr {img_snr} --diff_encodings {encoding_loc}
init_fibres {work_dir}/output/init.tct --acs {acs_mean} --base_intensity 1.0 --num_fibres 1 --img_dims "{dim} {dim} {dim}"    
time metropolis {work_dir}/output/noise.mif {work_dir}/output/init.tct {work_dir}/output/samples.tst --like_snr {like_snr} --diff_encodings {encoding_loc} --prior_acs {acs_stddev} {acs_mean}
""".format(work_dir=work_dir, dim=args.img_dims, acs_mean=args.prior_acs[1], img_snr=args.img_snr, like_snr=args.like_snr, encoding_loc=os.path.join(hpc.get_project_dir(), 'params', 'diffusion', 'encoding_1.b'), acs_stddev=args.prior_acs[0])
    # Submit job to que
    hpc.submit_job(SCRIPT_NAME, cmd_line, args.np, work_dir, output_dir, que_name=args.que_name)
