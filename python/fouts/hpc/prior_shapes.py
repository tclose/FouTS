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
parser.add_argument('--prior_acs', nargs=2, default=(1000, 0.1), type=float, help='The standard deviation and mean of the ACS (Apparent Connection Strength) prior.')
parser.add_argument('--img_snr', default=5, type=float, help='The snr to used in the noisy image')
parser.add_argument('--like_snr', default=20, type=float, help='The assumed snr to used in the likelihood function in the metropolis sampling')
parser.add_argument('--degree', default=7, type=int, help='The degree of the generated Fourier tracts')
parser.add_argument('--output_dir', default=None, type=str, help='The parent directory in which the output directory will be created (defaults to $HOME/Output)')
parser.add_argument('--num_runs', default=1, type=int, help='The number of runs to submit to the que')
parser.add_argument('--dry_run', action='store_true', help='Only perform a dry run (create jobscript then quit)')
parser.add_argument('--np', type=int, default=1, help='The the number of processes to use for the simulation (default: %(default)s)')
parser.add_argument('--que_name', type=str, default='short', help='The the que to submit the job to(default: %(default)s)')
parser.add_argument('--prior_freq', default=[15.0], type=float, nargs='+', help='The scaling of the frequency prior')
parser.add_argument('--prior_aux_freq', default=[60.0], type=float, nargs='+', help='The scaling of the frequency prior')
parser.add_argument('--prior_density_high', default=[1], type=float, nargs='+', help='The scaling of the density prior')
parser.add_argument('--prior_density_low', default=[0.01], type=float, nargs='+', help='The scaling of the density prior')
parser.add_argument('--prior_hook', default=[100000.0], type=float, nargs='+', help='The scaling of the density prior')
parser.add_argument('--combo', action='store_true', help='Instead of treating each ranging parameter sequence as the 1..N values for that parameter, all combinations of the provided parameters are tested.')
args = parser.parse_args()
# For the following parameters to this script, ensure that number of parameter values match, or if they are a singleton 
# list it is assumed to be constant and that value that value is replicated to match the number of other of other 
# parameters in the set. Otherwise if the '--combo' option is provided then loop through all combinations of the 
# provided parameters. 
ranging_param_names = ['prior_freq', 'prior_aux_freq', 'prior_density_low', 'prior_density_high', 'prior_hook']
ranging_params = hpc.combo_params(args, ranging_param_names, args.combo)
for i in xrange(args.num_runs):
    for prior_freq, prior_aux_freq, prior_density_low, prior_density_high, prior_hook in zip(*ranging_params):
        required_dirs = [os.path.join('params', 'fibre', 'tract', 'masks', 'mcmc', 'metropolis'), os.path.join('params', 'diffusion')]
        # Create work directory and get path for output directory
        work_dir, output_dir = hpc.create_work_dir(SCRIPT_NAME, args.output_dir, required_dirs=required_dirs)
        with open(os.path.join(work_dir, 'output', 'summary.txt'), 'w') as f:
            for par_name in ranging_param_names:
                f.write('{par}: {val}\n'.format(par=par_name, val=eval(par_name)))
        # Set up command to run the script
        cmd_line = """
# Select one of the 60 orientations and write it to a new file        
head -n $(( $RANDOM * 60 / 32767 )) {work_dir}/params/diffusion/encoding_60.b | tail -n 1 > {work_dir}/output/encoding.b \

generate_image dummy.tct {work_dir}/output/noise.mif --empty --img_dims "{dim} {dim} {dim}" --noise_ref_signal 1 \
--noise_snr {img_snr} --diff_encodings {work_dir}/output/encoding.b

init_fibres {work_dir}/output/init.tct --degree {degree} --acs {acs_mean} --base_intensity 1.0 --num_fibres 1 \
--img_dims "{dim} {dim} {dim}" --width_epsilon 0.01 --length_epsilon 0.01  

time metropolis {work_dir}/output/noise.mif {work_dir}/output/init.tct {work_dir}/output/samples.tst \
--like_snr {like_snr} --diff_encodings {work_dir}/output/encoding.b --prior_acs {acs_stddev} {acs_mean} \
--walk_step_location {work_dir}/params/fibre/tract/masks/mcmc/metropolis/default{degree}.tct \
--prior_freq {prior_freq} {prior_aux_freq} -prior_density {prior_density_high} {prior_density_low} 100 \
-prior_hook {prior_hook} 100 15 
    """.format(work_dir=work_dir, dim=args.img_dims, degree=args.degree, acs_mean=args.prior_acs[1], img_snr=args.img_snr,
            like_snr=args.like_snr, acs_stddev=args.prior_acs[0], prior_freq=prior_freq, prior_aux_freq=prior_aux_freq,
             prior_density_low=prior_density_low, prior_density_high=prior_density_high, prior_hook=prior_hook)
        # Submit job to que
        hpc.submit_job(SCRIPT_NAME, cmd_line, args.np, work_dir, output_dir, que_name=args.que_name,
                                                                                        copy_to_output=required_dirs)
