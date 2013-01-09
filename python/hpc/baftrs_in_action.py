#!/usr/bin/env python
"""
 This script submits scripts to generate random prior shapes from noisy images and submit them to the job que on a sun
 grid engine

 Author: Tom Close (tclose@oist.jp)
 Created: 6/8/2012
"""
# Required imports
import hpc
import argparse
import os.path
import time
#Name of the script for the output directory and submitted mpi job
SCRIPT_NAME = 'baftrs_in_action'
DATASETS = [os.path.join('donald', 'fornix.mif')]
#DATASETS = [os.path.join('donald', 'fornix.mif'), os.path.join('heath', 'fornix.mif'), os.path.join('lisa', 'fornix.mif')]
# Required dirs for the script to run
REQUIRED_DIRS = ['params/image/reference', 'params/diffusion', 'params/fibre/tract/masks/mcmc/metropolis']
# Arguments that can be given to the script
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--step_scale', default=0.015, type=float, help='The scale of the steps used for the metropolis sampling (default: %(default)s)')
parser.add_argument('--num_iterations', default=100000, type=int, help='The number of interations in the metropolis sampling (default: %(default)s)')
parser.add_argument('--sample_period', default=1000, type=int, help='The sample period of the metropolis sampling (default: %(default)s)')
parser.add_argument('--num_norms', type=int, default=3, help='The number of times the num_width_sections and num_length_sections are normalised %(default)s') 
parser.add_argument('--num_width_sections', default=4, help='The number of samples to use across a Fourier tracts cross-section (default: %(default)s)')
parser.add_argument('--interp_type', default='sinc', type=str, help='The type of interpolation used in the reference image (default: %(default)s)')
parser.add_argument('--interp_extent', default=1, type=int, help='The interpolation extent used in the reference image (default: %(default)s)')
parser.add_argument('--assumed_interp_extent', default=1, type=int, help='The interpolation type used in the likelihood images (default: %(default)s)')
parser.add_argument('--degree', default=5, type=int, help='The degree of the fibre used to sample from (default: %(default)s)')
parser.add_argument('--prior_freq', default=[15.0], type=float, nargs='+', help='The scaling of the frequency prior (default: %(default)s)')
parser.add_argument('--prior_aux_freq', default=[60.0], type=float, nargs='+', help='The scaling of the frequency prior (default: %(default)s)')
parser.add_argument('--prior_density_high', default=[1], type=float, nargs='+', help='The scaling of the density prior (default: %(default)s)')
parser.add_argument('--prior_density_low', default=[1], type=float, nargs='+', help='The scaling of the density prior (default: %(default)s)')
parser.add_argument('--prior_hook', default=[100000.0], type=float, nargs='+', help='The scaling of the density prior (default: %(default)s)')
parser.add_argument('--prior_thin', default=[0.0], type=float, nargs='+', help='The scaling of the density prior (default: %(default)s)')
parser.add_argument('--width_epsilon', default=0.15, type=float, help='The amount of width epsilon to use (default: %(default)s)')
parser.add_argument('--length_epsilon', default=0.15, type=float, help='The amount of length epsilon to use (default: %(default)s)')
parser.add_argument('--img_snr', default=20.0, type=float, help='The snr to used in the noisy image')
parser.add_argument('--like_snr', default=[20.0], type=float, nargs='+',
                    help='The assumed snr to used in the likelihood function in the metropolis sampling')
parser.add_argument('--init_offset', default=(115, 117, 75), nargs=3, type=float,
                    help='The offset of the ROI in which the inital tracts are generated')
parser.add_argument('--init_extent', default=(2, 2, 2), nargs=3, type=int,
                    help='The extent (in number of voxels) of the ROI in which the inital tracts are generated')
parser.add_argument('--seed', type=int, help="The random seed for the whole run")
parser.add_argument('--output_dir', default=None, type=str, help='The parent directory in which the output directory \
will be created (defaults to $HOME/Output)')
parser.add_argument('--num_runs', default=1, type=int, help='The number of runs to submit to the que %(default)s')
parser.add_argument('--dry_run', action='store_true', help='Only perform a dry run (create jobscript then quit)')
parser.add_argument('--np', type=int, default=1, help='The the number of processes to use for the simulation \
(default: %(default)s)')
parser.add_argument('--que_name', type=str, default='long', help='The the que to submit the job to \
(default: %(default)s)')
parser.add_argument('--combo', action='store_true', help="Instead of treating each ranging parameter sequence as the 1..N values " \
                                                         "for that parameter, all combinations of the provided parameters are tested.")
parser.add_argument('--estimate_response', action='store_true', help="Uses an estimated diffusion response function instead of the default tensor one")
args = parser.parse_args()
# For the following parameters to this script, ensure that number of parameter values match, or if they are a singleton 
# list it is assumed to be constant and that value that value is replicated to match the number of other of other 
# parameters in the set. Otherwise if the '--combo' option is provided then loop through all combinations of the 
# provided parameters. 
ranging_param_names = ['prior_freq', 'prior_aux_freq', 'prior_density_low', 'prior_density_high',
                            'prior_hook', 'prior_thin', 'like_snr']
ranging_params = hpc.combo_params(args, ranging_param_names, args.combo)
# Generate a random seed to seed the random number generators of the cmds
if not args.seed:
    seed = int(time.time() * 100)
    print "Using seed {}".format(seed)
else:
    seed = args.seed
for i in xrange(args.num_runs):
    for prior_freq, prior_aux_freq, prior_density_low, prior_density_high, prior_hook, prior_thin, like_snr in zip(*ranging_params):

        for dataset in DATASETS:
            # Create work directory and get path for output directory
            work_dir, output_dir = hpc.create_work_dir(SCRIPT_NAME, args.output_dir, required_dirs=REQUIRED_DIRS)
            with open(os.path.join(work_dir, 'summary.txt'), 'w') as f:
                f.write(dataset + '\n')
                for par_name in ranging_param_names:
                    f.write('{par}: {val}\n'.format(par=par_name, val=eval(par_name)))
            # Create a file in the output directory with just the dataset printed in it (usefulf for quickly 
            # determining what the dataset is
            with open(os.path.join(work_dir, 'output', 'dataset_name.txt'), 'w') as dataset_file:
                dataset_file.write(dataset + '\n')
            # Strip dataset of symbols for tract number and img dimension
            # Get the dataset path
            dataset_path = os.path.join(work_dir, 'params', 'image', 'reference', dataset)
            dataset_dir = os.path.dirname(dataset)
            if args.estimate_response:
                response_path = os.path.join(work_dir, 'params', 'image', 'reference', dataset_dir, 'response.txt')
                response_str = "-diff_response {}".format(response_path)
            else:
                response_str = ""
            response_b0_path = os.path.join(work_dir, 'params', 'image', 'reference', dataset_dir, 'response.b0.txt')
            cmd_line = """             
# Create initial_fibres of appropriate degree
init_fibres {work_dir}/output/init.tct -degree {args.degree} -num_fibres 1 -width_epsilon {args.width_epsilon} \
-length_epsilon {args.length_epsilon} -edge_buffer 0.0 -seed {init_seed} -base_intensity 1.0 \
-img_dims "{init_extent[0]} {init_extent[1]} {init_extent[2]}" \
-img_offset "{init_offset[0]} {init_offset[1]} {init_offset[2]}" 

# Run metropolis
metropolis {dataset_path} {work_dir}/output/init.tct {work_dir}/output/samples.tst -like_snr {like_snr} \
-exp_interp_extent {args.assumed_interp_extent} -walk_step_scale {args.step_scale} -num_iter {args.num_iterations} \
-sample_period {args.sample_period} -diff_encodings_location {work_dir}/params/diffusion/encoding_60.b \
-seed {seed} -prior_freq {prior_freq} {prior_aux_freq} -prior_density {prior_density_high} \
{prior_density_low} 100 -prior_hook {prior_hook} 100 15 -prior_thin {prior_thin} 2 \
-exp_num_width_sections {args.num_width_sections} -exp_type {args.interp_type} {response_str} \
-exp_b0 `cat {b0_path}` -diff_warn -walk_step_location \
{work_dir}/params/fibre/tract/masks/mcmc/metropolis/default{args.degree}.tct

    """.format(work_dir=work_dir, dataset_path=dataset_path, args=args,
               seed=seed, init_seed=seed + 1, prior_freq=prior_freq, prior_aux_freq=prior_aux_freq,
               prior_density_low=prior_density_low, prior_density_high=prior_density_high, prior_hook=prior_hook,
               prior_thin=prior_thin, like_snr=like_snr, response_str=response_str, b0_path=response_b0_path,
               init_extent=args.init_extent, init_offset=args.init_offset)
            # Submit job to que
            hpc.submit_job(SCRIPT_NAME, cmd_line, args.np, work_dir, output_dir, que_name=args.que_name,
                                                                dry_run=args.dry_run, copy_to_output=['summary.txt'])
        seed += 2
