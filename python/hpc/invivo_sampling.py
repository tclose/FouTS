#!/usr/bin/env python
"""
 This script submits scripts to generate random prior shapes from noisy images 
 and submit them to the job que on a sun grid engine

 Author: Tom Close (tclose@oist.jp)
 Created: 6/8/2012
"""
# Required imports
import hpc
import argparse
import os.path
import time

#Name of the script for the output directory and submitted mpi job
SCRIPT_NAME = 'invivo_sampling'
# Required dirs for the script to run
REQUIRED_DIRS = ['params/image/reference', 'params/diffusion',
                 'params/fibre/tract/masks/mcmc/metropolis']
# Arguments that can be given to the script
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--step_scale', default=0.015, type=float,
                    help="The scale of the steps used for the metropolis "
                         "sampling (default: %(default)s)")
parser.add_argument('--num_iterations', default=75000, type=int,
                    help="The number of interations in the metropolis sampling "
                         "(default: %(default)s)")
parser.add_argument('--sample_period', default=250, type=int,
                    help="The sample period of the metropolis sampling "
                          "(default: %(default)s)")
parser.add_argument('--num_norms', type=int, default=1,
                    help="The number of times the num_width_sections and "
                         "num_length_sections are normalised "
                         "(default: %(default)s)")
parser.add_argument('--samples_per_acs', type=float, default=2,
                    help="The number of strands to use per ACS when "
                         "renormalising the tractlets at the start of each "
                         "major iteration (default: %(default)s)")
parser.add_argument('--samples_per_length', type=float, default=5,
                    help="The number of length samples to use per ACS when "
                         "renormalising the tractlets at the start of each "
                         "major iteration (default: %(default)s)")
parser.add_argument('--min_length_samples', type=float, default=25,
                    help="The number of strands to use per ACS when "
                         "renormalising the tractlets at the start of each "
                         "major iteration (default: %(default)s)")
parser.add_argument('--min_width_samples', type=float, default=4,
                    help="The number of length samples to use per ACS when "
                         "renormalising the tractlets at the start of each "
                         "major iteration (default: %(default)s)")
parser.add_argument('--interp_type', default='sinc', type=str,
                    help="The type of interpolation used in the reference "
                         "image (default: %(default)s)")
parser.add_argument('--interp_extent', default=1, type=int,
                    help="The interpolation extent used in the reference image "
                         "(default: %(default)s)")
parser.add_argument('--assumed_interp_extent', default=1, type=int,
                    help="The interpolation type used in the likelihood images "
                         "(default: %(default)s)")
parser.add_argument('--degree', default=5, type=int,
                    help="The degree of the fibre used to sample from "
                         "(default: %(default)s)")
parser.add_argument('--prior_freq', default=[15.0], type=float, nargs='+',
                    help="The scaling of the frequency prior (default: "
                         "%(default)s)")
parser.add_argument('--prior_aux_freq', default=[60.0], type=float, nargs='+',
                    help="The scaling of the frequency prior (default: "
                         "%(default)s)")
parser.add_argument('--prior_density_high', default=[1], type=float, nargs='+',
                    help="The scaling of the density prior (default: "
                         "%(default)s)")
parser.add_argument('--prior_density_low', default=[1], type=float, nargs='+',
                    help="The scaling of the density prior (default: "
                         "%(default)s)")
parser.add_argument('--prior_hook', default=[100000.0], type=float, nargs='+',
                    help="The scaling of the density prior (default: "
                         "%(default)s)")
parser.add_argument('--prior_thin', default=[0.0], type=float, nargs='+',
                    help="The scaling of the density prior (default: "
                         "%(default)s)")
parser.add_argument('--width_epsilon', default=0.15, type=float,
                    help="The amount of width epsilon to use (default: "
                         "%(default)s)")
parser.add_argument('--length_epsilon', default=0.15, type=float,
                    help="The amount of length epsilon to use (default: "
                         "%(default)s)")
parser.add_argument('--img_snr', default=20.0, type=float,
                    help="The snr to used in the noisy image")
parser.add_argument('--like_snr', default=[20.0], type=float, nargs='+',
                    help="The assumed snr to used in the likelihood function "
                         "in the metropolis sampling")
parser.add_argument('--init_extent', default=(1, 1, 1), nargs=3, type=int,
                    help="The extent (in number of voxels) of the ROI in which "
                         "the inital tracts are generated")
parser.add_argument('--random_seed', type=int, help="The random seed for the whole "
                                             "run")
parser.add_argument('--output_dir', default=None, type=str,
                    help="The parent directory in which the output directory "
                         "will be created (defaults to $HOME/Output)")
parser.add_argument('--num_runs', default=1, type=int,
                    help="The number of runs to submit to the que")
parser.add_argument('--dry_run', action='store_true',
                    help="Only perform a dry run (create jobscript then quit)")
parser.add_argument('--np', type=int, default=1,
                    help="The the number of processes to use for the simulation"
                         "(default: %(default)s)")
parser.add_argument('--que_name', type=str, default='long',
                    help="The the que to submit the job to  (default: "
                         "%(default)s)")
parser.add_argument('--combo', action='store_true',
                    help="Instead of treating each ranging parameter sequence "
                         " as the 1..N values for that parameter, all "
                         "combinations of the provided parameters are tested.")
#parser.add_argument('--estimate_response', action='store_true',
#                    help="Uses an estimated diffusion response function instead"
#                         " of the default tensor one")
parser.add_argument('--dataset', type=str,
                    default=os.path.join('donald', 'images', 'corpus_callosum.150.mif'),
                    help="The dataset to use (default: %(default)s).")
parser.add_argument('--seed_positions', type=str,
                    default=os.path.join('donald', 'corpus_callosum.150.seed.txt'),
                    help="The seed locations of the tracts.")
parser.add_argument('--seed_pos_stddev', type=float, default=0.5,
                    help="The standard deviation about the seed positions")
parser.add_argument('--intensity_estimate', type=str,
                    default=os.path.join('donald', 'intensity.150.txt'),
                    help="The estimated base intensity of the image.")
args = parser.parse_args()
# For the following parameters to this script, ensure that number of parameter 
# values match, or if they are a singleton list it is assumed to be constant 
# and that value that value is replicated to match the number of other of other 
# parameters in the set. Otherwise if the '--combo' option is provided then loop
# through all combinations of the provided parameters. 

ranging_param_names = ['prior_freq', 'prior_aux_freq', 'prior_density_low',
                       'prior_density_high', 'prior_hook', 'prior_thin',
                       'like_snr']
ranging_params = hpc.combo_params(args, ranging_param_names, args.combo)
# Generate a random seed to seed the random number generators of the cmds
if not args.random_seed:
    random_seed = int(time.time() * 100)
    print "Using random_seed {}".format(random_seed)
else:
    random_seed = args.random_seed
for i in xrange(args.num_runs):
    for (prior_freq, prior_aux_freq, prior_density_low, prior_density_high,
                prior_hook, prior_thin, like_snr) in zip(*ranging_params):

        # Create work directory and get path for output directory
        work_dir, output_dir = hpc.create_work_dir(SCRIPT_NAME,
                                                   args.output_dir,
                                                   required_dirs=\
                                                   REQUIRED_DIRS)
        with open(os.path.join(work_dir, 'summary.txt'), 'w') as f:
            f.write(args.dataset + '\n')
            for par_name in ranging_param_names:
                f.write('{par}: {val}\n'.format(par=par_name,
                                                val=eval(par_name)))
        # Create a file in the output directory with just the dataset 
        # printed in it (useful for quickly determining what the dataset is)
        with open(os.path.join(work_dir, 'output', 'dataset_name.txt'),
                  'w') \
                as dataset_file:
            dataset_file.write(args.dataset + '\n')
        # Strip dataset of symbols for tract number and img dimension
        # Get the dataset path
        dataset_path = os.path.join(work_dir, 'params', 'image',
                                    'reference', args.dataset)
        dataset_dir = os.path.dirname(args.dataset)
        response_b0_path = os.path.join(work_dir, 'params', 'image',
                                        'reference', dataset_dir,
                                        'response.b0.txt')
        cmd_line = \
""" 
init_fibres  {work_dir}/output/init.tct -degree {args.degree} \
-width_epsilon {args.width_epsilon} \
-length_epsilon {args.length_epsilon} -random_seed {init_seed} \
-seed_positions {work_dir}/params/image/reference/{args.seed_positions} \
-centre_stddev {args.seed_pos_stddev} -base_intensity 1.0

# Normalise the density of the initial tracts
normalise_density {work_dir}/output/init.tct

# Calculate appropriate number of length and width samples                
calculate_num_samples --samples_per_acs {args.samples_per_acs} \
--samples_per_length {args.samples_per_length} --strategy max \
--min_length_sections {args.min_length_samples} \
--min_width_sections {args.min_width_samples} \
{work_dir}/output/init.tct {work_dir}/num_samples 

# Set required properties for the initial tracts
set_properties {work_dir}/output/init.tct \
-set length_epsilon {args.length_epsilon} \
-set width_epsilon {args.width_epsilon} -set base_intensity 1.0
                
# Run metropolis
metropolis {dataset_path} {work_dir}/output/init.tct \
{work_dir}/output/samples.tst -like_snr {like_snr} \
-exp_interp_extent {args.assumed_interp_extent} \
-walk_step_scale {args.step_scale} -num_iter {args.num_iterations} \
-sample_period {args.sample_period} -seed {metrop_seed} \
-diff_encodings_location {work_dir}/params/diffusion/encoding_60.b \
-prior_freq {prior_freq} {prior_aux_freq} \
-prior_density {prior_density_high} {prior_density_low} 100 \
-prior_hook {prior_hook} 100 15 -prior_thin {prior_thin} 2 \
-exp_num_length_sections `cat {work_dir}/num_samples.length.txt` \
-exp_num_width_sections `cat {work_dir}/num_samples.width.txt` \
-exp_type {args.interp_type} \
-exp_base_intensity `cat {work_dir}/params/image/reference/{args.intensity_estimate}` \
-diff_warn -walk_step_location \
{work_dir}/params/fibre/tract/masks/mcmc/metropolis/default{args.degree}.tct

# Get last sample
select_fibres {work_dir}/output/samples.tst \
{work_dir}/output/last.tct --include {last_sample}

""".format(work_dir=work_dir, dataset_path=dataset_path, args=args,
                     init_seed=random_seed, metrop_seed=random_seed + 1, prior_freq=prior_freq,
                     prior_aux_freq=prior_aux_freq,
                     prior_density_low=prior_density_low,
                     prior_density_high=prior_density_high,
                     prior_hook=prior_hook,
                     prior_thin=prior_thin, like_snr=like_snr,
                     b0_path=response_b0_path,
                     last_sample=(args.num_iterations // args.sample_period) - 1)
            # Submit job to que
        hpc.submit_job(SCRIPT_NAME, cmd_line, args.np, work_dir, output_dir,
                       que_name=args.que_name, dry_run=args.dry_run,
                       copy_to_output=['summary.txt'])

