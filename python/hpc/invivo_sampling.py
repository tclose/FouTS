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
SCRIPT_NAME = 'invio_sampling'
# Required dirs for the script to run
REQUIRED_DIRS = ['params/image/reference', 'params/diffusion',
                 'params/fibre/tract/masks/mcmc/metropolis']
# Arguments that can be given to the script
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--step_scale', default=0.015, type=float,
                    help="The scale of the steps used for the metropolis "
                         "sampling (default: %(default)s)")
parser.add_argument('--num_iterations', default=100000, type=int,
                    help="The number of interations in the metropolis sampling "
                         "(default: %(default)s)")
parser.add_argument('--sample_period', default=1000, type=int,
                    help="The sample period of the metropolis sampling "
                          "(default: %(default)s)")
parser.add_argument('--num_norms', type=int, default=3,
                    help="The number of times the num_width_sections and "
                         "num_length_sections are normalised "
                         "(default: %(default)s)")
parser.add_argument('--samples_per_acs', type=float, default=100,
                    help="The number of strands to use per ACS when "
                         "renormalising the tractlets at the start of each "
                         "major iteration (default: %(default)s)")
parser.add_argument('--samples_per_length', type=float, default=100,
                    help="The number of length samples to use per ACS when "
                         "renormalising the tractlets at the start of each "
                         "major iteration (default: %(default)s)")
parser.add_argument('--min_length_samples', type=float, default=15,
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
parser.add_argument('--seed', type=int, help="The random seed for the whole "
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
parser.add_argument('--estimate_response', action='store_true',
                    help="Uses an estimated diffusion response function instead"
                         " of the default tensor one")
parser.add_argument('--dataset', action='append', type=str,
                    default=[os.path.join('donald', 'fornix.mif')],
                    help="The datasets to use (default: %(default)s).")
parser.add_argument('--init_locations', action='append', type=str,
                    default=[((111, 119, 73), (111, 80, 75))],
                    help="The initial locations of the tracts (default: "
                         "%(default)s).")
args = parser.parse_args()
# For the following parameters to this script, ensure that number of parameter 
# values match, or if they are a singleton list it is assumed to be constant 
# and that value that value is replicated to match the number of other of other 
# parameters in the set. Otherwise if the '--combo' option is provided then loop
# through all combinations of the provided parameters. 
if len(args.dataset) != len(args.init_locations):
    raise Exception("The same amount of 'dataset' and 'init_locations' must be"
                    "supplied ({} and {} supplied respectively)"
                    .format(args.dataset, args.init_locations))
ranging_param_names = ['prior_freq', 'prior_aux_freq', 'prior_density_low',
                       'prior_density_high', 'prior_hook', 'prior_thin',
                       'like_snr']
ranging_params = hpc.combo_params(args, ranging_param_names, args.combo)
# Generate a random seed to seed the random number generators of the cmds
if not args.seed:
    seed = int(time.time() * 100)
    print "Using seed {}".format(seed)
else:
    seed = args.seed
for i in xrange(args.num_runs):
    for (prior_freq, prior_aux_freq, prior_density_low, prior_density_high,
                prior_hook, prior_thin, like_snr) in zip(*ranging_params):

        for dataset, init_locations in zip(args.dataset, args.init_locations):
            # Create work directory and get path for output directory
            work_dir, output_dir = hpc.create_work_dir(SCRIPT_NAME,
                                                       args.output_dir,
                                                       required_dirs=\
                                                       REQUIRED_DIRS)
            with open(os.path.join(work_dir, 'summary.txt'), 'w') as f:
                f.write(dataset + '\n')
                for par_name in ranging_param_names:
                    f.write('{par}: {val}\n'.format(par=par_name,
                                                    val=eval(par_name)))
            # Create a file in the output directory with just the dataset 
            # printed in it (useful for quickly determining what the dataset is)
            with open(os.path.join(work_dir, 'output', 'dataset_name.txt'),
                      'w') \
                    as dataset_file:
                dataset_file.write(dataset + '\n')
            # Strip dataset of symbols for tract number and img dimension
            # Get the dataset path
            dataset_path = os.path.join(work_dir, 'params', 'image',
                                        'reference', dataset)
            dataset_dir = os.path.dirname(dataset)
            if args.estimate_response:
                response_path = os.path.join(work_dir, 'params', 'image',
                                             'reference', dataset_dir,
                                             'response.txt')
                response_str = "-diff_response {}".format(response_path)
            else:
                response_str = ""
            response_b0_path = os.path.join(work_dir, 'params', 'image',
                                            'reference', dataset_dir,
                                            'response.b0.txt')
            cmd_line = \
"""
# Initialise empty set to hold the initial locations of the tractlets
new_fibres -set_size 0 -acs 1.0 {work_dir}/output/init_0.tct
            
""".format(work_dir=work_dir)
            if not len(init_locations):
                raise Exception("Initial locations list was empty")
            for i, init_location in enumerate(init_locations):
                cmd_line += \
"""             
# Create initial_fibres of appropriate degree
init_fibres {work_dir}/init_location_{i}.tct -degree {args.degree} \
-num_fibres 1 -width_epsilon {args.width_epsilon} \
-length_epsilon {args.length_epsilon} -edge_buffer 0.0 -seed {seed} \
-base_intensity 1.0 \
-img_dims "{args.init_extent[0]} {args.init_extent[1]} {args.init_extent[2]}" \
-img_offset "{init_location[0]} {init_location[1]} {init_location[2]}"

# Add new tract to initial set
combine_fibres {work_dir}/output/init_0.tct {work_dir}/init_location_{i}.tct \
{work_dir}/tmp.tct
mv {work_dir}/tmp.tct {work_dir}/output/init_0.tct
mv {work_dir}/tmp.tctx {work_dir}/output/init_0.tctx

""".format(work_dir=work_dir, i=i, init_location=init_location, args=args,
           seed=seed)
                seed += 1
            # Normalise the density of the initial 
            # tracts
            cmd_line += \
"""
# Normalise the density of the initial tracts
normalise_density {work_dir}/output/init_0.tct

""".format(work_dir=work_dir, args=args)

            for norm_count in range(args.num_norms):
                cmd_line += \
""" 
# Calculate appropriate number of length and width samples                
calculate_num_samples --samples_per_acs {args.samples_per_acs} \
--samples_per_length {args.samples_per_length} --strategy max \
--min_length_samples {args.min_length_samples} \
--min_width_samples {args.min_width_samples} \
{work_dir}/output/init_{norm_count}.tct {work_dir}/num_samples 

# Set required properties for the initial tracts
set_properties {work_dir}/output/init_{norm_count}.tct \
-set length_epsilon {args.length_epsilon} \
-set width_epsilon {args.width_epsilon} -set base_intensity 1.0
                
# Run metropolis
metropolis {dataset_path} {work_dir}/output/init_{norm_count}.tct \
{work_dir}/output/samples_{norm_count}.tst -like_snr {like_snr} \
-exp_interp_extent {args.assumed_interp_extent} \
-walk_step_scale {args.step_scale} -num_iter {args.num_iterations} \
-sample_period {args.sample_period} -seed {seed} \
-diff_encodings_location {work_dir}/params/diffusion/encoding_60.b \
-prior_freq {prior_freq} {prior_aux_freq} \
-prior_density {prior_density_high} {prior_density_low} 100 \
-prior_hook {prior_hook} 100 15 -prior_thin {prior_thin} 2 \
-exp_num_length_sections `cat {work_dir}/num_samples.length.txt` \
-exp_num_width_sections `cat {work_dir}/num_samples.width.txt` \
-exp_type {args.interp_type} {response_str} -exp_b0 `cat {b0_path}` \
-diff_warn -walk_step_location \
{work_dir}/params/fibre/tract/masks/mcmc/metropolis/default{args.degree}.tct

# Get last sample
select_fibres {work_dir}/output/samples_{norm_count}.tst \
{work_dir}/output/init_{next_norm_count}.tct --include {last_sample}

""".format(work_dir=work_dir, dataset_path=dataset_path, args=args,
                seed=seed, init_seed=seed + 1, prior_freq=prior_freq,
                prior_aux_freq=prior_aux_freq,
                prior_density_low=prior_density_low,
                prior_density_high=prior_density_high,
                prior_hook=prior_hook,
                prior_thin=prior_thin, like_snr=like_snr,
                response_str=response_str, b0_path=response_b0_path,
                norm_count=norm_count, next_norm_count=norm_count + 1,
                last_sample=(args.num_iterations // args.sample_period) - 1)
                # Increment the seed used for the random initialisation
                seed += 1
            cmd_line += \
"""
# Rename the final init_??.tct to last.tct
mv {work_dir}/output/init_{num_norms}.tct {work_dir}/output/last.tct
mv {work_dir}/output/init_{num_norms}.tctx {work_dir}/output/last.tctx
""".format(work_dir=work_dir, num_norms=args.num_norms)
            # Submit job to que
            hpc.submit_job(SCRIPT_NAME, cmd_line, args.np, work_dir, output_dir,
                           que_name=args.que_name, dry_run=args.dry_run,
                           copy_to_output=['summary.txt'])

