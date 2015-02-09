#!/usr/bin/env python
"""
 This script submits scripts to generate random prior shapes from noisy images
 and submit them to the job que on a sun grid engine

 Author: Tom Close (tclose@oist.jp)
 Created: 6/8/2012
"""
# Required imports
from fouts import hpc
import argparse
import os.path
import time

# Name of the script for the output directory and submitted mpi job
SCRIPT_NAME = 'phantom_sampling'
# Required dirs for the script to run
REQUIRED_DIRS = ['params/image/phantom/paper/from_tcks/one/trim',
                 'params/diffusion',
                 'params/fibre/tract/masks/mcmc/metropolis']


def sampling_cmd(args, work_dir, random_seed, phantom_index):
    dataset_path = os.path.join(
        work_dir, 'params', 'image', 'phantom', 'paper', 'from_tcks',
        args.true_interp_style, 'trim', str(args.voxel_res))
    mask_path = os.path.join(
        work_dir, 'params', 'fibre', 'tract', 'masks', 'mcmc', 'metropolis',
        'default.tct')
    cmd = ("""
        init_fibres {work_dir}/init.tct -degree {degree} -num_fibres  \\
        {num_tracts} -img_dims 3,3,3 -curve_stddev 0.001 \\
        -base_intensity 1 -width_epsilon {width_epsilon} \\
        -width_stddev {init_width_stddev} -length_stddev {init_length} \\
        -width_mean {init_width_mean}  -length_epsilon {length_epsilon} \\
        -img_vox_lengths {voxel_res},{voxel_res},{voxel_res} \\

        # run the metropolis algorithm
        metropolis {dataset_path}/{phantom_index}.mif {work_dir}/init.tct \\
        {work_dir}/samples.tst -like_snr {assumed_snr} \\
        -exp_type {interp_type} -exp_interp_extent {interp_extent} \\
        -walk_step_scale {step_scale} -seed {seed} \\
        -exp_num_width_sections {num_width_sections} \\
        -exp_base_intensity `cat {dataset_path}/maxb0.ratio.txt` \\
        -exp_num_length_sections {num_length_sections} \\
        -walk_step_location {mask_path} \\
        -num_iterations {num_iterations} -sample_period {sample_period} \\
        -prior_freq {prior_freq} {prior_aux_freq} \\
        -prior_density {prior_density_high} {prior_density_low} 100 \\
        -prior_hook {prior_hook} 100 15 -diff_warn -save_image
        """.format(
        work_dir=work_dir, phantom_index=phantom_index,
        step_scale=args.step_scale, voxel_res=args.voxel_res,
        num_iterations=args.num_iterations,
        sample_period=args.sample_period,
        num_length_sections=args.num_length_sections,
        dataset_path=dataset_path,
        num_width_sections=args.num_width_sections, num_tracts=args.num_tracts,
        interp_type=args.interp_type, interp_extent=args.interp_extent,
        init_perturb_stddev=args.init_perturb_stddev,
        mask_path=mask_path, degree=args.degree, prior_freq=args.prior_freq,
        prior_aux_freq=args.prior_aux_freq,
        prior_density_high=args.prior_density_high,
        prior_density_low=args.prior_density_low, prior_hook=args.prior_hook,
        assumed_snr=args.assumed_snr, init_acs=args.init_acs,
        init_length=args.init_length, init_width_mean=args.init_width_mean,
        init_width_stddev=args.init_width_stddev,
        width_epsilon=args.width_epsilon,
        length_epsilon=args.length_epsilon, seed=random_seed))
    return cmd

# Arguments that can be given to the script
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--voxel_res', default=0.15, type=float,
                    help="The resolution of the image to fit "
                         "(default: %(default)s)")
parser.add_argument('--step_scale', default=0.001, type=float,
                    help="The scale of the steps used for the metropolis "
                         "sampling (default: %(default)s)")
parser.add_argument('--num_iterations', default=50000, type=int,
                    help="The number of interations in the metropolis sampling"
                         " (default: %(default)s)")
parser.add_argument('--sample_period', default=1000, type=int,
                    help="The sample period of the metropolis sampling "
                         "(default: %(default)s)")
parser.add_argument('--num_length_sections', default=20, type=float)
parser.add_argument('--num_width_sections', default=4, type=float)
parser.add_argument('--num_tracts', default=5, type=float,
                    help="The number of tracts to fit (default: %(default)s")
parser.add_argument('--interp_type', default='sinc', type=str,
                    help="The type of interpolation used in the reference "
                         "image (default: %(default)s)")
parser.add_argument('--interp_extent', default=1, type=int,
                    help="The interpolation extent used in the reference image"
                         " (default: %(default)s)")
parser.add_argument('--true_interp_style', default='one', type=str,
                    help="The style of interpolation extent used in the "
                         "reference image, either 'full' or 'one'")
parser.add_argument('--init_perturb_stddev', default=0.2, type=float)
parser.add_argument('--degree', default=3, type=int,
                    help="The degree of the fibre used to sample from "
                         "(default: %(default)s)")
parser.add_argument('--prior_freq', default=15.0, type=float,
                    help="The scaling of the frequency prior (default: "
                         "%(default)s)")
parser.add_argument('--prior_aux_freq', default=60.0, type=float,
                    help="The scaling of the auxiliary frequency prior "
                         "(default: %(default)s)")
parser.add_argument('--prior_density_high', default=1, type=float,
                    help="The scaling of the density prior (default: "
                         "%(default)s)")
parser.add_argument('--prior_density_low', default=1, type=float,
                    help="The scaling of the density prior (default: "
                         "%(default)s)")
parser.add_argument('--prior_hook', default=100000.0, type=float,
                    help="The scaling of the density prior (default: "
                         "%(default)s)")
parser.add_argument('--assumed_snr', default=20, type=float)
parser.add_argument('--init_acs', default=0.5, type=float)
parser.add_argument('--init_length', default=0.01, type=float)
parser.add_argument('--init_width_mean', default=0.025, type=float)
parser.add_argument('--init_width_stddev', default=0.0025, type=float)
parser.add_argument('--width_epsilon', default=0.001, type=float)
parser.add_argument('--length_epsilon', default=0.001, type=float)
parser.add_argument('--random_seed', type=int, default=None,
                    help="The random seed for the whole run")
parser.add_argument('--output_dir', default=None, type=str,
                    help="The parent directory in which the output directory "
                         "will be created (defaults to $HOME/Output)")
parser.add_argument('--num_runs', default=1, type=int,
                    help="The number of runs to submit to the que "
                         "(default: %(default)s).")
parser.add_argument('--dry_run', action='store_true',
                    help="Only perform a dry run (create jobscript then quit)")
parser.add_argument('--np', type=int, default=1,
                    help="The the number of processes to use for the "
                         "simulation (default: %(default)s)")
parser.add_argument('--que_name', type=str, default='long',
                    help="The the que to submit the job to  (default: "
                         "%(default)s)")
args = parser.parse_args()

# Generate a random seed to seed the random number generators of the cmds
if args.random_seed is None:
    random_seed = int(time.time() * 100)
    print "Using random_seed {}".format(random_seed)
else:
    random_seed = args.random_seed
for run_i in xrange(args.num_runs):
    for phantom_i in xrange(1, 11):
        # Create work directory and get path for output directory
        work_dir, output_dir = hpc.create_work_dir(
            SCRIPT_NAME, args.output_dir, required_dirs=REQUIRED_DIRS)
        # Create a file in the output directory with just the dataset
        # printed in it (useful for quickly determining what the dataset is)
        with open(os.path.join(work_dir, 'output', 'dataset_name.txt'),
                  'w') as dataset_file:
            dataset_file.write('{}/{}\n'.format(args.voxel_res, phantom_i))
        # Create directory to save images for debugging.
        os.mkdir(os.path.join(work_dir, 'images'))
        cmd_line = sampling_cmd(args=args, work_dir=work_dir,
                                phantom_index=phantom_i,
                                random_seed=random_seed + 1)
        # Submit job to que
        hpc.submit_job(SCRIPT_NAME, cmd_line, args.np, work_dir, output_dir,
                       que_name=args.que_name, dry_run=args.dry_run,
                       copy_to_output=['summary.txt'])
        # Increment the random seed by two (one for the init and another for
        # the metropolis) so the next run is run with a different random seed
        random_seed += 2

# #select tractlets to tracks parameters,
# strand_per_vol=10000,
# true_strands_per_area=100000,
# include_range=50:99,
# true_num_length_sections=1000,
#
# tdi_degree=20,
# tdi_vox=0.0075,0.0075,0.0444,
# tdi_suffix='',
# trim_dist=0.375,
#
# #cleanup samples parameters,
# display_min_vol=0.0,
#
# true_ext='tck',
# num_tractlets=8,
# num_final_tractlets=8,
# save_strands_per_area=1000,
