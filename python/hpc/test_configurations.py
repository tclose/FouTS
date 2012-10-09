#!/usr/bin/env python
"""
 This script submits scripts to generate random prior shapes from noisy images and submit them to the job que on a sun
 grid engine

 Author: Tom Close (tclose@oist.jp)
 Created: 6/8/2012
"""
#Name of the script for the output directory and submitted mpi job
SCRIPT_NAME = 'test_configurations'
CONFIGURATIONS = ['layer-n5-d5', 'x-curve-z_y-curve--z', 'x-small', 'x_y', 'x-big', 'x-pos-yz', 'x',
                'yz-curve-x', 'x-curve-y_x-curve--y', 'x-rotate-big', 'x_xxy', 'x-curve-z',
                'x-rotate', 'x_xy', 'x-curve-z_y-curve--z']
REQUIRED_DIRS = ['params/fibre/tract/test_configurations', 'params/diffusion']
# Required imports
import hpc
import argparse
import os.path
import re
import subprocess as sp
import time
# Arguments that can be given to the script
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--img_dim', default=10, type=int, help='The size of the noisy image to fit against')
parser.add_argument('--step_scale', default=0.001, type=float, help='The scale of the steps used for the metropolis sampling')
parser.add_argument('--num_iterations', default=200000, type=int, help='The number of interations in the metropolis sampling')
parser.add_argument('--sample_period', default=1000, type=int, help='The sample period of the metropolis sampling')
parser.add_argument('--degree', default=3, type=int, help='The degree of the strands to fit')
parser.add_argument('--num_width_sections', default=15, help='The number of samples to use across a Fourier tracts cross-section')
parser.add_argument('--interp_type', default='sinc', type=str, help='The type of interpolation used in the reference image')
parser.add_argument('--interp_extent', default=1, type=int, help='The interpolation extent used in the reference image')
parser.add_argument('--assumed_interp_extent', default=1, type=int, help='The interpolation type used in the likelihood images')
parser.add_argument('--prior_freq', default=[25.0], type=float, action='append', help='The scaling of the frequency prior (currently unused)')
parser.add_argument('--prior_aux_freq', default=[75.0], type=float, action='append', help='The scaling of the frequency prior (currently unused)')
parser.add_argument('--prior_density_high', default=[1.0], type=float, action='append', help='The scaling of the density prior (currently unused)')
parser.add_argument('--prior_density_low', default=[1.0], type=float, action='append', help='The scaling of the density prior (currently unused)')
parser.add_argument('--prior_hook', default=[50.0], type=float, action='append', help='The scaling of the density prior (currently unused)')
parser.add_argument('--img_snr', default=5.0, type=float, help='The snr to used in the noisy image')
parser.add_argument('--like_snr', default=[20.0], type=float, action='append', help='The assumed snr to used in the likelihood function in \
the metropolis sampling')
parser.add_argument('--output_dir', default=None, type=str, help='The parent directory in which the output directory \
will be created (defaults to $HOME/Output)')
parser.add_argument('--num_runs', default=1, type=int, help='The number of runs to submit to the que')
parser.add_argument('--dry_run', action='store_true', help='Only perform a dry run (create jobscript then quit)')
parser.add_argument('--np', type=int, default=1, help='The the number of processes to use for the simulation \
(default: %(default)s)')
parser.add_argument('--que_name', type=str, default='short', help='The the que to submit the job to \
(default: %(default)s)')
args = parser.parse_args()
# For the following parameters to this script, ensure that number of parameter values match, or if they are a singleton 
# list it is assumed to be constant and that value that value is replicated to match the number of other of other 
# parameters in the set
num_param_sets = max((len(args.prior_freq), len(args.prior_aux_freq), len(args.prior_density_low),
                                            len(args.prior_density_high), len(args.prior_hook), len(args.like_snr)))
if len(args.prior_freq) == 1:
    args.prior_freq *= num_param_sets
elif len(args.prior_freq) != num_param_sets:
    raise Exception('Number of prior frequency params ({}) does not match number of params ({})'.
                                                                    format(len(args.prior_freq), num_param_sets))
if len(args.prior_aux_freq) == 1:
    args.prior_aux_freq *= num_param_sets
elif len(args.prior_aux_freq) != num_param_sets:
    raise Exception('Number of prior auxiliary frequency params ({}) does not match number of params ({})'.
                                                                    format(len(args.prior_aux_freq), num_param_sets))
if len(args.prior_density_high) == 1:
    args.prior_density_high *= num_param_sets
elif len(args.prior_density_high) != num_param_sets:
    raise Exception('Number of prior frequency params ({}) does not match number of params ({})'.
                                                                    format(len(args.prior_density_high), num_param_sets))
if len(args.prior_density_low) == 1:
    args.prior_density_low *= num_param_sets
elif len(args.prior_density_low) != num_param_sets:
    raise Exception('Number of prior frequency params ({}) does not match number of params ({})'.
                                                                    format(len(args.prior_density_low), num_param_sets))
if len(args.prior_hook) == 1:
    args.prior_hook *= num_param_sets
elif len(args.prior_hook) != num_param_sets:
    raise Exception('Number of prior frequency params ({}) does not match number of params ({})'.
                                                                    format(len(args.prior_hook), num_param_sets))
if len(args.like_snr) == 1:
    args.like_snr *= num_param_sets
elif len(args.like_snr) != num_param_sets:
    raise Exception('Number of prior frequency params ({}) does not match number of params ({})'.
                                                                    format(len(args.like_snr), num_param_sets))
# Get parameters directory
param_dir = os.path.join(hpc.get_project_dir(), 'params')
# Get reference signal to compare noise snr against (the maximum b0 reading from a single straight x tract
gen_img_cmd = "generate_image {param_dir}/fibre/tract/noise_ref.tct /tmp/noise_ref.mif -exp_type {args.interp_type} \
-exp_interp_extent {args.interp_extent} -clean -exp_num_width_sections {args.num_width_sections} \
-diff_encodings_location {param_dir}/diffusion/encoding_60.b".format(param_dir=param_dir, args=args)
try:
#    print sp.check_output("generate_image --help", shell=True, env=os.environ.copy())
    sp.check_call(gen_img_cmd, shell=True, env=os.environ.copy())
except Exception as e:
    raise Exception('Generate image command: ''{0}'' caused an error ''{1}'''.format(gen_img_cmd, e))
noise_ref_signal = sp.check_output('maxb0 /tmp/noise_ref.mif', shell=True, env=os.environ.copy())
# Generate a random seed to seed the random number generators of the cmds
seed = int(time.time() * 100)
for i in xrange(args.num_runs):
    for prior_freq, prior_aux_freq, prior_density_low, prior_density_high, prior_hook, like_snr in zip(args.prior_freq,
                   args.prior_aux_freq, args.prior_density_low, args.prior_density_low, args.prior_hook, args.like_snr):
        for config in CONFIGURATIONS:
            # Create work directory and get path for output directory
            work_dir, output_dir = hpc.create_work_dir(SCRIPT_NAME, args.output_dir, required_dirs=REQUIRED_DIRS)
            # Strip configuration of symbols for tract number and img dimension
            config_name = re.sub('\-n[0-9]+', '', config) # Strip tract number
            config_name = re.sub('\-d[0-9]+', '', config_name) # Strip dimension
            # Get the configuration path
            config_path = '{work_dir}/params/fibre/tract/test_configurations/{config_name}.tct'.format(work_dir=work_dir,
                                                                                                 config_name=config_name)
            # Try to extract number of tracts to use for sampling from configuration
            num_tracts = re.findall('\-n[0-9]+', config)
            # If number of tracts not given in configuration name, match that of the configuration file.
            if num_tracts:
                num_tracts = int(num_tracts[0].lstrip('-n'))
            else:

                try:
                    num_tracts = int(sp.check_output("fibre_info {config_path} | grep total_count | awk '{{print $2}}'".
                                                                            format(config_path=config_path), shell=True))
                except OSError as e:
                    raise Exception("Could not open configurations file ''{0}'' for reading (returned: {1})".format(config_path, e))
                except ValueError as e:
                    raise Exception("Could not read 'total_count' property from configuration file '{config_path}', could \
    be corrupted.\n(Error: {e})".format(config_path=config_path, e=e))
                else:
                    num_tracts = int(num_tracts)
            img_dim = re.findall('\-d[0-9]+', config)
            if img_dim:
                img_dim = float(img_dim[0].lstrip('-d'))
            else:
                img_dim = args.img_dim
            cmd_line = """      
    # Generate image        
    generate_image {config_path} {work_dir}/output/image.mif $noise \
    img_dims "{img_dim} {img_dim} {img_dim}" -exp_type {args.interp_type} \
    -exp_interp_extent {args.interp_extent} -noise_snr {args.img_snr} -noise_type gaussian -noise_ref_signal {noise_ref_signal} \
    -diff_encodings_location {work_dir}/params/diffusion/encoding_60.b
    
    # Initialise fibres
    init_fibres {work_dir}/output/init.tct -num_fibres {num_tracts} \
    -img_dims "{img_dim} {img_dim} {img_dim}" -degree {args.degree} -seed {init_seed} -acs 0.01 -base_intensity 1.0
    
    # Run metropolis
    metropolis {work_dir}/output/image.mif {work_dir}/output/init.tct {work_dir}/output/samples.tst -like_snr {like_snr} \
    -exp_interp_extent {args.assumed_interp_extent} -walk_step_scale {args.step_scale} -num_iter {args.num_iterations} \
    -sample_period {args.sample_period} -diff_encodings_location {work_dir}/params/diffusion/encoding_60.b \
    -seed {metropolis_seed} -prior_freq {prior_freq} {prior_aux_freq} -prior_density {prior_density_high} \
    {prior_density_low} 100 -prior_hook {prior_hook} 100
        
    # Run analysis
    stats_fibres {config_path} {work_dir}/output/samples.tst
    """.format(work_dir=work_dir, config_path=config_path, config=config, args=args, noise_ref_signal=noise_ref_signal,
               num_tracts=num_tracts, img_dim=img_dim, init_seed=seed, metropolis_seed=seed + 1, prior_freq=prior_freq,
               prior_aux_freq=prior_aux_freq, prior_density_low=prior_density_low, prior_density_high=prior_density_high,
               prior_hook=prior_hook, like_snr=like_snr)
            # Submit job to que
            hpc.submit_job(SCRIPT_NAME, cmd_line, args.np, work_dir, output_dir, que_name=args.que_name, dry_run=args.dry_run,
                                                                                             copy_to_output=REQUIRED_DIRS)
        seed += 2
