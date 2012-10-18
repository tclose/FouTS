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
CONFIGURATIONS = ['x-curve-y_x-curve--y']
CONFIGURATIONS = ['x', 'x-n2', 'x_xy', 'x_xy-n3']
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
parser.add_argument('--img_dim', default=3, type=int, help='The size of the noisy image to fit against')
parser.add_argument('--step_scale', default=0.001, type=float, help='The scale of the steps used for the metropolis sampling')
parser.add_argument('--num_iterations', default=100000, type=int, help='The number of interations in the metropolis sampling')
parser.add_argument('--sample_period', default=1000, type=int, help='The sample period of the metropolis sampling')
parser.add_argument('--degree', default=3, type=int, help='The degree of the strands to fit')
parser.add_argument('--num_width_sections', default=4, help='The number of samples to use across a Fourier tracts cross-section')
parser.add_argument('--interp_type', default='sinc', type=str, help='The type of interpolation used in the reference image')
parser.add_argument('--interp_extent', default=1, type=int, help='The interpolation extent used in the reference image')
parser.add_argument('--assumed_interp_extent', default=1, type=int, help='The interpolation type used in the likelihood images')
parser.add_argument('--prior_freq', default=[25.0], type=float, nargs='+', help='The scaling of the frequency prior')
parser.add_argument('--prior_aux_freq', default=[75.0], type=float, nargs='+', help='The scaling of the frequency prior')
parser.add_argument('--prior_density_high', default=[1.0], type=float, nargs='+', help='The scaling of the density prior')
parser.add_argument('--prior_density_low', default=[1.0], type=float, nargs='+', help='The scaling of the density prior')
parser.add_argument('--prior_hook', default=[50.0], type=float, nargs='+', help='The scaling of the density prior')
parser.add_argument('--prior_thin', default=[0.0], type=float, nargs='+', help='The scaling of the density prior')
parser.add_argument('--width_epsilon', default=[0.0], type=float, nargs='+', help='The amount of width epsilon to use')
parser.add_argument('--length_epsilon', default=[0.0], type=float, nargs='+', help='The amount of length epsilon to use')
parser.add_argument('--img_snr', default=20.0, type=float, help='The snr to used in the noisy image')
parser.add_argument('--like_snr', default=[20.0], type=float, nargs='+', help='The assumed snr to used in the likelihood function in \
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
ranging_params = [args.prior_freq, args.prior_aux_freq, args.prior_density_low,
                  args.prior_density_high, args.prior_hook, args.prior_thin, args.like_snr, args.width_epsilon,
                  args.length_epsilon]
num_param_sets = max([len(p) for p in ranging_params])
for par in ranging_params:
    if len(par) == 1:
        par *= num_param_sets
    elif len(par) != num_param_sets:
        raise Exception('Number of params ({}) does not match number of params ({})'.
                                                                    format(len(par), num_param_sets))
# Get parameters directory
param_dir = os.path.join(hpc.get_project_dir(), 'params')
output_parent_dir = os.path.realpath(os.path.join(os.environ['HOME'], 'output'))
# Get reference signal to compare noise snr against (the maximum b0 reading from a single straight x tract
gen_img_cmd = "generate_image {param_dir}/fibre/tract/noise_ref.tct {output_dir}/noise_ref.mif -exp_type {args.interp_type} \
-exp_interp_extent {args.interp_extent} -clean -exp_num_width_sections {args.num_width_sections} \
-diff_encodings_location {param_dir}/diffusion/encoding_60.b".format(param_dir=param_dir,
                                      output_dir=output_parent_dir, args=args)
try:
#    print sp.check_output("generate_image --help", shell=True, env=os.environ.copy())
    sp.check_call(gen_img_cmd, shell=True, env=os.environ.copy())
except Exception as e:
    raise Exception('Generate image command: ''{0}'' caused an error ''{1}'''.format(gen_img_cmd, e))
print "Generated noise reference image {}/noise_ref.mif', with command: {}".format(output_parent_dir, gen_img_cmd)
noise_ref_signal = sp.check_output('maxb0 {}/noise_ref.mif'.format(output_parent_dir), shell=True,
                                                                                                env=os.environ.copy())
# Generate a random seed to seed the random number generators of the cmds
seed = int(time.time() * 100)
for i in xrange(args.num_runs):
    for prior_freq, prior_aux_freq, prior_density_low, prior_density_high, prior_hook, prior_thin, like_snr, \
                                                                width_epsilon, length_epsilon in zip(*ranging_params):
        for config in CONFIGURATIONS:
            # Create work directory and get path for output directory
            work_dir, output_dir = hpc.create_work_dir(SCRIPT_NAME, args.output_dir, required_dirs=REQUIRED_DIRS)
            # Create a file in the output directory with just the configuration printed in it (usefulf for quickly 
            # determining what the configuration is
            with open(os.path.join(work_dir, 'output', 'config_name'), 'w') as config_name_file:
                config_name_file.write(config + '\n')
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
# Copy configuration to output directory for future reference
cp {config_path} {work_dir}/output/config.tct
if [ -f {config_path}x ]; then
    cp {config_path}x {work_dir}/output/config.tctx
fi

# Generate image
generate_image {work_dir}/output/config.tct {work_dir}/output/image.mif \
-img_dims "{img_dim} {img_dim} {img_dim}" -exp_type {args.interp_type} \
-exp_interp_extent {args.interp_extent} -noise_snr {args.img_snr} -noise_type gaussian -noise_ref_signal {noise_ref_signal} \
-diff_encodings_location {work_dir}/params/diffusion/encoding_60.b -exp_num_width_sections {args.num_width_sections}

# Initialise fibres
init_fibres {work_dir}/output/init.tct -num_fibres {num_tracts} \
-img_dims "{img_dim} {img_dim} {img_dim}" -degree {args.degree} -seed {init_seed} -base_intensity 1.0

# Run metropolis
metropolis {work_dir}/output/image.mif {work_dir}/output/init.tct {work_dir}/output/samples.tst -like_snr {like_snr} \
-exp_interp_extent {args.assumed_interp_extent} -walk_step_scale {args.step_scale} -num_iter {args.num_iterations} \
-sample_period {args.sample_period} -diff_encodings_location {work_dir}/params/diffusion/encoding_60.b \
-seed {metropolis_seed} -prior_freq {prior_freq} {prior_aux_freq} -prior_density {prior_density_high} \
{prior_density_low} 100 -prior_hook {prior_hook} 100 -prior_thin {prior_thin} 2 -exp_num_width_sections {args.num_width_sections} \
 -exp_type {args.interp_type} -exp_width_epsilon {width_epsilon} -exp_length_epsilon {length_epsilon}
    
# Run analysis
stats_fibres {config_path} {work_dir}/output/samples.tst
    """.format(work_dir=work_dir, config_path=config_path, config=config, args=args, noise_ref_signal=noise_ref_signal,
               num_tracts=num_tracts, img_dim=img_dim, init_seed=seed, metropolis_seed=seed + 1, prior_freq=prior_freq,
               prior_aux_freq=prior_aux_freq, prior_density_low=prior_density_low, prior_density_high=prior_density_high,
               prior_hook=prior_hook, prior_thin=prior_thin, like_snr=like_snr, width_epsilon=width_epsilon,
               length_epsilon=length_epsilon)
            # Submit job to que
            hpc.submit_job(SCRIPT_NAME, cmd_line, args.np, work_dir, output_dir, que_name=args.que_name,
                                                                    dry_run=args.dry_run, copy_to_output=REQUIRED_DIRS)
        seed += 2
