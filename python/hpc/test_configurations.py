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
import re
import subprocess as sp

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--img_dims', default=10, type=int, help='The size of the noisy image to fit against')
parser.add_argument('--step_scale', default=0.001, type=float, help='')
parser.add_argument('--num_iter', default=200000, type=int, help='')
parser.add_argument('--sample_period', default=1000, type=int, help='')
parser.add_argument('--degree', default=3, type=int, help='The degree of the strands to fit')
parser.add_argument('--num_width_sections', default=15, help='The number of samples to use across a Fourier tracts cross-section')
parser.add_argument('--interp_type', default='sinc', type=str, help='')
parser.add_argument('--interp_extent', default=3, type=int, help='')
parser.add_argument('--assumed_interp_extent', default=1, type=int, help='')
parser.add_argument('--prior_freq', default=10, type=int, help='')
parser.add_argument('--prior_density', default=10, type=int, help='')
parser.add_argument('--img_snr', default=5, type=float, help='The snr to used in the noisy image')
parser.add_argument('--like_snr', default=20, type=float, help='The assumed snr to used in the likelihood function in \
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

CONFIGURATIONS = ['layer.tct', 'x-curve-z_y-curve-z.tct', 'x-small.tct', 'x_y.tct', 'x-big.tct', 'x-pos-yz.tct', 'x.tct',
                'yz-curve-x.tct', 'x-curve-y_x-curve--y.tct', 'x-rotate-big.tct', 'x_xxy.tct', 'x-curve-z.tct',
                'x-rotate.tct', 'x_xy.tct', 'x-curve-z_y-curve--z.tct']

REQUIRED_DIRS = ['fibre/tract/test_configurations', 'diffusion']


param_dir=os.path.join(hpc.get_project_dir(), 'params')
gen_img_cmd="generate_image {param_dir}/fibre/tract/noise_ref.tct /tmp/noise_ref.mif -exp_type {args.interp_type} \
-exp_interp_extent {args.interp_extent} -clean -exp_num_width_sections {args.num_width_sections}".format(
                                                                                        param_dir=param_dir, args=args)
print gen_img_cmd
# Get reference signal to compare noise snr against (the maximum b0 reading from a single straight x tract
sp.check_call(gen_img_cmd)
noise_ref_signal=sp.check_output('maxb0 /tmp/noise_ref.mif')

for i in xrange(args.num_runs):
    for config in CONFIGURATIONS:
        # Create work directory and get path for output directory
        work_dir, output_dir = hpc.create_work_dir(SCRIPT_NAME, args.output_dir, required_dirs=REQUIRED_DIRS)
        # Set up command to run the script
        if config.endswith('clean'):
            noise_options = '--noise_snr {0}'.format(args.img.snr)
        else:
            noise_options = ''
        num_tracts_str = re.findall('[0-9]+', config)
        if not num_tracts_str:
            num_tracts_str = sp.check_output('fibre_info {work_dir}/fibre/tract/test_configurations/{config}.tct | \
grep total_count | awk ''{print $2}''')
        try:
            num_tracts = int(num_tracts_str)
        except ValueError:
            raise Exception ('Could not convert number of tracts str ''{0}'' to integer'.format(num_tracts_str))
        config_name = re.sub('[0-9]', '', config)
        cmd_line = """      
# Generate image        
generate_image {work_dir}/params/fibre/tract/test_configurations/{config_name}.tct {work_dir}/{config}.mif $noise \
img_dims "{args.img_dim} {args.img_dim} {args.img_dim}" -exp_type {args.interp_type} \
-exp_interp_extent {args.interp_extent} -noise_snr {args.img_snr} -noise_type gauss -noise_ref_signal {noise_ref_signal} \
-diff_encoding_location {work_dir}/diffusion_encoding_60.b

# Initialise fibres
init_cmd="init_fibres {work_dir}/{config}_init.tct -num_fibres {num_tracts} \
-img_dims "{args.img_dim} {args.img_dim} {args.img_dim}" -degree {args.degree}

# Run metropolis
metropolis {work_dir}/{config}.mif {work_dir}/{config}_init.tct {work_dir}/{config}.tst -like_snr {args.like_snr} \
-exp_interp_extent {args.assumed_interp} -walk_step_scale {args.step_scale} -num_iter {args.num_iterations} \
-sample_period {args.sample_period} -diff_encoding_location {work_dir}/diffusion_encoding_60.b

# Run analysis
stats_fibres {work_dir}/{config_name}.tct {work_dir}/{config}.tst
""".format(work_dir=work_dir, config_name=config_name, config=config, args=args, noise_ref_signal=noise_ref_signal,
           num_tracts=num_tracts)
        # Submit job to que
        hpc.submit_job(SCRIPT_NAME, cmd_line, args.np, work_dir, output_dir, que_name=args.que_name, dry_run=args.dry_run,
                                                                                         copy_to_output=REQUIRED_DIRS)
