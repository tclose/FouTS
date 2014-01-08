#!/usr/bin/env python
"""
 This script submits scripts to generate random prior shapes from noisy images and submit them to the job que on a sun
 grid engine

 Author: Tom Close (tclose@oist.jp)
 Created: 6/8/2012
"""

# Name of the script for the output directory and submitted mpi job
SCRIPT_NAME = 'strand_number'

import hpc
import argparse
import os.path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--img_dims', default=10, type=int, help='The size of the noisy image to fit against')
parser.add_argument('--like_snr', default=20, type=float, help='The assumed snr to used in the likelihood function in the metropolis sampling')
parser.add_argument('--output_dir', default=None, type=str, help='The parent directory in which the output directory will be created (defaults to $HOME/Output)')
parser.add_argument('--num_runs', default=1, type=int, help='The number of runs to submit to the que')
parser.add_argument('--dry_run', action='store_true', help='Only perform a dry run (create jobscript then quit)')
parser.add_argument('--np', type=int, default=1, help='The the number of processes to use for the simulation (default: %(default)s)')
parser.add_argument('--que_name', type=str, default='short', help='The the que to submit the job to(default: %(default)s)')
parser.add_argument('--num_iterations', type=int, default=100000, help='the number of iterations to perform')
parser.add_argument('--sample_period', type=int, default=1000, help='the number of iterations before sampling')
parser.add_argument('--num_length', type=int, default=30, help='The number of samples along the fibre length')
parser.add_argument('--perturb_scale', type=float, default=0.001, help='the perturbation applied to the true configuration')
parser.add_argument('--interp_extent', type=float, default=1.0, help='the extent of the interpolation kernel')
args = parser.parse_args()
# For the following parameters to this script, ensure that number of parameter values match, or if they are a singleton
# list it is assumed to be constant and that value that value is replicated to match the number of other of other
# parameters in the set. Otherwise if the '--combo' option is provided then loop through all combinations of the
# provided parameters.

num_strands_n_scales = [(1, 0.001), (3, 0.001), (5, 0.001), (7, 0.001)]
required_dirs = [os.path.join('params', 'fibre', 'tract', 'single'),
                 os.path.join('params', 'fibre', 'strand', 'masks'),
                 os.path.join('params', 'diffusion')]
for num_strands, step_scale in num_strands_n_scales:
    # Create work directory and get path for output directory
    work_dir, output_dir = hpc.create_work_dir(SCRIPT_NAME, args.output_dir,
                                               required_dirs=required_dirs)

    acs = 1.0 / float(num_strands)
    # Set up command to run the script
    cmd_line = """
select_fibres {work_dir}/params/fibre/tract/single/x.tct {work_dir}/noise_ref.str -num_width 1

set_properties {work_dir}/noise_ref.str -set_elem acs 1.0

generate_image {work_dir}/noise_ref.str {work_dir}/noise_ref.mif \
-exp_num_length {args.num_length} -img_dim 1,1,1 -diff_encodings {work_dir}/params/diffusion/encoding_60.b \
 -exp_base_intensity 1 -diff_isotropic -exp_type sinc -exp_interp {args.interp_extent} -clean

NOISE_REF=`maxb0 {work_dir}/noise_ref.mif`

select_fibres {work_dir}/params/fibre/tract/single/x.tct {work_dir}/true_{num_strands}.str -num_width {num_strands}

set_properties {work_dir}/true_{num_strands}.str -set_elem acs {acs}

#Generate image with appropriate SNR.
generate_image {work_dir}/true_{num_strands}.str {work_dir}/image_{num_strands}.mif \
-exp_num_length {args.num_length} -img_dim 3,3,3 -diff_encodings {work_dir}/params/diffusion/encoding_60.b \
-exp_base_intensity 1 -diff_isotropic -exp_type sinc -exp_interp {args.interp_extent} -clean

perturb_fibres {work_dir}/true_{num_strands}.str {work_dir}/init_{num_strands}.str -std {args.perturb_scale} \
-scales_location {work_dir}/params/fibre/strand/masks/no_intens.str

#Perform MH sampling
metropolis {work_dir}/image_{num_strands}.mif {work_dir}/init_{num_strands}.str \
{work_dir}/samples_{num_strands}.sst -exp_num_length {args.num_length} \
-diff_encodings {work_dir}/params/diffusion/encoding_60.b \
-exp_base_intensity 1 -like_snr {args.like_snr} -walk_step_scale {step_scale} -num_iter {args.num_iterations} \
-sample_period {args.sample_period} -diff_isotropic -exp_type sinc -exp_interp {args.interp_extent} \
-like_ref_signal $NOISE_REF -walk_step_location \
{work_dir}/params/fibre/strand/masks/mcmc/metropolis/default{args.degree}.tct

stats_fibres {work_dir}/true_{num_strands}.str {work_dir}/samples_{num_strands}.sst
    """.format(work_dir=work_dir, num_strands=num_strands, step_scale=step_scale, args=args,
               acs=acs)
        # Submit job to que
    hpc.submit_job(SCRIPT_NAME, cmd_line, args.np, work_dir, output_dir, que_name=args.que_name,
                   copy_to_output=required_dirs, dry_run=args.dry_run)
