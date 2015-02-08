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
parser.add_argument('--output_dir', default=None, type=str, help='The parent directory in which the output directory will be created (defaults to $HOME/output)')
parser.add_argument('--dry_run', action='store_true', help='Only perform a dry run (create jobscript then quit)')
parser.add_argument('--np', type=int, default=1, help='The the number of processes to use for the simulation (default: %(default)s)')
parser.add_argument('--que_name', type=str, default='long', help='The the que to submit the job to (default: %(default)s)')
parser.add_argument('--num_iterations', type=int, default=100000, help='the number of iterations to perform (default: %(default)s)')
parser.add_argument('--sample_period', type=int, default=1000, help='the number of iterations before sampling (default: %(default)s)')
parser.add_argument('--num_length', type=int, default=30, help='The number of samples along the fibre length (default: %(default)s)')
parser.add_argument('--perturb_scale', type=float, default=0.001, help='the perturbation applied to the true configuration (default: %(default)s)')
parser.add_argument('--like_snr', default=20, type=float, help='The assumed snr to used in the likelihood function in the metropolis sampling (default: %(default)s)')
parser.add_argument('--interp_extent', type=float, default=1.0, help='the extent of the interpolation kernel (default: %(default)s)')
parser.add_argument('--step_scale', type=float, default=0.002, help='scale the metropolis steps (default: %(default)s)')
parser.add_argument('--num_strands', type=int, default=[1, 3, 5, 7, 8, 9], nargs='+',
                    help='The number of strands to test (default: %(default)s)')
parser.add_argument('--end_on_sphere_scale', type=float, default=0, help="The scale of the end_on_sphere prior (default: %(default)s)")
parser.add_argument('--end_on_sphere_radius', type=float, default=0.225, help="The radius of the end_on_sphere prior (default: %(default)s)")
args = parser.parse_args()
# For the following parameters to this script, ensure that number of parameter values match, or if they are a singleton
# list it is assumed to be constant and that value that value is replicated to match the number of other of other
# parameters in the set. Otherwise if the '--combo' option is provided then loop through all combinations of the
# provided parameters.

required_dirs = [os.path.join('params', 'fibre', 'tract', 'single'),
                 os.path.join('params', 'fibre', 'strand', 'masks'),
                 os.path.join('params', 'diffusion')]
for num_strands in args.num_strands:
    # Create work directory and get path for output directory
    work_dir, output_dir = hpc.create_work_dir(SCRIPT_NAME, args.output_dir,
                                               required_dirs=required_dirs)
    with open(os.path.join(work_dir, 'output', 'num_strands'), 'w') as f:
        f.write('{}\n'.format(num_strands))
    # Set up command to run the script
    cmd_line = """
select_fibres {work_dir}/params/fibre/tract/single/x.tct {work_dir}/noise_ref.str -num_width 1

set_properties {work_dir}/noise_ref.str -set_elem acs 1.0

generate_image {work_dir}/noise_ref.str {work_dir}/noise_ref.mif \
-exp_num_length {args.num_length} -img_dim 1,1,1 -diff_encodings {work_dir}/params/diffusion/encoding_60.b \
 -exp_base_intensity 1 -diff_isotropic -exp_type sinc -exp_interp {args.interp_extent} -clean

NOISE_REF=`maxb0 {work_dir}/noise_ref.mif`

select_fibres {work_dir}/params/fibre/tract/single/x.tct {work_dir}/output/true.str -num_width {num_strands}

ACS=$(echo "1.0/`fibre_info {work_dir}/output/true.str | grep total_count | awk '{{print $2}}'`" | bc -l)

set_properties {work_dir}/output/true.str -set_elem acs $ACS

#Generate image with appropriate SNR.
generate_image {work_dir}/output/true.str {work_dir}/output/image.mif \
-exp_num_length {args.num_length} -img_dim 3,3,3 -diff_encodings {work_dir}/params/diffusion/encoding_60.b \
-exp_base_intensity 1 -diff_isotropic -exp_type sinc -exp_interp {args.interp_extent} -clean

perturb_fibres {work_dir}/output/true.str {work_dir}/output/init.str -std {args.perturb_scale} \
-scales_location {work_dir}/params/fibre/strand/masks/no_intens.str

#Perform MH sampling
metropolis {work_dir}/output/image.mif {work_dir}/output/init.str \
{work_dir}/output/samples.sst -exp_num_length {args.num_length} \
-diff_encodings {work_dir}/params/diffusion/encoding_60.b \
-exp_base_intensity 1 -like_snr {args.like_snr} -walk_step_scale {args.step_scale} \
-num_iter {args.num_iterations} -sample_period {args.sample_period} -diff_isotropic -exp_type sinc \
-exp_interp {args.interp_extent} -like_ref_signal $NOISE_REF -walk_step_location \
{work_dir}/params/fibre/strand/masks/mcmc/metropolis/default.str \
-prior_end_on_sphere {args.end_on_sphere_scale} {args.end_on_sphere_radius}

stats_fibres {work_dir}/output/true.str {work_dir}/output/samples.sst
    """.format(work_dir=work_dir, num_strands=num_strands, args=args)
        # Submit job to que
    hpc.submit_job(SCRIPT_NAME, cmd_line, args.np, work_dir, output_dir, que_name=args.que_name,
                   dry_run=args.dry_run)
