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

CONFIGURATIONS=['x', 'x-big', 'x-small', 'x-pos-yz', 'xy', 'xz', 'xyz', 'x-curve-z', 'yz-curve-x', 'x_y', 'x_xy1', 
                'x_xy', 'x_xy3', 'x_xxy', 'x-curve-z_y-curve-z', 'layer5', 'x-rotate', 'x-curve-y_x-curve--y', 
                'x-rotate-big']

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--img_dims', default=10, type=int, help='The size of the noisy image to fit against')
parser.add_argument('--step_scale', default=0.001, type=float, help='')
parser.add_argument('--num_iter', default=200000, type=int, help='')
parser.add_argument('--sample_period', default=1000, type=int, help='')
parser.add_argument('--interp', default=3, type=int, help='')
parser.add_argument('--assumed_interp', default=1, type=int, help='')
parser.add_argument('--prior_freq', default=10, type=int, help='')
parser.add_argument('--prior_density', default=10, type=int, help='')
parser.add_argument('--img_snr', default=5, type=float, help='The snr to used in the noisy image')
parser.add_argument('--like_snr', default=20, type=float, help='The assumed snr to used in the likelihood function in the metropolis sampling')
parser.add_argument('--output_dir', default=None, type=str, help='The parent directory in which the output directory will be created (defaults to $HOME/Output)')
parser.add_argument('--num_runs', default=1, type=int, help='The number of runs to submit to the que')
parser.add_argument('--dry_run', action='store_true', help='Only perform a dry run (create jobscript then quit)')
parser.add_argument('--np', type=int, default=1, help='The the number of processes to use for the simulation (default: %(default)s)')
parser.add_argument('--que_name', type=str, default='short', help='The the que to submit the job to(default: %(default)s)')
args = parser.parse_args()

for i in xrange(args.num_runs):
    for config in CONFIGURATIONS:
        # Create work directory and get path for output directory
        work_dir, output_dir = hpc.create_work_dir(SCRIPT_NAME, args.output_dir)
        # Set up command to run the script
        if config.endswith('clean'):
            noise_options = '--noise_snr {0}'.format(args.img.snr)
        else:
            noise_options = ''
#        #Either pick up the number of tractlets to use from the configuration name or use the same as in the reference set.
#        tractlet_count=`echo $config | grep -o "[0-9]" | tr -d '\n'`
#        if [ ${#tractlet_count} = 0 ]; then
#          refname=$config
#          tractlet_count=`fibre_info $WORK_DIR/ref/$config.tct -prop count`
#        else
#          refname=${config//[0-9]/}
#        fi            
        cmd_line = """      
# Generate image        
generate_image {work_dir}/ref/$refname.tct {work_dir}/$run_index/$config.mif $noise -exp_type sinc -exp_interp_extent {interp} -noise_snr {img_snr} -noise_type gauss -noise_ref_signal {noise_ref_signal}"

# Initialise fibres
init_cmd="init_fibres {work_dir}/$run_index/${config}_init.tct -num_fibres $tractlet_count $init_radius -degree 3 -length $LENGTH_STDDEV -curve_stddev $CURVE_STDDEV"

# Run metropolis
metropolis {work_dir}/$run_index/$config.mif {work_dir}/$run_index/${config}_init.tct {work_dir}/$run_index/$config.tst -like_snr $NOISE_SNR -exp_type $ASSUMED_INTERP_TYPE -exp_interp_extent $ASSUMED_INTERP_EXTENT -walk_step_scale $WALK_STEP_SCALE -dty_prob_scale $DENS_PROB_SCALE -prior_tract_sheer $TRACT_SHEER_SCALE -prior_end_on_sphere $ENDPOINTS_SCALE $ENDPOINTS_RADIUS -prior_length $PRIOR_LENGTH_SCALE -num_iter $NUM_ITERATIONS -sample_period $SAMPLE_PERIOD"

# Run analysis
stats_fibres {work_dir}/ref/$refname.tct {work_dir}/$run_index/$config.tst
        """.format(work_dir=work_dir, dim=args.img_dims, acs_mean=args.prior_acs[1], img_snr=args.img_snr, like_snr=args.like_snr, encoding_loc=os.path.join(hpc.get_project_dir(), 'params', 'diffusion', 'encoding_1.b'), acs_stddev=args.prior_acs[0])
        # Submit job to que
        hpc.submit_job(SCRIPT_NAME, cmd_line, args.np, work_dir, output_dir, que_name=args.que_name)
