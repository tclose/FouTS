#!/bin/bash

WORK_DIR='/home/tclose/Documents/Thesis/MCMC/data/simple_tract'
DATA_DIR='/home/tclose/Data/Tractography'
WALK_STEP_SCALE=0.001
NUM_ITERATIONS=200000
SAMPLE_PERIOD=1000
MASK_LOCATION='/home/tclose/Data/Tractography/fibre/tract/masks/mcmc/metropolis/default.tct'
NOISE_SNR=20
NOISE_TYPE='gaussian'
INTERP_TYPE='sinc'
INTERP_EXTENT=3
ASSUMED_INTERP_TYPE='sinc'
ASSUMED_INTERP_EXTENT=1

CURVE_STDDEV=0.025
LENGTH_STDDEV=0.1

DENS_PROB_SCALE=0
TRACT_SHEER_SCALE=0
#ENDPOINTS_SCALE=100
ENDPOINTS_SCALE=100
ENDPOINTS_RADIUS=0.225
#PRIOR_LENGTH_SCALE=500
PRIOR_LENGTH_SCALE=0
NUM_RUNS=10
START_RUN=1

GENERATE_IMAGE=1
INIT_FIBRES=1
METROPOLIS=1
STATS=1

#CONFIGURATIONS='x x-big x-small x-pos-yz xy xz xyz x-curve-z yz-curve-x x_y x_xy1 x_xy x_xy3 x_xxy x-curve-z_y-curve-z layer5 x-rotate x-curve-y_x-curve--y' 
CONFIGURATIONS='x-rotate-big'
#run_index=2

for run_index in `seq $START_RUN $NUM_RUNS`; do

  mkdir -p $WORK_DIR/$run_index
  
  base_img_cmd="generate_image $DATA_DIR/fibre/tract/single/x.tct $WORK_DIR/tmp1.mif -exp_type $INTERP_TYPE -exp_interp_extent $INTERP_EXTENT -clean"
  
  display_and_run "$base_img_cmd"
  
  noise_ref_signal=`maxb0 $WORK_DIR/tmp1.mif`
  
  echo
  echo "Noise reference signal: $noise_ref_signal"
  echo
  
  for config in $CONFIGURATIONS; do
  
    #Check to see if the keyword 'clean' is appended to the configuration name and set the noise parameters accordingly
    if [ `echo $config | grep -c "clean"` = 1 ]; then
      noise=''
    else
      noise=$NOISE
    fi
  
  
    #Either pick up the number of tractlets to use from the configuration name or use the same as in the reference set.
    tractlet_count=`echo $config | grep -o "[0-9]" | tr -d '\n'`
    if [ ${#tractlet_count} = 0 ]; then
      refname=$config
      tractlet_count=`fibre_info $WORK_DIR/ref/$config.tct -prop count`
    else
      refname=${config//[0-9]/}
    fi
  
  
    if [ $GENERATE_IMAGE = 1 ]; then
      #Generate the reference image from the reference set
      img_cmd="generate_image $WORK_DIR/ref/$refname.tct $WORK_DIR/$run_index/$config.mif $noise -exp_type $INTERP_TYPE -exp_interp_extent $INTERP_EXTENT -noise_snr $NOISE_SNR -noise_type $NOISE_TYPE -noise_ref_signal $noise_ref_signal"
  
      display_and_run "$img_cmd"
  
    fi
    
  
    if [ $INIT_FIBRES = 1 ]; then
      init_cmd="init_fibres $WORK_DIR/$run_index/${config}_init.tct -num_fibres $tractlet_count $init_radius -degree 3 -length $LENGTH_STDDEV -curve_stddev $CURVE_STDDEV"
  
      display_and_run "$init_cmd"
  
      #Randomly initialise the tractlets
      #new_fibres $WORK_DIR/tmp1.tct -degree 3 -set_size $tractlet_count
      #perturb_fibres $WORK_DIR/tmp1.tct $WORK_DIR/$run_index/${config}_init.tct -std $INIT_PERTURB_STDDEV -scales $MASK_LOCATION
    fi
  
  
    if [ $METROPOLIS = 1 ]; then
      #Run the metropolis algorithm
      cmd="metropolis $WORK_DIR/$run_index/$config.mif $WORK_DIR/$run_index/${config}_init.tct $WORK_DIR/$run_index/$config.tst -like_snr $NOISE_SNR -exp_type $ASSUMED_INTERP_TYPE -exp_interp_extent $ASSUMED_INTERP_EXTENT -walk_step_scale $WALK_STEP_SCALE -dty_prob_scale $DENS_PROB_SCALE -prior_tract_sheer $TRACT_SHEER_SCALE -prior_end_on_sphere $ENDPOINTS_SCALE $ENDPOINTS_RADIUS -prior_length $PRIOR_LENGTH_SCALE -num_iter $NUM_ITERATIONS -sample_period $SAMPLE_PERIOD"
  
      display_and_run "$cmd"
    fi
  
    if [ $STATS = 1 ]; then
  
      stats_cmd="stats_fibres $WORK_DIR/ref/$refname.tct $WORK_DIR/$run_index/$config.tst"
  
      display_and_run "$stats_cmd"
  
    fi
  
  done

done

#Cleanup temporary files.
rm $WORK_DIR/tmp*
