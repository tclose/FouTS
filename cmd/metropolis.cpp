/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

 This file is part of Bayesian Tractlet Sampling (BTS).

 BTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 BTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with BTS.  If not, see <http://www.gnu.org/licenses/>.

 */


extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}


#include <ctime>

#include "bts/cmd.h"

#include "math/matrix.h"
#include "progressbar.h"


#include "bts/common.h"
#include "bts/file.h"


#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/strand/set.h"


#include "bts/diffusion/model.h"
#include "bts/image/expected/trilinear/buffer.h"
#include "bts/image/expected/gaussian/buffer.h"
#include "bts/image/expected/quartic/buffer.h"

#include "bts/image/observed/buffer.h"

#include "bts/prob/uniform.h"
#include "bts/prob/prior.h"





#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood/one_sided_gaussian.h"
#include "bts/prob/likelihood/gaussian.h"

#include "bts/mcmc/proposal/walker.h"
#include "bts/mcmc/proposal/distribution.h"
#include "bts/mcmc/proposal/distribution/gaussian.h"

#include "bts/mcmc/metropolis.h"

#include "bts/fibre/strand/set/walker.h"
#include "bts/fibre/tractlet/set/walker.h"

#include "bts/math/common.h"


#include "bts/inline_functions.h"

using namespace BTS;



SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Runs Metropolis-Hastings sampling on a given image from a given starting configuration of tractlets or strands.",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("input_image", "The image the tractlets will be fit against.").type_image_in(),

  Argument ("inital_tractlets", "A starting point for the MCMC algorithm.").type_file (),

  Argument ("samples_location", "The location where the samples will be saved.").type_file (),

  Argument()
};


OPTIONS = {


  Option ("num_iterations", "The Number of MCMC iterations to take of the prob distribution over the tractlets")
   + Argument ("num_iterations", "").type_integer (1, MCMC::Metropolis::NUM_ITERATIONS_DEFAULT, LARGE_INT),

  Option ("sample_period", "The number of MCMC iterations that will be performed before a sample is saved")
   + Argument ("sample_period", "").type_integer (1, MCMC::Metropolis::SAMPLE_PERIOD_DEFAULT, LARGE_INT),

  Option ("burn_num_iterations", "The Number of MCMC iterations to take as part of the Burn-in.")
   + Argument ("burn_num_iterations", "").type_integer (0, MCMC::Metropolis::BURN_NUM_ITERATIONS_DEFAULT, LARGE_INT),

  Option ("burn_sample_period", "The number of MCMC iterations that will be performed before a sample is saved in the Burn-in.")
   + Argument ("burn_sample_period", "").type_integer (1, MCMC::Metropolis::BURN_SAMPLE_PERIOD_DEFAULT, LARGE_INT),

  Option ("burn_enforce_bounds", "Enforce the image bounds during the burn-in period. Otherwise the signal outside the image expected, given the fibre configuration will be recorded and its likelihood compared with empty (isotropic) intensities."),

  Option("anneal_frac_start", "Over the burn-in period annealing can be simulated via linear linear cooling specified by a scalar fraction that begins a defined value and linearly increases to 1.0.")
   + Argument("anneal_frac_start","").type_float(0.0, MCMC::Metropolis::ANNEAL_FRAC_START_DEFAULT, 1.0),

  Option ("seed", "The random seed that is passed to the random generator")
   + Argument ("seed", ""),

  Option ("prior_only", "Only sample from the prior distribution."),

  Option ("unverbose", "Turn off verbose output of sampling."),

  Option ("save_images", "Save both observed and expected images for debugging."),

  DIFFUSION_PARAMETERS,

  EXPECTED_IMAGE_PARAMETERS,

  LIKELIHOOD_PARAMETERS,

  PRIOR_PARAMETERS,
  
  PROPOSAL_WALKER_PARAMETERS,

  PROPOSAL_DISTRIBUTION_PARAMETERS,

  COMMON_PARAMETERS,

Option() };




EXECUTE {

//-----------------//
//  Load Arguments //
//-----------------//

  std::string obs_image_location = argument[0];
  std::string initial_location = argument[1];
  std::string samples_location = argument[2];
  std::string burn_samples_location = File::strip_extension (samples_location) + ".burn." + File::extension(samples_location);

  MR::Image::Header header (obs_image_location);

  if (header.ndim() != 4)
    throw Exception ("dwi image should contain 4 dimensions");



//----------------------------------//
//  Get and Set Optional Parameters //
//----------------------------------//

  size_t num_iterations =               MCMC::Metropolis::NUM_ITERATIONS_DEFAULT;
  size_t sample_period =                MCMC::Metropolis::SAMPLE_PERIOD_DEFAULT;
  size_t burn_num_iterations =          MCMC::Metropolis::BURN_NUM_ITERATIONS_DEFAULT;
  size_t burn_sample_period =           MCMC::Metropolis::BURN_SAMPLE_PERIOD_DEFAULT;
  bool burn_enforce_bounds =            false;
  double anneal_frac_start =            MCMC::Metropolis::ANNEAL_FRAC_START_DEFAULT;
  size_t seed =                         time(NULL);
  bool prior_only =                     false;
  bool verbose =                        true;
  bool save_images =                    false;



  Options opt = get_options("num_iterations");
  if (opt.size())
    num_iterations = opt[0][0];

  opt = get_options("sample_period");
  if (opt.size())
    sample_period = opt[0][0];

  opt = get_options("burn_num_iterations");
  if (opt.size())
    burn_num_iterations = opt[0][0];

  opt = get_options("burn_sample_period");
  if (opt.size())
    burn_sample_period = opt[0][0];

  opt = get_options("burn_enforce_bounds");
  if (opt.size())
    burn_enforce_bounds = true;
  
  opt = get_options("anneal_frac_start");
  if (opt.size())
    anneal_frac_start = opt[0][0];

  opt = get_options("seed");
  if (opt.size()) {
    std::string seed_string = opt[0][0];
    seed = to<size_t>(seed_string);
  } else
    std::cout << "No random seed supplied. Using timestamp: " << seed << std::endl;

  opt = get_options("prior_only");
  if (opt.size())
    prior_only = true;

  opt = get_options("unverbose");
  if (opt.size())
    verbose = false;

  opt = get_options("save_images");
  if (opt.size())
    save_images = true;

  // Loads parameters to construct Diffusion::Model ('diff_' prefix)
  SET_DIFFUSION_PARAMETERS;

  // Loads parameters to construct Image::Expected::*::Buffer ('img_' prefix)
  SET_EXPECTED_IMAGE_PARAMETERS;

  // Loads parameters to construct Prob::Likelihood ('like_' prefix)
  SET_LIKELIHOOD_PARAMETERS;
  
  // Loads parameters to construct Prob::Prior ('prior_' prefix)  
  SET_PRIOR_PARAMETERS;

  // Loads parameters to construct Proposal::Distribution ('walk_' prefix)
  SET_PROPOSAL_WALKER_PARAMETERS(initial_location);

  // Loads parameters to construct Proposal::Distribution ('walk_' prefix)
  SET_PROPOSAL_DISTRIBUTION_PARAMETERS;

  // Loads parameters that are common to all commands.
  SET_COMMON_PARAMETERS;


  //--------------------------------//
  //  Set up reference image buffer //
  //--------------------------------//

  Image::Observed::Buffer obs_image(obs_image_location, Diffusion::Encoding::Set(diff_encodings));

  //If gradient scheme is included in reference image header, use that instead of default (NB: Will override any gradients passed to '-diff_encodings' option).
  if (header.get_DW_scheme().rows()) {
    diff_encodings = header.get_DW_scheme();
    diff_encodings_location = "From observed image";
  }


  if ((header.count("diff_response_SH")) && (Math::matlab_str(diff_response_SH) != header["diff_response_SH"]))
    std::cout << std::endl << "Warning! Diffusion response function harmonics (" << Math::matlab_str(diff_response_SH)  << ") do not match reference image (" << header["diff_response_SH"] + ")!" << std::endl;


  //----------------------------//
  //  Initialize Expected Image //
  //----------------------------//



  Diffusion::Model diffusion_model = Diffusion::Model::factory (diff_encodings,
                                                                diff_response_SH,
                                                                diff_adc,
                                                                diff_fa,
                                                                diff_isotropic);

  Image::Expected::Buffer* exp_image = Image::Expected::Buffer::factory(exp_type,
                                                        obs_image,
                                                        diffusion_model,
                                                        exp_num_length_sections,
                                                        exp_num_width_sections,
                                                        exp_interp_extent,
                                                        exp_enforce_bounds,
                                                        exp_half_width);

  //----------------------------------------------------------------//
  // Auto-generate base acs initial value/b_intens_gauss_mean //
  //----------------------------------------------------------------//

  if (exp_base_intensity < 0)
    exp_base_intensity = exp_image->base_intensity_default(obs_image, initial_location);


  //-----------------------//
  // Initialize Likelihood //
  //-----------------------//

  Prob::Likelihood* likelihood = Prob::Likelihood::factory(like_type, obs_image, exp_image, like_snr, like_b0_include, like_outside_scale, like_ref_b0, like_ref_signal);



  //----------------------//
  // Initialize Proposals //
  //----------------------//

  gsl_rng* rand_gen = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rand_gen, seed);


  MCMC::Proposal::Distribution *proposal_distribution = MCMC::Proposal::Distribution::factory(prop_distr_type, rand_gen);

//
//  Proposal::Intensity* acs_proposer = Proposal::Intensity::factory(walk_intens_type, proposal_distribution, walk_intens_scale);
//
//  Proposal::BaseIntensity* base_intensity_proposer = Proposal::BaseIntensity::factory(walk_base_intens_type, proposal_distribution, walk_base_intens_scale);



  //-------------------------//
  //  Set Output Properties  //
  //-------------------------//

  std::map<std::string,std::string> run_properties;

  run_properties["Method"]                = "metropolis";
  run_properties["sample_period"]         = str(sample_period);
  run_properties["seed"]                  = str(seed);
  run_properties["burn_sample_period"]    = str(burn_sample_period);
  run_properties["burn_num_iterations"]   = str(burn_num_iterations);
  run_properties["obs_image"]             = obs_image_location;
  run_properties["initial_state"]         = Fibre::Base::Object::load_matlab_str(initial_location);
  run_properties["initial_state_location"]= initial_location;
  run_properties["prior_only"]            = str(prior_only);
  run_properties["num_iterations"]        = str(num_iterations);
  run_properties["burn_enforce_bounds"]   = str(burn_enforce_bounds);
  run_properties["anneal_frac_start"]     = str(anneal_frac_start);


  ADD_DIFFUSION_PROPERTIES(run_properties);

  ADD_LIKELIHOOD_PROPERTIES(run_properties);

  ADD_EXPECTED_IMAGE_PROPERTIES(run_properties);

  ADD_PRIOR_PROPERTIES(run_properties);

  ADD_PROPOSAL_WALKER_PROPERTIES(run_properties);

  ADD_PROPOSAL_DISTRIBUTION_PROPERTIES(run_properties);

  ADD_COMMON_PROPERTIES(run_properties);


  Prob::Prior prior (prior_scale,
                    prior_freq_scale,
                    prior_freq_aux_scale,
                    prior_hook_scale,
                    prior_hook_num_points,
                    prior_density_high_scale,
                    prior_density_low_scale,
                    prior_density_num_points,
                    prior_acs_scale,
                    prior_acs_mean,
                    prior_length_scale,
                    prior_length_mean);


  //-------------------------//
  //  Sampling from Strands  //
  //-------------------------//


  if (File::has_extension<Fibre::Strand>(initial_location)) {


    //------------------------//
    //  Load Initial Strands  //
    //------------------------//

    Fibre::Strand::Set strands (initial_location);

    strands.set_base_intensity(exp_base_intensity);

    //---------------------//
    // Initialize Proposer //
    //---------------------//

    Fibre::Strand::Set::Walker* walker = Fibre::Strand::Set::Walker::factory(strands,
                                                                            walk_type,
                                                                            walk_step_scale,
                                                                            walk_step_location,
                                                                            proposal_distribution,
                                                                            walk_base_intens_scale);

    //------------------//
    // Perform sampling //
    //------------------//


    Fibre::Strand::Set burnt_strands;

    if (burn_num_iterations) {


      if (burn_enforce_bounds != exp_enforce_bounds)
        likelihood->set_enforce_bounds(burn_enforce_bounds);

      burnt_strands = MCMC::metropolis <Fibre::Strand::Set,
                                        Prob::Likelihood,
                                        Prob::Prior>

                                      ( strands,
                                        *likelihood,
                                        prior,
                                        *walker,
                                        burn_samples_location,
                                        run_properties,
                                        burn_num_iterations,
                                        burn_sample_period,
                                        rand_gen,
                                        anneal_frac_start,
                                        prior_only,
                                        verbose,
                                        save_images);


      if (burn_enforce_bounds != exp_enforce_bounds)
        likelihood->set_enforce_bounds(exp_enforce_bounds);


    } else
      burnt_strands = strands;


    MCMC::metropolis <Fibre::Strand::Set,
                      Prob::Likelihood,
                      Prob::Prior>

                     (burnt_strands,
                      *likelihood,
                      prior,
                      *walker,
                      samples_location,
                      run_properties,
                      num_iterations,
                      sample_period,
                      rand_gen,
                      1.0,
                      prior_only,
                      verbose,
                      save_images);



  //------------------------//
  //  Sampling from Tractlets  //
  //------------------------//

  } else if (File::has_extension<Fibre::Tractlet>(initial_location)) {


    //-------------//
    // Load Tractlets //
    //-------------//

    Fibre::Tractlet::Set tractlets (initial_location);

    tractlets.set_base_intensity(exp_base_intensity);


    //---------------------//
    // Initialize Proposer //
    //---------------------//

    Fibre::Tractlet::Set::Walker* walker = Fibre::Tractlet::Set::Walker::factory ( tractlets,
                                                                            walk_type,
                                                                            walk_step_scale,
                                                                            walk_step_location,
                                                                            proposal_distribution,
                                                                            walk_base_intens_scale);

    
    if (walk_type == "manifold") {

      Fibre::Tractlet::Set dummy (tractlets);

      //Ensures that tractlets lie on the manifold to begin with
      walker->step(tractlets, dummy);

      tractlets = dummy;

    }


    
    //------------------//
    // Perform sampling //
    //------------------//

    Fibre::Tractlet::Set burnt_tractlets;

    if (burn_num_iterations) {

      if (burn_enforce_bounds != exp_enforce_bounds)
        likelihood->set_enforce_bounds(burn_enforce_bounds);

      burnt_tractlets = MCMC::metropolis < Fibre::Tractlet::Set,
                                        Prob::Likelihood,
                                        Prob::Prior> (

                                  tractlets,
                                  *likelihood,
                                  prior,
                                  *walker,
                                  burn_samples_location,
                                  run_properties,
                                  burn_num_iterations,
                                  burn_sample_period,
                                  rand_gen,
                                  anneal_frac_start,
                                  prior_only,
                                  verbose,
                                  save_images
                                );

      if (burn_enforce_bounds != exp_enforce_bounds)
        likelihood->set_enforce_bounds(exp_enforce_bounds);


    } else
      burnt_tractlets = tractlets;

    MCMC::metropolis <Fibre::Tractlet::Set,
                      Prob::Likelihood,
                      Prob::Prior> (

                burnt_tractlets, 
                *likelihood,
                prior,
                *walker,
                samples_location,
                run_properties,
                num_iterations,
                sample_period,
                rand_gen,
                1.0,
                prior_only,
                verbose,
                save_images
              );

  } 


  delete exp_image;
  delete proposal_distribution;
  delete likelihood;

  gsl_rng_free(rand_gen);



}


