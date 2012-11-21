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

#include "bts/fibre/strand/set/momentum.h"
#include "bts/fibre/tractlet/set/momentum.h"
#include "bts/mcmc/proposal/distribution.h"
#include "bts/mcmc/proposal/distribution/gaussian.h"
#include "bts/mcmc/proposal/momentum/weighted.h"

#include "bts/math/common.h"

#include "bts/mcmc/hamiltonian.h"

#include "bts/file.h"



#include "bts/inline_functions.h"

using namespace BTS;


SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Runs Hamiltonian MCMC sampling on a given image from a given starting configuration of tractlets or strands.",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("input_image", "The image the tractlets will be fit against.").type_image_in(),

  Argument ("inital tractlets", "A starting point for the MCMC algorithm.").type_file (),

  Argument ("samples location", "The location where the samples will be saved.").type_file (),

  Argument()
};


OPTIONS = {


  Option ("num_samples", "The Number of MCMC iterations to take of the prob distribution over the tractlets")
   + Argument ("num_samples", "").type_integer (1, MCMC::Hamiltonian::NUM_SAMPLES_DEFAULT, LARGE_INT),

  Option ("num_leapfrog_steps", "The number of MCMC iterations that will be performed before a sample is saved")
   + Argument ("num_leapfrog_steps", "").type_integer (1, MCMC::Hamiltonian::NUM_LEAPFROG_STEPS_DEFAULT, LARGE_INT),

  Option ("burn_num_samples", "The Number of MCMC iterations to take as part of the Burn-in.")
   + Argument ("burn_num_samples", "").type_integer (0, MCMC::Hamiltonian::BURN_NUM_SAMPLES_DEFAULT, LARGE_INT),

  Option ("burn_num_leapfrog_steps", "The number of MCMC iterations that will be performed before a sample is saved in the Burn-in.")
   + Argument ("burn_num_leapfrog_steps", "").type_integer (1, MCMC::Hamiltonian::BURN_NUM_LEAPFROG_STEPS_DEFAULT, LARGE_INT),

  Option ("burn_snr", "The Signal-to-noise ratio assumed for the burn-in period.")
   + Argument ("burn_snr", "").type_float (SMALL_FLOAT, MCMC::Hamiltonian::BURN_SNR_DEFAULT, LARGE_FLOAT),

  Option ("seed", "The random seed that is passed to the random generator")
   + Argument ("seed", ""),

  Option ("prior_only", "Only sample from the prior distribution."),

  Option ("save_iterations", "Save all iterations not just samples."),

  DIFFUSION_PARAMETERS,

  EXPECTED_IMAGE_PARAMETERS,

  LIKELIHOOD_PARAMETERS,

  PRIOR_PARAMETERS,
  
  PROPOSAL_DISTRIBUTION_PARAMETERS,

  PROPOSAL_MOMENTUM_PARAMETERS,

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

  MR::Image::Header header (argument[0]);

  if (header.ndim() != 4)
    throw Exception ("dwi image should contain 4 dimensions");



//----------------------------------//
//  Get and Set Optional Parameters //
//----------------------------------//


  size_t num_samples = MCMC::Hamiltonian::NUM_SAMPLES_DEFAULT;
  size_t num_leapfrog_steps = MCMC::Hamiltonian::NUM_LEAPFROG_STEPS_DEFAULT;
  size_t burn_num_samples = MCMC::Hamiltonian::BURN_NUM_SAMPLES_DEFAULT;
  size_t burn_num_leapfrog_steps = MCMC::Hamiltonian::BURN_NUM_LEAPFROG_STEPS_DEFAULT;
  double burn_snr = MCMC::Hamiltonian::BURN_SNR_DEFAULT;
  size_t seed = time(NULL);
  bool prior_only = false;
  bool save_iterations = true;


  Options opt = get_options("num_samples");
  if (opt.size())
    num_samples = opt[0][0];

  opt = get_options("num_leapfrog_steps");
  if (opt.size())
    num_leapfrog_steps = opt[0][0];

  opt = get_options("burn_num_samples");
  if (opt.size())
    burn_num_samples = opt[0][0];

  opt = get_options("burn_num_leapfrog_steps");
  if (opt.size())
    burn_num_leapfrog_steps = opt[0][0];
    
  opt = get_options("burn_snr");
  if (opt.size())
    burn_snr = opt[0][0];

  opt = get_options("seed");
  if (opt.size()) {
    std::string seed_string = opt[0][0];
    seed = to<size_t>(seed_string);
  } else
    std::cout << "No random seed supplied. Using timestamp: " << seed << std::endl;

  opt = get_options("prior_only");
  if (opt.size())
    prior_only = true;

  opt = get_options("save_iterations");
  if (opt.size())
    save_iterations = true;

  // Loads parameters to construct Diffusion::Model ('diff_' prefix)
  SET_DIFFUSION_PARAMETERS;

  // Loads parameters to construct Image::Expected::*::Buffer ('img_' prefix)
  SET_EXPECTED_IMAGE_PARAMETERS;
  
  // Loads parameters to construct Prob::Likelihood ('like_' prefix)
  SET_LIKELIHOOD_PARAMETERS;
  
  // Loads parameters to construct Prob::Prior ('prior_' prefix)  
  SET_PRIOR_PARAMETERS;

  // Loads parameters to construct Proposal::Distribution (any variables with 'prop_' prefix)
  SET_PROPOSAL_DISTRIBUTION_PARAMETERS;

  // Loads parameters to construct Proposal::Distribution ('prop_' prefix)
  SET_PROPOSAL_MOMENTUM_PARAMETERS(initial_location);

  // Loads parameters that are common to all commands.
  SET_COMMON_PARAMETERS;


  //--------------------------------//
  //  Set up reference image buffer //
  //--------------------------------//


  Image::Observed::Buffer obs_image(obs_image_location, Diffusion::Encoding::Set(diff_encodings));


  //If gradient scheme is included in reference image header, use that instead of default (NB: Will override any gradients passed to '-diff_encodings' option).
  if (header.get_DW_scheme().rows()) {
    diff_encodings = header.get_DW_scheme();
    diff_encodings_location = "from observed image";
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
                                                                diff_isotropic,
                                                                diff_warn_b_mismatch);

  Image::Expected::Buffer* exp_image = Image::Expected::Buffer::factory(exp_type,
                                                        obs_image,
                                                        diffusion_model,
                                                        exp_num_length_sections,
                                                        exp_num_width_sections,
                                                        exp_interp_extent,
                                                        exp_enforce_bounds,
                                                        exp_half_width);




  //----------------------------------------------------------------//
  // Auto-generate base intensity initial value/b_intens_gauss_mean //
  //----------------------------------------------------------------//

  double base_intensity = exp_image->get_base_intensity(exp_b0);


  //-----------------------//
  // Initialize Likelihood //
  //-----------------------//

  Prob::Likelihood* likelihood    = Prob::Likelihood::factory(like_type, obs_image, exp_image, like_snr, like_b0_include, like_outside_scale, like_ref_b0, like_ref_signal);



  //----------------------------------//
  // Initialize Proposal Distribution //
  //----------------------------------//

  gsl_rng* rand_gen = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rand_gen, seed);

  MCMC::Proposal::Distribution* proposal_distribution = MCMC::Proposal::Distribution::factory(prop_distr_type, rand_gen);



//  Proposal::Intensity* intensity_proposer = Proposal::Intensity::factory(prop_intens_type, proposal_distribution, prop_intens_scale);
//
//  Proposal::BaseIntensity* base_intensity_proposer = Proposal::BaseIntensity::factory(prop_base_intens_type, proposal_distribution, prop_base_intens_scale);



  //-------------------------//
  //  Set Output Properties  //
  //-------------------------//

  std::map<std::string,std::string> run_properties;

  run_properties["Method"]                = "hamiltonian";
  run_properties["num_leapfrog_steps"]         = str(num_leapfrog_steps);
  run_properties["seed"]           = str(seed);
  run_properties["num_samples"]           = str(num_samples);
  run_properties["burn_num_leapfrog_steps"]    = str(burn_num_leapfrog_steps);
  run_properties["burn_num_samples"]      = str(burn_num_samples);
  run_properties["burn_snr"]              = str(burn_snr);
  run_properties["obs_image"]             = obs_image_location;
  run_properties["initial_state"]         = Fibre::Base::Object::load_matlab_str(initial_location);
  run_properties["initial_state_location"]= initial_location;

  run_properties["prior_only"]            = str(prior_only);

  ADD_DIFFUSION_PROPERTIES(run_properties);

  ADD_LIKELIHOOD_PROPERTIES(run_properties);

  ADD_EXPECTED_IMAGE_PROPERTIES(run_properties);

  ADD_PRIOR_PROPERTIES(run_properties);

  ADD_PROPOSAL_DISTRIBUTION_PROPERTIES(run_properties);

  ADD_PROPOSAL_MOMENTUM_PROPERTIES(run_properties);

  ADD_COMMON_PROPERTIES(run_properties);


  //-------------------------//
  //  Sampling from Strands  //
  //-------------------------//

  //------------------//
  // Initialize Prior //
  //------------------//

  Prob::Prior prior (prior_scale,
                    prior_freq_scale,
                    prior_freq_aux_scale,
                    prior_hook_scale,
                    prior_hook_num_points,
                    prior_hook_num_width_sections,
                    prior_density_high_scale,
                    prior_density_low_scale,
                    prior_density_num_points,
                    prior_acs_scale,
                    prior_acs_mean,
                    prior_length_scale,
                    prior_length_mean,
                    prior_thinness_scale,
                    prior_thinness_power);


  if (File::has_extension<Fibre::Strand>(initial_location)) {


    //------------------------//
    //  Load Initial Strands  //
    //------------------------//

    Fibre::Strand::Set strands (initial_location);

    if (base_intensity)
      strands.set_base_intensity(base_intensity);


    //---------------------//
    // Initialize Momentum //
    //---------------------//

    MCMC::Proposal::Momentum momentum = MCMC::Proposal::Momentum::factory ( strands,
                                                                            momen_step_scale,
                                                                            momen_step_location,
                                                                            momen_weight_location,
                                                                            proposal_distribution);


    //------------------//
    // Perform sampling //
    //------------------//

    Fibre::Strand::Set burnt_strands;



    if (burn_num_samples) {

      likelihood->set_assumed_snr(burn_snr, like_ref_b0, like_ref_signal);

      burnt_strands = MCMC::hamiltonian < Fibre::Strand::Set,
                                          Prob::Likelihood,
                                          Prob::Prior> (
                                            strands,
                                            *likelihood,
                                            prior,
                                            momentum,
                                            burn_samples_location,
                                            run_properties,
                                            burn_num_samples,
                                            burn_num_leapfrog_steps,
                                            rand_gen,
                                            prior_only,
                                            save_iterations
                                          );

      likelihood->set_assumed_snr(like_snr, like_ref_b0, like_ref_signal);

    } else
      burnt_strands = strands;


    MCMC::hamiltonian < Fibre::Strand::Set,
                        Prob::Likelihood,
                        Prob::Prior> (
      burnt_strands,
      *likelihood,
      prior,
      momentum,
      samples_location,
      run_properties,
      num_samples,
      num_leapfrog_steps,
      rand_gen,
      prior_only,
      save_iterations
    );



  //------------------------//
  //  Sampling from Tractlets  //
  //------------------------//


  } else if (File::has_extension<Fibre::Tractlet>(initial_location)) {


    //-------------//
    // Load Tractlets //
    //-------------//

    Fibre::Tractlet::Set tractlets (initial_location);

    if (base_intensity)
      tractlets.set_base_intensity(base_intensity);


    //---------------------//
    // Initialize Proposer //
    //---------------------//


    MCMC::Proposal::Momentum momentum = MCMC::Proposal::Momentum::factory ( tractlets,
                                                                              momen_step_scale,
                                                                              momen_step_location,
                                                                              momen_weight_location,
                                                                              proposal_distribution);



    
    //------------------//
    // Perform sampling //
    //------------------//

    Fibre::Tractlet::Set burnt_tractlets;
    
    if (burn_num_samples) {

      likelihood->set_assumed_snr(burn_snr, like_ref_b0, like_ref_signal);

      burnt_tractlets = MCMC::hamiltonian <Fibre::Tractlet::Set,
                                        Prob::Likelihood,
                                        Prob::Prior> (

                                      tractlets,
                                      *likelihood,
                                      prior,
                                      momentum,
                                      burn_samples_location,
                                      run_properties,
                                      burn_num_samples,
                                      burn_num_leapfrog_steps,
                                      rand_gen,
                                      prior_only,
                                      save_iterations
                                    );

      likelihood->set_assumed_snr(like_snr, like_ref_b0, like_ref_signal);

    } else
      burnt_tractlets = tractlets;



    MCMC::hamiltonian <Fibre::Tractlet::Set,
                       Prob::Likelihood,
                       Prob::Prior> (

                burnt_tractlets, 
                *likelihood,
                prior,
                momentum,
                samples_location,
                run_properties,
                num_samples,
                num_leapfrog_steps,
                rand_gen,
                prior_only,
                save_iterations
              );


  } else
    throw Exception ("Unrecognised extension of initial state location, '" + initial_location + "'.");


  std::cout << "MATLAB plot commands:\n\nplot_fibres " << samples_location << "\nscan.py " << samples_location << std::endl;


  gsl_rng_free(rand_gen);

  delete exp_image;
  delete proposal_distribution;
  delete likelihood;


}













