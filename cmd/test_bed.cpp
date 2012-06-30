/*
    Copyright 2009 Brain Research Institute, Melbourne, Australia

    Written by Thomas G. Close, 04/03/2009.

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


#include "bts/cmd.h"

#include "bts/common.h"

#include "bts/math/common.h"

#include "bts/prob/test/gaussian.h"
#include "bts/prob/test/landscape.h"
#include "bts/prob/test/landscape/peak.h"
#include "bts/prob/test/bayes_log_regression.h"
#include "bts/prob/uniform.h"
#include "bts/mcmc/proposal/distribution.h"
#include "bts/mcmc/proposal/distribution/gaussian.h"
#include "bts/mcmc/proposal/walker.h"
#include "bts/mcmc/proposal/momentum.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"

#include "bts/mcmc/metropolis.h"
#include "bts/mcmc/hamiltonian.h"
#include "bts/mcmc/riemannian.h"

#include "bts/triple.h"


#include "bts/mcmc/state.h"
//#include "bts/mcmc/state/writer.h"
#include "bts/mcmc/state/tensor/writer.h"
#include "bts/utilities/reader.h"

#include "bts/fibre/strand/set/momentum.h"
#include "bts/fibre/tractlet/set/momentum.h"

#include "bts/inline_functions.h"

#include "bts/diffusion/response.h"



using namespace BTS;

const char* STEP_LOCATION_DEFAULT = "/home/tclose/Code/Tractography/params/test_bed/scaling.tes";

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "dummy program",
  "",
  NULL
};

ARGUMENTS = {
    Argument ("initial location", "The location of the initial state.").type_file (),
    Argument ("samples location", "The location where the samples will be saved.").type_file (),
  Argument()
};


const char* ALGORITHM_TYPE_DEFAULT = "metropolis";
const char* PROB_TYPE_DEFAULT = "Prob::Test::Gaussian";
const double STEP_SCALE_DEFAULT = 0.1;
const size_t NUM_ITERATIONS_DEFAULT = 10;
const size_t SAMPLE_PERIOD_DEFAULT = 1e2;
const size_t NUM_LEAPFROG_STEPS_DEFAULT = 30;
const size_t NUM_NEWTON_STEPS_DEFAULT = 6;
const char* WALKER_STEP_LOCATION_DEFAULT = "/home/tclose/data/mcmc/metropolis/params/walker/scaling.sta";
const char* MOMENTUM_STEP_LOCATION_DEFAULT = "/home/tclose/data/mcmc/hamiltonian/params/momentum/scaling.sta";
const char* AXIS_SCALES_DEFAULT_LOCATION_DEFAULT = "/home/tclose/data/prob/test/gaussian/axis_scales.sta";

OPTIONS = {

  Option ("algorithm", "Select the algorithm to use (either 'metropolis', 'hamiltonian' at this stage).")
   + Argument ("algorithm", "").type_text(),

  Option ("prob_type", "Select the probability object (either 'Prob::Test::Gaussian' or 'Prob::Test::Peaks' at this stage.")
   + Argument ("prob_type", ""),

  Option ("num_samples", "The Number of MCMC iterations to take of the prob distribution over the tractlets")
   + Argument ("num_samples", "").type_integer (1, NUM_ITERATIONS_DEFAULT, LARGE_INT),

  Option ("sample_period", "The number of MCMC iterations that will be performed before a sample is saved")
   + Argument ("sample_period", "").type_integer (1, SAMPLE_PERIOD_DEFAULT, LARGE_INT),

  Option ("num_leapfrog_steps", "The number of leapfrog steps to be performed in the 'hamiltonian' and 'riemannian' algorithms.")
   + Argument ("num_leapfrog_steps", "").type_integer (1, NUM_LEAPFROG_STEPS_DEFAULT, LARGE_INT),

  Option ("num_newton_steps", "The number of Newton steps that will be performed in calculating the non-separable Hamiltonian in the 'riemannian' algorithm.")
   + Argument ("num_newton_steps", "").type_integer (1, NUM_NEWTON_STEPS_DEFAULT, LARGE_INT),

  Option ("seed", "The random seed that is passed to the random generator")
   + Argument ("seed", ""),
  
  Option ("axis_scales_location", "Location of the file containing the scales of the Gaussian axes")
   + Argument ("axis_scales_location", ""),

  Option ("blr_location", "Location of the Bayesian log regression data.")
   + Argument ("blr_location", "").type_text (Prob::Test::BayesLogRegression::DATA_LOCATION_DEFAULT),

  Option ("blr_poly_order", "Order of the polynomials used in the Bayesian Log Regression.")
   + Argument ("blr_poly_order", "").type_integer (1, Prob::Test::BayesLogRegression::POLY_ORDER_DEFAULT, LARGE_INT),

  Option ("blr_prior_variance", "Variance of the prior used in the Bayesian Log Regression")
   + Argument ("blr_prior_variance", "").type_float (SMALL_FLOAT, Prob::Test::BayesLogRegression::PRIOR_VARIANCE_DEFAULT, LARGE_FLOAT),

  Option ("step_scale", "Overall step scale.")
   + Argument ("step_scale", "").type_float (SMALL_FLOAT, STEP_SCALE_DEFAULT, LARGE_FLOAT),

  Option ("step_location", "Name of the file containing the relative scales of the step sizes.")
   + Argument ("step_location", ""),
  
  Option ("save_iterations", "Save every iteration in seperate file (with '.iter' inserted before the extension)."),
  
  Option ("suppress_stdout", "Don't print probability to standard out."),

  PROPOSAL_DISTRIBUTION_PARAMETERS,

  PROPOSAL_MOMENTUM_PARAMETERS,

  COMMON_PARAMETERS,

Option() };



EXECUTE {

  std::string initial_location = argument[0];
  std::string samples_location = argument[1];

  Options opt;
  
  std::string algorithm_type = ALGORITHM_TYPE_DEFAULT;
  std::string prob_type = PROB_TYPE_DEFAULT;
  size_t num_samples = NUM_ITERATIONS_DEFAULT;
  size_t sample_period = SAMPLE_PERIOD_DEFAULT;
  size_t num_leapfrog_steps = NUM_LEAPFROG_STEPS_DEFAULT;
  size_t num_newton_steps = NUM_NEWTON_STEPS_DEFAULT;
  size_t seed = time(NULL);
  std::string axis_scales_location = AXIS_SCALES_DEFAULT_LOCATION_DEFAULT;
  std::string blr_location = Prob::Test::BayesLogRegression::DATA_LOCATION_DEFAULT;
  size_t blr_poly_order = Prob::Test::BayesLogRegression::POLY_ORDER_DEFAULT;
  float blr_prior_variance = Prob::Test::BayesLogRegression::PRIOR_VARIANCE_DEFAULT;
  double step_scale = STEP_SCALE_DEFAULT;
  std::string step_location;

  bool save_iterations = true;  
  bool suppress_print = false;  

  
  opt = get_options("algorithm");
  if (opt.size())
    algorithm_type = opt[0][0].c_str(); //.get_string()

  opt = get_options("prob_type");
  if (opt.size())
    prob_type = opt[0][0].c_str();
  
  opt = get_options("num_samples");
  if (opt.size())
    num_samples = opt[0][0];

  opt = get_options("sample_period");
  if (opt.size())
    sample_period = opt[0][0];

  opt = get_options("num_leapfrog_steps");
  if (opt.size())
    num_leapfrog_steps = opt[0][0];

  opt = get_options("num_newton_steps");
  if (opt.size())
    num_newton_steps = opt[0][0];

  opt = get_options("seed");
  if (opt.size()) {
    std::string seed_string = opt[0][0];
    seed = to<size_t>(seed_string);
  } else
    std::cout << "No random seed supplied. Using timestamp: " << seed << std::endl;

  opt = get_options("axis_scales_location");
  if (opt.size())
    axis_scales_location = opt[0][0].c_str();


  opt = get_options("blr_location");
  if (opt.size())
    blr_location = opt[0][0].c_str();

  opt = get_options("blr_poly_order");
  if (opt.size())
    blr_poly_order = opt[0][0];

  opt = get_options("blr_prior_variance");
  if (opt.size())
    blr_prior_variance = opt[0][0];

  opt = get_options("step_scale");
  if (opt.size())
    step_scale = opt[0][0];

  opt = get_options("save_iterations");
  if (opt.size())
    save_iterations = true;
    
  opt = get_options("suppress_print");
  if (opt.size())
    suppress_print = true;

  if (algorithm_type == "metropolis")
    step_location = WALKER_STEP_LOCATION_DEFAULT;

  else
    step_location = MOMENTUM_STEP_LOCATION_DEFAULT;


  opt = get_options("step_location");
  if (opt.size())
    step_location = opt[0][0].c_str();


  // Loads parameters to construct Proposal::Distribution (any variables with 'prop_' prefix)
  SET_PROPOSAL_DISTRIBUTION_PARAMETERS;

  // Loads parameters to construct Proposal::Distribution ('prop_' prefix)
  SET_PROPOSAL_MOMENTUM_PARAMETERS(initial_location);

  // Loads parameters to construct Prob::Test::Peaks (any variables with 'test_gen_' prefix)
  SET_TEST_LANDSCAPE_PARAMETERS;

  // Loads parameters that are common to all commands.
  SET_COMMON_PARAMETERS;


  //-------------------------------------//
  //  Initialize random number generator //
  //-------------------------------------//


  gsl_rng* rand_gen = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rand_gen, seed);

 
  //---------------------------//
  // Initialise the test state //
  //---------------------------//

  MCMC::State state(initial_location);

  //-----------------------------------//
  // Initialize Proposal::Distribution //
  //-----------------------------------//

  MCMC::Proposal::Distribution *proposal_distribution = MCMC::Proposal::Distribution::factory(prop_distr_type, rand_gen);

  MCMC::State step (step_location);

  step *= step_scale;

  MCMC::Proposal::Walker walker (proposal_distribution, step);
  MCMC::Proposal::Momentum momentum (proposal_distribution, step);

  //----------------------------------------------//
  //  Initialize SampleWriter and its Properties  //
  //----------------------------------------------//


  // Set the 'run-wide' properties.
  std::map<std::string,std::string> run_properties;

  run_properties["algorithm"]                 = algorithm_type;
  run_properties["prob_type"]                 = prob_type;
  run_properties["sample_period"]             = str(sample_period);
  run_properties["num_leapfrog_steps"]        = str(num_leapfrog_steps);
  run_properties["num_newton_steps"]          = str(num_newton_steps);
  run_properties["num_samples"]               = str(num_samples);
  run_properties["seed"]                      = str(seed);
  run_properties["step"]                      = Math::matlab_str(step);
  run_properties["step_scale"]                = str(step_scale);
  run_properties["step_location"]             = step_location;


  ADD_PROPOSAL_DISTRIBUTION_PROPERTIES(run_properties);

  ADD_TEST_LANDSCAPE_PROPERTIES(run_properties);
  
  ADD_COMMON_PROPERTIES(run_properties);

  //-------------------//
  // Run the algorithm //
  //-------------------//

  if (prob_type == "Prob::Test::Landscape") {

    Prob::Test::Landscape test_peaks(lnd_location);

    if (algorithm_type == "metropolis") {

//        MCMC::metropolis <MCMC::State,Prob::Likelihood,Prob::Prior> (
//                    state,
//                    test_peaks,
//                    test_peaks,
//                    walker,
//                    samples_location,
//                    run_properties,
//                    num_samples * sample_period,
//                    sample_period,
//                    rand_gen,
//                    true,
//                    suppress_print
//                  );

      throw Exception ("removed for debugging.");

      } else if (!strcmp(algorithm_type.c_str(), "hamiltonian")) {

        MCMC::hamiltonian <MCMC::State,Prob::Test::Landscape,Prob::Test::Landscape> (
                  state,
                  test_peaks,
                  test_peaks,
                  momentum,
                  samples_location,
                  run_properties,
                  num_samples,
                  num_leapfrog_steps,
                  rand_gen,
                  true,
                  save_iterations,
                  suppress_print
                );

      } else
        throw Exception ("Unrecognised value for '-algorithm' option \"" + algorithm_type + "\", can be either \"metropolis\" or \"hamiltonian\".");

      if (!lnd_location.size())
        test_peaks.save(samples_location + ".lnd");

  } else if (prob_type == "Prob::Test::Gaussian") {

    run_properties["axis_scales_location"]            = str(axis_scales_location);

    MCMC::State axis_scale (axis_scales_location);
    Prob::Test::Gaussian test_gaussian(axis_scale);

    if (algorithm_type == "metropolis") {
      
//      MCMC::metropolis <MCMC::State, Prob::Likelihood, Prob::Prior> (
//                  state,
//                  test_gaussian,
//                  test_gaussian,
//                  walker,
//                  samples_location,
//                  run_properties,
//                  num_samples * sample_period,
//                  sample_period,
//                  rand_gen,
//                  true,
//                  suppress_print
//                );

      throw Exception("");
                
      } else if (!strcmp(algorithm_type.c_str(), "hamiltonian")) {                  

        MCMC::hamiltonian <MCMC::State, Prob::Test::Gaussian, Prob::Test::Gaussian> (
                  state, 
                  test_gaussian,
                  test_gaussian,
                  momentum,
                  samples_location,
                  run_properties,
                  num_samples,
                  num_leapfrog_steps,
                  rand_gen,
                  true,
                  save_iterations,                
                  suppress_print
                );
                
      } else if (!strcmp(algorithm_type.c_str(), "riemannian")) {


        MCMC::Proposal::Momentum::Weighted::NonSeparable nonseparable_momentum  = MCMC::Proposal::Momentum::Weighted::NonSeparable::factory( state,
                                                                                                                                momen_step_scale,
                                                                                                                                momen_step_location,
                                                                                                                                proposal_distribution,
                                                                                                                                num_newton_steps);

        Prob::Uniform uniform;

        MCMC::riemannian <MCMC::State, Prob::Test::Gaussian, Prob::Uniform> (
                  state,
                  test_gaussian,
                  uniform,
                  nonseparable_momentum,
                  samples_location,
                  run_properties,
                  num_samples,
                  num_leapfrog_steps,
                  num_newton_steps,
                  rand_gen,
                  true,
                  save_iterations,
                  suppress_print
                );

      } else
        throw Exception ("Unrecognised value for '-algorithm' option \"" + algorithm_type + "\", can be either 'metropolis', 'hamiltonian' or 'riemannian'.");

  } else if (prob_type == "Prob::Test::BayesLogRegression") {

    Prob::Test::BayesLogRegression bayes_log_regression (blr_location, blr_poly_order, blr_prior_variance);

    if (!strcmp(algorithm_type.c_str(), "riemannian")) {


      MCMC::Proposal::Momentum::Weighted::NonSeparable nonseparable_momentum  = MCMC::Proposal::Momentum::Weighted::NonSeparable::factory( state,
                                                                                                                             momen_step_scale,
                                                                                                                             momen_step_location,
                                                                                                                             proposal_distribution,
                                                                                                                             num_newton_steps);
      Prob::Uniform uniform;

      MCMC::riemannian <MCMC::State, Prob::Test::BayesLogRegression, Prob::Uniform> (
               state,
               bayes_log_regression,
               uniform,
               nonseparable_momentum,
               samples_location,
               run_properties,
               num_samples,
               num_leapfrog_steps,
               num_newton_steps,
               rand_gen,
               true,
               save_iterations,
               suppress_print
             );

    } else
      throw Exception ("Unrecognised value for '-algorithm' option \"" + algorithm_type + "\", can only be 'riemannian'.");

  } else
    throw Exception ("Unrecognised probability type " + prob_type + "', can be either 'Prob::Test::Landscape' or 'Prob::Test::Gaussian'.");

  
  std::ofstream prop_file;

  prop_file.open((samples_location + ".prp").c_str());

  for (std::map<std::string, std::string>::iterator prop_it = run_properties.begin(); prop_it != run_properties.end(); ++prop_it)
    prop_file << prop_it->first << ": " << prop_it->second << std::endl;

  prop_file.close();


  if (prob_type == "Prob::Test::Landscape") {
    std::cout << "\n\nMATLAB plot command:\n\nplot_test_scan " << lnd_location << " -state_filename " << samples_location;

    if (algorithm_type == "metropolis")
      std::cout << " -mark 0";

  } else {

    std::cout << "\n\nMATLAB plot command:\n\nplot_test_bed " << samples_location;

    std::cout << std::endl << std::endl;
  }



}


