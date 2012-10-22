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


#ifndef __bts_mcmc_metropolis_h__
#define __bts_mcmc_metropolis_h__

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

#include <map>

#include "progressbar.h"


#include "bts/mcmc/common.h"

#include "bts/mcmc/annealer.h"

#include "bts/image/expected/buffer.h"

#include "bts/common.h"



namespace BTS {

  namespace MCMC {
    

    namespace Metropolis {

      const size_t    NUM_ITERATIONS_DEFAULT = 1e5;
      const size_t    SAMPLE_PERIOD_DEFAULT = 1e3;
      const size_t    BURN_NUM_ITERATIONS_DEFAULT = 0; //1e4;
      const size_t    BURN_SAMPLE_PERIOD_DEFAULT = 1e2;
      const double  BURN_SNR_DEFAULT = 20;
      const double  ANNEAL_FRAC_START_DEFAULT = 1.0; //0.05;

    }

    template <typename State, typename Likelihood, typename Prior> State metropolis ( State& initial_x,
                                                                                      Likelihood& likelihood,
                                                                                      Prior& prior,
                                                                                      typename State::Walker& walker,
                                                                                      const std::string& samples_location,
                                                                                      const std::map<std::string,std::string>& run_properties,
                                                                                      size_t num_iterations,
                                                                                      size_t sample_period,
                                                                                      gsl_rng* rand_gen,
                                                                                      double anneal_frac_start = 1.0,
                                                                                      bool prior_only = false,
                                                                                      bool verbose = true,
                                                                                      bool save_images = false) {


      if (save_images)
        likelihood.get_observed_image().save(samples_location + ".obs.mif");

      std::vector<std::string> sample_header;

      sample_header.push_back(LOG_PROB_PROP);
      sample_header.push_back(ANNEAL_LOG_PROB_PROP);
      sample_header.push_back(ACCEPTANCE_RATIO_PROP);
      sample_header.push_back(ELAPSED_TIME_PROP);
      sample_header.push_back("densities");

#ifndef TEST_BED
//      sample_header.push_back("total_signal");
#endif

      sample_header.push_back("likelihood");
      sample_header.push_back("prior");

      std::vector<std::string> elem_header;

//      State::append_characteristic_keys(elem_header);

      std::vector<std::string> components_list = prior.list_components();

      sample_header.insert(sample_header.end(), components_list.begin(), components_list.end());

      typename State::Writer samples (samples_location, initial_x, sample_header, elem_header, run_properties);
//      typename State::Writer iterations (File::strip_extension(samples_location) + ".iter."  + File::extension(samples_location), run_properties, sample_header);

      State x = initial_x;
      State prop_x = x;

      double prior_px = prior.log_prob(x);


      double likelihood_px;

      MCMC::Annealer annealer(num_iterations, anneal_frac_start);

      if (prior_only)
        likelihood_px = 0;
      else
        likelihood_px = likelihood.log_prob(x);


      double px = likelihood_px * annealer.factor() + prior_px;


    //-------------------------//
    //  Take the MCMC samples  //
    //-------------------------//

      size_t num_samples = num_iterations / sample_period;

      // Initialise the progress bar
      MR::ProgressBar progress_bar  ("Generating " + str(num_samples) + " Metropolis-Hastings MCMC samples ...", num_samples);


      for (size_t sample_i = 0; sample_i < num_samples; sample_i++) {

        size_t accepted = 0;

        clock_t sample_starttime = clock();

        for (size_t iteration_i = 0; iteration_i < sample_period; iteration_i++) {

          // Draw a new sample from the proposal distribution centred on the current state.
          walker.step(x, prop_x, 1.0 / MR::Math::sqrt(annealer.factor()));

          //Calculate the unnormalised probability of the stepd state.
          double prop_prior_px = prior.log_prob(prop_x);

          double prop_likelihood_px;

          if (prior_only)
            prop_likelihood_px = 0;
          else
            prop_likelihood_px = likelihood.log_prob(prop_x);

          double prop_px = prop_likelihood_px * annealer.factor() + prop_prior_px;


//          prop_x.set_extend_prop(LOG_PROB_PROP] = str(prop_px);
//          prop_x.set_extend_prop("likelihood"]     = str(likelihood_px);
//          prop_x.set_extend_prop("prior"]          = str(prior_px);
//
//          std::map<std::string, double> component_values = prior.get_component_values(prop_x);
//
//          prop_x.properties.insert(component_values.begin(), component_values.end());
//
//          iterations.append(prop_x);

          // Divide the unnormalised probability of the stepd state by the unnormalised probability of the current state
          // (remembering that we are dealing with log probability so it is implemented as a subtraction).
          double a = prop_px - px;

//          std::cout << "prop_px: " << prop_px << " px: " << px << " a: " << a << std::endl;

//          prop_x.save("/home/tclose/data/tmp.tct");

          //If the ratio is greater than a uniform value between 0 and 1 then accept the stepd step.
          if ( (a > 0) || log(gsl_ran_flat(rand_gen, 0.0, 1.0)) <= a ) {

            // Accept the stepd tractlets
            x = prop_x;
            px = prop_px;

            //debugging only
            likelihood_px = prop_likelihood_px;
            prior_px = prop_prior_px;

            accepted++;

          }

          annealer.increment();

          px = likelihood_px * annealer.factor() + prior_px;

        }

        
        // Calculate stats about the current sample
        double acceptance_ratio = ((double)accepted)/(double)sample_period;
        double elapsed_time = (double)(clock() - sample_starttime) / (double)CLOCKS_PER_SEC;

        
        std::vector<double> areas = x[0].cross_sectional_areas(100);
        std::vector<double> densities(100);
        for (size_t i = 0; i < 100; ++i)
          densities[i] = areas[i] / x[0].acs();

        // Record sample stats.
        x.set_extend_prop(ANNEAL_LOG_PROB_PROP,str(likelihood_px + prior_px));
        x.set_extend_prop(LOG_PROB_PROP,str(px));
        x.set_extend_prop(ACCEPTANCE_RATIO_PROP,str(acceptance_ratio));
        x.set_extend_prop(ELAPSED_TIME_PROP,str(elapsed_time));
        x.set_extend_prop("densities", str(densities));

#ifndef TEST_BED
        //-------- Debugging ---------//
//        Image::Expected::Buffer& exp_image = *likelihood.get_expected_image().clone();
//
//        Image::Expected::Buffer& blank_image = *exp_image.clone();
//        Image::Expected::Buffer& dummy_image_gradient = *exp_image.clone();
//        blank_image.zero();

//        double total_signal = exp_image.sum_strong_difference(blank_image, dummy_image_gradient,0);
#endif


        State dummy_gradient = x;

        x.set_extend_prop("likelihood", str(likelihood_px));
        x.set_extend_prop("prior", str(prior_px));

#ifndef TEST_BED
//        x.set_extend_prop("total_signal"]     = str(total_signal);
#endif

        std::map<std::string, double> component_values = prior.get_component_values(x);

        for (std::map<std::string, double>::iterator comp_it = component_values.begin(); comp_it != component_values.end(); ++comp_it)
          x.set_extend_prop(comp_it->first, str(comp_it->second));

//        x.properties.insert(component_values.begin(), component_values.end());

#ifndef TEST_BED
//        likelihood.get_expected_image().save("/home/tclose/git/BaFTrS/output/signal/iter_" + str(sample_i) + ".mif");

//        delete &exp_image;
//        delete &blank_image;
//        delete &dummy_image_gradient;
#endif

        //-------- End Debugging ------//

        x.set_characteristics();

        // Save sample.
        samples.append(x);


        if (save_images) {
          likelihood.get_expected_image().save(samples_location + ".exp.mif");
        }


        // Print out sample properties.
        if (verbose) {
          std::cout << std::endl;
          std::cout << "Iteration: " << (sample_i+1) * sample_period << "/" << num_iterations <<  ", ";
          std::cout << "log[px]: " << px << ", ";
          std::cout << "acceptance ratio: " << acceptance_ratio << ", ";
          std::cout << "elapsed time: " << elapsed_time;
          std::cout << std::endl;
        }

        progress_bar++;

      }

      //MR::ProgressBar::done();

      State x_return(x);

      return x_return;
            
    }
    


  }

}

#endif
