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


#ifndef __bts_mcmc_hamiltonian_h__
#define __bts_mcmc_hamiltonian_h__

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

#include <map>

#include "progressbar.h"

#include "bts/fibre/properties/extended.h"
#include "bts/mcmc/common.h"

#include "bts/common.h"



namespace BTS {

  namespace MCMC {
    
    namespace Hamiltonian {

      const size_t NUM_SAMPLES_DEFAULT = 20;
      const size_t NUM_LEAPFROG_STEPS_DEFAULT = 30;
      const size_t BURN_NUM_SAMPLES_DEFAULT = 0;
      const size_t BURN_NUM_LEAPFROG_STEPS_DEFAULT = 30;
      const size_t BURN_SNR_DEFAULT = 1;

    }


    template <typename State, typename Likelihood, typename Prior> State hamiltonian (State& initial_x,
                                                                                      Likelihood& likelihood,
                                                                                      Prior& prior,
                                                                                      MCMC::Proposal::Momentum& momentum,
                                                                                      const std::string& samples_location,
                                                                                      const std::map<std::string,std::string>& run_properties,
                                                                                      size_t num_samples,
                                                                                      size_t num_leapfrog_steps,
                                                                                      gsl_rng* rand_gen,
                                                                                      bool prior_only = false,
                                                                                      bool save_iterations = false,
                                                                                      bool suppress_print = false ) {


      std::vector<std::string> sample_header;

      sample_header.push_back(LOG_PROB_PROP);
//      sample_header.push_back(ELAPSED_TIME_PROP);
      sample_header.push_back(H_PROP);
      sample_header.push_back(PROPOSED_H_PROP);
//      sample_header.push_back("total_signal");
//      sample_header.push_back("likelihood");
//      sample_header.push_back("prior");

      std::vector<std::string> components_list = prior.list_components();

      sample_header.insert(sample_header.end(), components_list.begin(), components_list.end());
   
      typename State::Writer samples (samples_location, initial_x, sample_header, run_properties);

      typename State::Writer iterations;
//      typename State::Writer iteration_gradients;
//      typename State::Momentum::Writer iteration_momentums;
//      typename State::Writer iteration_prior_gradients;
//      typename State::Writer iteration_curv_gradients;
//      typename State::Writer iteration_end_gradients;
//      typename State::Writer iteration_all_gradients;
//      typename State::Writer iteration_likelihood_gradients;

      SamplePropertyHeader iteration_property_header;

      iteration_property_header.push_back("log_px");
      iteration_property_header.push_back("pred_d_log_px");
      iteration_property_header.push_back("act_d_log_px");
      iteration_property_header.push_back("grad_norm2");
      iteration_property_header.push_back("log_kinetic_energy");
      
//      iteration_property_header.push_back("curv_pred_d_log_px");
//      iteration_property_header.push_back("curv_act_d_log_px");
//
//      iteration_property_header.push_back("end_pred_d_log_px");
//      iteration_property_header.push_back("end_act_d_log_px");
//
//      iteration_property_header.push_back("all_pred_d_log_px");
//      iteration_property_header.push_back("all_act_d_log_px");


      if (save_iterations) {
        iterations.create(File::strip_extension(samples_location) + ".iter." + File::extension(samples_location), initial_x,  iteration_property_header, run_properties);

//        iteration_gradients.create(File::strip_extension(samples_location) + ".gradient." + File::extension(samples_location), run_properties, std::vector<std::string>());
//        iteration_momentums.create(File::strip_extension(samples_location) + ".momentum." + File::extension(samples_location), run_properties, std::vector<std::string>());
//        iteration_likelihood_gradients.create(File::strip_extension(samples_location) + ".like_gradient." + File::extension(samples_location), run_properties, std::vector<std::string>());
//        iteration_prior_gradients.create(File::strip_extension(samples_location) + ".prior_gradient." + File::extension(samples_location), run_properties, std::vector<std::string>());
//
//        iteration_curv_gradients.create(File::strip_extension(samples_location) + ".curv_gradient." + File::extension(samples_location), run_properties, std::vector<std::string>());
//        iteration_end_gradients.create(File::strip_extension(samples_location) + ".end_gradient." + File::extension(samples_location), run_properties, std::vector<std::string>());
//        iteration_all_gradients.create(File::strip_extension(samples_location) + ".all_gradient." + File::extension(samples_location), run_properties, std::vector<std::string>());
      }
        


      State x = initial_x;

      State zero = x;
      zero.zero();

      State gradient (zero);
      State prior_gradient (zero);
      State likelihood_gradient (zero);
      
//      double curv_log_prob = 0.0;
//      double end_log_prob = 0.0;
//      double all_log_prob = 0.0;

      State curv_gradient (zero);
      State end_gradient (zero);
      State all_gradient (zero);

      double start_log_kinetic_energy = NAN;

      double px = prior.log_prob(x, prior_gradient);

      if (!prior_only)
        px += likelihood.log_prob(x, likelihood_gradient);

      gradient = prior_gradient + likelihood_gradient;  

      //-------------------------//
      //  Take the MCMC samples  //
      //-------------------------//

      size_t iteration_count = 0;

      // Initialise the progress bar
      MR::ProgressBar progress_bar  ("Generating " + str(num_samples) + " Hamiltonian MCMC samples ...", num_samples);

      for (size_t sample_i = 0; sample_i < num_samples; sample_i++) {

        momentum.randomize();
        
        double log_kinetic_energy = momentum.log_kinetic_energy(); // momentum.norm2() / 2.0;
        
        if (!suppress_print)
          start_log_kinetic_energy = log_kinetic_energy;
          
        // NB: Since we want to find maxima not minima the probability of x is inverted when compared from the classical algorithm 
        // (as implemented in Mackay, D. 2003. Information Theory, Inference and Learning Algorithms, Cambridge University Press, pp388.)        
        
        double H = log_kinetic_energy - px;

        clock_t sample_starttime = clock();
        
        State prop_x = x;
        State prop_gradient = gradient;
        double prop_px = px;

        if (save_iterations) {
//          curv_log_prob = prior.component_log_prob("curv_mag", prop_x, curv_gradient);
//          end_log_prob = prior.component_log_prob("end_in_sphere", prop_x, end_gradient);
//          all_log_prob = prior.component_log_prob("all_in_cube", prop_x, all_gradient);
        }


        // Make 'num_leapfrog_steps' leapfrog steps.
        

        for (size_t leapfrog_i = 0; leapfrog_i < num_leapfrog_steps; leapfrog_i++) {

          // NB: Since we want to find maxima not minima the gradient of x is inverted when compared from the classical algorithm. 
          momentum.half_update_momentum(prop_gradient); // += prop_gradient * (leapfrog_step_size / 2.0);
          
          if (save_iterations) {

            prop_x.set_extend_prop("pred_d_log_px", str(momentum.predicted_change(prop_gradient)));
//            prop_x.set_extend_prop("curv_pred_d_log_px"] = str(momentum.predicted_change(curv_gradient));
//            prop_x.set_extend_prop("end_pred_d_log_px"] = str(momentum.predicted_change(end_gradient));
//            prop_x.set_extend_prop("all_pred_d_log_px"] = str(momentum.predicted_change(all_gradient));

          }

          // Debug purposes only.
          double prev_prop_px = prop_px;

          momentum.update_state(prop_x); //prop_x += leapfrog_step_size * momentum;

          prop_px = prior.log_prob(prop_x, prior_gradient);

          if (!prior_only)
            prop_px += likelihood.log_prob(prop_x, likelihood_gradient);
          
          prop_gradient = prior_gradient + likelihood_gradient;

          // NB: Since we want to find maxima not minima the gradient of x is inverted when compared from the classical algorithm.
          momentum.half_update_momentum(prop_gradient);  // momentum += prop_gradient * (leapfrog_step_size / 2.0);
          
          if (save_iterations) {

            prop_x.set_extend_prop("log_px", str(prop_px));
            prop_x.set_extend_prop("act_d_log_px", str(prop_px - prev_prop_px));
            prop_x.set_extend_prop("grad_norm2", str(prop_gradient.norm2()));
            prop_x.set_extend_prop("log_kinetic_energy", str(momentum.log_kinetic_energy()));
//            iteration_gradients.append(prop_gradient);
//            iteration_likelihood_gradients.append(likelihood_gradient);
//            iteration_prior_gradients.append(prior_gradient);
//            iteration_momentums.append(momentum);


//            double new_curv_log_prob = prior.component_log_prob("curv_mag", prop_x, curv_gradient);
//            double new_all_log_prob = prior.component_log_prob("end_in_sphere", prop_x, all_gradient);
//            double new_end_log_prob = prior.component_log_prob("all_in_cube", prop_x, end_gradient);
//
//            prop_x.set_extend_prop("curv_act_d_log_px"] = str(new_curv_log_prob - curv_log_prob);
//            prop_x.set_extend_prop("end_act_d_log_px"] = str(new_end_log_prob - end_log_prob);
//            prop_x.set_extend_prop("all_act_d_log_px"] = str(new_all_log_prob - all_log_prob);
//
//            iteration_curv_gradients.append(curv_gradient);
//            iteration_end_gradients.append(end_gradient);
//            iteration_all_gradients.append(all_gradient);
//
//            curv_log_prob = new_curv_log_prob;
//            end_log_prob = new_end_log_prob;
//            all_log_prob = new_all_log_prob;

            std::map<std::string, std::string> prop_component_values = prior.get_component_values(prop_x);

            for (std::map<std::string, std::string>::iterator comp_it = prop_component_values.begin(); comp_it != prop_component_values.end(); ++comp_it)
              prop_x.set_extend_prop(comp_it->first, comp_it->second);

            iterations.append(prop_x);

          }
          ++iteration_count;
        }
        
        breakpoint();
        
        // Decide whether to accept 'step' (or rather chain of leapfrog steps).
        
        log_kinetic_energy = momentum.log_kinetic_energy(); //momentum.norm2() / 2.0;
        
        // NB: Since we want to find maxima not minima the gradient of x is inverted when compared from the classical algorithm.
        double prop_H = log_kinetic_energy - prop_px;
        

        // NB: Since we want to find maxima not minima the order of this subtraction is reversed.
        double dH =  prop_H - H;



        double elapsed_time = (double)(clock() - sample_starttime) / (double)CLOCKS_PER_SEC;


        // Print out sample properties.
        if (!suppress_print) {
          std::cout << std::endl;
          std::cout << "Sample: " << sample_i << ", ";
//          std::cout << "state: " << x << ", ";
          std::cout << "log[px]: " << px << ", ";
          std::cout << "Start kinetic energy: " << start_log_kinetic_energy << ", ";
          std::cout << "H: " << H << ", ";
//          std::cout << "proposed state: " << prop_x << ", ";
          std::cout << "Proposed log[px]: " << prop_px << ", ";
          std::cout << "Final kinetic energy: " << log_kinetic_energy << ", ";          
          std::cout << "Proposed H: " << prop_H << ", ";
          std::cout << "Elapsed time: " << elapsed_time;
         
        }


        // If Hamiltonian decreases (i.e. there is less energy in the system) or the increase is below log[Y] (where Y is a uniformly distributed random variable between 0 and 1), accept the proposer.
        
        if ((dH < 0) || (log(gsl_ran_flat(rand_gen, 0,1)) < -dH)) {
          x = prop_x;
          gradient = prop_gradient;
          px = prop_px;
          if (!suppress_print)
            std::cout << ", Accepted!." << std::endl;
        } else if (!suppress_print)
            std::cout << ", Not accepted." << std::endl;

        // Record sample stats.
        
        x.set_extend_prop(LOG_PROB_PROP, str(px));
//        x.set_extend_prop(ELAPSED_TIME_PROP] = str(elapsed_time);
        x.set_extend_prop(H_PROP, str(H));
        x.set_extend_prop(PROPOSED_H_PROP, str(prop_H));
//        x.set_extend_prop("likelihood"]          = str(likelihood_log_px);
//        x.set_extend_prop("prior"]               = str(prior_px);
//        x.set_extend_prop("total_signal"]        = str(total_signal);

        std::map<std::string, std::string> component_values = prior.get_component_values(x);

        for (std::map<std::string, std::string>::iterator comp_it =component_values.begin(); comp_it != component_values.end(); ++comp_it)
          prop_x.set_extend_prop(comp_it->first, comp_it->second);

        // Save sample.
        samples.append(x);


        progress_bar++;

      }


      State out_x = x;

      return out_x;
            
    }
    


  }

}

#endif
