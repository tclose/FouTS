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


#ifndef __bts_mcmc_riemannian_h__
#define __bts_mcmc_riemannian_h__

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

#include <map>

#include "progressbar.h"
#include "math/cholesky.h"

#include "bts/fibre/strand/set/tensor.h"
#include "bts/fibre/tractlet/set/tensor.h"
#include "bts/mcmc/common.h"

#include "bts/mcmc/proposal/momentum/weighted/non_separable.h"
#include "bts/mcmc/naninf_exception.h"


#include "bts/common.h"




namespace BTS {

  namespace MCMC {

    namespace Riemannian {

      const size_t NUM_SAMPLES_DEFAULT = 100;
      const size_t NUM_LEAPFROG_STEPS_DEFAULT = 30;
      const size_t NUM_NEWTON_STEPS_DEFAULT = 6;
      const size_t BURN_NUM_SAMPLES_DEFAULT = 0;
      const size_t BURN_NUM_LEAPFROG_STEPS_DEFAULT = 30;
      const size_t BURN_NUM_NEWTON_STEPS_DEFAULT = 6;
      const size_t BURN_SNR_DEFAULT = 1;
      const float PRECONDITION_DEFAULT = 1000;

    }

    template <typename State_T, typename Prior_T, typename Likelihood_T> class Posterior {

      //Public static variables, nested classes and typedefs
      public:


      //Protected member variables
      protected:

        Prior_T& prior;
        Likelihood_T& likelihood;
        double precondition;

        size_t dimension;

        State_T prior_gradient, likelihood_gradient;
        typename State_T::Tensor prior_fisher, likelihood_fisher;
        std::vector<typename State_T::Tensor> prior_fisher_gradient, likelihood_fisher_gradient;

        typename State_T::Tensor::Writer fisher_writer;

      //Public member functions
      public:

        Posterior(const State_T& state, Prior_T& prior, Likelihood_T& likelihood, double precondition = 0.0)
          : prior(prior),
            likelihood(likelihood),
            precondition(precondition),
            dimension(state.vsize()),
            prior_gradient(state),
            likelihood_gradient(state),
            prior_fisher(state),
            likelihood_fisher(state) {

          prior_gradient.zero();
          likelihood_gradient.zero();

          prior_fisher.zero();
          likelihood_fisher.zero();

          for (size_t elem_i = 0; elem_i < dimension; ++elem_i) {
            prior_fisher_gradient.push_back(prior_fisher);
            likelihood_fisher_gradient.push_back(likelihood_fisher);
          }

          fisher_writer.create("/home/tclose/data/mcmc/riemannian/fishers.sta.tnr", state);

        }


        double      log_prob_and_fisher (const State_T& state, State_T& gradient, typename State_T::Tensor& fisher, std::vector<typename State_T::Tensor>& fisher_gradient) {

          double px = prior.log_prob_and_fisher(state, prior_gradient, prior_fisher, prior_fisher_gradient);
          px += likelihood.log_prob_and_fisher(state, likelihood_gradient, likelihood_fisher, likelihood_fisher_gradient);

          gradient = prior_gradient + likelihood_gradient;
          fisher = prior_fisher;
          fisher += likelihood_fisher;

          if (precondition) {

            for (size_t elem_i = 0; elem_i < dimension; ++elem_i)
              fisher(elem_i,elem_i) += precondition;

          }

          for (size_t elem_i = 0; elem_i < dimension; ++elem_i) {
            fisher_gradient[elem_i] = prior_fisher_gradient[elem_i];
            fisher_gradient[elem_i] += likelihood_fisher_gradient[elem_i];
          }

          fisher_writer.append(fisher);

          return px;

        }


        double    log_prob_and_fisher (const State_T& state, State_T& gradient, typename State_T::Tensor& fisher) {

          double px = prior.log_prob_and_fisher(state, prior_gradient, prior_fisher);
          px += likelihood.log_prob_and_fisher(state, likelihood_gradient, likelihood_fisher);

          gradient = prior_gradient;
          gradient += likelihood_gradient;

          fisher = prior_fisher;
          fisher += likelihood_fisher;

          if (precondition) {

            for (size_t elem_i = 0; elem_i < dimension; ++elem_i)
              fisher(elem_i,elem_i) += precondition;

          }

          fisher_writer.append(fisher);

          return px;

        }

      //Protected member functions
      protected:


    };


    template <typename State_T, typename Likelihood_T, typename Prior_T> State_T riemannian (State_T& initial_x,
                                                                                      Likelihood_T& likelihood,
                                                                                      Prior_T& prior,
                                                                                      MCMC::Proposal::Momentum::Weighted::NonSeparable& momentum,
                                                                                      const std::string& samples_location,
                                                                                      const std::map<std::string,std::string>& run_properties,
                                                                                      size_t num_samples,
                                                                                      size_t num_leapfrog_steps,
                                                                                      size_t num_newton_steps,
                                                                                      gsl_rng* rand_gen,
                                                                                      double precondition = 0.0,
                                                                                      bool prior_only = false,
                                                                                      bool save_iterations = false,
                                                                                      bool suppress_print = false ) {



      Posterior<State_T,Prior_T,Likelihood_T> posterior(initial_x, prior, likelihood, precondition);

      std::vector<std::string> sample_header;

      sample_header.push_back(LOG_PROB_PROP);
      sample_header.push_back(H_PROP);
      sample_header.push_back(PROPOSED_H_PROP);


      std::vector<std::string> components_list = prior.list_components();

      sample_header.insert(sample_header.end(), components_list.begin(), components_list.end());

      typename State_T::Writer samples (samples_location, initial_x, sample_header, run_properties);

      typename State_T::Writer iterations, gradient_iterations;


      SamplePropertyHeader iteration_property_header;

      iteration_property_header.push_back("log_px");
      iteration_property_header.push_back("pred_d_log_px");
      iteration_property_header.push_back("act_d_log_px");
      iteration_property_header.push_back("grad_norm2");
      iteration_property_header.push_back("log_kinetic_energy");



      if (save_iterations) {
        iterations.create(File::strip_extension(samples_location) + ".iter." + File::extension(samples_location), initial_x, iteration_property_header, run_properties);
        gradient_iterations.create(File::strip_extension(samples_location) + ".grad_iter." + File::extension(samples_location), initial_x, SamplePropertyHeader(), run_properties);
      }



      State_T x = initial_x;

      size_t dimension = x.vsize();

      State_T gradient(x);
      gradient.zero();

      typename State_T::Tensor fisher(x);
      fisher.zero();

      std::vector<typename State_T::Tensor> fisher_gradient;
      for (size_t elem_i = 0; elem_i < dimension; ++elem_i)
        fisher_gradient.push_back(fisher);


      typename State_T::Tensor fisher_chol (fisher);


      double px = posterior.log_prob_and_fisher(x, gradient, fisher, fisher_gradient);

      fisher_chol = fisher;

      MR::Math::Cholesky::decomp(fisher_chol);


      //-------------------------//
      //  Take the MCMC samples  //
      //-------------------------//

      size_t iteration_count = 0;

      // Initialise the progress bar
      MR::ProgressBar progress_bar  ("Generating " + str(num_samples) + " Riemannian Hamiltonian MCMC samples ...", num_samples);

      for (size_t sample_i = 0; sample_i < num_samples; sample_i++) {

        momentum.randomize(fisher_chol);

        double log_kinetic_energy = momentum.log_kinetic_energy(fisher_chol); //(LOG_2 + LOG_PI * dimension + fisher_det) / 2.0

        // NB: Since we want to find maxima not minima the probability of x is inverted when compared from the classical algorithm
        // (as implemented in Mackay, D. 2003. Information Theory, Inference and Learning Algorithms, Cambridge University Press, pp388.)

        double H = log_kinetic_energy - px;

        clock_t sample_starttime = clock();

        State_T prop_x = x;
        State_T prop_gradient = gradient;
        typename State_T::Tensor prop_fisher = fisher;
        typename State_T::Tensor prop_fisher_chol = fisher_chol;
        double prop_px = px;


        //Randomly select to evolve forwards or backwards in time.
        double time_direction;

        if (gsl_rng_uniform(rand_gen) > 0.5)
          time_direction = 1.0;
        else
          time_direction = -1.0;

//        MR::Math::Vector<double> time_direction_matrix ("/home/tclose/data/mcmc/state/time_direction.txt");
//
//        time_direction = time_direction_matrix[0];

        try {

          // Make 'num_leapfrog_steps' leapfrog steps.
          for (size_t leapfrog_i = 0; leapfrog_i < num_leapfrog_steps; leapfrog_i++) {


              if (save_iterations) {

                prop_x.set_extend_prop("pred_d_log_px",str(momentum.predicted_change(prop_gradient, prop_fisher_chol, time_direction)));

              }

              // NB: Since we want to find maxima not minima the gradient of x is inverted when compared from the classical algorithm.
              momentum.half_update_momentum<State_T>(prop_gradient, prop_fisher_chol, fisher_gradient, time_direction, num_newton_steps);

              // Debug purposes only.
              double prev_prop_px = prop_px;

              MR::Math::Vector<double> fishinv_momen(dimension), tmp_fishinv_momen(dimension);
              MR::Math::Vector<double>& prop_x_vector = prop_x;
              State_T tmp_x (prop_x);
              MR::Math::Vector<double>& tmp_x_vector = tmp_x;

              MR::Math::Cholesky::solve(fishinv_momen, fisher_chol, momentum.momentum());

              for (size_t newton_i = 0; newton_i < num_newton_steps; ++newton_i) {

                posterior.log_prob_and_fisher(tmp_x, prop_gradient, prop_fisher);

                MR::Math::Cholesky::decomp(prop_fisher_chol = prop_fisher);
                MR::Math::Cholesky::solve(tmp_fishinv_momen, prop_fisher_chol, momentum.momentum());

                for (size_t elem_i = 0; elem_i < dimension; ++elem_i) {
                  tmp_x_vector[elem_i] = prop_x_vector[elem_i] + (fishinv_momen[elem_i] + tmp_fishinv_momen[elem_i]) * momentum.step_sizes()[elem_i] * time_direction * 0.5;
                  if (isnan(tmp_x_vector[elem_i]) || isinf(tmp_x_vector[elem_i]))
                    throw NanInfException();
    //                throw Exception ("State_T has NaN values in " + str(newton_i) + " newton step of " + str(leapfrog_i) + " leapfrog step.\n\n" + str(tmp_x_vector) + "\n\n");
                }


              }

              prop_x = tmp_x;

              prop_px = posterior.log_prob_and_fisher(prop_x, prop_gradient, prop_fisher, fisher_gradient);

              MR::Math::Cholesky::decomp(prop_fisher_chol = prop_fisher);

              // NB: Since we want to find maxima not minima the gradient of x is inverted when compared from the classical algorithm.
              momentum.half_update_momentum<State_T>(prop_gradient, prop_fisher_chol, fisher_gradient, time_direction);

              if (save_iterations) {

                prop_x.set_extend_prop("log_px",str(prop_px));
                prop_x.set_extend_prop("act_d_log_px",str(prop_px - prev_prop_px));
                prop_x.set_extend_prop("grad_norm2",str(prop_gradient.norm2()));
                prop_x.set_extend_prop("log_kinetic_energy",str(momentum.log_kinetic_energy(prop_fisher_chol)));

                std::map<std::string, std::string> prop_component_values = prior.get_component_values(prop_x);

                for (std::map<std::string, std::string>::iterator comp_it = prop_component_values.begin(); comp_it != prop_component_values.end(); ++comp_it)
                  prop_x.set_extend_prop(comp_it->first, comp_it->second);

                iterations.append(prop_x);
                gradient_iterations.append(prop_gradient);

              }


              ++iteration_count;

          }

          // Decide whether to accept 'step' (or rather chain of leapfrog steps).
          double prop_log_kinetic_energy = momentum.log_kinetic_energy(prop_fisher_chol);

          // NB: Since we want to find maxima not minima the gradient of x is inverted when compared from the classical algorithm.
          double prop_H = prop_log_kinetic_energy - prop_px;

          // NB: Since we want to find maxima not minima the order of this subtraction is reversed.
          double dH =  prop_H - H;

          double elapsed_time = (double)(clock() - sample_starttime) / (double)CLOCKS_PER_SEC;

          // Print out sample properties.
          if (!suppress_print) {
            std::cout << std::endl;
            std::cout << "Sample: " << sample_i << ", ";
            std::cout << "log[px]: " << px << ", ";
            std::cout << "Start kinetic energy: " << log_kinetic_energy << ", ";
            std::cout << "H: " << H << ", ";
            std::cout << "Proposed log[px]: " << prop_px << ", ";
            std::cout << "Final kinetic energy: " << prop_log_kinetic_energy << ", ";
            std::cout << "Proposed H: " << prop_H << ", ";
            std::cout << "Elapsed time: " << elapsed_time;

          }


          // If Hamiltonian decreases (i.e. there is less energy in the system) or the increase is below log[Y] (where Y is a uniformly distributed random variable between 0 and 1), accept the proposer.
          if ((dH < 0) || (log(gsl_ran_flat(rand_gen, 0,1)) < -dH)) {
            x = prop_x;
            gradient = prop_gradient;
            fisher = prop_fisher;
            fisher_chol = prop_fisher_chol;
            px = prop_px;

            if (!suppress_print)
              std::cout << ", Accepted!." << std::endl;

          } else if (!suppress_print)
              std::cout << ", Not accepted." << std::endl;

          // Record sample stats.

          x.set_extend_prop(LOG_PROB_PROP, str(px));
          x.set_extend_prop(H_PROP, str(H));
          x.set_extend_prop(PROPOSED_H_PROP,  str(prop_H));


          std::map<std::string, std::string> component_values = prior.get_component_values(x);

          for (std::map<std::string,std::string>::iterator comp_it = component_values.begin(); comp_it != component_values.end(); ++comp_it)
            x.set_extend_prop(comp_it->first, comp_it->second);


        } catch (const NanInfException& e) {

          std::cout << std::endl;
          std::cout << "Sample: " << sample_i << " - ";
          std::cout << "Not accepted due to NaN or Inf values found" << std::endl;

        }

        // Save sample.
        samples.append(x);


        progress_bar++;

      }


      State_T out_x = x;

      return out_x;

    }



  }

}

#endif
