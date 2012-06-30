/*
 *  test_gaussian.cpp
 *  tractography
 *
 *  Created by Tom Close on 16/12/09.
 *  Copyright 2009 National ICT Australia (NICTA). All rights reserved.
 *
 */

#include "bts/prob/test/gaussian.h"
#include "bts/mcmc/state/tensor.h"

#include "bts/prob/inline_functions.h"

#include "bts/image/buffer.cpp.h"

namespace BTS {

  namespace Prob {

    namespace Test {

      const char*            Gaussian::AXIS_SCALES_LOCATION_DEFAULT = "/home/tclose/data/prob/test/gaussian/axis_scales.sta";

  
      double                 Gaussian::log_prob(const MCMC::State& state) {

        double log_prob = 0.0;

        for (size_t i = 0; i < state.size(); i++)
          log_prob += -MR::Math::pow2(state[i]) * axis_scales[i];

        return log_prob;

      }


      double                 Gaussian::log_prob(const MCMC::State& state, MCMC::State& gradient) {

        for (size_t i = 0; i < state.size(); i++)
          gradient[i] = -2.0 * state[i] * axis_scales[i];

        return log_prob(state);

      }
    

      double                  Gaussian::log_prob(const MCMC::State& state, MCMC::State& gradient, MCMC::State::Tensor& hessian) {

        hessian.zero();

        for (size_t i = 0; i < state.size(); i++)
          hessian(i,i) = 2.0 * axis_scales[i];

        return log_prob(state, gradient);

      }

      double                  Gaussian::log_prob_and_fisher(const MCMC::State& state, MCMC::State& gradient, MCMC::State::Tensor& fisher, std::vector<MCMC::State::Tensor>& fisher_gradient) {

        fisher_gradient.resize(state.size());

        for (size_t i = 0; i < state.size(); i++) {
          fisher_gradient[i].resize(state.size(), state.size());
          fisher_gradient[i] = 0.0;
        }

        return log_prob(state, gradient, fisher);

      }

    }
    
  }
  
}
