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
 
 
#ifndef __bts_prob_testgaussian__
#define __bts_prob_testgaussian__



#include "bts/prob/likelihood.h"
#include "bts/prob/prior.h"

//TODO: remove
#include "bts/image/expected/buffer.h"
#include "bts/mcmc/state.h"

namespace BTS {

  namespace Prob {

    namespace Test {

      class Gaussian {

        public:

        const static char* AXIS_SCALES_LOCATION_DEFAULT;


        protected:

          MR::Math::Vector<double>& axis_scales;

        public:

          //TODO: remove
          BTS::Image::Expected::Buffer* exp_image;

          Gaussian(MR::Math::Vector<double>& axis_scales) : axis_scales(axis_scales) {}

          double                         log_prob(const MCMC::State& test_state);

          double                         log_prob(const MCMC::State& test_state, MCMC::State& gradient);

          double                         log_prob(const MCMC::State& test_state, MCMC::State& gradient, MCMC::State::Tensor& hessian);

          double                         log_prob_and_fisher(const MCMC::State& test_state, MCMC::State& gradient, MCMC::State::Tensor& fisher)
            { return log_prob(test_state,gradient,fisher); }


          double                         log_prob_and_fisher(const MCMC::State& test_state, MCMC::State& gradient, MCMC::State::Tensor& fisher, std::vector<MCMC::State::Tensor>& fisher_gradient);

          void                           set_assumed_snr(double assumed_snr, const std::string& ref_b0, double ref_signal) {}

          void                           set_enforce_bounds(bool flag) {}

          std::vector<std::string> list_components() { return std::vector<std::string>(); }

          template <typename T> std::map<std::string, std::string> get_component_values(T& t) { return std::map<std::string, std::string>(); }


        friend std::ostream&             operator<< (std::ostream& stream, const Gaussian& test);

      };
      
      inline std::ostream&               operator<< (std::ostream& stream, const Gaussian& test) {

      
        stream << "Axis scales: " << test.axis_scales << std::endl;

        return stream;
      
      }
      
    }
    
  }
  
}

#endif
