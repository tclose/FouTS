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
 
 
#ifndef __bts_prob_uniform__
#define __bts_prob_uniform__

#include "bts/prob/likelihood.h"
#include "bts/prob/prior.h"
#include "bts/prob/uniform.h"

#include "bts/mcmc/state.h"
#include "bts/mcmc/state/tensor.h"


namespace BTS {

  namespace Prob {

    //Used as a filler instead of adding any particular prior.
    class Uniform {

      public:

        const static std::string NAME;
    
      public:
      

        double               get_scalar()
          { return 1.0; }

        Uniform*    clone() const
          { return new Uniform(*this); }

        double       log_prob(const Fibre::Strand& strand)
          { return 0.0; }

        double       log_prob(const Fibre::Tractlet& tractlet)
          { return 0.0; }

        double       log_prob(const Fibre::Strand::Set& strand)
          { return 0.0; }

        double       log_prob(const Fibre::Tractlet::Set& tractlet)
          { return 0.0; }

        double       log_prob(const Fibre::Strand& strand, Fibre::Strand& gradient)
          { gradient.zero(); return 0.0; }

        double       log_prob(const Fibre::Tractlet& tractlet, Fibre::Tractlet& gradient)
          { gradient.zero(); return 0.0; }

        double       log_prob(const Fibre::Strand::Set& strand, Fibre::Strand::Set& gradient)
          { gradient.zero(); return 0.0; }

        double       log_prob(const Fibre::Tractlet::Set& strand, Fibre::Tractlet::Set& gradient)
          { gradient.zero(); return 0.0; }    

        double       log_prob(const MCMC::State& state, MCMC::State& gradient)
          { gradient.zero(); return 0.0; }

        double       log_prob(const double& f, double& gradient)
          { gradient = 0.0; return 0.0; }

        double       log_prob_and_fisher(const MCMC::State& state, MCMC::State& gradient, MCMC::State::Tensor& fisher)
          { gradient.zero(); fisher.zero(); return 0.0; }

        double       log_prob_and_fisher(const MCMC::State& state, MCMC::State& gradient, MCMC::State::Tensor& fisher, std::vector<MCMC::State::Tensor>& fisher_gradient)
          { gradient.zero(); fisher.zero(); for (size_t i =0; i < fisher_gradient.size(); ++i) fisher_gradient[i] = 0.0; return 0.0; }

        double       component_log_prob(const std::string&, const BTS::Fibre::Strand& strand, BTS::Fibre::Strand& gradient)
          { gradient.zero(); return 0.0; }

        double       component_log_prob(const std::string&, const BTS::Fibre::Strand::Set& strands, BTS::Fibre::Strand::Set& gradient)
        { gradient.zero(); return 0.0; }

        double       component_log_prob(const std::string&, const BTS::Fibre::Tractlet& tractlet, BTS::Fibre::Tractlet& gradient)
                  { gradient.zero(); return 0.0; }

        double       component_log_prob(const std::string&, const BTS::Fibre::Tractlet::Set& tractlets, BTS::Fibre::Tractlet::Set& gradient)
          { gradient.zero(); return 0.0; }

        double                  log_prob(double expected, double observed) { return 0; }


        double                  b0_log_prob(double expected, double observed) { return 0; }


        double                  log_prob(double expected, double observed, double& d_lprob) { return 0; }


        double                  log_prob(double expected, double observed, double& d_lprob, double& d2_lprob2) { return 0; }



        void         set_enforce_bounds(bool flag) {}

        std::vector<std::string> list_components() { return std::vector<std::string>(); }

        template <typename T> std::map<std::string, std::string> get_component_values(T& t) { return std::map<std::string, std::string>(); }

        const std::string& get_name()
          { return NAME; }

        void      set_assumed_snr(double assumed_snr, const std::string& ref_b0, double ref_signal) {}

    };
    
  }

}
 
#endif
