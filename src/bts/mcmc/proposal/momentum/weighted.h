/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 15/09/2010.

    This file is part of MRtrix.

    MRtrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRtrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef __bts_mcmc_proposal_momentum_weighted_h__
#define __bts_mcmc_proposal_momentum_weighted_h__

#include "bts/mcmc/proposal/momentum.h"

namespace BTS {

	namespace MCMC {

	  namespace Proposal {

      class Momentum::Weighted : public Momentum {

        //Public nested classes
        public:

          class NonSeparable;

        public:


          template <typename T> static Weighted  factory(const T& state,
                                                        double step_scale,
                                                        const std::string& step_location,
                                                        MCMC::Proposal::Distribution* const proposal_distribution);

        protected:

          MR::Math::Vector<double> momen_step;

        //Public static functions
        public:

        //Public member functions
        public:

          Weighted(Distribution* const prop_distr = 0) : Momentum(prop_distr) {}


          Weighted(Distribution* const proposal_distribution, const MR::Math::Vector<double>& relative_step_sizes);


          ~Weighted() {}


          Weighted(const Weighted& W)
            : Momentum(W), momen_step(W.size()) {}


          Weighted&                               operator=(const Weighted& W)
            { Momentum::operator=(W); momen_step.resize(W.size()); return *this; }


          virtual void                            update_state(MR::Math::Vector<double>& state, const MR::Math::Matrix<double>& weights_chol, double time_direction = 1.0) const;


          template <typename T> void              update_state(T& state, const MR::Math::Matrix<double>& weights_chol, double time_direction = 1.0) const
            { MR::Math::Vector<double>& state_vector = state; update_state(state_vector, weights_chol, time_direction);  }


          void                                    randomize(const MR::Math::Matrix<double>& cholesky);


          double                                  log_kinetic_energy(const MR::Math::Matrix<double>& cholesky) const;

          //For debugging purposes.
          double                                  predicted_change(const MR::Math::Vector<double>& gradient, const MR::Math::Matrix<double>& weights_chol, double time_direction = 1.0) const;


          template <typename T> double            predicted_change(const T& gradient, const MR::Math::Matrix<double>& weights_chol, double time_direction = 1.0) const
            { return predicted_change(gradient, weights_chol, time_direction); }

        protected:

          double                                  log_normalization(double determinant) const;




      };

	  }

	}

}

#endif /* __bts_mcmc_proposal_momentum_weighted_h__ */
