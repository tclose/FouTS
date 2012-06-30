/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, Aug 2, 2010.

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

#ifndef __bts_fibre_tractlet_set_momentum_h__
#define __bts_fibre_tractlet_set_momentum_h__

#include "bts/fibre/tractlet/set.h"
#include "bts/mcmc/proposal/distribution.h"
#include "bts/mcmc/proposal/momentum.h"

namespace BTS {

	namespace Fibre {

		class Tractlet::Set::Momentum {

			//Public static variables, nested classes and typedefs
			public:

        static const double  STEP_SCALE_DEFAULT;
        static const char*   STEP_LOCATION_DEFAULT;

			protected:

        MCMC::Proposal::Momentum mcmc_momentum;

			//Public static functions.
			public:

        static Momentum      factory( const Set& state,
                                      double step_scale,
                                      const std::string& step_location,
                                      MCMC::Proposal::Distribution* const proposal_distribution);



			//Public member functions
			public:

      Momentum(MCMC::Proposal::Distribution* const prop_distr = 0)
        : mcmc_momentum(prop_distr) {}

      Momentum(MCMC::Proposal::Distribution* const proposal_distribution,
               const Set& relative_step_sizes)

        : mcmc_momentum(proposal_distribution, relative_step_sizes) {}


      void              half_update_momentum(const Set& gradient)
        { mcmc_momentum.half_update_momentum(gradient); }

      void              update_state(Set& state) const
        { MR::Math::Vector<double>& state_vector = state; mcmc_momentum.update_state(state_vector); }


      double            predicted_change(const Set& gradient) const
        { return mcmc_momentum.predicted_change(gradient); }

      double            log_kinetic_energy() const
        { return mcmc_momentum.log_kinetic_energy(); }

      void            randomize()
        { mcmc_momentum.randomize(); }

		};

	}

}

#endif /* __bts_fibre_tractlet_set_momentum_h__ */
