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

#include "bts/fibre/tractlet/set/momentum.h"




namespace BTS {

	namespace Fibre{

    const double                Tractlet::Set::Momentum::STEP_SCALE_DEFAULT = 1.0;
    const char*                 Tractlet::Set::Momentum::STEP_LOCATION_DEFAULT = "/home/tclose/data/mcmc/hamiltonian/params/stepsizes/default.tct";

    Tractlet::Set::Momentum       Tractlet::Set::Momentum::factory( const Set& state,
                                                                double step_scale,
                                                                const std::string& step_location,
                                                                MCMC::Proposal::Distribution* const proposal_distribution) {
      Fibre::Tractlet::Set step;

      if (step_location.size()) {
        step.load(step_location);

        //If the loaded step only contains one strand use this as the template for the whole set
        if (step.size() == 1) {
          Fibre::Tractlet::Set full_step(state.size(), state.degree());
          full_step.copy_props(step);
          full_step.copy_elem_props(step);

          for (size_t elem_i = 0; elem_i < state.size(); ++elem_i)
            full_step[elem_i] = step[0];

          step = full_step;

        }

        else if (step.size() != state.size())
          throw Exception ("Loaded step sizes' size (" + str(step.size()) + ") does not match state (" + str(state.size()) + ").");

      } else {
        step = state;
        step.set(step_scale);
      }


      return Fibre::Tractlet::Set::Momentum (proposal_distribution, step);

    }

	}

}
