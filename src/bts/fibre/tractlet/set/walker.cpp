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

#include "bts/fibre/tractlet/set/walker.h"

namespace BTS {

	namespace Fibre {


    const char*        Tractlet::Set::Walker::TYPE_DEFAULT            = "standard";
    const double       Tractlet::Set::Walker::STEP_SCALE_DEFAULT      = 0.001;
    const std::string  Tractlet::Set::Walker::STEP_LOCATION_DEFAULT   = str(getenv("HOME")) + str("/git/BaFTrS/params/fibre/tract/masks/mcmc/metropolis/default.tct");


    Tractlet::Set::Walker*    Tractlet::Set::Walker::factory(const Tractlet::Set& state,
                                                      const std::string& type,
                                                      double step_scale,
                                                      const std::string& step_location,
                                                      MCMC::Proposal::Distribution* distribution,
                                                      double base_intens_scale) {



      Walker* walker;

      if (type == "standard") {
        Fibre::Tractlet::Set step;

        if (step_location.size()) {
          step.load(step_location);
          if (!step.props_match(state))
            throw Exception("Properties do not match between loaded properties ( " + str(step.prop_keys()) + "," + str(step.elem_prop_keys()) + ") and state properties (" + str(state.prop_keys()) + ", " + str(state.elem_prop_keys()) + ").");

          if (step.size() == 1) {

            Fibre::Tractlet step_template = step[0];

            step.resize(state.size());

            for (size_t tract_i = 0; tract_i < state.size(); ++tract_i) {

              step[tract_i] = step_template;

              if (state[tract_i].degree() > step_template.degree())
                throw Exception ("Tractlet (" + str(tract_i) + ") is larger than step template.");
              else if (state[tract_i].degree() < step_template.degree())
                step[tract_i].redegree(state[tract_i].degree());

            }

          } else if (step.size() != state.size())
            throw Exception ("Loaded step sizes' size (" + str(step.size()) + ") does not match state (" + str(state.size()) + ").");

          step *= step_scale;

        } else {
          step = state;
          step.zero();
          step.set(step_scale);
        }

        walker = new Walker (distribution, step);

      } else

        throw Exception ("Unrecognised Fibre::Tractlet::Set::Walker type '" + type + "'.");

      return walker;



    }


    void                      Tractlet::Set::Walker::step(const Tractlet::Set& current, Tractlet::Set& proposed, double scalar) const {

      MCMC::Proposal::Walker::step(current, proposed, scalar);

    }


	}


}
