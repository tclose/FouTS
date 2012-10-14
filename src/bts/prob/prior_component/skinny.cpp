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

#include "math/math.h"
#include "bts/prob/prior_component/skinny.h"


namespace BTS {

  namespace Prob {

    namespace PriorComponent {

      const double            Skinny::SCALE_DEFAULT     = 0.0; //10
      const size_t            Skinny::POWER_DEFAULT     = 2;
      const std::string       Skinny::NAME              = "Skinny";

      double          Skinny::log_prob(const Fibre::Tractlet tract) {
        if (power != 2)
          throw Exception("Power can only be 2 at this point");
        return - scale * MR::Math::pow2(tract[1].norm() / tract.acs());
      }


      double          Skinny::log_prob(const Fibre::Tractlet tract, Fibre::Tractlet gradient) {
        if (power != 2)
          throw Exception("Power can only be 2 at this point");
        assert(!gradient.is_owner());
        // If gradient hasn't been initialised, initialise it to the size of the tractlet, otherwise check its degree.
        if (!gradient.degree()) {
          gradient = tract;
          gradient.zero();
        } else if (gradient.degree() != tract.degree())
          throw Exception ("Tractlet degree (" + str(tract.degree()) + ") does not match that of the supplied gradient (" +
                                                                                        str(gradient.degree()) + ").");

        double lprob = -scale * MR::Math::pow2(tract[1].norm() / tract.acs());

        gradient.set_acs(2.0 * scale * MR::Math::pow2(tract[1].norm()) / MR::Math::pow3(tract.acs()));
        gradient[1][X] = -2.0 * scale * tract[1][X] * tract[1].norm() / MR::Math::pow2(tract.acs());
        gradient[1][Y] = -2.0 * scale * tract[1][Y] * tract[1].norm() / MR::Math::pow2(tract.acs());
        gradient[1][Z] = -2.0 * scale * tract[1][Z] * tract[1].norm() / MR::Math::pow2(tract.acs());

        return lprob;
      }

    }

  }

}
