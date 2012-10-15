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
#include "bts/prob/prior_component/length.h"


namespace BTS {

  namespace Prob {

    namespace PriorComponent {

      const double            Length::SCALE_DEFAULT     = 0.0; //10
      const double            Length::MEAN_DEFAULT     = 0.05;
      const std::string       Length::NAME              = "length";

      double          Length::log_prob(const Fibre::Strand strand) {

        return - scale * MR::Math::pow2(strand[1].norm() - mean);

      }


      double          Length::log_prob(const Fibre::Strand strand, Fibre::Strand gradient) {

        assert(!gradient.is_owner());

        // If gradient hasn't been initialised, initialise it to the size of the tractlet, otherwise check its degree.
        if (!gradient.degree()) {
          gradient = strand;
          gradient.zero();
        } else if (gradient.degree() != strand.degree())
          throw Exception ("Tractlet degree (" + str(strand.degree()) + ") does not match that of the supplied gradient (" +
                                                                                        str(gradient.degree()) + ").");

        double lprob = -scale * MR::Math::pow2(strand[1].norm() - mean);

        gradient[1][X] = -2.0 * scale * strand[1][X] * (strand[1].norm() - mean) / strand[1].norm();
        gradient[1][Y] = -2.0 * scale * strand[1][Y] * (strand[1].norm() - mean) / strand[1].norm();
        gradient[1][Z] = -2.0 * scale * strand[1][Z] * (strand[1].norm() - mean) / strand[1].norm();

        return lprob;
      }

    }

  }

}
