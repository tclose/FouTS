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


#include "bts/fibre/track/set.h"
#include "bts/prob/prior_component/density.h"
#include "bts/prob/inline_functions.h"

namespace BTS {

  namespace Prob {

    namespace PriorComponent {

      const double            Density::HIGH_SCALE_DEFAULT  = 1.0;
      const double            Density::LOW_SCALE_DEFAULT   = 1.0;
      const std::string       Density::NAME                = "density";
      const size_t            Density::NUM_POINTS_DEFAULT  = 100;

      double             Density::log_prob(const Fibre::Tractlet tractlet, Fibre::Tractlet gradient) {

        double log_prob = 0.0;
        gradient.invalidate();

        std::vector<double> areas = tractlet.cross_sectional_areas(num_points);

        for (size_t point_i = 0; point_i < num_points; ++point_i) {

          // Divide the Apparrent Connection Strength (ACS) of the tractlet by the perpendicular component of the area.
          double density = tractlet.acs() / areas[point_i];

          // Arbitrary function that makes the probability of the density being either very large or very low improbable.
          log_prob -= high_scale* MR::Math::pow2(1.0 - density);
          log_prob -= low_scale * MR::Math::pow2(1.0 - 1.0 / density);

        }

        return log_prob / (double)num_points;


      }

    }

  }

}
