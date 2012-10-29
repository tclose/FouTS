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


#ifndef __bts_prob_prior_strand_hook_h__
#define __bts_prob_prior_strand_hook_h__

#include "math/math.h"

#include "bts/fibre/strand/set.h"

namespace BTS {

  namespace Prob {

    namespace PriorComponent {

      class Hook {

        public:

          const static double             SCALE_DEFAULT;
          const static std::string        NAME;
          const static bool    ACTIVE_DEFAULT;
          const static size_t NUM_POINTS_DEFAULT;
          const static size_t NUM_WIDTH_SECTIONS_DEFAULT;


        //Protected member variables.
        protected:

          double scale;
          size_t num_points;
          size_t num_width_sections;

        public:

          Hook (double scale, size_t num_points, size_t num_width_sections)
          : scale(scale), num_points(num_points), num_width_sections(num_width_sections) {}

          Hook(const Hook& h)
            : scale(h.scale), num_points(h.num_points), num_width_sections(h.num_width_sections) {}

          Hook&  operator=(const Hook& h)
            { scale = h.scale; num_points = h.num_points; num_width_sections = h.num_width_sections; return *this; }

          ~Hook() {}

          Hook*       clone() const
            { return new Hook(*this); }

          double          log_prob(const Fibre::Strand strand, Fibre::Strand gradient);

          double          log_prob(const Fibre::Strand& strand, Fibre::Strand& gradient, Fibre::Strand::Tensor& hessian);

          double          log_prob(const Fibre::Tractlet strand, Fibre::Tractlet gradient);

          double          log_prob(const Fibre::Tractlet& strand, Fibre::Tractlet& gradient, Fibre::Tractlet::Tensor& hessian);

          const std::string&  get_name()
            { return NAME; }
      };

    }

  }

}



#endif
