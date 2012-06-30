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


#ifndef __bts_prob_priorcomponent_acs_h__
#define __bts_prob_priorcomponent_acs_h__

#include "math/math.h"

#include "bts/fibre/strand/set.h"

namespace BTS {

  namespace Prob {

    namespace PriorComponent {

      class ACS {

        public:

          const static double             SCALE_DEFAULT;
          const static double             MEAN_DEFAULT;
          const static std::string        NAME;

        //Protected member variables.
        protected:

          double scale;
          double mean;

        public:

          ACS (double scale, double mean)
          : scale(scale), mean(mean) {}

          ACS(const ACS& m)
            : scale(m.scale), mean(m.mean) {}

          ACS&  operator=(const ACS& m)
            { scale = m.scale; mean = m.mean; return *this; }

          ~ACS() {}

          ACS*       clone() const
            { return new ACS(*this); }

          double          log_prob(const Fibre::Tractlet tractlet);

          double          log_prob(const Fibre::Tractlet tractlet, Fibre::Tractlet gradient);

          const std::string&  get_name()
            { return NAME; }
      };

    }

  }

}



#endif
