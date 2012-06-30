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


#ifndef __bts_prob_prior_strand_curvature_h__
#define __bts_prob_prior_strand_curvature_h__

#include "math/math.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/strand/tensor.h"
#include "bts/fibre/strand/set/tensor.h"

namespace BTS {

  namespace Prob {

    namespace PriorComponent {

      class Magnitude {

        public:

          const static double             SCALE_DEFAULT;
          const static double             AUX_SCALE_DEFAULT;
          const static std::string        NAME;

        //Protected member variables.
        protected:

          double scale;
          double aux_scale;

        public:

          Magnitude (double scale, double aux_scale)
            : scale(scale), aux_scale(aux_scale) {}

          Magnitude(const Magnitude& m)
            : scale(m.scale), aux_scale(aux_scale) {}

          Magnitude&  operator=(const Magnitude& m)
            { scale = m.scale; aux_scale = m.aux_scale; return *this; }

          ~Magnitude() {}

          Magnitude*      clone() const
            { return new Magnitude(*this); }

          double          log_prob(const Fibre::Strand strand, Fibre::Strand gradient);

          double          log_prob(const Fibre::Tractlet, Fibre::Tractlet gradient);

          double          log_prob(const Fibre::Strand strand, Fibre::Strand gradient, Fibre::Strand::Tensor hessian)  { throw Exception ("Not implemented yet."); }

          double          log_prob_and_fisher(const Fibre::Strand strand, Fibre::Strand gradient, Fibre::Strand::Tensor fisher) { throw Exception ("Not implemented yet."); }


          double          log_prob_and_fisher(const Fibre::Strand strand, Fibre::Strand gradient, Fibre::Strand::Tensor fisher, std::vector<Fibre::Strand::Tensor>& fisher_gradient) { throw Exception ("Not implemented yet."); }

          const std::string&  get_name()
            { return NAME; }
      };

    }

  }

}



#endif
