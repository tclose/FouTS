/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

 This file is part of Fourier Tract Sampling (FouTS).

 FouTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FouTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FTS.  If not, see <http://www.gnu.org/licenses/>.

 */


#ifndef __bts_prob_prior_strand_endonsphere_h__
#define __bts_prob_prior_strand_endonsphere_h__

#include "bts/fibre/strand/set.h"
#include "bts/image/observed/buffer.h"

namespace FTS {

  namespace Prob {

      namespace PriorComponent {

        class EndOnSphere {

          public:

            const static std::string      NAME;
            const static double           SCALE_DEFAULT;
            const static double           RADIUS_DEFAULT;

          protected:

            double scale;
            double radius;

          public:

            static double           radius_default(size_t dim, double vox_length);

            static double           radius_default(const Image::Observed::Buffer& obs_image);

          public:

            EndOnSphere(double scale, double radius)
              : scale(scale), radius(radius) {}


            EndOnSphere(const EndOnSphere& s)
              : scale(s.scale), radius(s.radius) {}

            EndOnSphere&       operator= (const EndOnSphere& s)
              { scale = s.scale; radius = s.radius; return *this; }

            EndOnSphere*       clone() const
              { return new EndOnSphere(*this); }


            double       log_prob(const Fibre::Strand& strand, Fibre::Strand& gradient);

//            double       log_prob(const Fibre::Strand& strand, Fibre::Strand& gradient, Fibre::Strand::Tensor& hessian);

//            double       log_prob_and_fisher(const Fibre::Strand& strand, Fibre::Strand& gradient, Fibre::Strand::Tensor& fisher);
//
//            double       log_prob_and_fisher(const Fibre::Strand& strand, Fibre::Strand& gradient, Fibre::Strand::Tensor& fisher, std::vector<Fibre::Strand::Tensor>& fisher_gradient);

            const std::string&  get_name()
              { return NAME; }

        };

     }


  }

}

#endif /* __bts_prob_prior_strand_endonsphere_h__ */
