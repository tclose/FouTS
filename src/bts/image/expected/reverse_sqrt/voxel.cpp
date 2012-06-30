/*
    Copyright 2010 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close on May 19, 2010.

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

#include "bts/image/expected/reverse_sqrt/voxel.h"
#include "bts/image/expected/reverse_sqrt/buffer.h"

#include "bts/image/inline_functions.h"
#include "bts/diffusion/inline_functions.h"

namespace BTS {

  namespace Image {

    namespace Expected {

      namespace ReverseSqrt {

        Voxel::Voxel(Buffer& buffer, const Index& coord)
          : Expected::Voxel(buffer, coord, buffer.diffusion_model), image(&buffer) {}


        double         Expected::ReverseSqrt::Voxel::interpolate(const Coord& pos) {

          Coord disp = pos - centre();

          double scalar = 1.0;

          if (disp.upper_bounded(1.0) && disp.lower_bounded(-1.0)) {

            for (size_t dim_i = 0; dim_i < 3; dim_i++)
              scalar *= MR::Math::sqrt(1.0 - MR::Math::abs(disp[dim_i]));// * 3.0 / 2.0;

          } else
            scalar = 0.0;

          return scalar;

        }


        double         Expected::ReverseSqrt::Voxel::interpolate(const Coord& pos, Coord& pos_gradient) {

          throw Exception ("Expected::ReverseSqrt::Voxel::interpolate(const Coord&, Coord&) not implemented yet.");

        }

      }

    }

  }

}

