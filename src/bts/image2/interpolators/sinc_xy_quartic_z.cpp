/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 10/09/2012.

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

#include <math>

#include "bts/image2/interpolators/sinc_xy_quartic_z.h"
#include "bts/math/vector.h"
#include "bts/math/matrix.h"

namespace BTS {

	namespace Image2 {

	  namespace Interpolators {

	    MR::Math::Vector<double> SincXYQuarticZ::interpolations(Math::Matrix positions, Coord voxel_centre,
                                                                Math::Vector output, Math::Matrix work) {

	      // Calculate displacement of positions from voxel centre
	      work = positions;
	      for (size_t dim_i = 0; dim_i < 2; ++dim_i)
	        work.column(dim_i) -= voxel_centre[dim_i];

	      // Loop through each position and calculate the
        for (size_t pos_i = 0; pos_i < positions.rows(); ++pos_i) {
          Coord interp = work.row(pos_i);
          for (size_t dim_i = 0; dim_i < 2; ++dim_i) {
            if (interp[dim_i] == 0.0)
              interp[dim_i] = 1.0;
            else if (interp[dim_i < sinc_ext)
              interp[dim_i] = sin(M_PI * interp[dim_i) / (M_PI * interp[dim_i);
            else
              interp[dim_i] = 0.0;
          }

          interp[Z] = UPTO_HERE;

          output = work

	      }

	      return output;
	    }

	  }


	}

}
