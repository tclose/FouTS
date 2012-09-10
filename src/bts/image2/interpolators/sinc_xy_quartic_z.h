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

#ifndef __bts_image2_interpolators_sincxyquarticz_h__
#define __bts_image2_interpolators_sincxyquarticz_h__

#include "bts/image2/interpolator.h"
#include "bts/math/matrix.h"
#include "bts/math/vector.h"

namespace BTS {

	namespace Image2 {

	  namespace Interpolators {

      class SincXYQuarticZ : public BTS::Image2::Interpolator {

        //Public static variables, nested classes and typedefs
        public:


        //Protected member variables
        protected:

          double sinc_ext;

        //Public member functions
        public:

          SincXYQuarticZ(size_t sinc_extent)
            : sinc_ext((double)sinc_extent) {}

          MR::Math::Vector<double> interpolations(Math::Matrix positions, Coord voxel_centre, Math::Vector output,
                                                                                                Math::Matrix work);

        //Protected member functions
        protected:


      };

	  }

	}

}

#endif /* __bts_image2_interpolators_sincxyquarticz_h__ */
