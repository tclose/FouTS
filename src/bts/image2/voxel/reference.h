/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 22/07/2012.

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

#ifndef __bts_image2_voxel_reference_h__
#define __bts_image2_voxel_reference_h__

namespace BTS {

  namespace Image2 {

    namespace Voxel {

      class Reference {

        //Public static variables, nested classes and typedefs
        public:


        //Protected member variables
        protected:

          MR::Math::Vector<double> values;

        //Public member functions
        public:

          size_t size() const
            { return values.size(); }


          void initialise(size_t num_encodings) {
            values.resize(num_encodings);
          }


        //Protected member functions
        protected:


      };

    }

	}

}

#endif /* __bts_image2_reference_h__ */
