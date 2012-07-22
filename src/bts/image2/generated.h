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

#ifndef __bts_image2_generated_h__
#define __bts_image2_generated_h__

#include "bts/image2/reference.h"
#include "bts/image2/voxel/generated.h"

namespace BTS {

	namespace Image2 {

		class Generated : public Reference_tpl<Voxel::Generated> {

			//Public static variables, nested classes and typedefs
			public:

		    typedef Voxel::Generated Voxel;

			protected:

		    template <typename T> class Lookup : public std::vector<int> {

		      public:
		        int offset;

		        T& operator[](size_t index) {
		          T* value;
		          // If Lookup hasn't been initialised, append a new row and set the offset to the new index.
		          if (!size()) {
		            offset = index;
		            value = &push_back(T());
		          // Else if the index is within range, return the
		          } else if (index >= offset && (index - offset) < size())
		            value = &this->std::vector<int>::operator[](index - offset);
		          else {
		            int index_diff = offset - index;
		            resize(size() + abs(index_diff));
		            // If the new index is below the current one, shift all the indices up by one (could be really expensive
		            // I suppose), but it is meant to be expensive write, fast read as this case hopefully shouldn't happen
		            // very often.
		            if (index_diff < 0)
		              for (size_t i = size(); i > size() + index_diff; --i)
		                this->std::vector<int>::operator[](index) = this->std::vector<int>::operator[](index+index_diff);
		            value = &this->std::vector<int>::operator[](index - offset);
		          }
		          return *value;
		        }

		    };

			//Protected member variables
			protected:

		    MR::Math::Matrix<double> position_matrix;
        MR::Math::Matrix<double> tangent_matrix;

		    Triple<size_t> neigh_size;
		    Lookup< Lookup< Lookup< std::vector<Voxel::Generated*> > > > neigh_lookup;

			//Public member functions
			public:


		
			//Protected member functions
			protected:

		    std::vector<Voxel*>& find_neighbourhood(Coord pos) {
		      Index neigh_index((int)floor(pos[X]-0.5),(int)floor(pos[Y]-0.5),(int)floor(pos[Z]-0.5));
		      std::vector<Voxel*>& neigh = neigh_lookup[neigh_index[X]][neigh_index[Y]][neigh_index[Z]];
		      if (!neigh.size()) {
	          neigh.resize(neigh_size[X] * neigh_size[Y] * neigh_size[Z]);
	          for (size_t x = 0; x < neigh_size[X]; ++x)
	            for (size_t y = 0; y < neigh_size[Y]; ++y)
	              for (size_t z = 0; z < neigh_size[Z]; ++z)
	                neigh[i] = operator[]((int)floor(pos[X]) )
		      }
		      return neigh;
		    }

		};

	}

}

#endif /* __bts_image2_generated_h__ */
