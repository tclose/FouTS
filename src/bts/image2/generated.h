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

		      public:

		        T& operator[](size_t index) {
              int index_diff;
              //If Lookup hasn't been initialised yet, manually set the index diff to one.
              if (size())
                index_diff = offset - index;
              else {
                index_diff = 1;
                offset = index;
		          }
		          // If Lookup isn't big enough, append new rows, shifting the current elements if necessary (i.e. if index
              // < offset), and set the offset to the new index.
		          if (index < offset || index_diff >= size()) {
		            resize(size() + abs(index_diff));
		            // If the new index is below the current one, shift all the indices up by one (could be really expensive
		            // I suppose), but it is meant to be expensive write, fast read as this case hopefully shouldn't happen
		            // very often.
		            if (index_diff < 0) {
		              for (int i = size()-1; i >= size() + index_diff; --i)
		                this->std::vector<int>::operator[](index) = this->std::vector<int>::operator[](index+index_diff);
		              offset = index;
		            }
		          }
              return this->std::vector<int>::operator[](index - offset);
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

		    Voxel&         operator[](int x, int y, int z) {
		      return operator[](Index(x,y,z));
        }

		    Voxel&         operator[](const Index& index) {
		      return voxels[index];
	      }
		
			//Protected member functions
			protected:

		    std::vector<Voxel*>& find_neighbourhood(Coord pos) {
		      Index neigh_index((int)(pos[X]-0.5),(int)(pos[Y]-0.5),(int)(pos[Z]-0.5));
		      std::vector<Voxel*>& neigh = neigh_lookup[neigh_index[X]][neigh_index[Y]][neigh_index[Z]];
		      if (!neigh.size()) {
	          neigh.resize(neigh_size[X] * neigh_size[Y] * neigh_size[Z]);
	          for (int x = -neigh_size[X]; x < neigh_size[X]; ++x)
	            for (int y = -neigh_size[Y]; y < neigh_size[Y]; ++y)
	              for (int z = -neigh_size[Z]; z < neigh_size[Z]; ++z)
	                neigh[x * neigh_size[Y] * neigh_size[Z] + y * neigh_size[Z] t] =
	                                               operator[](neigh_index[X] + x, neigh_index[Y] + y, neigh_index[Z] + z);
		      }
		      return neigh;
		    }

		};

	}

}

#endif /* __bts_image2_generated_h__ */
