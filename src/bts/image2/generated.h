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
#include "bts/image2/index.h"

namespace BTS {
    
    namespace Image2 {
        
        class Generated: public Reference_tpl<Voxel::Generated> {
                
                //Public static variables, nested classes and typedefs
            public:
                
                typedef Voxel::Generated Voxel;

            protected:
                
                template<typename T> class Lookup: public std::vector<T> {
                        
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
                                // If the new index is below the current one, shift all the indices up by one (could be expensive
                                // I suppose, but it is meant to be expensive write, fast read as this case hopefully shouldn't happen
                                // very often).
                                if (index_diff < 0) {
                                    for (int i = size() - 1; i >= size() + index_diff; --i)
                                        this->std::vector<T>::operator[](index) = this->std::vector<
                                                T>::operator[](index + index_diff);
                                    offset = index;
                                }
                            }
                            return this->std::vector<T>::operator[](index - offset);
                        }
                        
                };

                //Protected member variables
            protected:
//        size_t num_encodings;
//        Triple<size_t> dims;
//        Triple<double> voxel_size;
//        Triple<double> offset;
//        std::map<Index,T> voxels;
//        T empty_voxel;
                
                Reference* reference;
                Diffusion::Model diff_model;
                Interpolator* interpolator;
                MR::Math::Matrix<double> position_matrix;
                MR::Math::Matrix<double> tangent_matrix;

                Triple<size_t> neigh_size;
                Lookup<Lookup<Lookup<std::vector<Voxel::Generated*> > > > neigh_lookup;

                //Public member functions
            public:
                
                Voxel& operator[](int x, int y, int z) {
                    return operator[](Index(x, y, z));
                }
                
                Voxel& operator[](const Index& index) {
                    Voxel& voxel = voxels[index];
                    if (!voxel.size())
                        initialise_voxel(voxel, index);
                    return voxel;
                }
                
                bool has_index(int x, int y, int z) {
                    return has_index(Index(x, y, z));
                }
                
                bool has_index(const Index& index) {
                    return (bool) voxels.count(index);
                }
                
                //Protected member functions
            protected:
                
                std::vector<Voxel*>& find_neighbourhood(Coord pos) {
                    size_t x_index = (int) (pos[X] - 0.5);
                    size_t y_index = (int) (pos[Y] - 0.5);
                    size_t z_index = (int) (pos[Z] - 0.5);
                    std::vector<Voxel*>& neigh = neigh_lookup[x_index][y_index][z_index];
                    if (!neigh.size())
                        initialise_neighbourhood(neigh, Index(x_index, y_index, z_index));
                    return neigh;
                }
                
                std::vector<Voxel*>& initialise_neighbourhood(std::vector<Voxel*>& neigh,
                                                              const Index& index) {
                    neigh.resize(neigh_size[X] * neigh_size[Y] * neigh_size[Z]);
                    size_t count = 0;
                    for (int x = -neigh_size[X]; x < neigh_size[X]; ++x)
                        for (int y = -neigh_size[Y]; y < neigh_size[Y]; ++y)
                            for (int z = -neigh_size[Z]; z < neigh_size[Z]; ++z)
                                neigh[count++] = operator[](index[X] + x, index[Y] + y,
                                        index[Z] + z);
                    return neigh;
                }
                
                Voxel& initialise_voxel(Voxel& voxel, const Index& index) {
                    Voxel* voxel_reference = 0;
                    if (reference)
                        voxel_reference = &(*reference)[index];
                    voxel.initialise(num_encodings, diff_model, interpolator, voxel_reference)
                    return voxel;
                }
                
        };
    
    }

}

#endif /* __bts_image2_generated_h__ */
