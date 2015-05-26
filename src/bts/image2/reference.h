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

#ifndef __bts_image2_reference_h__
#define __bts_image2_reference_h__

#include <map>

#include "bts/image2/index.h"
#include "bts/image2/voxel/reference.h"

namespace FTS {
    
    namespace Image2 {
        
        template<typename T> class Reference_tpl {
                
                //Public static variables, nested classes and typedefs
            public:
                
                //Protected member variables
            protected:
                
                size_t num_encodings;
                Triple<size_t> dims;
                Triple<double> voxel_size;
                Triple<double> offset;
                std::map<Index, T> voxels;
                T empty_voxel;

                //Public member functions
            public:
                
                Reference_tpl(size_t num_encodings, const Triple<size_t>& dims,
                              const Triple<double>& voxel_size, const Triple<double>& offset)
                        : num_encodings(num_encodings), dims(dims), voxel_size(voxel_size), offset(
                                  offset), empty_voxel(num_encodings, 0.0) {
                }
                
                T& operator[](int x, int y, int z) const {
                    return operator[](Index(x, y, z));
                }
                
                T& operator[](const Index& index) const {
                    typename std::map<Index, T>::iterator voxel_it = voxels.find(index);
                    if (voxel_it != voxels.end())
                        return *voxel_it;
                    else
                        return empty_voxel;
                }
                
                //Protected member functions
            protected:
                
        };
        
        class Reference: public Reference_tpl<Voxel::Reference> {
                
                //Public static variables, nested classes and typedefs
            public:
                
                //Protected member variables
            protected:
                
                //Public member functions
            public:
                
                //Protected member functions
            protected:
                
        };
    
    }

}

#endif /* __bts_image2_reference_h__ */
