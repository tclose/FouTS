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

#ifndef __bts_image_container_buffer_h__
#define __bts_image_container_buffer_h__

namespace FTS {
    
    namespace Image {
        
        namespace Container {
            
            template<typename T> class Buffer;
        
        }
    
    }

}

#include <set>
#include "bts/image/container/voxel.h"
#include "bts/fibre/strand/section.h"
#include "bts/fibre/tractlet/section.h"
#include "bts/fibre/tractlet/tensor.h"
#include "bts/fibre/strand/tensor.h"

#define LOOP(op) \
for (size_t z = 0; z < this->dim(Z); z++) { \
  for (size_t y = 0; y < this->dim(Y); y++) { \
    for (size_t x = 0; x < this->dim(X); x++) { \
      op \
    } \
  } \
}

namespace FTS {
    
    namespace Image {
        
        namespace Container {
            
            template<typename T> class Buffer: public Image::Buffer_tpl<Container::Voxel<T> > {
                    
                public:
                    
                    typedef std::vector<Buffer> Set;

                protected:
                    
                    size_t number_directions;

                public:
                    
                    Buffer()
                            : Image::Buffer_tpl<Voxel<T> >(false), number_directions(0) {
                    }
                    
                    Buffer(Triple<size_t> dims, size_t num_encodings)
                            : Image::Buffer_tpl<Voxel<T> >(dims, false), number_directions(
                                      num_encodings) {
                    }
                    
                    Buffer(const Buffer& B)
                            : Image::Buffer_tpl<Voxel<T> >(B), number_directions(
                                      B.number_directions) {
                    }
                    
                    ~Buffer() {
                    }
                    
                    Buffer* clone() {
                        return new Buffer(*this);
                    }
                    
                    Buffer& operator=(const Buffer& B) {
                        Image::Buffer_tpl<Voxel<T> >::operator=(B);
                        number_directions = B.number_directions;
                        return *this;
                    }
                    
                    size_t num_encodings() const {
                        return number_directions;
                    }
                    
                    Buffer& clear() {
                        this->Image::Buffer_tpl<Voxel<T> >::clear();
                        return *this;
                    }
                    
                    Buffer& zero() {
                        this->Image::Buffer_tpl<Voxel<T> >::zero();
                        return *this;
                    }
                    
                    void reset(const Triple<size_t>& dimensions, size_t num_encodings) {
                        Image::Buffer_tpl<Voxel<T> >::reset(dimensions);
                        number_directions = num_encodings;
                    }
                    
                    void reset(const Triple<size_t>& dimensions, size_t num_encodings,
                               bool enforce_bounds) {
                        Image::Buffer_tpl<Voxel<T> >::reset(dimensions, enforce_bounds);
                        number_directions = num_encodings;
                    }
                    
                    void quick_clear() {
                        for (typename Image::Buffer_tpl<Voxel<T> >::iterator vox_it = this->begin();
                                vox_it != this->end(); ++vox_it)
                            vox_it->second.quick_empty();
                    }
                    
                    inline Voxel<T>& operator()(int x, int y, int z) {
                        return operator()(Index(x, y, z));
                    }
                    
                    inline Voxel<T>& operator()(const Index& c) {
                        Voxel<T>& vox = Image::Buffer_tpl<Voxel<T> >::operator()(c);
                        vox.undo_empty();
                        return vox;
                    }
                    
                    bool is_empty(const Index& c) const {
                        typename Image::Buffer_tpl<Voxel<T> >::const_iterator vox_it =
                                this->voxels.find(c);
                        return (vox_it == this->end()) ? true : vox_it->second.is_empty();
                    }
                    
                    bool is_empty(size_t x, size_t y, size_t z) const {
                        return is_empty(Index(x, y, z));
                    }
                    
                protected:
                    
                    Container::Voxel<T> new_voxel(const Index& coord) {
                        return Container::Voxel<T>(number_directions);
                    }
                    
            };
        
        }
    
    }

}

#undef LOOP

#endif //__bts_image_container_buffer_h__
