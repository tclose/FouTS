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

#ifndef __bts_image_container_voxel_h__
#define __bts_image_container_voxel_h__

namespace FTS {
    
    namespace Image {
        
        namespace Container {
            
            template<typename T> class Voxel;
        
        }
    
    }

}

#include "bts/image/voxel.h"

#define LOOP(op) \
for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++) {\
  op \
}

namespace FTS {
    
    namespace Image {
        
        namespace Container {
            
            template<typename T> class Voxel: public Image::Voxel<T> {
                    
                protected:
                    
                    bool fake_empty;

                public:
                    
                    Voxel()
                            : fake_empty(false) {
                    }
                    
                    Voxel(size_t num_encodings)
                            : Image::Voxel<T>(num_encodings), fake_empty(false) {
                    }
                    
                    Voxel(size_t num_encodings, const T& default_value)
                            : Image::Voxel<T>(num_encodings, default_value), fake_empty(false) {
                    }
                    
                    Voxel(const Voxel<T>& v)
                            : Image::Voxel<T>(v), fake_empty(false) {
                    }
                    
                    Voxel(const Image::Voxel<T>& v)
                            : Image::Voxel<T>(v), fake_empty(false) {
                    }
                    
                    Voxel& operator=(const Voxel<T>& v) {
                        this->Image::Voxel<T>::operator=(v);
                        fake_empty = v.fake_empty;
                        return *this;
                    }
                    
                    Voxel& operator*=(double M) {
                        assert(!fake_empty);
                        LOOP(this->operator[](encode_i) *= M;);
                        return *this;
                    }
                    
                    Voxel& operator/=(double M) {
                        assert(!fake_empty);
                        LOOP(this->operator[](encode_i) /= M;);
                        return *this;
                    }
                    
                    Voxel& zero() {
                        assert(!fake_empty);
                        this->Image::Voxel<T>::zero();
                        return *this;
                    }
                    
                    void quick_empty() {
                        fake_empty = true;
                    }
                    
                    bool is_empty() const {
                        return fake_empty;
                    }
                    
                    void undo_empty() {
                        fake_empty = false;
                    }
                    
                    void reset(const T& template_value) {
                        
                        assert(this->num_encodings());
                        
                        if (this->operator[](0).degree() != template_value.degree()) {
                            T zero_value(template_value);
                            zero_value.zero();
                            LOOP(this->operator[](encode_i) = zero_value;)
                        } else {
                            LOOP(this->operator[](encode_i).zero();)
                        }
                        
                    }
                    
            };
        
        }
    
    }

}

#undef LOOP

#endif
