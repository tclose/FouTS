/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

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

#ifndef __bts_image_reference_buffer_h__
#define __bts_image_reference_buffer_h__

namespace BTS {
    
    namespace Image {
        
        namespace Reference {
            
            template<typename T> class Buffer;
        
        }
    
    }

}

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

namespace BTS {
    
    namespace Image {
        
        namespace Reference {
            
            template<typename T> class Buffer: public Image::Buffer_tpl<std::vector<T*> > {
                    
                public:
                    
                    typedef std::map<size_t, Buffer<T> > Set;

//          class Set : public std::map<Index, Buffer<T> > {
//
//            Set() : std::map<Index, Buffer<T> >() {}
//            Set(const std::map<Index,Buffer<T> >& m) : std::map<Index,Buffer<T> >(m) {}
//            ~Set() {};
//
//          };
                    
                public:
                    
                    Buffer()
                            : Image::Buffer_tpl<std::vector<T*> >(false) {
                    }
                    
                    Buffer(const Buffer& B)
                            : Image::Buffer_tpl<std::vector<T*> >(B) {
                    }
                    
                    ~Buffer() {
                    }
                    
                    Buffer* clone() {
                        return new Buffer(*this);
                    }
                    
                    Buffer& operator=(const Buffer& B) {
                        Image::Buffer_tpl<std::vector<T*> >::operator=(B);
                        return *this;
                    }
                    
                    void clear_references() {
                        for (typename Image::Buffer_tpl<std::vector<T*> >::iterator vox_it =
                                this->begin(); vox_it != this->end(); ++vox_it)
                            vox_it->second.clear();
                    }
                    
            };
        
        }
    
    }

}

#undef LOOP

#endif //__bts_image_reference_buffer_h__
