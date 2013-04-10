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

#include "bts/image/expected/reverse_sqrt/buffer.h"
#include "bts/image/expected/reverse_sqrt/voxel.h"

#include "bts/fibre/strand/tensor.h"
#include "bts/fibre/tractlet/tensor.h"

#include "bts/image/inline_functions.h"
#include "bts/diffusion/inline_functions.h"

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
        
        namespace Expected {
            
            namespace ReverseSqrt {
                
                const std::string Buffer::SHORT_NAME = "reverse_sqrt";
                
                Buffer::Buffer(const Triple<size_t>& dimensions, const Triple<double>& voxel_sizes,
                               const Diffusion::Model& diffusion_model, size_t num_sections,
                               size_t num_strands, double extent,
                               const Triple<double>& corner_offset, bool enforce_bounds)
                        : Buffer_tpl<Voxel>(dimensions, voxel_sizes, diffusion_model, num_sections,
                                num_strands, extent, corner_offset, enforce_bounds) {
                    name_init();
                }
                
                Buffer::Buffer(const Buffer& buffer)
                        : Buffer_tpl<Voxel>(buffer) {
                    name_init();
                }
                
                std::ostream& operator<<(std::ostream& stream, const Buffer& buffer) {
                    return buffer.to_stream(stream);
                }
            
            }
        
        }
    
    }

}

#undef LOOP
