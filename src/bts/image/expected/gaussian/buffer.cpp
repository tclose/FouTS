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

#include "bts/image/expected/gaussian/buffer.h"
#include "bts/image/expected/gaussian/voxel.h"

#include "bts/fibre/strand/tensor.h"
#include "bts/fibre/tractlet/tensor.h"

#include "bts/image/inline_functions.h"
#include "bts/diffusion/inline_functions.h"

namespace FTS {
    
    namespace Image {
        
        namespace Expected {
            
            namespace Gaussian {
                
                const std::string Buffer::SHORT_NAME = "gaussian";
                
                Buffer::Buffer(const Triple<size_t>& dimensions, const Triple<double>& voxel_sizes,
                               const Diffusion::Model& diffusion_model, size_t num_sections,
                               size_t num_strands, double extent, double half_width,
                               const Triple<double>& corner_offset, bool enforce_bounds)
                        
                        : Buffer_tpl<Voxel>(dimensions, voxel_sizes, diffusion_model, num_sections,
                                  num_strands, extent, corner_offset, enforce_bounds), variance(
                                  half_width / (SQRT_LOG_2 * SQRT_2)) {
                    name_init();
                }
                
                Buffer::Buffer(const Buffer& buffer)
                        : Buffer_tpl<Voxel>(buffer), variance(buffer.variance) {
                    name_init();
                }
                
                std::ostream& operator<<(std::ostream& stream, const Buffer& buffer) {
                    return buffer.to_stream(stream);
                }
            
            }
        
        }
    
    }

}

