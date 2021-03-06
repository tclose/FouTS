/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Aug 30, 2010.

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

#include "image/header.h"

#include "bts/common.h"

#include "bts/image/expected/buffer.h"
#include "bts/image/expected/trilinear/buffer.h"
#include "bts/image/expected/gaussian/buffer.h"
#include "bts/image/expected/quartic/buffer.h"
#include "bts/image/expected/sinc/buffer.h"
#include "bts/image/expected/top_hat/buffer.h"
#include "bts/image/expected/reverse_sqrt/buffer.h"
#include "bts/image/expected/realistic/buffer.h"

#include "bts/image/voxel.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"

#include "bts/image/inline_functions.h"

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
        
        namespace Expected {
            
            const Triple<size_t> Buffer::DIMS_DEFAULT = Triple<size_t>(3, 3, 3);
            const Triple<double> Buffer::VOX_LENGTHS_DEFAULT = Triple<double>(0.15, 0.15, 0.15);
            const Triple<double> Buffer::OFFSETS_DEFAULT = Triple<double>::Zeros;
            const std::string Buffer::FILE_EXTENSION = "mif";
            
            const double Buffer::INTERP_EXTENT_DEFAULT = 1.0;
            const double Buffer::HALF_WIDTH_DEFAULT = 0.392470007505158;
            const char* Buffer::TYPE_DEFAULT = "sinc";
            const size_t Buffer::NUM_LENGTH_SECTIONS_DEFAULT = 15;
            const size_t Buffer::NUM_WIDTH_SECTIONS_DEFAULT = 4;
            const bool Buffer::ENFORCE_BOUNDS_DEFAULT = false;
            
            Buffer* Buffer::factory(const std::string& type, const Triple<size_t>& dims,
                                    const Triple<double>& vox_lengths,
                                    const Diffusion::Model& diffusion_model,
                                    size_t num_length_sections, size_t num_width_sections,
                                    double interp_extent, const Triple<double>& offsets,
                                    bool enforce_bounds, double gaussian_half_width) {
                
                Buffer* image;
                
                if (type == Trilinear::Buffer::SHORT_NAME)
                    image = new Trilinear::Buffer(dims, vox_lengths, diffusion_model,
                            num_length_sections, num_width_sections, interp_extent, offsets,
                            enforce_bounds);
                
                else if (type == Gaussian::Buffer::SHORT_NAME)
                    image = new Gaussian::Buffer(dims, vox_lengths, diffusion_model,
                            num_length_sections, num_width_sections, interp_extent,
                            gaussian_half_width, offsets, enforce_bounds);
                
                else if (type == Quartic::Buffer::SHORT_NAME)
                    
                    image = new Quartic::Buffer(dims, vox_lengths, diffusion_model,
                            num_length_sections, num_width_sections, interp_extent, offsets,
                            enforce_bounds);
                
                else if (type == Sinc::Buffer::SHORT_NAME)
                    
                    image = new Sinc::Buffer(dims, vox_lengths, diffusion_model,
                            num_length_sections, num_width_sections, interp_extent, offsets,
                            enforce_bounds);
                
                else if (type == TopHat::Buffer::SHORT_NAME)
                    
                    image = new TopHat::Buffer(dims, vox_lengths, diffusion_model,
                            num_length_sections, num_width_sections, interp_extent, offsets,
                            enforce_bounds);
                
                else if (type == ReverseSqrt::Buffer::SHORT_NAME)
                    
                    image = new ReverseSqrt::Buffer(dims, vox_lengths, diffusion_model,
                            num_length_sections, num_width_sections, interp_extent, offsets,
                            enforce_bounds);
                
                else if (type == Realistic::Buffer::SHORT_NAME)
                    
                    image = new Realistic::Buffer(dims, vox_lengths, diffusion_model,
                            num_length_sections, num_width_sections, interp_extent, offsets,
                            enforce_bounds);
                
                else
                    throw Exception(
                            "Unrecognised interpolation type '" + type
                            + "' passed to option '-exp_type'.");
                
                return image;
                
            }
            
            Image::Buffer Buffer::clean_buffer() {
                
                Image::Buffer clean_buffer(dims(), num_encodings());
                
                LOOP(clean_buffer(x,y,z) = operator()(x,y,z);)
                
                return clean_buffer;
                
            }
            
            bool Buffer::dims_match(const Observed::Buffer& reference) {
                
                return (reference.dims() == this->dims());
                
            }
        
            std::map<Index, Voxel>::iterator Buffer::begin() {
                throw Exception("not implemented");
            }

            std::map<Index, Voxel>::iterator Buffer::end() {
                throw Exception("not implemented");
            }
        }
    
    }

}
