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

#include "bts/image/expected/gaussian/voxel.h"

#ifndef __bts_image_expected_gaussian_buffer_h__
#define __bts_image_expected_gaussian_buffer_h__

#include "bts/image/expected/buffer.h"
#include "bts/image/expected/buffer_tpl.h"

namespace BTS {
    
    namespace Image {
        
        namespace Expected {
            
            namespace Gaussian {
                
                class Buffer: public Expected::Buffer, public Expected::Buffer_tpl<Gaussian::Voxel> {
                        
                        //Public static constants
                    public:
                        
                        const static std::string SHORT_NAME;

                        //Protected member variables
                    protected:
                        
                        double variance;

                    public:
                        
                        Buffer()
                                : Buffer_tpl<Gaussian::Voxel>(this) {
                        }
                        
                        Buffer(const Triple<size_t>& dimensions, const Triple<double>& voxel_sizes,
                               const Diffusion::Model& diffusion_model, size_t num_sections =
                                       Expected::Buffer::NUM_LENGTH_SECTIONS_DEFAULT,
                               size_t num_strands = Expected::Buffer::NUM_WIDTH_SECTIONS_DEFAULT,
                               double interp_extent = Expected::Buffer::INTERP_EXTENT_DEFAULT,
                               double half_width = Expected::Buffer::HALF_WIDTH_DEFAULT,
                               const Triple<double>& corner_offset = Triple<double>(),
                               bool enforce_bounds = true);

                        Buffer(const Buffer& buffer);

                        ~Buffer() {
                        }
                        
                        double gauss_exp_var() {
                            return variance;
                        }
                        
                        EXPECTED_BUFFER_FUNCTIONS
                        ;
                        
                };
                
                std::ostream& operator<<(std::ostream& stream, const Buffer& buffer);
            
            }
        
        }
    
    }

}

#endif
