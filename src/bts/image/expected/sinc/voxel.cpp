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

#include "bts/image/expected/sinc/voxel.h"
#include "bts/image/expected/sinc/buffer.h"

#include "bts/fibre/strand/basic_section/tensor.h"
#include "bts/fibre/tractlet/section/tensor.h"

#include "bts/image/inline_functions.h"
#include "bts/diffusion/inline_functions.h"

namespace FTS {
    
    namespace Image {
        
        namespace Expected {
            
            namespace Sinc {
                
                Voxel::Voxel(Buffer& buffer, const Index& coord)
                        : Expected::Voxel(buffer, coord, buffer.diffusion_model), image(&buffer) {
                }
                
                double Voxel::interpolate(const Coord& pos) {
                    
                    Coord disp = pos - centre();
                    
                    Coord interpolate;
                    
                    double interpolation;
                    
                    //Truncate the sinc function at a consitent distance from the voxel centre.
                    if (disp.lower_bounded(-image->get_extent()) && disp.upper_bounded(
                                image->get_extent())) {
                        
                        for (size_t dim_i = 0; dim_i < 3; ++dim_i) {
                            
                            if (disp[dim_i] == 0.0)
                                interpolate[dim_i] = 1.0;
                            else
                                interpolate[dim_i] = MR::Math::sin(M_PI * disp[dim_i])
                                        / (M_PI * disp[dim_i]);
                        }
                        
                        interpolation = interpolate[X] * interpolate[Y] * interpolate[Z];
                        
                    } else
                        interpolation = 0.0;
                    
                    return interpolation;
                    
                }
            
            }
        
        }
    
    }

}
