/*
 Copyright 2010 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close on May 19, 2010.

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

#include "bts/image/expected/realistic/voxel.h"
#include "bts/image/expected/realistic/buffer.h"

#include "bts/image/inline_functions.h"
#include "bts/diffusion/inline_functions.h"

namespace FTS {
    
    namespace Image {
        
        namespace Expected {
            
            namespace Realistic {
                
                Voxel::Voxel(Buffer& buffer, const Index& coord)
                        : Expected::Voxel(buffer, coord, buffer.diffusion_model), image(&buffer) {
                }
                
                double Expected::Realistic::Voxel::interpolate(const Coord& pos) {
                    
                    Coord disp = pos - this->centre();
                    
                    Coord interpolate;
                    
                    double interpolation;
                    
                    //Truncate the sinc function at a consitent distance from the voxel centre.
                    if (disp.lower_bounded(-image->get_extent()) && disp.upper_bounded(
                                image->get_extent())) {
                        
                        //Do X and Y as a Sinc
                        for (size_t dim_i = 0; dim_i < 2; ++dim_i) {
                            
                            if (disp[dim_i] == 0.0)
                                interpolate[dim_i] = 1.0;
                            else
                                interpolate[dim_i] = MR::Math::sin(M_PI * disp[dim_i])
                                        / (M_PI * disp[dim_i]);
                        }
                        
                        if ((disp[Z] < 1) && (disp[Z] > -1))
                            interpolate[Z] = (MR::Math::pow4(disp[Z])
                                    - 2.0 * MR::Math::pow2(disp[Z])
                                              + 1);
                        else
                            interpolate[Z] = 0.0;
                        
                        interpolation = interpolate[X] * interpolate[Y] * interpolate[Z];
                        
                    } else
                        interpolation = 0.0;
                    
                    return interpolation;
                    
                }
                
                double Expected::Realistic::Voxel::interpolate(const Coord& pos,
                                                               Coord& pos_gradient) {
                    
                    throw Exception(
                            "Expected::Realistic::Voxel::interpolate(const Coord&, Coord&) not implemented yet.");
                    
                }
            
            }
        
        }
    
    }

}

