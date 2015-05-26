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

#include "bts/image/expected/trilinear/voxel.h"
#include "bts/image/expected/trilinear/buffer.h"

#include "bts/image/inline_functions.h"
#include "bts/diffusion/inline_functions.h"

namespace FTS {
    
    namespace Image {
        
        namespace Expected {
            
            namespace Trilinear {
                
                Voxel::Voxel(Buffer& buffer, const Index& coord)
                        : Expected::Voxel(buffer, coord, buffer.diffusion_model), image(&buffer) {
                }
                
                double Expected::Trilinear::Voxel::interpolate(const Coord& pos) {
                    
                    Coord disp = pos - centre();
                    Coord interpolate = Coord::Ones - disp.abs();
                    
                    double scalar;
                    
                    if (interpolate.non_negative()) {
                        
                        scalar = interpolate[X] * interpolate[Y] * interpolate[Z];
                        
                    } else {
                        
                        scalar = 0.0;
                    }
                    
                    return scalar;
                    
                }
                
                double Expected::Trilinear::Voxel::interpolate(const Coord& pos,
                                                               Coord& pos_gradient) {
                    
                    Coord disp = pos - centre();
                    
                    Coord interpolate = Coord::Ones - disp.abs();
                    
#ifndef GRADIENT_NOT_REQUIRED
                    Coord disp_sign = disp.sign();
#endif
                    
                    double scalar;
                    
                    if (interpolate.non_negative()) {
                        
                        scalar = interpolate[X] * interpolate[Y] * interpolate[Z];
                        
#ifndef GRADIENT_NOT_REQUIRED
                        pos_gradient[X] = -disp_sign[X] * interpolate[Y] * interpolate[Z];
                        pos_gradient[Y] = -disp_sign[Y] * interpolate[X] * interpolate[Z];
                        pos_gradient[Z] = -disp_sign[Z] * interpolate[X] * interpolate[Y];
#endif
                        
                    } else {
                        
                        scalar = 0.0;
                        
#ifndef GRADIENT_NOT_REQUIRED
                        pos_gradient.zero();
#endif
                        
                    }
                    
                    return scalar;
                    
                }
            
            }
        
        }
    
    }

}

