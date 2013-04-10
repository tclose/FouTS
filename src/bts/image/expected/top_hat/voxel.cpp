/*
 Copyright 2010 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close on May 19, 2010.

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

#include "bts/image/expected/top_hat/voxel.h"
#include "bts/image/expected/top_hat/buffer.h"

#include "bts/image/inline_functions.h"
#include "bts/diffusion/inline_functions.h"

namespace BTS {
    
    namespace Image {
        
        namespace Expected {
            
            namespace TopHat {
                
                Voxel::Voxel(Buffer& buffer, const Index& index)
                        : Expected::Voxel(buffer, index, buffer.diffusion_model), image(&buffer) {
                }
                
                double Voxel::interpolate(const Coord& pos) {
                    
                    Coord disp = pos - centre();
                    
                    double scalar;
                    
                    if (disp.upper_bounded(0.5) && disp.lower_bounded(-0.5))
                        scalar = 1.0;
                    else
                        scalar = 0.0;
                    
                    return scalar;
                    
                }
            
            }
        
        }
    
    }

}

