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

#include "bts/image/expected/gaussian/voxel.h"
#include "bts/image/expected/gaussian/buffer.h"

#include "bts/image/inline_functions.h"
#include "bts/diffusion/inline_functions.h"

namespace FTS {
    
    namespace Image {
        
        namespace Expected {
            
            namespace Gaussian {
                
                Voxel::Voxel(Buffer& buffer, const Index& coord)
                        : Expected::Voxel(buffer, coord, buffer.diffusion_model), image(&buffer) {
                }
                
                double Voxel::interpolate(const Coord& pos) {
                    
                    Coord disp = pos - centre();
                    
                    Coord interpolate(
                            exp(-MR::Math::pow2(disp[X]) / (2.0 * image->gauss_exp_var())),
                            exp(-MR::Math::pow2(disp[Y]) / (2.0 * image->gauss_exp_var())),
                            exp(-MR::Math::pow2(disp[Z]) / (2.0 * image->gauss_exp_var())));
                    
                    interpolate /= sqrt(2.0 * M_PI * image->gauss_exp_var());
                    
                    return interpolate[X] * interpolate[Y] * interpolate[Z];
                    
                }
                
                double Voxel::interpolate(const Coord& pos, Coord& pos_gradient) {
                    
                    Coord disp = pos - centre();
                    
                    Coord interpolate(
                            exp(-MR::Math::pow2(disp[X]) / (2.0 * image->gauss_exp_var())),
                            exp(-MR::Math::pow2(disp[Y]) / (2.0 * image->gauss_exp_var())),
                            exp(-MR::Math::pow2(disp[Z]) / (2.0 * image->gauss_exp_var())));
                    
                    interpolate /= sqrt(2.0 * M_PI * image->gauss_exp_var());
                    
                    double scalar;
                    
                    scalar = interpolate[X] * interpolate[Y] * interpolate[Z];
                    
                    pos_gradient[X] = -disp[X] * scalar / image->gauss_exp_var();
                    pos_gradient[Y] = -disp[Y] * scalar / image->gauss_exp_var();
                    pos_gradient[Z] = -disp[Z] * scalar / image->gauss_exp_var();
                    
                    return scalar;
                    
                }
            
            }
        
        }
    
    }

}

