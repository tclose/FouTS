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

#ifndef __bts_image_expected_trilinear_voxel_h__
#define __bts_image_expected_trilinear_voxel_h__

#include "bts/common.h"

#include "bts/triple.h"

#include "bts/image/expected/buffer.h"
#include "bts/image/observed/buffer.h"
#include "bts/image/expected/voxel.h"
#include "bts/image/expected/direction.h"

#include "bts/diffusion/model.h"

namespace BTS {
    
    namespace Image {
        
        namespace Expected {
            
            namespace Trilinear {
                
                class Voxel: public Expected::Voxel {
                        
                    public:
                        
                        typedef Trilinear::Buffer Buffer;

                    protected:
                        
                        Buffer* image;

                    public:
                        
                        Voxel() {
                        }
                        
                        Voxel(Buffer& buffer, const Index& coord);

                        Voxel(const Voxel& voxel)
                                : Expected::Voxel(voxel), image(voxel.image) {
                        }
                        
                        Voxel& operator=(const Voxel& v) {
                            Expected::Voxel::operator=(v);
                            image = v.image;
                            return *this;
                        }
                        
                        ~Voxel() {
                        }
                        
                        //      protected:
                        
                        double interpolate(const Coord& pos);

                        double interpolate(const Coord& pos, Coord& gradient);

                        double interpolate(const Coord& triple, Coord& gradient,
                                           Coord::Tensor& hessian) {
                            throw Exception("Gradient is not defined for Top-hat interpolation.");
                        }
                        
                };
            
            }
        
        }
    
    }

}

#endif
