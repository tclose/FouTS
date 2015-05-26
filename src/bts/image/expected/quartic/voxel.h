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

#ifndef __bts_image_expected_quartic_voxel_h__
#define __bts_image_expected_quartic_voxel_h__

#include "bts/common.h"

#include "bts/triple.h"

#include "bts/image/expected/buffer.h"
#include "bts/image/observed/buffer.h"
#include "bts/image/expected/voxel.h"
#include "bts/image/expected/direction.h"

#include "bts/diffusion/model.h"

namespace FTS {
    
    namespace Image {
        
        namespace Expected {
            
            namespace Quartic {
                
                class Voxel: public Expected::Voxel {
                        
                    public:
                        
                        typedef Quartic::Buffer Buffer;

                    protected:
                        
                        Buffer* image;

                    public:
                        
                        Voxel() : image(NULL) {
                        }
                        
                        Voxel(Buffer& buffer, const Index& coord);

                        //! doesn't copy across the parent image.
                        Voxel(const Voxel& v)
                                : Expected::Voxel(v), image(v.image) {
                        }
                        
                        //! doesn't copy across the parent image.
                        Voxel& operator=(const Voxel& v) {
                            Expected::Voxel::operator=(v);
                            image = v.image;
                            return *this;
                        }
                        
                        ~Voxel() {
                        }
                        
                        void set_image(Buffer* image) {
                            this->image = image;
                        }
                        
                        double interpolate(const Fibre::Strand::BasicSection& section);

                        /*!The Strand::Section and Tractlet::Section interpolation methods are not strictly necessary as they copy virtual functions from the Expected::Voxel base class.
                         However they will allow the interpolate functions to be inlined and since this is an inner loop it could be critical.*/
                        double interpolate(const Fibre::Strand::BasicSection& section,
                                           Fibre::Strand::BasicSection& gradient);

                        double interpolate(const Fibre::Strand::BasicSection& section,
                                           Fibre::Strand::BasicSection& gradient,
                                           Fibre::Strand::BasicSection::Tensor& hessian);

                        double interpolate(const Fibre::Tractlet::Section& section,
                                           Fibre::Strand::BasicSection& gradient) {
                            return interpolate((const Fibre::Strand::BasicSection&) section,
                                    (Fibre::Strand::BasicSection&) gradient);
                        }
                        
                        double interpolate(const Fibre::Tractlet::Section& section,
                                           Fibre::Strand::BasicSection& gradient,
                                           Fibre::Tractlet::Section::Tensor& hessian);

                        //      protected:
                        
                        double interpolate(const Coord& pos);

                        double interpolate(const Coord& pos, Coord& gradient);

                        double interpolate(const Coord& triple, Coord& gradient,
                                           Coord::Tensor& hessian);
                        
                };
            
            }
        
        }
    
    }

}

#endif
