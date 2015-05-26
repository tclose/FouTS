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

#include "bts/image/expected/direction.h"

#ifndef __bts_image_expected_voxel_h__
#define __bts_image_expected_voxel_h__

namespace FTS {
    
    namespace Image {
        
        namespace Expected {
            
            class Voxel;
        
        }
    
    }

}

#include "bts/common.h"

#include "bts/triple.h"
#include "bts/coord.h"

#include "bts/image/expected/buffer.h"
#include "bts/image/observed/buffer.h"
#include "bts/image/observed/voxel.h"

#include "bts/image/expected/direction.h"

#include "bts/fibre/tractlet/section.h"
#include "bts/fibre/strand/basic_section/tensor.h"
#include "bts/fibre/strand/section.h"

#include "bts/diffusion/model.h"
#include "bts/image/expected/buffer.h"

namespace FTS {
    
    namespace Image {
        
        namespace Expected {
            
            class Voxel: public Observed::Voxel {
                    
                protected:
                    
                    std::vector<Direction> directions;

                    double precalc_interpolation;
                    Fibre::Tractlet::Section precalc_gradient;
                    Fibre::Tractlet::Section::Tensor precalc_hessian;

                public:
                    
                    Voxel() : precalc_interpolation(0.0) {
                    }
                    
                    Voxel(Expected::Buffer& exp_image, const Index& index,
                          Diffusion::Model& diffusion_model);

                    Voxel(const Voxel& v);

                    Voxel& operator=(const Voxel& v);

                    virtual ~Voxel() {
                    }
                    
                    const Diffusion::Encoding& encoding(size_t encode_index) const;

                    Voxel& zero() {
                        Image::Voxel<double>::zero();
                        precalc_interpolation = NAN;
                        precalc_gradient.invalidate(), precalc_hessian.invalidate();
                        return *this;
                    }
                    
                    Voxel& negate() {
                        Image::Voxel<double>::negate();
                        precalc_interpolation = NAN;
                        precalc_gradient.invalidate(), precalc_hessian.invalidate();
                        return *this;
                    }
                    
                    double precalculated_interpolation() const {
                        return precalc_interpolation;
                    }
                    
                    const Fibre::Tractlet::Section& precalculated_gradient() const {
                        return precalc_gradient;
                    }
                    
                    const Fibre::Tractlet::Section::Tensor& precalculated_hessian() const {
                        return precalc_hessian;
                    }
                    
                    void precalculate_interpolation(Fibre::Strand::Section& section) {
                        section.precalc_interpolation = interpolate(section);
                    }
                    
                    void precalculate_interpolation_gradient(Fibre::Strand::Section& section) {
                        section.precalc_interpolation = interpolate(section,
                                section.precalc_interp_gradient);
                    }
                    
                    void precalculate_interpolation_gradient_and_hessian(
                            Fibre::Strand::Section& section) {
                        section.precalc_interpolation = interpolate(section,
                                section.precalc_interp_gradient, section.precalc_interp_hessian);
                    }
                    
                    virtual double interpolate(const Fibre::Strand::BasicSection& section);

                    virtual double interpolate(const Fibre::Strand::BasicSection& section,
                                               Fibre::Strand::BasicSection& gradient);

                    virtual double interpolate(const Fibre::Strand::BasicSection& section,
                                               Fibre::Strand::BasicSection& gradient,
                                               Fibre::Strand::BasicSection::Tensor& hessian);

                    //Default to the standard Strand::Section interpolation.
                    virtual double interpolate(const Fibre::Tractlet::Section& section,
                                               Fibre::Strand::BasicSection& gradient) {
                        return interpolate((const Fibre::Strand::BasicSection&) section,
                                (Fibre::Strand::BasicSection&) gradient);
                    }
                    
                    //Default to the standard Strand::Section interpolation.
                    virtual double interpolate(const Fibre::Tractlet::Section& section,
                                               Fibre::Strand::BasicSection& gradient,
                                               Fibre::Tractlet::Section::Tensor& hessian);

                    Direction& direction(size_t index) {
                        return directions[index];
                    }
                    
                    Voxel operator-() const {
                        Voxel answer(*this);
                        for (size_t encode_i = 0; encode_i < num_encodings(); encode_i++)
                            answer[encode_i] = -operator[](encode_i);
                        return answer;
                    }
                    
                    //      protected:
                    
                    virtual double interpolate(const Coord& triple) {
                        throw Exception(
                                "If interpolate(const Fibre::Strand::BasicSection&,) is not overidden in derived class then interpolate(const Coord&) must be provided.");
                    }
                    
                    virtual double interpolate(const Coord& triple, Coord& gradient) {
                        throw Exception(
                                "If interpolate(const Fibre::Strand::BasicSection&, Fibre::Strand::BasicSection&) is not overidden in derived class then interpolate(const Coord&, Coord&) must be provided.");
                    }
                    
                    virtual double interpolate(const Coord& triple, Coord& gradient,
                                               Coord::Tensor& hessian) {
                        throw Exception(
                                "If interpolate(const Fibre::Strand::BasicSection&, Fibre::Strand::BasicSection&, Fibre::Strand::BasicSection::Tensor&) is not overidden in derived class then interpolate(const Coord&, Coord&, Coord::Tensor&) must be provided.");
                    }
                    
            };
        
        }
    
    }

}

#endif
