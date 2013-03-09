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

#include "bts/image/expected/voxel.h"
#include "bts/fibre/strand/basic_section/tensor.h"
#include "bts/fibre/tractlet/section/tensor.h"
#include "bts/coord/tensor.h"

#include "bts/image/inline_functions.h"

namespace BTS {
    
    namespace Image {
        
        namespace Expected {
            
            Voxel::Voxel(Expected::Buffer& exp_image, const Index& coord,
                         Diffusion::Model& diffusion_model)
                    : Observed::Voxel(diffusion_model.num_encodings(), coord), precalc_interpolation(
                              NAN) {
                
                for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++)
                    directions.push_back(Direction(*this, diffusion_model[encode_i]));
                
            }
            
            Voxel::Voxel(const Voxel& v)
                    : Observed::Voxel(v), precalc_interpolation(v.precalc_interpolation), precalc_gradient(
                              v.precalc_gradient) {
                
                for (size_t encode_i = 0; encode_i < v.directions.size(); encode_i++)
                    directions.push_back(Direction(*this, *v.directions[encode_i].response));
                
            }
            
            Voxel& Voxel::operator=(const Voxel& v) {
                
                Observed::Voxel::operator=(v);
                precalc_interpolation = v.precalc_interpolation;
                precalc_gradient = v.precalc_gradient;
                
                for (size_t encode_i = 0; encode_i < v.directions.size(); encode_i++)
                    directions.push_back(Direction(*this, *v.directions[encode_i].response));
                
                return *this;
                
            }
            
            const Diffusion::Encoding& Voxel::encoding(size_t encode_index) const {
                
                return directions[encode_index].encoding();
                
            }
            
            double Voxel::interpolate(const Fibre::Strand::BasicSection& section) {
                
                double pos_interp = interpolate(section.position());
                
                double tan_norm = section.tangent().norm();
                
                return pos_interp * tan_norm;
                
            }
            
            double Voxel::interpolate(const Fibre::Strand::BasicSection& section,
                                      Fibre::Strand::BasicSection& gradient) {
                
                Coord pos_interp_gradient;
                Coord tan_norm_gradient;
                
                double pos_interp = interpolate(section.position(), pos_interp_gradient);
                double tan_norm = section.tangent().norm(tan_norm_gradient);
                
                gradient.zero();
                
                gradient.position() = tan_norm * pos_interp_gradient;
                gradient.tangent() = pos_interp * tan_norm_gradient;
                
                return pos_interp * tan_norm;
                
            }
            
            double Voxel::interpolate(const Fibre::Strand::BasicSection& section,
                                      Fibre::Strand::BasicSection& gradient,
                                      Fibre::Strand::BasicSection::Tensor& hessian) {
                
                Coord pos_interp_gradient;
                Coord tan_norm_gradient;
                Coord::Tensor interp_hessian;
                Coord::Tensor norm_hessian;
                
                double pos_interp = interpolate(section.position(), pos_interp_gradient,
                        interp_hessian);
                double tan_norm = section.tangent().norm(tan_norm_gradient, norm_hessian);
                
                gradient.zero();
                
                gradient.position() = tan_norm * pos_interp_gradient;
                gradient.tangent() = pos_interp * tan_norm_gradient;
                
                hessian.zero();
                
                hessian.position(X).position() = tan_norm * interp_hessian[X];
                hessian.position(Y).position() = tan_norm * interp_hessian[Y];
                hessian.position(Z).position() = tan_norm * interp_hessian[Z];
                
                hessian.position(X).tangent() = pos_interp_gradient[X] * tan_norm_gradient;
                hessian.position(Y).tangent() = pos_interp_gradient[Y] * tan_norm_gradient;
                hessian.position(Z).tangent() = pos_interp_gradient[Z] * tan_norm_gradient;
                
                hessian.tangent(X).tangent() = pos_interp * norm_hessian[X];
                hessian.tangent(Y).tangent() = pos_interp * norm_hessian[Y];
                hessian.tangent(Z).tangent() = pos_interp * norm_hessian[Z];
                
                hessian.tangent(X).position() = pos_interp_gradient * tan_norm_gradient[X];
                hessian.tangent(Y).position() = pos_interp_gradient * tan_norm_gradient[Y];
                hessian.tangent(Z).position() = pos_interp_gradient * tan_norm_gradient[Z];
                
                return pos_interp * tan_norm;
                
            }
            
            double Voxel::interpolate(const Fibre::Tractlet::Section& section,
                                      Fibre::Strand::BasicSection& gradient,
                                      Fibre::Tractlet::Section::Tensor& hessian) {
                
                Fibre::Strand::BasicSection::Tensor strand_hessian;
                
                double interp = interpolate((const Fibre::Strand::BasicSection&) section,
                        (Fibre::Strand::BasicSection&) gradient, strand_hessian);
                
                hessian = strand_hessian;
                
                return interp;
                
            }
        
        }
    
    }

}
