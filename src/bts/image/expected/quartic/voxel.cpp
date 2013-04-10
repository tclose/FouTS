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

#include "bts/image/expected/quartic/voxel.h"
#include "bts/image/expected/quartic/buffer.h"

#include "bts/fibre/strand/basic_section/tensor.h"
#include "bts/fibre/tractlet/section/tensor.h"

#include "bts/image/inline_functions.h"
#include "bts/diffusion/inline_functions.h"

namespace BTS {
    
    namespace Image {
        
        namespace Expected {
            
            namespace Quartic {
                
                Voxel::Voxel(Buffer& buffer, const Index& coord)
                        : Expected::Voxel(buffer, coord, buffer.diffusion_model), image(&buffer) {
                }
                
//----------------------//
// Non Hessian versions //
//----------------------//
                
                double Voxel::interpolate(const Coord& pos) {
                    
                    Coord disp = pos - this->centre();
                    
                    Coord interpolate;
                    
                    double interpolation;
                    
                    if (disp.lower_bounded(-1.0) && disp.upper_bounded(1.0)) {
                        
                        for (size_t dim_i = 0; dim_i < 3; ++dim_i)
                            interpolate[dim_i] = MR::Math::pow4(disp[dim_i])
                                    - 2.0 * MR::Math::pow2(disp[dim_i])
                                                 + 1;
                        
                        interpolation = interpolate[X] * interpolate[Y] * interpolate[Z];
                        
                    } else {
                        interpolation = 0.0;
                    }
                    
                    return interpolation;
                    
                }
                
                double Voxel::interpolate(const Coord& pos, Coord& pos_gradient) {
                    
                    Coord disp = pos - this->centre();
                    
                    Coord interpolate;
                    
                    double interpolation;
                    
                    if (disp.lower_bounded(-1.0) && disp.upper_bounded(1.0)) {
                        
                        for (size_t dim_i = 0; dim_i < 3; ++dim_i)
                            interpolate[dim_i] = MR::Math::pow4(disp[dim_i])
                                    - 2.0 * MR::Math::pow2(disp[dim_i])
                                                 + 1;
                        
                        interpolation = interpolate[X] * interpolate[Y] * interpolate[Z];
                        
#ifndef GRADIENT_NOT_REQUIRED
                        pos_gradient[X] = (4.0 * MR::Math::pow3(disp[X]) - 4.0 * disp[X])
                                * interpolate[Y] * interpolate[Z];
                        pos_gradient[Y] = (4.0 * MR::Math::pow3(disp[Y]) - 4.0 * disp[Y])
                                * interpolate[X] * interpolate[Z];
                        pos_gradient[Z] = (4.0 * MR::Math::pow3(disp[Z]) - 4.0 * disp[Z])
                                * interpolate[X] * interpolate[Y];
#endif //GRADIENT_NOT_REQUIRED not defined
                    } else {
                        interpolation = 0.0;
                        pos_gradient = Coord::Zeros;
                    }
                    
                    return interpolation;
                    
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
                
//------------------//
// Hessian versions //
//------------------//
                
                double Voxel::interpolate(const Coord& pos, Coord& gradient,
                                          Coord::Tensor& hessian) {
                    
                    Coord disp = pos - this->centre();
                    
                    Coord interpolate;
                    
                    double interpolation;
                    
                    if (disp.lower_bounded(-1.0) && disp.upper_bounded(1.0)) {
                        
                        for (size_t dim_i = 0; dim_i < 3; ++dim_i)
                            interpolate[dim_i] = MR::Math::pow4(disp[dim_i])
                                    - 2.0 * MR::Math::pow2(disp[dim_i])
                                                 + 1;
                        
                        interpolation = interpolate[X] * interpolate[Y] * interpolate[Z];
                        
                        Coord partial;
                        
                        partial[X] = (4.0 * MR::Math::pow3(disp[X]) - 4.0 * disp[X]);
                        partial[Y] = (4.0 * MR::Math::pow3(disp[Y]) - 4.0 * disp[Y]);
                        partial[Z] = (4.0 * MR::Math::pow3(disp[Z]) - 4.0 * disp[Z]);
                        
                        gradient[X] = partial[X] * interpolate[Y] * interpolate[Z];
                        gradient[Y] = interpolate[X] * partial[Y] * interpolate[Z];
                        gradient[Z] = interpolate[X] * interpolate[Y] * partial[Z];
                        
                        hessian[X][X] = (12.0 * MR::Math::pow2(disp[X]) - 4.0) * interpolate[Y]
                                        * interpolate[Z];
                        hessian[X][Y] = hessian[Y][X] = partial[X] * partial[Y] * interpolate[Z];
                        hessian[X][Z] = hessian[Z][X] = partial[X] * interpolate[Y] * partial[Z];
                        
                        hessian[Y][Y] = interpolate[X] * (12.0 * MR::Math::pow2(disp[Y]) - 4.0)
                                        * interpolate[Z];
                        hessian[Y][Z] = hessian[Z][Y] = interpolate[X] * partial[Y] * partial[Z];
                        
                        hessian[Z][Z] = interpolate[X] * interpolate[Y]
                                        * (12.0 * MR::Math::pow2(disp[Z]) - 4.0);
                        
                    } else {
                        interpolation = 0.0;
                        gradient.zero();
                        hessian.zero();
                    }
                    
                    return interpolation;
                    
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
                                          Fibre::Tractlet::Section& gradient,
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

}
