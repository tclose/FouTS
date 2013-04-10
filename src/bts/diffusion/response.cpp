/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Aug 12, 2010.

 This file is part of MRtrix.

 MRtrix is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 MRtrix is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef __bts_diffusion_response_cpp_h__
#define __bts_diffusion_response_cpp_h__

#include "bts/diffusion/response.h"
#include "bts/coord/tensor.h"

namespace BTS {
    
    namespace Diffusion {
        
        //TODO: Encapsulate logic into inline functions to avoid replication between weighting functions.
        
        double Response::weighting(const Coord& tangent) const {
            
            //Initalize the weighting with the 0th order coefficient.
            double weight = coeffs[0];
            
            double tan_dot_orient = Encoding::orient.dot(tangent);
            double tan_dot_tan = tangent.norm2();
            
            //Otherwise the total added value = 0 anyway.
            if (tan_dot_orient != 0.0 && tan_dot_tan != 0.0) {
                
                double tan_dot_orient2 = MR::Math::pow2(tan_dot_orient);
                
                double cos2 = tan_dot_orient2 / tan_dot_tan;
                
                //Will hold the even powers of the cosine of the angle between DW and tangent.
                double cos2p = 1.0;
                
                for (size_t p_int = 1; p_int < coeffs.size(); ++p_int) {
                    
                    const double& coeff_p = coeffs[p_int];
                    
                    cos2p *= cos2;
                    
                    weight += coeff_p * cos2p;
                    
                }
                
            }
            
            return weight;
            
        }
        
        double Response::weighting(const Coord& tangent, Coord& gradient) const {
            
            //Initalize the weighting with the 0th order coefficient.
            double weight = coeffs[0];
            
#ifndef GRADIENT_NOT_REQUIRED
            gradient.zero();    //coeffs[0]; Not sure why this was coeffs[0].
#endif
            
            double tan_dot_orient = Encoding::orient.dot(tangent);
            double tan_dot_tan = tangent.norm2();
            
            //Otherwise the total added value = 0 anyway.
            if (tan_dot_orient != 0.0 && tan_dot_tan != 0.0) {
                
                double tan_dot_orient2 = MR::Math::pow2(tan_dot_orient);
                
                double cos2 = tan_dot_orient2 / tan_dot_tan;
                
                //Will hold the even powers of the cosine of the angle between DW and tangent.
                double cos2p = 1.0;
                
                for (size_t p_int = 1; p_int < coeffs.size(); ++p_int) {
                    
                    const double& coeff_p = coeffs[p_int];
                    
                    double p = (double) p_int;
                    
                    cos2p *= cos2;
                    
                    weight += coeff_p * cos2p;
                    
#ifndef GRADIENT_NOT_REQUIRED
                    gradient += coeff_p * 2 * p * cos2p
                                * (Encoding::orient / tan_dot_orient - tangent / tan_dot_tan);
#endif
                    
                }
                
            }
            
            return weight;
            
        }
        
        double Response::weighting(const Coord& tangent, Coord& gradient,
                                   Coord::Tensor& hessian) const {
            
            //Initalize the weighting with the 0th order coefficient.
            double weight = coeffs[0];
            
#ifndef GRADIENT_NOT_REQUIRED
            gradient.zero();    //coeffs[0]; Not sure why this was coeffs[0].
            
#ifndef HESSIAN_NOT_REQUIRED
            hessian.zero();
#endif
            
#endif
            
            double tan_dot_orient = Encoding::orient.dot(tangent);
            double tan_dot_tan = tangent.norm2();
            
            //Otherwise the total added value = 0 anyway.
            if (tan_dot_orient != 0.0 && tan_dot_tan != 0.0) {
                
                double tan_dot_orient2 = MR::Math::pow2(tan_dot_orient);
                
#ifndef HESSIAN_NOT_REQUIRED
                double tan_dot_tan2 = MR::Math::pow2(tan_dot_tan);
#endif
                
                double cos2 = tan_dot_orient2 / tan_dot_tan;
                
                //Will hold the even powers of the cosine of the angle between DW and tangent.
                double cos2p = 1.0;
                
                for (size_t p_int = 1; p_int < coeffs.size(); ++p_int) {
                    
                    const double& coeff_p = coeffs[p_int];
                    
                    double p = (double) p_int;
                    
                    cos2p *= cos2;
                    
                    weight += coeff_p * cos2p;
                    
#ifndef GRADIENT_NOT_REQUIRED
                    
                    gradient += coeff_p * 2 * p * cos2p
                                * (Encoding::orient / tan_dot_orient - tangent / tan_dot_tan);
                    
#ifndef HESSIAN_NOT_REQUIRED
                    
                    //Note that the Hessian calculation could be slightly more efficient by only calculating the upper triangle
                    //of the matrix but has been implemented this way for mathematical clarity.
                    hessian +=
                            coeff_p * 2
                            * p
                            * cos2p
                            * (((2 * p - 1) / tan_dot_orient2) * orient.outer(orient) - (2 * p
                                    / (tan_dot_orient * tan_dot_tan))
                                                                                        * (tangent.outer(
                                                                                                   orient)
                                                                                           + orient.outer(
                                                                                                   tangent))
                               + ((2 * p + 2) / tan_dot_tan2) * tangent.outer(tangent)
                               - Coord::Tensor::create_diagonal(1 / tan_dot_tan));
                    
#endif
                    
#endif
                    
                }
                
            }
            
            return weight;
            
        }
        
        std::ostream& operator<<(std::ostream& stream, const Diffusion::Response& response) {
            
            stream << "b=" << response.b;
            
            stream << ", orient=" << response.orient;
            
            stream << ", coeffs=" << response.coeffs;
            
            stream << std::endl;
            
            return stream;
            
        }
    
    }

}

#endif /* __bts_diffusion_response_cpp_h__ */
