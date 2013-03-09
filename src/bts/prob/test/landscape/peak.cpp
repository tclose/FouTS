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

#include "bts/prob/test/landscape/peak.h"

#include "bts/prob/inline_functions.h"

#include "bts/math/common.h"
#include "math/matrix.h"
#include "bts/image/buffer.cpp.h"

namespace BTS {
    
    namespace Prob {
        
        namespace Test {
            
            Landscape::Peak::Peak(size_t ndims, const MCMC::State& centre,
                                  const MCMC::State::Tensor& widths, double height, size_t type)
                    : type(type), ndims(ndims), height(height), centre(centre), widths(widths)

            {
                
                if (ndims != widths.rows())
                    throw Exception(
                            "Size of centre (" + str(ndims) + ") and length of widths ("
                            + str(widths.rows()) + ") do not match.");
                
            }
            
            double Landscape::Peak::log_prob(const MCMC::State& point, MCMC::State& gradient) {
                
                double log_prob;
                
                if (type == GAUSSIAN)
                    log_prob = gaussian_log_prob(point, gradient);
                else if (type == PYRAMID)
                    log_prob = pyramid_log_prob(point, gradient);
                else
                    throw Exception("Unrecognised peak type, " + str(type) + ".");
                
                return log_prob;
                
            }
            
            double Landscape::Peak::log_prob(const MCMC::State& point, MCMC::State& gradient,
                                             MCMC::State::Tensor& hessian) {
                
                double log_prob;
                
                if (type == GAUSSIAN)
                    log_prob = gaussian_log_prob(point, gradient, hessian);
                else if (type == PYRAMID)
                    throw Exception("Hessian is not useful for pyramid type peaks.");
                else
                    throw Exception("Unrecognised peak type, " + str(type) + ".");
                
                return log_prob;
                
            }
            
            double Landscape::Peak::log_prob(const MCMC::State& point, MCMC::State& gradient,
                                             MCMC::State::Tensor& hessian,
                                             std::vector<MCMC::State::Tensor>& rank3_hessian) {
                
                double log_prob;
                
                if (type == GAUSSIAN)
                    log_prob = gaussian_log_prob(point, gradient, hessian, rank3_hessian);
                else if (type == PYRAMID)
                    throw Exception(
                            "Hessian or a 3rd order Hessian is not useful for pyramid type peaks.");
                else
                    throw Exception("Unrecognised peak type, " + str(type) + ".");
                
                return log_prob;
                
            }
            
            double Landscape::Peak::gaussian_log_prob(const MCMC::State& point,
                                                      MCMC::State& gradient,
                                                      MCMC::State::Tensor& hessian) {
                
                MCMC::State disp = point - this->centre;
                
                MCMC::State widths_mult_disp(ndims);
                
                MR::Math::mult(widths_mult_disp, widths, disp);
                
                double exponent = -MR::Math::dot(disp, widths_mult_disp);
                
                if (exponent > 0)
                    throw Exception("Tensor was not a positive definite matrix.");
                
                double log_prob = height * exp(exponent);
                
                gradient = widths_mult_disp;
                gradient *= -2.0 * log_prob;
                
                hessian = Math::outer(widths_mult_disp, widths_mult_disp);
                hessian *= 4.0 * log_prob;
                
                MR::Math::Matrix<double> hess_part = widths;
                hess_part *= -2.0 * log_prob;
                hessian += hess_part;
                
                return log_prob;
                
            }
            
            double Landscape::Peak::gaussian_log_prob(
                    const MCMC::State& point, MCMC::State& gradient, MCMC::State::Tensor& hessian,
                    std::vector<MCMC::State::Tensor>& rank3_hessian) {
                
                MCMC::State disp = point - this->centre;
                
                MCMC::State widths_mult_disp(ndims);
                
                MR::Math::mult(widths_mult_disp, widths, disp);
                
                double exponent = -MR::Math::dot(disp, widths_mult_disp);
                
                if (exponent > 0)
                    throw Exception("Tensor was not a positive definite matrix.");
                
                double log_prob = height * exp(exponent);
                
                gradient = widths_mult_disp;
                gradient *= -2.0 * log_prob;
                
                MR::Math::Matrix<double> outer_hess_part = Math::outer(widths_mult_disp,
                        widths_mult_disp);
                outer_hess_part *= 4.0 * log_prob;
                
                MR::Math::Matrix<double> width_hess_part = widths;
                width_hess_part *= -2.0 * log_prob;
                
                hessian = outer_hess_part;
                hessian += width_hess_part;
                
                for (size_t dim_i = 0; dim_i < ndims; dim_i++) {
                    
                    rank3_hessian.push_back(MCMC::State::Tensor(point.size()));
                    
                    MCMC::State::Tensor& d_hessian = rank3_hessian[dim_i];
                    
                    MR::Math::Matrix<double> outer_hess_part1 = outer_hess_part;
                    outer_hess_part1 *= -2.0 * widths_mult_disp[dim_i];
                    
                    MR::Math::Matrix<double> outer_hess_part2 = Math::outer(widths.column(dim_i),
                            widths_mult_disp);
                    outer_hess_part2 += Math::outer(widths_mult_disp, widths.column(dim_i));
                    outer_hess_part2 *= 4.0 * log_prob;
                    
                    MR::Math::Matrix<double> width_hess_part1 = width_hess_part;
                    width_hess_part1 *= -2.0 * widths_mult_disp[dim_i];
                    
                    d_hessian = outer_hess_part1;
                    d_hessian += outer_hess_part2;
                    d_hessian += width_hess_part1;
                    
                }
                
                return log_prob;
                
            }
            
            double Landscape::Peak::pyramid_log_prob(const MCMC::State& point,
                                                     MCMC::State& gradient) {
                
                MCMC::State disp = point - this->centre;
                
                gradient.resize(ndims);
                
                double log_prob = MR::Math::abs(height);
                
                for (size_t dim_i = 0; dim_i < ndims; dim_i++)
                    log_prob -= MR::Math::abs(disp[dim_i]) * widths(dim_i, dim_i);
                
                if (log_prob > 0) {
                    
                    for (size_t dim_i = 0; dim_i < disp.size(); dim_i++)
                        gradient[dim_i] = -widths(dim_i, dim_i) * Math::sign(disp[dim_i]);
                    
                    gradient *= Math::sign(height);
                    log_prob *= Math::sign(height);
                    
                } else {
                    log_prob = 0.0;
                    gradient.zero();
                }
                
                return log_prob;
            }
            
            std::ostream& operator<<(std::ostream& stream, const Landscape::Peak& peak) {
                
                stream << "centre: " << peak.centre << ", height: " << peak.height << " type: "
                << peak.type << ", widths: " << peak.widths << std::endl;
                
                return stream;
                
            }
        
        }
    
    }

}
