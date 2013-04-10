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

#include "bts/prob/prior_component/frequency.h"

#include "bts/prob/inline_functions.h"

namespace BTS {
    
    namespace Prob {
        
        namespace PriorComponent {
            
            const double Frequency::SCALE_DEFAULT = 5.0;
            const double Frequency::AUX_SCALE_DEFAULT = 10.0;
            const std::string Frequency::NAME = "frequency";
            
            double Frequency::log_prob(const Fibre::Strand strand, Fibre::Strand gradient) {
                
                double lprob = 0.0;
                gradient.zero();
                
                const Triple<double>& v1 = strand[1];
                
                double v1_norm2 = v1.norm2();
                
                for (size_t degree_i = 2; degree_i < strand.degree(); degree_i++) {
                    
                    const Triple<double>& vm = strand[degree_i];
                    double vm_norm2 = vm.norm2();
                    
                    if (v1_norm2) {
                        
                        lprob -= scale * MR::Math::pow2(degree_i) * vm_norm2 / v1_norm2;
                        
#ifndef GRADIENT_NOT_REQUIRED
                        //The gradient of the curve magnitude w.r.t. Vm.
                        gradient[degree_i] -= vm * 2.0 * MR::Math::pow2(degree_i) * scale
                                / v1_norm2;
                        
                        //The gradient of the curve magnitude w.r.t. V1.
                        gradient[1] += v1 * 2.0 * MR::Math::pow2(degree_i) * scale * vm_norm2
                                / (v1_norm2 * v1_norm2);
#endif
                        
                    } else if (vm_norm2)    // If v1 is zero and vm is non zero set the probability to be very low
                        lprob = -LARGE_FLOAT;
                }
                
                return lprob;
                
            }
            
//      double Frequency::log_prob(const Fibre::Strand& strand, Fibre::Strand& gradient, Fibre::Strand::Tensor& hessian) {
//
//        double lprob = 0.0;
//        gradient.zero();
//
//        hessian.reset(strand);
//
//        const Triple<double>& v1 = strand[1];
//
//        double v1_norm2 = v1.norm2();
//
//        for (size_t degree_i = 2; degree_i < strand.degree(); degree_i++) {
//
//          const Triple<double>& vm = strand[degree_i];
//          double vm_norm2 = vm.norm2();
//
//          lprob -= scale * degree_i * vm_norm2 / v1_norm2;
//
//  #ifndef GRADIENT_NOT_REQUIRED
//          //The gradient of the curve magnitude w.r.t. Vm.
//          gradient[degree_i] -= vm * 2.0 * degree_i * scale / v1_norm2;
//
//          //The gradient of the curve magnitude w.r.t. V1.
//          gradient[1] +=  v1 * 2.0 * scale * degree_i * vm_norm2/(v1_norm2 * v1_norm2);
//
//  #ifndef HESSIAN_NOT_REQUIRED
//
//          for (size_t dim_i1 = 0; dim_i1 < 3; dim_i1++) {
//
//            hessian(degree_i,dim_i1,degree_i,dim_i1) -= 2.0 * degree_i * scale / v1_norm2;
//            hessian(1, dim_i1, 1, dim_i1) += 2.0 * scale * degree_i * vm_norm2 / MR::Math::pow2(v1_norm2);
//
//            for (size_t dim_i2 = 0; dim_i2 < 3; dim_i2++) {
//
//              hessian(degree_i, dim_i1, 1, dim_i2) += 4.0 * scale * degree_i * vm[dim_i1] * v1[dim_i2] / MR::Math::pow2(v1_norm2);
//              hessian(1, dim_i1, degree_i, dim_i2) += 4.0 * scale * degree_i * v1[dim_i1] * vm[dim_i2] / MR::Math::pow2(v1_norm2);
//
//              hessian(1, dim_i1, 1, dim_i2) -= 8.0 * scale * degree_i * v1[dim_i1] * v1[dim_i2] * vm_norm2 / MR::Math::pow3(v1_norm2);
//            }
//
//
//          }
//
//
//  #endif
//
//  #endif
//
//        }
//
//
//        return lprob;
//
//      }
            
            double Frequency::log_prob(const Fibre::Tractlet tractlet, Fibre::Tractlet gradient) {
                
                gradient.zero();
                double lprob = log_prob(tractlet[0], gradient[0]);
                
                const Coord v1 = tractlet(0, 1);
                
                double v1_norm2 = v1.norm2();
                
                for (size_t ax_i = 1; ax_i < 3; ++ax_i) {
                    for (size_t degree_i = 1; degree_i < tractlet.degree(); degree_i++) {
                        
                        const Coord vm = tractlet(ax_i, degree_i);
                        double vm_norm2 = vm.norm2();
                        
                        if (v1_norm2) {
                            
                            lprob -= scale * MR::Math::pow2(degree_i) * vm_norm2 / v1_norm2;
                            
#ifndef GRADIENT_NOT_REQUIRED
                            //The gradient of the curve magnitude w.r.t. Vm.
                            gradient(ax_i, degree_i) -= vm * 2.0 * MR::Math::pow2(degree_i) * scale
                                    / v1_norm2;
                            
                            //The gradient of the curve magnitude w.r.t. V1.
                            gradient(0, 1) += v1 * 2.0 * MR::Math::pow2(degree_i) * scale * vm_norm2
                                    / (v1_norm2 * v1_norm2);
#endif
                            
                        } else if (vm_norm2)    // If v1 is zero and vm is non zero set the probability to be very low
                            lprob = -LARGE_FLOAT;
                    }
                }
                
                return lprob;
                
            }
        
//
//      double Frequency::log_prob_and_fisher(const Fibre::Strand& strand, Fibre::Strand& gradient, Fibre::Strand::Tensor& fisher) {
//
//        double log_prob = 0.0;
//        gradient.zero();
//
//        fisher.reset(strand);
//
//        const Coord v1 = strand[1];
//
//        double v1_norm2 = v1.norm2();
//
//        for (size_t degree_i = 2; degree_i < strand.degree(); degree_i++) {
//
//          Coord vm = strand[degree_i];
//          double vm_norm2 = vm.norm2();
//
//          log_prob -= degree_i * scale * vm_norm2 / v1_norm2;
//
//  #ifndef GRADIENT_NOT_REQUIRED
//          //The gradient of the curve magnitude w.r.t. Vm.
//          gradient[degree_i] -= vm * 2.0 * degree_i * scale / v1_norm2;
//
//          //The gradient of the curve magnitude w.r.t. V1.
//          gradient[1] +=  v1 * 2.0 * degree_i * scale * vm_norm2/(v1_norm2 * v1_norm2);
//
//          for (size_t dim_i1 = 0; dim_i1 < 3; dim_i1++)
//            fisher(degree_i, dim_i1, degree_i, dim_i1) += 2.0 * degree_i * scale / v1_norm2;
//
//  #endif
//
//        }
//
//
//        return log_prob;
//
//      }
//
//
//      double Frequency::log_prob_and_fisher(const Fibre::Strand& strand, Fibre::Strand& gradient, Fibre::Strand::Tensor& fisher, std::vector<Fibre::Strand::Tensor>& fisher_gradient) {
//
//        double lprob = 0.0;
//        gradient.zero();
//
//        size_t vector_size = strand.vsize();
//
//        fisher.reset(strand);
//
//        //TODO: Make reset function to do this
//        fisher_gradient.resize(vector_size);
//        for (size_t elem_i = 0; elem_i < vector_size; ++elem_i)
//          fisher_gradient[elem_i].reset(strand);
//
//        const Coord v1 = strand[1];
//
//        double v1_norm2 = v1.norm2();
//
//        for (size_t degree_i = 2; degree_i < strand.degree(); degree_i++) {
//
//          const Triple<double>& vm = strand[degree_i];
//          double vm_norm2 = vm.norm2();
//
//          lprob -= degree_i * scale * vm_norm2 / v1_norm2;
//
//  #ifndef GRADIENT_NOT_REQUIRED
//          //The gradient of the curve magnitude w.r.t. Vm.
//          gradient[degree_i] -= vm * 2.0 * degree_i * scale / v1_norm2;
//
//          //The gradient of the curve magnitude w.r.t. V1.
//          gradient[1] +=  v1 * 2.0 * degree_i * scale * vm_norm2/(v1_norm2 * v1_norm2);
//
//          for (size_t dim_i1 = 0; dim_i1 < 3; dim_i1++) {
//
//            fisher(degree_i, dim_i1, degree_i, dim_i1) += 2.0 * degree_i * scale / v1_norm2;
//
//            //TODO: Make this independent of whether the strand has an acs;
//            for (size_t dim_i2 = 0; dim_i2 < 3; dim_i2++)
//              fisher_gradient[dim_i2 + 3 + strand.has_var_acs()](degree_i, dim_i1, degree_i, dim_i1) -= 4.0 * degree_i * scale * v1[dim_i2] / MR::Math::pow2(v1_norm2);
//
//          }
//
//
//  #endif
//
//        }
//
//
//        return lprob;
//
//      }
        
        }
    
    }

}
