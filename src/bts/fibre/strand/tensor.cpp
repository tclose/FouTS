/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 17/08/2010.

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

#include "bts/fibre/strand/basic_section.h"
#include "bts/fibre/strand/section.h"
#include "bts/fibre/strand/tensor.h"
#include "bts/fibre/strand/basic_section/tensor.h"

namespace FTS {
    
    namespace Fibre {
        
        void Strand::Tensor::add_section_hessian(const Strand& strand, const Section& section,
                                                 const BasicSection& gradient,
                                                 const BasicSection::Tensor& hessian) {
            
            Fibre::Strand acs_inc = (section.position_coeffs.size(), 0.0);
            
            acs_inc.add_section_gradient(strand, section, gradient);
            
            if (has_var_acs()) {
                acs() += acs_inc;
                acs().set_acs(0.0);
            }
            
            for (size_t row_degree_i = 0; row_degree_i < degree(); row_degree_i++)
                for (size_t row_dim_i = 0; row_dim_i < 3; row_dim_i++) {
                    for (size_t col_degree_i = 0; col_degree_i < degree(); col_degree_i++)
                        for (size_t col_dim_i = 0; col_dim_i < 3; col_dim_i++) {
                            double& elem = operator()(row_degree_i, row_dim_i, col_degree_i,
                                    col_dim_i);
                            elem += section.position_coeffs[row_degree_i] * section.position_coeffs[col_degree_i]
                                    * hessian.position(row_dim_i).position(col_dim_i)
                                    + section.position_coeffs[row_degree_i] * section.tangent_coeffs[col_degree_i]
                                      * hessian.position(row_dim_i).tangent(col_dim_i)
                                      * section.length_fraction
                                    + section.tangent_coeffs[row_degree_i] * section.position_coeffs[col_degree_i]
                                      * hessian.tangent(row_dim_i).position(col_dim_i)
                                      * section.length_fraction
                                    + section.tangent_coeffs[row_degree_i] * section.tangent_coeffs[col_degree_i]
                                      * hessian.tangent(row_dim_i).tangent(col_dim_i)
                                      * MR::Math::pow2(section.length_fraction);
                            
                        }
                    
                    if (has_var_acs())
                        row(row_degree_i * 3 + row_dim_i).var_acs() +=
                                acs_inc[row_degree_i][row_dim_i];
                }
            
        }
    
    }

}

