/*
 Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

 Written by Thomas G Close, 5/05/09.

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

#ifndef __bts_math_svd_h__
#define __bts_math_svd_h__

extern "C" {
#include <gsl/gsl_linalg.h>
}

#include "math/matrix.h"
#include "math/vector.h"

#include "bts/common.h"

namespace BTS {
    
    namespace Math {
        
        class USV {
                
            public:
                
                MR::Math::Matrix<double> u;
                MR::Math::Matrix<double> v;
                MR::Math::Vector<double> sigma;

                USV(size_t rows, size_t columns)
                        : u(rows, columns), v(columns, columns), sigma(columns) {
                }
                ~USV() {
                }
                
        };
        
        // Protects M from being copied over by U in the gsl SVD computation.
        inline USV svd(const MR::Math::Matrix<double>& matrix) {
            
            USV usv(matrix.rows(), matrix.columns());
            
            usv.u = matrix;
            
            MR::Math::Vector<double> work(usv.sigma.size());
            
            gsl_linalg_SV_decomp(usv.u.gsl(), usv.v.gsl(), usv.sigma.gsl(), work.gsl());
            
            return usv;
            
        }
        
        // 'matrix' will be replaced by 'u' in the SVD computation.  This saves copying although the input matrix is destroyed.
        inline MR::Math::Vector<double>& svd(MR::Math::Matrix<double>& matrix,
                                             MR::Math::Matrix<double>& v,
                                             MR::Math::Vector<double>& sigma) {
            
            MR::Math::Vector<double> work(sigma.size());
            
            gsl_linalg_SV_decomp(matrix.gsl(), v.gsl(), sigma.gsl(), work.gsl());
            
            return sigma;
            
        }
    
    }

}

#endif
