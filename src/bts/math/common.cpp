/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Sep 13, 2010.

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

#include "bts/math/common.h"

namespace BTS {
    
    namespace Math {
        
        MR::Math::Matrix<double> random_rotation(gsl_rng* rand_gen, size_t ndims) {
            
            MR::Math::Matrix<double> rotate_matrix(ndims, ndims);
            
            for (size_t col_i = 0; col_i < ndims; col_i++) {
                
                MR::Math::Vector<double> ith_vector = rotate_matrix.column(col_i);
                
                for (size_t dim_i = 0; dim_i < ndims; dim_i++)
                    ith_vector[dim_i] = gsl_ran_gaussian(rand_gen, 1.0);
                
                MR::Math::normalise(ith_vector);
                
                for (size_t prev_col_i = 0; prev_col_i < col_i; prev_col_i++) {
                    
                    MR::Math::Vector<double> proj = rotate_matrix.column(prev_col_i);
                    
                    double proj_scalar = MR::Math::dot(ith_vector, proj);
                    
                    proj *= proj_scalar;
                    
                    ith_vector -= proj;
                    
                    MR::Math::normalise(ith_vector);
                    
                }
                
            }
            
            return rotate_matrix;
            
        }
        
        std::vector<bool> binary_string(size_t size, size_t index) {
            
            std::vector<bool> string(size);
            
            size_t mask = 1;
            
            for (size_t i = 0; i < size; ++i) {
                string[i] = index & mask;
                mask = mask << 1;
            }
            
            return string;
            
        }
        
        std::string matlab_str(const MR::Math::Vector<double>& v) {
            
            std::ostringstream stream;
            
            stream << "[";
            
            for (size_t i = 0; i < v.size() - 1; i++)
                stream << v[i] << " ";
            
            stream << v[v.size() - 1];
            
            stream << "]";
            
            return stream.str();
            
        }
        
        std::string matlab_str(const MR::Math::Matrix<double>& m) {
            
            std::ostringstream stream;
            stream << "[";
            
            for (size_t row_i = 0; row_i < m.rows(); row_i++) {
                for (size_t col_i = 0; col_i < m.columns(); col_i++) {
                    stream << m(row_i, col_i);
                    if (col_i != m.columns() - 1)
                        stream << ", ";
                }
                if (row_i != m.rows() - 1)
                    stream << ";";
            }
            
            stream << "]";
            
            return stream.str();
            
        }
        
        double median(const std::vector<double>& vector) {
            
            std::vector<double> v = vector;
            double median_value;
            size_t halfway_index = v.size() / 2;
            if (v.size() % 2 == 0) {
                std::nth_element(v.begin(), v.begin() + halfway_index, v.end());
                median_value = v[halfway_index];
            } else {
                std::nth_element(v.begin(), v.begin() + halfway_index, v.end());
                std::nth_element(v.begin(), v.begin() + halfway_index + 1, v.end());
                median_value = (v[halfway_index] + v[halfway_index + 1]) / 2.0;
            }
            return median_value;
        }
        
        double trace(const MR::Math::Matrix<double>& m) {
            
            assert(m.rows() == m.columns());
            
            double sum = 0.0;
            
            for (size_t i = 0; i < m.rows(); ++i)
                sum += m(i, i);
            
            return sum;
            
        }
    
    }

}
