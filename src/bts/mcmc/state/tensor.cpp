/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 11/08/2010.

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

#include "bts/mcmc/state/tensor.h"

namespace FTS {
    
    namespace MCMC {
        
        State::Tensor::Tensor(const State& template_value)
                : MR::Math::Matrix<double>(template_value.size(), template_value.size()) {
            
            for (size_t row_i = 0; row_i < template_value.size(); ++row_i)
                this->row(row_i) = template_value;
            
        }
        
        State::Tensor& State::Tensor::invalidate() {
            
            MR::Math::Matrix<double>::operator=(NAN);
            
            return *this;
            
        }
        
        State::Tensor& State::Tensor::zero() {
            
            MR::Math::Matrix<double>::operator=(0.0);
            
            return *this;
            
        }
    
    }

}
