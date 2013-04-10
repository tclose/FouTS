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

#ifndef __bts_mcmc_state_tensor_h__
#define __bts_mcmc_state_tensor_h__

#include "bts/mcmc/state.h"
#include "math/matrix.h"

namespace BTS {
    
    namespace MCMC {
        
        class State::Tensor: public MR::Math::Matrix<double> {
                
                //Public static variables, nested classes and typedefs
            public:
                
                class Writer;

                //Public member functions
            public:
                
                Tensor() {
                }
                
                Tensor(size_t size)
                        : MR::Math::Matrix<double>(size, size) {
                }
                
                Tensor(const State& value);

                Tensor(const MR::Math::Matrix<double>& tensor)
                        : MR::Math::Matrix<double>(tensor) {
                }
                
                ~Tensor() {
                }
                
                Tensor& operator=(const MR::Math::Matrix<double>& tensor) {
                    MR::Math::Matrix<double>::operator=(tensor);
                    return *this;
                }
                
                size_t length() const {
                    return MR::Math::Matrix<double>::rows();
                }
                
                Tensor& invalidate();

                Tensor& zero();

                MR::Math::Matrix<double> matrix() const {
                    return *this;
                }
                
                void from_matrix(const MR::Math::Matrix<double>& matrix) {
                    MR::Math::Matrix<double>::operator=(matrix);
                }
                
        };
    
    }

}

#endif /* __bts_mcmc_state_tensor_h__ */
