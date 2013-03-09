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

#include "bts/mcmc/state/tensor/writer.h"

namespace BTS {
    
    namespace MCMC {
        
        const std::string State::Tensor::Writer::HEADER_PREAMBLE = "%%% Keys %%%:";
        
        void State::Tensor::Writer::create(const std::string& location, size_t length) {
            
            this->length = length;
            
            out.open(location.c_str());
            
            out << HEADER_PREAMBLE;
            
            for (size_t i = 0; i < length; i++)
                out << "elem" << i << " ";
            
            out << std::endl;
            
        }
        
        void State::Tensor::Writer::append(const MCMC::State::Tensor& tensor) {
            
            if (length != tensor.length())
                throw Exception(
                        "Length of tensor matrix (" + str(tensor.length())
                        + ") did not match size of template value vector (" + str(length) + ").");
            
            out << tensor;
            
        }
    
    }

}
