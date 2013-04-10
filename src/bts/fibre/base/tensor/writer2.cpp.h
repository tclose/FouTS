/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 10/08/2010.

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

#ifndef __bts_tensor_writer2_cpp_h__
#define __bts_tensor_writer2_cpp_h__

#include "bts/fibre/tensor/writer2.h"

namespace BTS {
    
    namespace Fibre {
        
        namespace Tensor {
            
            template<typename T> const std::string Writer2<T>::HEADER_PREAMBLE = "%%% Keys %%%:";
            
            template<typename T> void Writer2<T>::create(const std::string& location,
                                                         const T& template_value) {
                
                header = Unzipper::header(template_value);
                template_length = Unzipper::size(template_value);
                
                out.open(location.c_str());
                
                out << HEADER_PREAMBLE;
                
                for (std::vector<std::string>::iterator key_it = header.begin();
                        key_it != header.end(); ++key_it)
                    out << *key_it << " ";
                
                out << std::endl;
                
            }
            
            template<typename T> void Writer2<T>::append(const MR::Math::Matrix<double>& tensor) {
                
                if (template_length != tensor.rows())
                    throw Exception(
                            "Length of tensor matrix (" + str(tensor.rows())
                            + ") did not match size of template value vector ("
                            + str(template_length) + ").");
                
                out << tensor;
                
                out.flush();
                
            }
        
        }
    
    }

}

#endif /* __bts_tensor_writer2_cpp_h__ */
