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

#ifndef __bts_fibre_base_tensorwriter_cpp_h__
#define __bts_fibre_base_tensorwriter_cpp_h__

#include "bts/fibre/base/tensor_writer.h"

namespace BTS {
    
    namespace Fibre {
        
        namespace Base {
            
            template<typename T> const std::string TensorWriter<T>::HEADER_PREAMBLE =
                    "%%% Keys %%%:";
            
            template<typename T> void TensorWriter<T>::create(
                    const std::string& location, const typename T::RowType& template_value) {
                
                hdr = header(template_value);
                template_length = template_value.vsize();
                
                out.open(location.c_str());
                
                out << HEADER_PREAMBLE;
                
                for (std::vector<std::string>::iterator key_it = hdr.begin(); key_it != hdr.end();
                        ++key_it)
                    out << *key_it << " ";
                
                out << std::endl;
                
            }
            
            template<typename T> void TensorWriter<T>::append(const T& tensor) {
                
                if (template_length != tensor.rows())
                    throw Exception(
                            "Number of rows in tensor matrix (" + str(tensor.rows())
                            + ") did not match size of template (" + str(template_length) + ").");
                
                if (template_length != tensor.columns())
                    throw Exception(
                            "Number of columns in tensor matrix (" + str(tensor.columns())
                            + ") did not match size of template (" + str(template_length) + ").");
                
                out << (MR::Math::Matrix<double>&) tensor;
                
            }
            
            template<typename T> std::vector<std::string> TensorWriter<T>::header(
                    const typename T::RowType& tmplate) {
                
                std::vector<std::string> headr;
                
                header_element(tmplate, headr, "");
                
                return headr;
                
            }
            
            template<typename T> void TensorWriter<T>::header_element(
                    const Strand& strand, std::vector<std::string>& header,
                    const std::string& prefix) {
                
                for (size_t prop_i = 0; prop_i < strand.num_props(); ++prop_i)
                    header.push_back(prefix + strand.prop_key(prop_i) + ">");
                
                for (uint degree_i = 0; degree_i < strand.degree(); degree_i++)
                    header_element(strand[degree_i], header,
                            prefix + "degree" + str(degree_i) + ">");
                
            }
            
            template<typename T> void TensorWriter<T>::header_element(
                    const Tractlet& tractlet, std::vector<std::string>& header,
                    const std::string& prefix) {
                
                for (size_t prop_i = 0; prop_i < tractlet.num_props(); ++prop_i)
                    header.push_back(prefix + tractlet.prop_key(prop_i) + ">");
                
                for (uint ax_i = 0; ax_i < 3; ax_i++)
                    header_element(tractlet[ax_i], header, prefix + "axis" + str(ax_i) + ">");
                
            }
            
            template<typename T> void TensorWriter<T>::header_element(
                    const Tractlet::Set& set, std::vector<std::string>& header,
                    const std::string& prefix) {
                
                for (size_t prop_i = 0; prop_i < set.num_props(); ++prop_i)
                    header.push_back(prefix + set.prop_key(prop_i) + ">");
                
                for (uint elem_i = 0; elem_i < set.size(); elem_i++)
                    header_element(set[elem_i], header, prefix + "tractlet" + str(elem_i) + ">");
                
            }
            
            template<typename T> void TensorWriter<T>::header_element(
                    const Strand::Set& set, std::vector<std::string>& header,
                    const std::string& prefix) {
                
                for (size_t prop_i = 0; prop_i < set.num_props(); ++prop_i)
                    header.push_back(prefix + set.prop_key(prop_i) + ">");
                
                for (uint elem_i = 0; elem_i < set.size(); elem_i++)
                    header_element(set[elem_i], header, prefix + "strand" + str(elem_i) + ">");
                
            }
            
            template<typename T> void TensorWriter<T>::header_element(
                    const Coord& element, std::vector<std::string>& header,
                    const std::string& prefix) {
                
                for (uint dim_i = 0; dim_i < 3; dim_i++)
                    header.push_back(prefix + "dim" + str(dim_i) + ">");
                
            }
        
        }
    
    }

}

#endif /* __bts_fibre_base_tensorwriter_cpp_h__ */
