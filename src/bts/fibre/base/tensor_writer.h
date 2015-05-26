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

#ifndef __bts_fibre_base_tensorwriter_h__
#define __bts_fibre_base_tensorwriter_h__

#include "bts/fibre/strand.h"
#include "bts/fibre/tractlet.h"

namespace FTS {
    
    namespace Fibre {
        
        namespace Base {
            
            template<typename T> class TensorWriter {
                    
                    //Public static variables, nested classes and typedefs
                public:
                    
                    static const std::string HEADER_PREAMBLE;

                    //Protected member variables
                protected:
                    
                    std::vector<std::string> hdr;
                    size_t template_length;
                    std::ofstream out;

                public:
                    
                    static std::vector<std::string> header(const typename T::RowType& tmplate);

                    static void header_element(const Strand& tmplate,
                                               std::vector<std::string>& header,
                                               const std::string& prefix);
                    static void header_element(const Tractlet& tmplate,
                                               std::vector<std::string>& header,
                                               const std::string& prefix);
                    static void header_element(const Strand::Set& tmplate,
                                               std::vector<std::string>& header,
                                               const std::string& prefix);
                    static void header_element(const Tractlet::Set& tmplate,
                                               std::vector<std::string>& header,
                                               const std::string& prefix);
                    static void header_element(const Coord& tmplate,
                                               std::vector<std::string>& header,
                                               const std::string& prefix);

                    //Public member functions
                public:
                    
                    TensorWriter() {
                    }
                    
                    TensorWriter(const std::string& location,
                                 const typename T::RowType& template_value =
                                         typename T::RowType()) {
                        create(location, template_value);
                    }
                    
                    TensorWriter(const std::string& location, size_t size) {
                        create(location, size);
                    }
                    
                    ~TensorWriter() {
                        close();
                    }
                    
                    void close() {
                        out.close();
                    }
                    
                    void create(const std::string& location, size_t size) {
                        create(location, T(size));
                    }
                    
                    void create(const std::string& location,
                                const typename T::RowType& template_value);

                    void append(const T& tensor);
                    
            };
        
        }
    
    }

}

#endif /* __bts_fibre_base_tensorwriter_h__ */
