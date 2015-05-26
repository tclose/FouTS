/*
 Copyright 2010 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close on 16/06/2010.

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

#ifndef __bts_utilities_reader_h__
#define __bts_utilities_reader_h__

#include "bts/utilities/common.h"

namespace FTS {
    
    namespace Utilities {
        
        template<typename T> class Reader {
                
            public:
                
                std::ifstream fin;

                Reader() {
                }
                ;

                Reader(const std::string& location, Properties& properties, Header& key_list) {
                    open(location);
                }
                
                Reader(const std::string& location) {
                    open(location);
                }
                
                ~Reader() {
                    close();
                }
                
                void open(const std::string& location);

                void open(const std::string& location, const Properties& dummy_properties,
                          const Header& dummy_header) {
                    open(location);
                }
                
                void close() {
                    fin.close();
                }
                
                bool next(T& state) {
                    Properties dummy;
                    return next(state, dummy);
                }
                
                bool next(T& state, Properties& properties_row);
                
        };
    
    }

}

#endif /* __bts_utilities_reader_h__ */
