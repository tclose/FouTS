/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Jul 26, 2010.

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

#ifndef __bts_fibre_setwriter_cpp_h__
#define __bts_fibre_setwriter_cpp_h__

#include "bts/fibre/base/set.h"

namespace BTS {
    
    namespace Fibre {
        
        namespace Base {
            
            template<typename T> void SetWriter<T>::create(
                    const std::string& location, const std::vector<const char*>& props,
                    const std::vector<const char*>& elem_props,
                    const std::vector<std::string> extend_prop_keys,
                    const std::vector<std::string> extend_elem_prop_keys,
                    const std::map<std::string, std::string>& file_props) {
                
                std::vector<std::string> elem_header;
                
                //Add the elem props into the elem header first (doesn't really need to be first but I like it that way when I read the files manually)
                for (size_t prop_i = 0; prop_i < elem_props.size(); ++prop_i)
                    elem_header.push_back(elem_props[prop_i]);
                
                //Insert the extended properties into the set and elem headers.
                elem_header.insert(elem_header.end(), extend_elem_prop_keys.begin(),
                        extend_elem_prop_keys.end());
                
                Writer<T>::create(location, props, extend_prop_keys, file_props);
                
                elem_prop_hdr = elem_header;
                
                std::string set_properties_location = location + "xx";
                
                if (elem_prop_hdr.size()) {
                    
                    ext_elem_out.open((set_properties_location).c_str());
                    
                    Writer<T>::write_props_header(ext_elem_out, elem_prop_hdr);
                    
                } else if (File::exists(set_properties_location))
                    File::remove(set_properties_location);
                
            }
            
            template<typename T> void SetWriter<T>::append(const T& set) {
                
                Base::Writer<T>::append(set, set.get_extend_props());
                
                for (size_t elem_i = 0; elem_i < set.size(); elem_i++) {
                    
                    std::map<std::string, std::string> properties_row =
                            set.get_extend_elem_prop_row(elem_i);    //
                                    
                    set[elem_i].insert_props(properties_row);
                    
                    write_properties_row(this->elem_prop_hdr, properties_row, this->ext_elem_out);
                    
                }
                
                this->ext_elem_out << row_seperator(this->count) << std::endl;
                
            }
            
            template<typename T> void SetWriter<T>::append(
                    const T& set, std::map<std::string, std::string>& properties) {
                
                T new_set(set);
                
                std::map<std::string, std::string> extend_props = new_set.get_extend_props();
                
                for (std::map<std::string, std::string>::iterator prop_it = extend_props.begin();
                        prop_it != extend_props.end(); ++prop_it)
                    new_set.set_extend_prop(prop_it->first, prop_it->second);
                
                append(new_set);
                
            }
            
            template<typename T> void SetTextWriter<T>::write(const Coord& coord) {
                
                this->out << "\n" << coord[X] << " " << coord[Y] << " " << coord[Z];
                this->out.flush();
                
            }
        
        }
    
    }

}

#endif /* __bts_fibre_set_writer_cpp_h__ */
