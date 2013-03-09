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

#ifndef __bts_fibre_setreader_cpp_h__
#define __bts_fibre_setreader_cpp_h__

namespace BTS {
    
    namespace Fibre {
        
        namespace Base {
            
            template<typename T> void SetReader<T>::open(const std::string& location) {
                
                Reader<T>::open(location);
                
                std::string elem_properties_location = location + "xx";
                
                if (File::is_file(elem_properties_location)) {
                    
                    elem_in.open(elem_properties_location.c_str());
                    
                    elem_prop_hdr = Reader<T>::read_props_header(elem_in);
                    ext_elem_prp_keys = elem_prop_hdr;
                    elem_prp_keys = Object::extract_props<typename T::Element>(ext_elem_prp_keys);
                    
                    elem_begin = elem_in.tellg();
                    
                } else {
                    
                    elem_prop_hdr.clear();
                    ext_elem_prp_keys.clear();
                    elem_prp_keys.clear();
                    
                    elem_begin = 0;
                    
                }
                
            }
            
            template<typename T> void SetReader<T>::rewind() {
                
                Reader<T>::rewind();
                
                if (elem_begin)
                    elem_in.seekg(elem_begin);
                
            }
            
            template<typename T> bool SetReader<T>::next(T& set) {
                
                std::map<std::string, std::string> properties;
                
                set.reset(this->prp_keys, elem_prp_keys);
                
                if (Reader<T>::next_basic(set)) {
                    
                    //Set 'set-wide' properties
                    if (this->prop_hdr.size()) {
                        std::string line;
                        
                        if (!std::getline(this->prop_in, line))
                            throw Exception(
                                    "No corresponding line found in extended properties file.");
                        
                        properties = parse_props_line(this->prop_hdr, line);
                        
                        set.extract_and_set_props(properties);
                        set.set_extend_props(properties);
                        
                    }
                    
                    //Set element properties.
                    set.clear_extend_elem_props();
                    
                    if (elem_prop_hdr.size()) {
                        
                        for (std::vector<std::string>::iterator prop_it = ext_elem_prp_keys.begin();
                                prop_it != ext_elem_prp_keys.end(); ++prop_it)
                            set.add_extend_elem_prop(*prop_it, "");
                        
                        std::string line;
                        
                        if (!std::getline(elem_in, line))
                            throw Exception("No corresponding row found in set properties file.");
                        
                        size_t elem_index = 0;
                        
                        while (!is_row_seperator(line)) {
                            
                            if (elem_index >= set.size())
                                throw Exception(
                                        "Loaded properties contain more rows than elements in the set");
                            
                            std::map<std::string, std::string> elem_prop_row =
                                    Reader<T>::parse_props_line(elem_prop_hdr, line);
                            
                            //FIXME: Needs to be implemented
                            std::map<std::string, std::string> ext_elem_prop_row =
                                    set[elem_index].extract_and_set_props(elem_prop_row);
                            
                            set.set_extend_elem_prop_row(ext_elem_prop_row, elem_index);
                            
                            ++elem_index;
                            
                            if (!std::getline(elem_in, line))
                                throw Exception("Set properties row is incomplete.");
                        }
                        
                    }
                    
                } else
                    return false;
                
                return true;
                
            }
            
            template<typename T> bool SetReader<T>::is_row_seperator(const std::string& line) {
                
                bool is_seperator = false;
                
                if (line.size() > 12)
                    is_seperator = !strcmp(line.substr(0, 8).c_str(), "--- END ")
                            && !strcmp(line.substr(line.size() - 4, 4).c_str(), " ---");
                
                return is_seperator;
                
            }
            
            template<typename T> Coord SetTextReader<T>::read_coord() {
                
                Coord coord;
                
                if (!this->in.good())
                    throw Exception("Corrupted fibre object file '" + this->file_location + "'.");
                
                std::string next_value;
                
                for (int dim_i = 0; dim_i < 3; dim_i++) {
                    
                    this->in >> next_value;
                    
                    if (next_value == "nan")
                        coord[dim_i] = NAN;
                    else if (next_value == "inf")
                        coord[dim_i] = INFINITY;
                    else if (next_value == "-inf")
                        coord[dim_i] = -INFINITY;
                    else if (strspn(next_value.c_str(), "-e0.123456789") == next_value.size())
                        coord[dim_i] = to<double>(next_value);
                    else
                        dim_i--;    //If next_value is not a number empty skip. (this is mainly for the start of the file which has some funny bits).
                        
                }
                
                return coord;
            }
        
        }
    
    }

}

#endif /* __bts_fibre_set_reader_cpp_h__ */
