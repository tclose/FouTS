/*
 Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

 Written by Thomas G Close, 5/05/09.

 This file is part of Bayesian Tractlet Sampling (BTS).

 BTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 BTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with BTS.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef __bts_fibre_setreader_h__
#define __bts_fibre_setreader_h__

#include <vector>
#include <string>
#include <iostream>

#include "bts/fibre/base/set.h"
#include "bts/fibre/base/reader.h"

namespace BTS {
    
    namespace Fibre {
        
        namespace Base {
            
            template<typename T> class SetReader: public Reader<T> {
                    
                    // Protected member fields.
                protected:
                    
                    std::vector<std::string> elem_prop_hdr;
                    std::vector<const char*> elem_prp_keys;
                    std::vector<std::string> ext_elem_prp_keys;

                    std::ifstream elem_in;
                    off64_t elem_begin;

                    // Public methods and constructors.
                public:
                    
                    SetReader() {
                    }
                    
                    SetReader(const std::string& location) {
                        open(location);
                    }
                    
                    ~SetReader() {
                        close();
                    }
                    
                    void open(const std::string& location);

                    void rewind();

                    bool next(T& fibre);

                    bool next(T& fibre, std::map<std::string, std::string>& properties) {
                        if (next(fibre)) {
                            properties = fibre.get_extend_props();
                            return true;
                        } else
                            return false;
                    }
                    
                    void close() {
                        Base::Reader<T>::close();
                        elem_in.close();
                        elem_prop_hdr.clear();
                    }
                    
                    std::vector<std::string> extend_elem_prop_keys() const {
                        return ext_elem_prp_keys;
                    }
                    
                    std::vector<const char*> elem_prop_keys() const {
                        return elem_prp_keys;
                    }
                    
                    std::vector<std::string> elem_header() const {
                        return elem_prop_hdr;
                    }
                    
                protected:
                    
                    static bool is_row_seperator(const std::string& line);
                    
            };
            
            template<typename T> class SetTextReader: public SetReader<T> {
                    
                    //Public static variables, nested classes and typedefs
                public:
                    
                    SetTextReader() {
                    }
                    
                    SetTextReader(const std::string& location)
                            : Base::SetReader<T>(location) {
                    }
                    
                    Coord read_coord();
                    
            };
        
        }
    
    }

}

#include "bts/fibre/base/set_reader.cpp.h"

#endif
