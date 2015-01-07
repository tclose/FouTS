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

#ifndef __bts_fibre_base_reader_forward_h__
#define __bts_fibre_base_reader_forward_h__

namespace BTS {
    
    namespace Fibre {
        
        namespace Base {
            
            template<typename T> class Reader;
        
        }
    
    }

}
#include "data_type.h"
#include "get_set.h"
#include "bts/coord.h"

#include "bts/common.h"

#include "bts/file.h"

namespace BTS {
    
    namespace Fibre {
        
        namespace Base {
            
            template<typename T> class Reader {
                    
                public:
                    
                    static const Coord END_OF_DATA;
                    const static std::string PROP_FILE_PREAMBLE;

                public:
                    
                    static std::map<std::string, std::string> parse_props_line(
                            const std::vector<std::string>& header, const std::string& line);

                    static std::vector<std::string> read_props_header(std::ifstream& line);

                    static void read_props_preamble(std::ifstream& file_in);

                    static void check_good(std::ifstream& in) {
                        if (!in.good())
                            throw Exception("Could not read from properties file.");
                    }
                    
                protected:
                    
                    std::ifstream in;
                    MR::DataType dtype;
                    size_t count;
                    std::string file_location;
                    std::map<std::string, std::string> ext_props;
                    std::vector<std::string> prop_hdr;
                    std::vector<const char*> prp_keys;
                    std::vector<std::string> ext_prp_keys;
                    std::ifstream prop_in;
                    int64_t prop_begin;
                    int64_t begin;

                public:
                    
                    Reader() {
                    }
                    
                    Reader(const std::string& location) {
                        open(location);
                    }
                    
                    virtual ~Reader() {
                        close();
                    }
                    
                    bool next(T& fibre);

                    bool next(T& fibre, std::map<std::string, std::string>& properties_row);

                    void close() {
                        in.close();
                        prop_in.close();
                    }
                    
                    void open(const std::string& location);

                    //TODO: Refactor to extend_elem_prop_keys() or somethinnessg
                    std::vector<std::string> extend_prop_keys() const {
                        return ext_prp_keys;
                    }
                    
                    std::vector<const char*> prop_keys() const {
                        return prp_keys;
                    }
                    
                    std::map<std::string, std::string> get_extend_props() const {
                        return ext_props;
                    }
                    
                    void rewind();

                protected:
                    
                    std::vector<std::string> prop_header() const {
                        return prop_hdr;
                    }
                    
                    bool next_basic(T& fibre);

#ifdef FAST_READWRITE // Means that Text::Writer will no longer work but avoids a virtual function call here (maybe a bit OTT but an option if you feel like it)
                    Coord read_coord ();
#else
                    virtual Coord read_coord();
#endif
                    
                    void read(Coord& coord, const Coord& first_t);

                    template<typename U> void read(U& fibre_object, const Coord& first_t);
                    
            };
            
            template<typename T> class TextReader: public Fibre::Base::Reader<T> {
                    
                    //Public static variables, nested classes and typedefs
                public:
                    
                    TextReader() {
                    }
                    
                    TextReader(const std::string& location)
                            : Fibre::Base::Reader<T>(location) {
                    }
                    
                    Coord read_coord();
                    
            };
        
        }
    
    }

}

#include "bts/fibre/base/reader.cpp.h"

#endif

