/*
 Copyright 2010 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close on Jun 3, 2010.

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

#ifndef __bts_fibre_base_reader_cpp_h__
#define __bts_fibre_base_reader_cpp_h__

#include "mrtrix.h"
#include "file/key_value.h"

#include "bts/fibre/base/object.h"

namespace BTS {
    
    namespace Fibre {
        
        namespace Base {
            
            template<typename T> const Coord Reader<T>::END_OF_DATA(INFINITY, INFINITY, INFINITY);
            
            template<typename T> const std::string Reader<T>::PROP_FILE_PREAMBLE =
                    "%%% Extended Properties File %%% - keys: ";
            
            template<typename T> void Reader<T>::open(const std::string& location) {
                
                ext_props.clear();
                
                MR::File::KeyValue kv(location, "mrtrix tracks");
                std::string data_file;
                
                while (kv.next()) {
                    std::string key = MR::lowercase(kv.key());
                    if (key == "file")
                        data_file = kv.value();
                    else if (key == "datatype")
                        dtype.parse(kv.value());
                    else if (key != "roi")
                        ext_props[key] = kv.value();
                }
                
                if (dtype == MR::DataType::Undefined)
                    throw Exception(
                            "no datatype specified in fibre object file '" + location + "'");
                if (dtype != MR::DataType::Float32LE && dtype != MR::DataType::Float32BE
                    && dtype != MR::DataType::Float64LE && dtype != MR::DataType::Float64BE)
                    throw Exception(
                            "only supported datatype for tracks location are Float32LE, Float32BE, Float64LE & Float64BE (in tracks location \"" + location
                            + "\")");
                
                if (data_file.empty())
                    throw Exception(
                            "missing 'file' specification in fibre object file '" + location + "'");
                
                std::istringstream locations_stream(data_file);
                std::string fname;
                locations_stream >> fname;
                int64_t offset = 0;
                if (locations_stream.good()) {
                    try {
                        locations_stream >> offset;
                    } catch (...) {
                        throw Exception(
                                "invalid data offset specified, \"" + fname
                                + "\", in fibre object file '" + location + "'");
                    }
                }
                
                if (fname != ".")
                    fname = MR::Path::join(MR::Path::dirname(location), fname);
                else
                    fname = location;
                
                in.open(fname.c_str(), std::ios::in | std::ios::binary);
                if (!in)
                    throw Exception(
                            "error opening fibre object file \"" + fname + "\": "
                            + strerror(errno));
                in.seekg(offset);
                
                begin = in.tellg();
                
                //-------------------------
                // Load element_properties
                //--------------------------
                
                std::string elem_props_location = location + "x";
                
                if (File::is_file(elem_props_location)) {
                    
                    prop_in.open(elem_props_location.c_str());
                    
                    prop_hdr = read_props_header(prop_in);
                    ext_prp_keys = prop_hdr;
                    prp_keys = Object::extract_props<T>(ext_prp_keys);
                    
                    prop_begin = prop_in.tellg();
                    
                } else {
                    
                    prop_hdr.clear();
                    prp_keys.clear();
                    ext_prp_keys.clear();
                    
                    prop_begin = 0;
                    
                }
                
            }
            
            template<typename T> void Reader<T>::rewind() {
                
                in.seekg(begin);
                
                if (prop_begin)
                    prop_in.seekg(prop_begin);
                
            }
            
            template<typename T> bool Reader<T>::next_basic(T& fibre_object) {
                
                if (!in.good())
                    return false;
                
                Coord first_coord = read_coord();
                
                if ((Triple<double> ) first_coord == (Triple<double> ) END_OF_DATA)
                    return false;
                
                read(fibre_object, first_coord);
                
                return true;
                
            }
            
            template<typename T> bool Reader<T>::next(
                    T& fibre_object, std::map<std::string, std::string>& properties_row) {
                
                fibre_object.clear(prp_keys);
                properties_row.clear();
                
                if (next_basic(fibre_object)) {
                    
                    if (prop_hdr.size()) {
                        std::string line;
                        
                        if (!std::getline(prop_in, line))
                            throw Exception(
                                    "No corresponding line found in extended properties file.");
                        
                        properties_row = parse_props_line(prop_hdr, line);
                        
                        fibre_object.extract_and_set_props(properties_row);
                        
                    }
                    
                } else
                    return false;
                
                return true;
                
            }
            
            template<typename T> bool Reader<T>::next(T& fibre_object) {
                
                std::map<std::string, std::string> properties_row;
                
                return next(fibre_object, properties_row);
                
            }
            
            template<typename T> Coord Reader<T>::read_coord() {
                
                if (!in.good())
                    throw Exception("Corrupted fibre object file '" + file_location + "'.");
                
                Coord t;
                
                if (dtype == MR::DataType::Float64LE) {
                    
                    in.read((char*) &t[0], sizeof(double) * 3);
                    
                    t[X] = MR::ByteOrder::LE(t[X]);
                    t[Y] = MR::ByteOrder::LE(t[Y]);
                    t[Z] = MR::ByteOrder::LE(t[Z]);
                    
                } else if (dtype == MR::DataType::Float64BE) {
                    
                    in.read((char*) &t[0], sizeof(double) * 3);
                    
                    t[X] = MR::ByteOrder::BE(t[X]);
                    t[Y] = MR::ByteOrder::BE(t[Y]);
                    t[Z] = MR::ByteOrder::BE(t[Z]);
                    
                } else if (dtype == MR::DataType::Float32LE) {
                    
                    Triple<float> t_float;
                    
                    in.read((char*) &t_float[0], sizeof(float) * 3);
                    
                    t[X] = (double) MR::ByteOrder::LE(t_float[X]);
                    t[Y] = (double) MR::ByteOrder::LE(t_float[Y]);
                    t[Z] = (double) MR::ByteOrder::LE(t_float[Z]);
                    
                } else if (dtype == MR::DataType::Float32BE) {
                    
                    Triple<float> t_float;
                    
                    in.read((char*) &t_float[0], sizeof(float) * 3);
                    
                    t[X] = (double) MR::ByteOrder::BE(t_float[X]);
                    t[Y] = (double) MR::ByteOrder::BE(t_float[Y]);
                    t[Z] = (double) MR::ByteOrder::BE(t_float[Z]);
                    
                } else
                    throw Exception("Unsupported dataype '" + str(dtype.specifier()) + "'.");
                
                return t;
                
            }
            
            template<typename T> void Reader<T>::read(Coord& t, const Coord& first_t) {
                
                t = first_t;
                
            }
            
            template<typename T> template<typename U> inline void Reader<T>::read(
                    U& fibre_object, const Coord& first_t) {
                
                Coord first_coord = first_t;
                
                while ((Triple<double> ) first_coord != (Triple<double> ) U::FILE_SEPARATOR) {
                    
                    typename U::Element elem;
                    
                    elem.clear();
                    
                    read(elem, first_coord);
                    
                    fibre_object.push_back(elem, false);    // Push back ignoring properties.
                            
                    first_coord = read_coord();
                    
                }
                
            }
            
            template<typename T> std::map<std::string, std::string> Reader<T>::parse_props_line(
                    const std::vector<std::string>& header, const std::string& line) {
                
                std::istringstream line_stream(line);
                
                std::string value;
                std::map<std::string, std::string> prop_row;
                
                size_t count = 0;
                while (std::getline(line_stream, value, '\t')) {
                    if (count > header.size())
                        throw Exception(
                                "More elements found on line  (" + str(count)
                                + ") than exist in header (" + str(header.size())
                                + "). The line was {" + line + "} and the header is " + str(header)
                                + ".");
                    prop_row[header[count]] = value;
                    ++count;
                }
                
                return prop_row;
                
            }
            
            template<typename T> void Reader<T>::read_props_preamble(std::ifstream& file_in) {
                
                check_good(file_in);
                
                char file_intro[PROP_FILE_PREAMBLE.size() + 2];
                
                file_in.get(file_intro, PROP_FILE_PREAMBLE.size() + 1);
                
                if (strcmp(file_intro, PROP_FILE_PREAMBLE.c_str()))
                    throw Exception(
                            "Required file intro '" + PROP_FILE_PREAMBLE
                            + "' was not found at start of extended properties file (found '"
                            + std::string(file_intro) + "' instead).");
                
            }
            
            template<typename T> std::vector<std::string> Reader<T>::read_props_header(
                    std::ifstream& file_in) {
                
                read_props_preamble(file_in);
                
                check_good(file_in);
                
                std::vector<std::string> header;
                
                std::string key_line;
                
                if (!std::getline(file_in, key_line))
                    throw Exception("No header found in extended properties file.");
                
                std::istringstream line_stream(key_line);
                
                std::string key;
                
                header.clear();
                
                while (std::getline(line_stream, key, '\t'))
                    header.push_back(key);
                
                return header;
                
            }
            
            template<typename T> Coord TextReader<T>::read_coord() {
                
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

#endif /* __bts_fibre_reader_icpp__ */
