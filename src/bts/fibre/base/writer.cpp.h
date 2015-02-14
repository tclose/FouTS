/*
 Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

 Written by Thomas G Close, 04/06/2010.

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

#ifndef __bts_fibre_base_writer_cpp_h__
#define __bts_fibre_base_writer_cpp_h__

#include "mrtrix.h"
#include "bts/version.h"

namespace BTS {
    
    namespace Fibre {
        
        namespace Base {
            
            template<typename T> void Writer<T>::create(
                    const std::string& location, const std::vector<const char*>& prop_keys,
                    const std::vector<std::string> extend_prop_keys,
                    std::map<std::string, std::string> file_props) {

                //Combine prop keys and extended prop keys into one header vector.
                std::vector<std::string> prop_header;
                
                for (size_t prop_i = 0; prop_i < prop_keys.size(); ++prop_i)
                    prop_header.push_back(prop_keys[prop_i]);
                
                prop_header.insert(prop_header.end(), extend_prop_keys.begin(),
                        extend_prop_keys.end());

                file_props["software version"] = version_number_string();
                file_props["datetime"] = current_datetime();
                
                out.open(location.c_str(), std::ios::out | std::ios::binary);
                if (!out)
                    throw Exception(
                            "error creating tracks file \"" + location + "\": " + strerror(errno));
                
                out << "mrtrix tracks\nEND\n";
                for (std::map<std::string, std::string>::const_iterator i = file_props.begin();
                        i != file_props.end(); ++i)
                    out << i->first << ": " << i->second << "\n";
                
                out << "datatype: " << dtype.specifier() << "\n";
                int64_t data_offset = int64_t(out.tellp()) + 65;
                out << "file: . " << data_offset << "\n";
                out << "count: ";
                count_offset = out.tellp();
                out << "\nEND\n";
                out.seekp(0);
                out << "mrtrix tracks    ";
                out.seekp(data_offset);
                write(Reader<T>::END_OF_DATA);
                
                prop_hdr = prop_header;
                
                std::string props_row_location = location + "x";
                
                if (prop_hdr.size()) {
                    
                    ext_out.open((props_row_location).c_str());
                    
                    write_props_header(ext_out, prop_hdr);
                    
                } else if (File::exists(props_row_location))
                    
                    File::remove(props_row_location);
                
                this->count = this->total_count = 0;
                
            }
            
            template<typename T> void Writer<T>::append(
                    const T& fibre_object, std::map<std::string, std::string> properties_row) {
                
                // Record position of current END_OF_DATA marker, so it can be replaced as the last step, thus ensuring
                // that the subsequent fibre object is completely written before the END_OF_DATA marker is updated.
                int64_t start_pos = this->out.tellp();
                
                // Holds the Triple that will replace the END_OF_DATA marker in the final step
                Coord first_coord;
                
                if (fibre_object.size()) {
                    
                    // Writes the first fibre object minus the first coord, which is returned to replace the END_OF_DATA marker in the final step.
                    first_coord = write_first(fibre_object[0]);
                    
                    // Write the remaining fibre objects to the file
                    for (size_t fibre_i = 1; fibre_i < fibre_object.size(); fibre_i++)
                        write(fibre_object[fibre_i]);
                    
                    // Write the seperator sequence which denotes the end of the fibre object.
                    write(T::FILE_SEPARATOR);
                    
                } else
                    // If the fibre_object is of size 0, the first coord will be simply the file seperator.
                    first_coord = T::FILE_SEPARATOR;
                
                // Write the new END_OF_DATA marker at the new end of the file.
                write(Reader<T>::END_OF_DATA);
                
                this->out.flush();
                
                // Now that the rest of the object is completely written replace the END_OF_DATA marker with the first element to complete write.
                replace_end_of_file(first_coord, start_pos);
                
                // Add extended row properties to accompanying text file.
                append_extend_elem_props_row(fibre_object, properties_row);
                
                // Update the fibre object count.
                increment_count();
                
            }
            
            // 'write' and 'write_first' work by recursion.  To use this fibre the fibre objects are only required to implement the 'operator[]
            // ()' and have a const static Coord field called FILE_SEPARATOR.  The recursion stops when the operator[]() returns a Coord.
            template<typename T> template<typename U> inline void Writer<T>::write(
                    const U& fibre_object) {
                
                for (size_t fibre_i = 0; fibre_i < fibre_object.size(); fibre_i++)
                    write(fibre_object[fibre_i]);
                
                write(U::FILE_SEPARATOR);
                
            }
            
            template<typename T> template<typename U> inline Coord Writer<T>::write_first(
                    const U& fibre_object) {
                
                Coord first_coord;
                
                if (fibre_object.size()) {
                    
                    first_coord = write_first(fibre_object[0]);
                    
                    for (size_t fibre_i = 1; fibre_i < fibre_object.size(); fibre_i++)
                        write(fibre_object[fibre_i]);
                    
                    write(U::FILE_SEPARATOR);
                    
                } else
                    
                    first_coord = U::FILE_SEPARATOR;
                
                return first_coord;
                
            }
            
            template<typename T> void Writer<T>::write(const Coord& t) {
                
                using namespace MR::ByteOrder;
                
                if (this->dtype == MR::DataType::Float64LE) {
                    
                    double x[3];
                    
                    x[0] = MR::ByteOrder::LE(t[0]);
                    x[1] = MR::ByteOrder::LE(t[1]);
                    x[2] = MR::ByteOrder::LE(t[2]);
                    
                    this->out.write((const char*) &x[0], 3 * sizeof(double));
                    
                } else if (this->dtype == MR::DataType::Float64BE) {
                    
                    double x[3];
                    
                    x[0] = MR::ByteOrder::BE(t[0]);
                    x[1] = MR::ByteOrder::BE(t[1]);
                    x[2] = MR::ByteOrder::BE(t[2]);
                    
                    this->out.write((const char*) &x[0], 3 * sizeof(double));
                    
                } else if (this->dtype == MR::DataType::Float32LE) {
                    
                    float x[3];
                    
                    x[0] = MR::ByteOrder::LE(t[0]);
                    x[1] = MR::ByteOrder::LE(t[1]);
                    x[2] = MR::ByteOrder::LE(t[2]);
                    
                    this->out.write((const char*) &x[0], 3 * sizeof(float));
                    
                } else if (this->dtype == MR::DataType::Float32BE) {
                    
                    float x[3];
                    
                    x[0] = MR::ByteOrder::BE(t[0]);
                    x[1] = MR::ByteOrder::BE(t[1]);
                    x[2] = MR::ByteOrder::BE(t[2]);
                    
                    this->out.write((const char*) &x[0], 3 * sizeof(float));
                    
                } else
                    throw Exception("Unsupported dataype '" + str(this->dtype.specifier()) + "'.");
                
            }
            
            //NB:: that the first coord is saved until last before overwriting the previous END_OF_DATA marker.
            // This specialization of the template function write_first is the  termination of a recursion
            template<typename T> Coord Writer<T>::write_first(const Coord& t) {
                
                return t;
                
            }
            
            template<typename T> void Writer<T>::update_count() {
                
                int64_t current_pos = this->out.tellp();
                
                this->out.seekp(this->count_offset);
                this->out << this->count << "\ntotal_count: " << this->total_count << "\nEND\n";
                
                this->out.seekp(current_pos);
                
            }
            
            template<typename T> void Writer<T>::replace_end_of_file(const Coord& t,
                                                                     int64_t start_pos) {
                
                int64_t current = this->out.tellp();
                this->out.seekp(start_pos - sizeof_coord);
                write(t);
                this->out.seekp(current);
                
            }
            
            template<typename T> void Writer<T>::check_properties_row(
                    const std::map<std::string, std::string>& properties_row) {
                
                if (properties_row.size() != this->prop_hdr.size()) {
                    
                    //Generate list of header for error message.
                    std::vector<std::string> header;
                    for (std::map<std::string, std::string>::const_iterator key_it =
                            properties_row.begin(); key_it != properties_row.end(); key_it++)
                        header.push_back(key_it->first);
                    
                    //Throw Exception
                    throw Exception(
                            "Supplied row properties (" + str(header)
                            + ") does not match the number of supplied header (" + str(prop_hdr)
                            + "\").");
                }
                
            }
            
            template<typename T> void Writer<T>::append_extend_elem_props_row(
                    const T& fibre_object, std::map<std::string, std::string>& properties_row) {
                fibre_object.insert_props(properties_row);
                
                write_properties_row(this->prop_hdr, properties_row, this->ext_out);
                
            }
            
            template<typename T> void Writer<T>::append_extend_elem_props_row(
                    std::map<std::string, std::string>& properties_row) {
                
                write_properties_row(this->prop_hdr, properties_row, this->ext_out);
                
            }
            
            template<typename T> void Writer<T>::write_properties_row(
                    const std::vector<std::string>& properties_header,
                    std::map<std::string, std::string>& properties_row, std::ofstream& fout) {
                
                for (std::vector<std::string>::const_iterator hdr_it = properties_header.begin();
                        hdr_it != properties_header.end(); hdr_it++) {
                    
                    if (!properties_row.count(*hdr_it))
                        throw Exception(
                                "Did not find value corresponding to property '" + *hdr_it + "'.");
                    
                    fout << properties_row[*hdr_it] << '\t';
                    
                }
                
                fout << std::endl;
                
            }
            
            template<typename T> void Writer<T>::init() {
                
                if (this->dtype == MR::DataType::Float64LE || this->dtype
                        == MR::DataType::Float64BE)
                    sizeof_coord = 3 * sizeof(double);
                else if (this->dtype == MR::DataType::Float32LE
                        || this->dtype == MR::DataType::Float32BE)
                    sizeof_coord = 3 * sizeof(float);
                else
                    throw Exception("Unsupported dataype '" + str(this->dtype.specifier()) + "'.");
                
            }
            
            template<typename T> void Writer<T>::write_props_header(
                    std::ofstream& out, const std::vector<std::string>& header) {
                
                out << Reader<T>::PROP_FILE_PREAMBLE;
                
                for (std::vector<std::string>::const_iterator key_it = header.begin();
                        key_it != header.end(); key_it++)
                    out << *key_it << '\t';
                
                out << std::endl;
                
            }
            
            template<typename T> void TextWriter<T>::write(const Coord& coord) {
                
                this->out << "\n" << coord[X] << " " << coord[Y] << " " << coord[Z];
                this->out.flush();
                
            }
        
        }
    
    }

}

#endif /* __bts_fibre_writer_icpp__ */
