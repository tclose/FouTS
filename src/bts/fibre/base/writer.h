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

#ifndef __bts_fibre_base_writer_forward_h__
#define __bts_fibre_base_writer_forward_h__

namespace BTS {
    
    namespace Fibre {
        
        namespace Base {
            
            template<typename T> class Writer;
        
        }
    
    }

}
#include "data_type.h"
#include "get_set.h"

#include "bts/coord.h"

#include "bts/file.h"

#include "bts/fibre/base/reader.h"
#include "bts/fibre/properties.h"

namespace BTS {
    
    namespace Fibre {
        
        namespace Base {
            
            template<typename T> class Writer {
                    
                public:
                    
                    static void write_props_header(std::ofstream& out,
                                                   const std::vector<std::string>& header);

                public:
                    
                    size_t count, total_count;

                    // Protected member fields.
                protected:
                    
                    std::ofstream out;
                    MR::DataType dtype;
                    int64_t count_offset;

                    std::vector<std::string> prop_hdr;
                    std::ofstream ext_out;

                    size_t sizeof_coord;

                    // Public methods and constructors
                public:
                    
                    Writer() {
                        init();
                    }
                    
                    template<typename U> Writer(const std::string& location,
                                                const U& template_or_reader,
                                                std::map<std::string, std::string> file_props =
                                                        Properties()) {
                        init();
                        create(location, template_or_reader, file_props);
                    }
                    
                    template<typename U> Writer(const std::string& location,
                                                const U& template_or_reader,
                                                const std::vector<std::string> extend_prop_keys,
                                                std::map<std::string, std::string> file_props =
                                                        Properties()) {
                        init();
                        create(location, template_or_reader, extend_prop_keys, file_props);
                    }
                    
                    template<typename U> Writer(const std::string& location,
                                                const std::vector<const char*>& prop_keys,
                                                const std::vector<std::string> extend_prop_keys,
                                                std::map<std::string, std::string> file_props =
                                                        Properties()) {
                        init();
                        create(location, prop_keys, extend_prop_keys, file_props);
                    }
                    
                    virtual ~Writer() {
                        close();
                    }
                    
                    template<typename U> void create(const std::string& location,
                                                     const U& template_or_reader,
                                                     std::map<std::string, std::string> file_props =
                                                             Properties()) {
                        create(location, Object::select_props<T>(template_or_reader.prop_keys()),
                                std::vector<std::string>(), file_props);
                    }
                    
                    template<typename U> void create(
                            const std::string& location, const U& template_or_reader,
                            const std::vector<std::string> extend_prop_keys,
                            std::map<std::string, std::string> file_props = Properties()) {
                        create(location, Object::select_props<T>(template_or_reader.prop_keys()),
                                std::vector<std::string>(), file_props);
                    }
                    
                    void append(const T& fibre_object,
                                std::map<std::string, std::string> properties_row = Properties());

                    void close() {
                        out.seekp(count_offset);
                        out << count << "\ntotal_count: " << total_count << "\nEND\n";
                        out.close();
                        ext_out.close();
                    }
                    
                    std::ofstream& get_out_stream() {
                        return this->out;
                    }
                    
                    int64_t get_count_offset() {
                        return this->count_offset;
                    }
                    
                    size_t get_count() {
                        return this->count;
                    }
                    
                    size_t get_total_count() {
                        return this->total_count;
                    }
                    
                    void override_count(size_t count) {
                        this->count = count;
                    }
                    
                    void override_total_count(size_t count) {
                        this->total_count = count;
                    }
                    
                    // Protected methods and constructors
                protected:
                    
#ifdef FAST_READWRITE //Will inline writing file to avoid overhead call.  However, this will remove the ability to write to text files (usually only used for debugging though).
                    void write (const Coord& coord);
#else
                    virtual void write(const Coord& coord);
#endif
                    
#ifdef FAST_READWRITE //Will inline writing file to avoid overhead call.  However, this will remove the ability to write to text files (usually only used for debugging though).
                    Coord write_first (const Coord& coord);
#else
                    virtual Coord write_first(const Coord& coord);
#endif
                    
                    template<typename U> void write(const U& fibre_object);

                    template<typename U> Coord write_first(const U& fibre_object);

                    void increment_count() {
                        ++this->count;
                        ++this->total_count;
                        update_count();
                    }
                    
                    void append_file_seperator();

                    void update_count();

                    void check_properties_row(
                            const std::map<std::string, std::string>& properties_row);

                    void append_extend_elem_props_row(
                            std::map<std::string, std::string>& properties_row);

                    void append_extend_elem_props_row(
                            const T& fibre_object,
                            std::map<std::string, std::string>& properties_row);

                    static void write_properties_row(
                            const std::vector<std::string>& key_list,
                            std::map<std::string, std::string>& properties_row,
                            std::ofstream& fout);

#ifdef FAST_WRITE //Will inline writing file to avoid overhead call.  However, this will remove the ability to write to text files.
                    void replace_end_of_file(const Coord& t, int64_t start_pos);
#else
                    virtual void replace_end_of_file(const Coord& t, int64_t start_pos);
#endif
                    
                    void init();

                    void create(const std::string& location,
                                const std::vector<const char*>& prop_keys,
                                const std::vector<std::string> extend_prop_keys,
                                std::map<std::string, std::string> file_props);

                    friend class Set<T> ;
                    
            };
            
            template<typename T> class TextWriter: public Fibre::Base::Writer<T> {
                    
                    //Public static variables, nested classes and typedefs
                public:
                    
                    TextWriter() {
                    }
                    
                    template<typename U> TextWriter(const std::string& location,
                                                    const U& template_or_reader,
                                                    std::map<std::string, std::string> file_props =
                                                            Properties())
                            : Writer<T>(location, template_or_reader, file_props) {
                    }
                    
                    template<typename U> TextWriter(const std::string& location,
                                                    const U& template_or_reader,
                                                    const std::vector<std::string> extend_prop_keys,
                                                    std::map<std::string, std::string> file_props =
                                                            Properties())
                            : Writer<T>(location, template_or_reader, extend_prop_keys, file_props) {
                    }
                    
                    void write(const Coord& coord);

                    void replace_end_of_file(const Coord& c, int64_t start_pos) {
                    }
                    
                    Coord write_first(const Coord& c) {
                        this->out.seekp((int) this->out.tellp() - 12);
                        write(c);
                        return c;
                    }
                    
                    using Fibre::Base::Writer<T>::write;
                    
            };
        
        }
    
    }

}

#include "bts/fibre/base/writer.cpp.h"

#endif /* __bts_fibre_base_writer_h__ */

