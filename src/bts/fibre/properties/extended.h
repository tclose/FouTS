/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

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

//
//#ifndef __bts_fibre_properties_extended_h__
//#define __bts_fibre_properties_extended_h__
//
//#include <vector>
//#include <string>
//#include <map>
//#include <fstream>
//
//#include "bts/common.h"
//
//#include "bts/fibre/properties.h"
//#include "bts/fibre/properties/header.h"
//
//namespace BTS {
//
//  namespace Fibre {
//
//    class Properties::Extended {
//
//      protected:
//
//        Header                                  header;
//        std::vector<std::string>                default_values;
//        std::vector< std::vector<std::string> > values;
//
//        //Used to locate the index of the property.
//        std::map<std::string, size_t> key_index;
//
//        size_t nrows;
//
//        //Used to temporarily store the file location for the printing of error messages.
//        std::string file_location;
//
//      public:
//
//        const static std::string  FILE_PREAMBLE;
//
//        Extended()
//          : nrows(0) {}
//
//
//        Extended(const Extended& ep)
//          : header(ep.header), default_values(ep.default_values), values(ep.values), key_index(ep.key_index), nrows(ep.nrows)   {}
//
//
//        Extended&                 operator= (const Extended& ep)
//          { header = ep.header; default_values = ep.default_values; values = ep.values; key_index = ep.key_index; nrows = ep.nrows; return *this; }
//
//
//        Extended&                 append(const Extended& ep) {
//
//          if (!num_extend_elem_props())
//            this->operator=(ep);
//
//          else {
//
//            size_t num_props = num_extend_elem_props();
//
//            header.insert(header.end(), ep.header.begin(), ep.header.end());
//            default_values.insert(default_values.end(), ep.default_values.begin(), ep.default_values.end());
//            values.insert(values.end(), ep.values.begin(), ep.values.end());
//
//            for (size_t key_i = 0; key_i < ep.num_extend_elem_props(); key_i++)
//              key_index[ep.header[key_i]] = num_props + key_i;
//
//          }
//
//          return *this;
//
//        }
//
//        size_t                              num_rows() const
//          { return nrows; }
//
//
//        size_t                              num_extend_elem_props() const
//          { return header.size(); }
//
//
//        void                                clear()
//          { header.clear(); default_values.clear(); values.clear(); key_index.clear(); nrows = 0; }
//
//
//        Header            	                get_header() const
//          { return header; }
//
//
//        std::vector<std::string>            get_default_values () const
//          { return default_values; }
//
//
//        void                                add_extend_elem_prop (const std::string& key, const std::string& default_value);
//
//
//        void                                add_extend_elem_props (const Header& header, const std::vector<std::string>& default_values);
//
//
//        void                                remove_extend_elem_prop(const std::string& key);
//
//
//        std::string                         get_extend_elem_prop(std::string key, size_t row_index) const;
//
//
//        template <typename U> U             get_extend_elem_prop(std::string key, size_t row_index) const
//          { return to<U>(get_extend_elem_prop(key,row_index)); }
//
//
//        Properties                          get_row(size_t row_index) const;
//
//
//        void                                set_extend_elem_prop(std::string key, std::string value, size_t row_index);
//
//
//        void                                set_row(const Properties& properties, size_t row_index);
//
//
//        bool                                has_extend_elem_prop(const std::string& key) const
//          { return key_index.count(key); }
//
//
//        void                                add_row()
//          { for (size_t key_i = 0; key_i < num_extend_elem_props(); key_i++) values[key_i].push_back(default_values[key_i]); nrows++; }
//
//        void                                erase_row(size_t index)
//          { for (size_t key_i = 0; key_i < num_extend_elem_props(); key_i++) values[key_i].erase(values[key_i].begin() + index); }
//
//
//        void                                add_row(const Properties& properties);
//
//
//        void                                extend_rows(size_t num_rows);
//
//
//        void                                load(std::string location);
//
//
//        void                                save(std::string location);
//
//
//        static std::vector<std::string>     parse_line(const std::string& line);
//
//
//        static std::string                  make_line(const std::vector<std::string>& values);
//
//
//        static Header                       read_props_header(std::ifstream& line);
//
//
//        static void                         write_header(std::ofstream& out, const Header& header);
//
//
//
//      protected:
//
//        void                                check_key (const std::string& key) const
//          { if (!key_index.count(key)) throw Exception ("Property '" + key + "' was not found."); }
//
//
//        void                                check_row_index (size_t row_index) const
//          { if (row_index >= num_rows()) throw Exception ("Properties row index '" + str(row_index) + "' out of range (" + str(num_rows()) + ")."); }
//
//
//        void                                check_row_size(const std::vector<std::string>& row)
//          { if (row.size() >= num_extend_elem_props()) throw Exception ("Size of new row (" + str(row.size()) + ") does not match number of properties (" + str(num_extend_elem_props()) +")."); }
//
//
//        static void                         check_good(std::ifstream& in)
//          { if (!in.good()) throw Exception ("Could not read from properties file."); }
//
//
//        static void                         read_props_preamble(std::ifstream& in);
//
//
//        void                                add_row(const std::vector<std::string>& values_row);
//
//
//
//    };
//
//
//
//  }
//}
//
//#include "bts/fibre/properties/header.h"
//
//#endif /* __bts_fibre_properties_extended_h__ */
