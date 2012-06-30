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


//#include "bts/fibre/properties/extended.h"
//
//namespace BTS {
//
//	namespace Fibre {
//
//    const std::string                        Properties::Extended::FILE_PREAMBLE           = "%%% Extended Properties File %%% - keys: ";
//
//
//    void                                     Properties::Extended::add_extend_elem_prop (const std::string& key, const std::string& default_value)                {
//
//      //If key is not already present.
//      if (!key_index.count(key)) {
//
//        key_index[key] = header.size();
//        header.push_back(key);
//        default_values.push_back(default_value);
//
//        values.push_back ( std::vector<std::string>(num_rows(), default_value) );
//
//      } else
//        default_values[key_index[key]] = default_value;
//
//    }
//
//
//    void                                     Properties::Extended::add_extend_elem_props (const std::vector<std::string>& header, const std::vector<std::string>& default_values)  {
//
//      if (header.size() != default_values.size())
//        throw Exception ("Number of header (" + str(header.size()) + ") and default values (" + str(default_values.size()) + ") do not match");
//
//      for (size_t key_i = 0; key_i < header.size(); key_i++)
//        add_extend_elem_prop(header[key_i], default_values[key_i]);
//
//    }
//
//
//    void                                     Properties::Extended::remove_extend_elem_prop (const std::string& key)                {
//
//      if (key_index.count(key)) {
//
//        size_t index = key_index[key];
//
//        header.erase(header.begin() + index);
//        default_values.erase(default_values.begin() + index);
//        values.erase(values.begin() + index);
//        key_index.erase(key);
//
//      }
//
//    }
//
//
//    std::string                                        Properties::Extended::get_extend_elem_prop(std::string key, size_t row_index) const                 {
//
//      check_key(key);
//      check_row_index(row_index);
//
//      size_t key_i = key_index.find(key)->second;
//
//      //Can't use operator[] for key_index since 'this' is const. NB: Key is guaranteed to exist in 'key_index' from the 'check_key()' function.
//      std::string value = values[key_i][row_index];
//
//      //FIXME: Doesn't need to be here just for debugging.
//      std::string value_out = value + "";
//
//      return value_out;
//
//    }
//
//
//    void                                     Properties::Extended::set_extend_elem_prop(std::string key, std::string value, size_t row_index)              {
//
//      check_key(key);
//      check_row_index(row_index);
//
//      values[key_index[key]][row_index] = value;
//
//    }
//
//
//    std::vector<std::string>                 Properties::Extended::parse_line(const std::string& line) {
//
//      std::istringstream line_stream(line);
//
//      std::string value;
//      std::vector<std::string> values_row;
//
//      while(std::getline(line_stream, value, '\t'))
//        values_row.push_back(value);
//
//      return values_row;
//
//    }
//
//
//    std::string                              Properties::Extended::make_line(const std::vector<std::string>& values) {
//
//      std::string line;
//
//      std::ostringstream line_stream;
//
//      for (std::vector<std::string>::const_iterator it = values.begin(); it != values.end(); ++it)
//        line_stream << it->c_str() << '\t';
//
//      line_stream << std::endl;
//
//      return line_stream.str();
//
//    }
//
//
//
//    void                                     Properties::Extended::read_props_preamble(std::ifstream& file_in) {
//
//      check_good(file_in);
//
//      char file_intro[FILE_PREAMBLE.size()+2];
//
//      file_in.get(file_intro, FILE_PREAMBLE.size()+1);
//
//      if (strcmp(file_intro, FILE_PREAMBLE.c_str()))
//        throw Exception ("Required file intro '" + FILE_PREAMBLE + "' was not found at start of extended properties file (found '" + std::string(file_intro) +"' instead).");
//
//    }
//
//
//    std::vector<std::string>  Properties::Extended::read_props_header(std::ifstream& file_in) {
//
//      read_props_preamble(file_in);
//
//      check_good(file_in);
//
//      std::vector<std::string> header;
//
//      std::string key_line;
//
//      if (!std::getline(file_in, key_line))
//        throw Exception ("No header found in extended properties file.");
//
//      std::istringstream line_stream(key_line);
//
//      std::string key;
//
//      header.clear();
//
//      while(std::getline(line_stream, key, '\t'))
//        header.push_back(key);
//
//      return header;
//
//    }
//
//
//    void                                     Properties::Extended::write_header(std::ofstream& out, const Header& header) {
//
//      out << FILE_PREAMBLE;
//
//      for (std::vector<std::string>::const_iterator key_it = header.begin(); key_it != header.end(); key_it++)
//        out << *key_it << '\t';
//
//      out << std::endl;
//
//    }
//
//
//    void                                     Properties::Extended::load(std::string location) {
//
//      clear();
//
//      std::ifstream in;
//      in.open(location.c_str());
//      this->file_location = location;
//
//      header = read_props_header(in);
//
//      for (size_t key_i = 0; key_i < header.size(); key_i++)
//        key_index[header[key_i]] = key_i;
//
//      std::string line;
//
//      while (std::getline(in, line))
//        add_row(parse_line(line));
//
//
//      file_location.clear();
//
//
//    }
//
//
//    void                                     Properties::Extended::save(std::string location) {
//
//      std::ofstream out;
//
//      out.open(location.c_str());
//
//      write_header(out, header);
//
//      for (size_t row_i = 0; row_i < num_rows(); row_i++) {
//
//        for (size_t key_i = 0; key_i < num_extend_elem_props(); key_i++)
//
//          out << values[key_i][row_i] << '\t';
//
//        out << std::endl;
//
//      }
//
//    }
//
//
//    Properties                               Properties::Extended::get_row(size_t row_index) const {
//
//      Properties properties;
//
//      for (size_t key_i = 0; key_i < header.size(); key_i++) {
//
//        properties[header[key_i]] = values[key_i][row_index];
//
//      }
//
//      return properties;
//
//    }
//
//
//    void                                     Properties::Extended::add_row(const std::vector<std::string>& values_row) {
//
//      check_row_size(values_row);
//
//      for (size_t key_i = 0; key_i < num_extend_elem_props(); key_i++)
//        this->values[key_i].push_back(values_row[key_i]);
//
//      nrows++;
//
//    }
//
//
//    void                                     Properties::Extended::add_row(const Properties& properties) {
//
//      std::vector<std::string> row(header.size());
//
//      for (size_t key_i = 0; key_i < header.size(); key_i++) {
//
//        Properties::const_iterator key_it = properties.find(header[key_i]);
//
//        if (key_it == properties.end())
//          throw Exception ("Property '" + header[key_i] + "' not found in added properties row.");
//
//        //Can't use operator[] for properties since properties are const.
//        row[key_i] = key_it->second;
//
//      }
//
//
//      for (size_t key_i = 0; key_i < num_extend_elem_props(); key_i++)
//        this->values[key_i].push_back(row[key_i]);
//
//      nrows++;
//
//    }
//
//
//    void                                     Properties::Extended::extend_rows(size_t new_num_rows) {
//
//      for (size_t key_i = 0; key_i < num_extend_elem_props(); key_i++) {
//        for (size_t row_count = num_rows(); row_count <= new_num_rows; row_count++)
//          values[key_i].push_back(default_values[key_i]);
//      }
//
//    }
//
//
//    void                                     Properties::Extended::set_row(const Properties& properties, size_t row_index) {
//
//      std::vector<std::string> row(header.size());
//
//      for (size_t key_i = 0; key_i < header.size(); key_i++) {
//
//        Properties::const_iterator key_it = properties.find(header[key_i]);
//
//        if (key_it == properties.end())
//          throw Exception ("Property '" + header[key_i] + "' not found in added properties row.");
//
//        //Can't use operator[] for properties since it is constant.
//        row[key_i] = key_it->second;
//
//      }
//
//      for (size_t key_i = 0; key_i < num_extend_elem_props(); key_i++)
//        this->values[key_i][row_index] = row[key_i];
//
//    }
//
//	}
//
//}

