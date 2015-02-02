/*
   Copyright 2008 Brain Research Institute, Melbourne, Australia

   Written by J-Donald Tournier, 27/06/08.

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

#ifndef __image_name_parser_h__
#define __image_name_parser_h__

#include "mrtrix.h"
#include "ptr.h"
#include "file/path.h"

namespace MR
{
  namespace Image
  {

    //! a class to interpret numbered filenames
    class NameParser
    {
      public:
        class Item
        {
          public:
            Item () : seq_length (0) { }

            void set_str (const std::string& s) {
              clear ();
              str = s;
            }

            void set_seq (const std::string& s) {
              clear ();
              if (s.size()) seq = parse_ints (s);
              seq_length = 1;
            }

            void clear () {
              str.clear();
              seq.clear();
              seq_length = 0;
            }

            std::string string () const {
              return (str);
            }

            const std::vector<int>& sequence () const {
              return (seq);
            }

            std::vector<int>& sequence () {
              return (seq);
            }

            bool is_string () const {
              return (seq_length == 0);
            }

            bool is_sequence () const {
              return (seq_length != 0);
            }

            size_t size () const {
              return (seq_length ? seq_length : str.size());
            }

            void calc_padding (size_t maxval = 0);

            friend std::ostream& operator<< (std::ostream& stream, const Item& item);

          protected:
            size_t seq_length;
            std::string str;
            std::vector<int> seq;
        };


        void parse (const std::string& imagename, size_t max_num_sequences = std::numeric_limits<size_t>::max());

        size_t num () const {
          return (array.size());
        }

        std::string spec () const {
          return (specification);
        }

        const Item& operator[] (size_t i) const {
          return (array[i]);
        }

        const std::vector<int>& sequence (size_t index) const {
          return (array[seq_index[index]].sequence());
        }

        size_t ndim () const {
          return (seq_index.size());
        }

        size_t index_of_sequence (size_t number = 0) const {
          return (seq_index[number]);
        }

        bool match (const std::string& file_name, std::vector<int>& indices) const;
        void calculate_padding (const std::vector<int>& maxvals);
        std::string name (const std::vector<int>& indices);
        std::string get_next_match (std::vector<int>& indices, bool return_seq_index = false);

        friend std::ostream& operator<< (std::ostream& stream, const NameParser& parser);

      private:
        std::vector<Item> array;
        std::vector<size_t> seq_index;
        std::string folder_name, specification, current_name;
        Ptr<Path::Dir> folder;

        void insert_str (const std::string& str) {
          Item item;
          item.set_str (str);
          array.insert (array.begin(), item);
        }

        void insert_seq (const std::string& str) {
          Item item;
          item.set_seq (str);
          array.insert (array.begin(), item);
          seq_index.push_back (array.size()-1);
        }
    };








    //! a class to hold a parsed image filename
    class ParsedName
    {
      public:
        ParsedName (const std::string& name, const std::vector<int>& index) : indices (index), filename (name) { }

        //! a class to hold a set of parsed image filenames
        class List { //: public std::vector< RefPtr<ParsedName> > {
          public:
            std::vector<int> parse_scan_check (const std::string& specifier, 
                size_t max_num_sequences = std::numeric_limits<size_t>::max());

            void scan (NameParser& parser);

            std::vector<int> count () const;

            size_t biggest_filename_size () const {
              return max_name_size;
            }

            size_t size () const { return list.size(); }

            const ParsedName& operator[] (size_t index) const { return *list[index]; }

          protected:
            std::vector< RefPtr<ParsedName> > list;
            void count_dim (std::vector<int>& dim, size_t& current_entry, size_t current_dim) const;
            size_t max_name_size;
        };



        std::string name () const {
          return filename;
        }
        size_t ndim () const {
          return indices.size();
        }
        int index (size_t num) const {
          return indices[num];
        }

        bool operator< (const ParsedName& pn) const;
        friend std::ostream& operator<< (std::ostream& stream, const ParsedName& pin);

      protected:
        std::vector<int>    indices;
        std::string         filename;

    };

  }
}

#endif

