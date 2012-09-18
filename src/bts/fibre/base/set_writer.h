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


#ifndef __bts_fibre_setwriter_h__
#define __bts_fibre_setwriter_h__

#include <iostream>

#include "bts/fibre/base/set.h"
#include "bts/fibre/base/writer.h"
#include "bts/fibre/base/set_reader.h"
#include "bts/fibre/properties.h"


namespace BTS { 

  namespace Fibre {

    namespace Base {

      template <typename T> class SetWriter : public Fibre::Base::Writer<T> {

        // Nested classes and typedefs
        public:

        // Protected member fields.
        protected:

          std::vector<std::string>    elem_prop_hdr;

          std::ofstream               ext_elem_out;

        // Public methods and constructors.
        public:

          SetWriter() {}


          template <typename U> SetWriter (const std::string& location,
                                           const U& set_or_reader,
                                           const std::map<std::string,std::string>& file_props = Properties())
           { create (location, set_or_reader, file_props); }


          template <typename U> SetWriter (const std::string& location,
                                             const U& set_or_reader,
                                             const std::vector<std::string>& extend_prop_keys,
                                             const std::map<std::string,std::string>& file_props = Properties())
            { create( location, set_or_reader, extend_prop_keys, file_props); }


          template <typename U> SetWriter( const std::string& location,
                                           const U& set_or_reader,
                                           const std::vector<std::string>& extend_prop_keys,
                                           const std::vector<std::string>& extend_elem_prop_keys,
                                           const std::map<std::string,std::string>& file_props = Properties())
            { create( location, set_or_reader, extend_prop_keys, extend_elem_prop_keys, file_props); }


          virtual ~SetWriter()
            { close(); }


          template <typename U> void create (const std::string& location,
                                             const U& set_or_reader,
                                             const std::map<std::string,std::string>& file_props = Properties()) {
            create( location,
                    Object::select_props<T>(set_or_reader.prop_keys()),
                    Object::select_props<typename T::Element>(set_or_reader.elem_prop_keys()),
                    set_or_reader.extend_prop_keys(),
                    set_or_reader.extend_elem_prop_keys(),
                    file_props);
          }


          template <typename U> void create (const std::string& location,
                                             const U& set_or_reader,
                                             const std::vector<std::string>& extend_prop_keys,
                                             const std::map<std::string,std::string>& file_props = Properties()) {
            create( location,
                    Object::select_props<T>(set_or_reader.prop_keys()),
                    Object::select_props<typename T::Element>(set_or_reader.elem_prop_keys()),
                    extend_prop_keys,
                    set_or_reader.extend_elem_prop_keys(),
                    file_props);
          }


          template <typename U> void create (const std::string& location,
                                             const U& set_or_reader,
                                             const std::vector<std::string>& extend_prop_keys,
                                             const std::vector<std::string>& extend_elem_prop_keys,
                                             const std::map<std::string,std::string>& file_props = Properties()) {
            create( location,
                    Object::select_props<T>(set_or_reader.prop_keys()),
                    Object::select_props<typename T::Element>(set_or_reader.elem_prop_keys()),
                    extend_prop_keys,
                    extend_elem_prop_keys,
                    file_props);
          }



          void append (const T& set);

          //Not recommended as it involves a copy only included to make a template function work.
          void append (const T& set, std::map<std::string,std::string>& properties);
  
          void close ()
            { Fibre::Base::Writer<T>::close(); ext_elem_out.close(); }
    
        // Protected methods and constructors
        protected:


          void create ( const std::string& location,
                        const std::vector<const char*>& props,
                        const std::vector<const char*>& elem_props,
                        const std::vector<std::string> extend_prop_keys = std::vector<std::string>(),
                        const std::vector<std::string> extend_elem_prop_keys = std::vector<std::string>(),
                        const std::map<std::string,std::string>& file_props = Properties());

          static std::string    row_seperator(size_t count)
          { return "--- END " + str(count) + " ---"; }

      };



      //TODO: Make use of DWI::Tractography::Writers template to Remove this text writer.
      template <typename T> class SetTextWriter : public SetWriter<T> {

        //Public static variables, nested classes and typedefs

        public:

          SetTextWriter()                                                                    {}


          template <typename U> SetTextWriter (const std::string& location,
                                           const U& set_or_reader,
                                           const std::map<std::string,std::string>& file_props = Properties())
                           : SetWriter<T>(location, set_or_reader, file_props) {}


          template <typename U> SetTextWriter (const std::string& location,
                                             const U& set_or_reader,
                                             const std::vector<std::string>& extend_prop_keys,
                                             const std::map<std::string,std::string>& file_props = Properties())
                           : SetWriter<T>(location, set_or_reader, extend_prop_keys, file_props) {}


          template <typename U> SetTextWriter( const std::string& location,
                                           const U& set_or_reader,
                                           const std::vector<std::string>& extend_prop_keys,
                                           const std::vector<std::string>& extend_elem_prop_keys,
                                           const std::map<std::string,std::string>& file_props = Properties())
                           : SetWriter<T>(location, set_or_reader, extend_prop_keys, extend_elem_prop_keys, file_props) {}


          void          write (const Coord& coord);

          void          replace_end_of_file(const Coord& t, off64_t start_pos)      {}

          Coord         write_first(const Coord& t)
            {  this->out.seekp((int)this->out.tellp() - 12); write (t); return t; }

        using Writer<T>::write;


      };

      
    }

  }

}

#include "bts/fibre/base/set_writer.cpp.h"


#endif
