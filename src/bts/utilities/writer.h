/*
    Copyright 2010 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close on 16/06/2010.

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

#ifndef __bts_utilities_writer_h__
#define __bts_utilities_writer_h__

#include "bts/utilities/common.h"

namespace BTS {

  namespace Utilities {

    template <typename T> class Writer {

      public:
    
        std::ofstream fout;

        Writer() {};

        Writer(const std::string& location)
        { create(location); }

        template <typename U> Writer( const std::string& location,
                                      const U& template_or_reader,
                                      const std::map<std::string,std::string>& properties = std::map<std::string,std::string>())
          { create(location); }

        template <typename U> Writer( const std::string& location,
                                      const U& template_or_reader,
                                      const Header& key_list,
                                      const std::map<std::string,std::string>& properties = std::map<std::string,std::string>())
          { create(location); }


        template <typename U> Writer( const std::string& location,
                                      const U& template_or_reader,
                                      const Header& key_list,
                                      const Header& elem_key_list,
                                      const std::map<std::string,std::string>& properties = std::map<std::string,std::string>())
          { create(location); }

        ~Writer()
          { close(); }

        //!Used in template functions, ignores properties.
        template <typename U> void create( const std::string& location,
                                      const U& template_or_reader,
                                      const std::map<std::string,std::string>& properties = std::map<std::string,std::string>())
          { create(location); }

        template <typename U> void create( const std::string& location,
                                      const U& template_or_reader,
                                      const Header& key_list,
                                      const std::map<std::string,std::string>& properties = std::map<std::string,std::string>())
          { create(location); }


        template <typename U> void create( const std::string& location,
                                      const U& template_or_reader,
                                      const Header& key_list,
                                      const Header& elem_key_list,
                                      const std::map<std::string,std::string>& properties = std::map<std::string,std::string>())
          { create(location); }



        void     create(const std::string& location);

        void     close()
          { fout.close(); }

        //!Used in template functions, ignores properties.
        void     append(const T& state, std::map<std::string,std::string>& properties_row)
          { append(state); }

        void     append(const T& state);


    };
    


	}

}


#endif /* __bts_utilities_writer_h__ */
