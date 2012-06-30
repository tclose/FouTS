/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 10/08/2010.

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

#ifndef __bts_fibre_tensor_writer2_h__
#define __bts_fibre_tensor_writer2_h__

namespace BTS {

	namespace Fibre {

	  namespace Tensor {

      template <typename T> class Writer2 {

        //Public static variables, nested classes and typedefs
        public:

          static const std::string     HEADER_PREAMBLE;

        //Protected member variables
        protected:

          std::vector<std::string>    header;
          size_t                  template_length;
          std::ofstream         out;

        //Public member functions
        public:

          Writer2() {}

          Writer2(const std::string& location, const T& template_value = T())
            { create(location, template_value); }

          Writer2(const std::string& location, size_t size)
            { create(location, size); }

          ~Writer2()
            { close(); }

          void close()
            { out.close(); }

          void create (const std::string& location, size_t size)
            { create(location, T(size)); }

          void create (const std::string& location, const T& template_value);

          void append(const MR::Math::Matrix<double>& tensor);


      };

	  }

	}

}

#endif /* __bts_fibre_tensor_writer2_h__ */
