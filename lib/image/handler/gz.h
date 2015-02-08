/*
   Copyright 2009 Brain Research Institute, Melbourne, Australia

   Written by J-Donald Tournier, 21/08/09.

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

#ifndef __image_handler_gz_h__
#define __image_handler_gz_h__

#include "image/handler/base.h"
#include "file/mmap.h"

namespace MR
{
  namespace Image
  {

    namespace Handler
    {

      class GZ : public Base
      {
        public:
          GZ (Header& header, size_t file_header_size, bool image_is_new) :
            Base (header, image_is_new), lead_in_size (file_header_size) {
            lead_in = file_header_size ? new uint8_t [file_header_size] : NULL ;
          }

          virtual ~GZ ();
          virtual void execute ();

          uint8_t* header () {
            return (lead_in);
          }

        protected:
          int64_t  bytes_per_segment;
          uint8_t* lead_in;
          size_t   lead_in_size;
      };

    }
  }
}

#endif


