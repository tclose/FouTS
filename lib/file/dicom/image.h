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

#ifndef __file_dicom_image_h__
#define __file_dicom_image_h__

#include "ptr.h"
#include "data_type.h"
#include "math/vector.h"
#include "file/dicom/element.h"

namespace MR
{
  namespace File
  {
    namespace Dicom
    {

      class Series;
      class Element;

      class Image
      {

        public:
          Image (Series* parent = NULL) :
            series (parent) {
            acq_dim[0] = acq_dim[1] = dim[0] = dim[1] = instance = acq = sequence = std::numeric_limits<size_t>::max();
            position_vector[0] = position_vector[1] = position_vector[2] = NAN;
            orientation_x[0] = orientation_x[1] = orientation_x[2] = NAN;
            orientation_y[0] = orientation_y[1] = orientation_y[2] = NAN;
            orientation_z[0] = orientation_z[1] = orientation_z[2] = NAN;
            distance = NAN;
            pixel_size[0] = pixel_size[1] = slice_thickness = slice_spacing = NAN;
            scale_intercept = 0.0;
            scale_slope = 1.0;
            bvalue = G[0] = G[1] = G[2] = NAN;
            data = bits_alloc = images_in_mosaic = data_size = 0;
            is_BE = false;
          }

          std::string  filename;
          std::string  sequence_name;
          std::string  manufacturer;
          Series*      series;

          size_t   acq_dim[2], dim[2], instance, acq, sequence;
          float    position_vector[3], orientation_x[3], orientation_y[3], orientation_z[3], distance;
          float    pixel_size[2], slice_thickness, slice_spacing, scale_slope, scale_intercept;
          float    bvalue, G[3];
          size_t   data, bits_alloc, images_in_mosaic, data_size;
          DataType data_type;
          bool     is_BE;

          void read ();
          void parse_item (Element& item, const std::string& dirname = "");
          void decode_csa (const uint8_t* start, const uint8_t* end);

          bool operator< (const Image& ima) const {
            if (acq != ima.acq) 
              return acq < ima.acq;
            assert (!isnan (distance));
            assert (!isnan (ima.distance));
            if (distance != ima.distance) 
              return distance < ima.distance;
            if (sequence != ima.sequence) 
              return sequence < ima.sequence;
            if (instance != ima.instance) 
              return instance < ima.instance;
            return false;
          }

          void print_fields (bool dcm, bool csa) const;

          void calc_distance () {
            if (images_in_mosaic) {
              float xinc = pixel_size[0] * float (dim[0] - acq_dim[0]) / 2.0;
              float yinc = pixel_size[1] * float (dim[1] - acq_dim[1]) / 2.0;
              for (size_t i = 0; i < 3; i++)
                position_vector[i] += xinc * orientation_x[i] + yinc * orientation_y[i];

              float normal[3];
              Math::cross (normal, orientation_x, orientation_y);
              if (Math::dot (normal, orientation_z) < 0.0) {
                orientation_z[0] = -normal[0];
                orientation_z[1] = -normal[1];
                orientation_z[2] = -normal[2];
              }
              else {
                orientation_z[0] = normal[0];
                orientation_z[1] = normal[1];
                orientation_z[2] = normal[2];
              }

            }
            else Math::cross (orientation_z, orientation_x, orientation_y);

            Math::normalise (orientation_z);
            distance = Math::dot (orientation_z, position_vector);
          }


      };

      std::ostream& operator<< (std::ostream& stream, const Image& item);


    }
  }
}


#endif


