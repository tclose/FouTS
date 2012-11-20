/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 16/07/2010.

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

#ifndef __bts_observed_buffer_cpp_h__
#define __bts_observed_buffer_cpp_h__

#include "dataset/loop.h"

namespace BTS {

  namespace Image {

    namespace Observed {

      template <typename T> Buffer_tpl<T>::Buffer_tpl(Triple<size_t> dimensions,
                                                  Triple<double> voxel_lengths,
                                                  Triple<double> corner_offsets,
                                                  bool enforce_bounds)

        : Image::Buffer_tpl<T>(dimensions, enforce_bounds),
          voxel_lengths(voxel_lengths),
          corner_offsets(corner_offsets) {

        if (!corner_offsets)
          this->corner_offsets =  default_corner_offset(dimensions, voxel_lengths);

      }


      template <typename T> void                        Buffer_tpl<T>::reset(const Triple<size_t>& dimensions,
                                                                                      const Triple<double>& voxel_lengths,
                                                                                      const Triple<double>& corner_offsets) {



        this->Image::Buffer_tpl<T>::reset(dimensions);
        this->voxel_lengths = voxel_lengths;
        this->corner_offsets = corner_offsets;

        if (!this->corner_offsets)
          this->corner_offsets =  default_corner_offset(dimensions, voxel_lengths);


      }




      template <typename T> std::ostream&         Buffer_tpl<T>::to_stream (std::ostream& stream) const {

        stream << "Dims: " << this->dims() << std::endl;

        stream << "Vox lengths: " << voxel_lengths << std::endl;

        stream << "Corner offsets: " << corner_offsets << std::endl;


        for (size_t z = 0; z < this->dim(Z); z++) {

          stream  << std::endl  << std::endl << "----------" << std::endl
                  << "   z=" << z << std::endl
                  << "----------" << std::endl;

          for (size_t y = 0; y < this->dim(Y); y++) {

            stream << std::endl << "y=" << y << std::endl << "----" << std::endl;

            for (size_t x = 0; x < this->dim(X); x++) {

              stream << "x=" << x << ": " << this->operator()(x,y,z) << std::endl;
            }
          }
        }


        return stream;

      }


      template <typename T> std::ostream& operator<< (std::ostream& stream, const Buffer_tpl<T>& buffer) {

        return buffer.to_stream(stream);

      }


      template <typename T> void                            Buffer_tpl<T>::save(const std::string& location) const {

        MR::Image::Header header;

        header.set_ndim(4);

        header.set_vox(X, this->vox_length(X));
        header.set_vox(Y, this->vox_length(Y));
        header.set_vox(Z, this->vox_length(Z));

        header.set_dim(X, this->dim(X));
        header.set_dim(Y, this->dim(Y));
        header.set_dim(Z, this->dim(Z));
        header.set_dim(DW, num_encodings());

        MR::Math::Matrix<float> transform(4,4);

        transform.identity();

        transform(0,3) = offset(X) + 0.5 * vox_length(X);
        transform(1,3) = offset(Y) + 0.5 * vox_length(Y);
        transform(2,3) = offset(Z) + 0.5 * vox_length(Z);

        header.set_transform(transform);

        MR::Math::Matrix<float> encodings_matrix (num_encodings(), 4);

        for (size_t encode_i = 0; encode_i < num_encodings(); encode_i++) {
          const Diffusion::Encoding& encode = this->encoding(encode_i);
          encodings_matrix(encode_i,X)  = encode[X];
          encodings_matrix(encode_i,Y)  = encode[Y];
          encodings_matrix(encode_i,Z)  = encode[Z];
          encodings_matrix(encode_i,DW) = encode.b_value();
        }

        header.set_DW_scheme(encodings_matrix);

        header.insert(properties().begin(), properties().end());

        File::clear_path(location);

        header.create(location);

        //Copy data from buffer to image.
        MR::Image::Voxel<double> datapoint (header);

        std::vector<size_t> loop_order(4);

        loop_order[0] = DW;
        loop_order[1] = X;
        loop_order[2] = Y;
        loop_order[3] = Z;

        MR::DataSet::LoopInOrder loop (loop_order);

        loop.start(datapoint);

        for (size_t z = 0; z < this->dim(Z); z++)
          for (size_t y = 0; y < this->dim(Y); y++)
            for (size_t x = 0; x < this->dim(X); x++) {

              Image::Voxel<double> vox;

              if (!this->is_empty(x,y,z)) {
                vox = this->operator()(x,y,z);
              }else {
                vox = Image::Voxel<double> (num_encodings());
              }

              for (size_t encode_i = 0; encode_i < num_encodings(); encode_i++) {
                datapoint.value() = vox[encode_i];
                loop.next(datapoint);
              }

            }

      }




      template <typename T> double                            Buffer_tpl<T>::rms(bool include_b0s) const {

        double sum = 0.0;

        for (typename Image::Buffer_tpl<T>::const_iterator vox_it = this->begin(); vox_it != this->end(); ++vox_it) {

          if (!include_b0s || vox_it->second.b0())
            sum += vox_it->second.rms();

        }

        return sum / this->num_voxels_in_bounds();

      }


      template <typename T> double                            Buffer_tpl<T>::max_b0() const {

        double max_b0 = 0;

        for (typename Image::Buffer_tpl<T>::const_iterator vox_it = this->begin(); vox_it != this->end(); ++vox_it)

          if (vox_it->second.b0() > max_b0)
            max_b0 = vox_it->second.b0();

        return max_b0;

      }


    }

  }

}

#endif /* __observed_buffer_tpl_cpp_h__ */
