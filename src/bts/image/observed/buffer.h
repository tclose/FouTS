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

#ifndef __bts_image_observed_buffer_h__
#define __bts_image_observed_buffer_h__

namespace BTS {
    
    namespace Image {
        
        namespace Observed {
            
            class Buffer;
        
        }
    
    }

}

#include "image/header.h"
#include "image/voxel.h"

#include "bts/triple.h"
#include "bts/image/buffer.h"
#include "bts/diffusion/encoding/set.h"
#include "bts/image/observed/buffer.h"

#include "bts/image/properties.h"

#include "bts/common.h"

#define LOOP(op) \
for (size_t z = 0; z < this->dim(Z); z++) { \
  for (size_t y = 0; y < this->dim(Y); y++) { \
    for (size_t x = 0; x < this->dim(X); x++) { \
      op \
    } \
  } \
}

namespace BTS {
    
    namespace Image {
        
        namespace Observed {
            
            template<typename T> class Buffer_tpl: public Image::Buffer_tpl<T> {
                    
                    //Public static functions.
                public:
                    
                    static Triple<double> default_corner_offset(
                            const Triple<size_t>& dimensions, const Triple<double>& voxel_lengths) {
                        return -voxel_lengths * dimensions * 0.5;
                    }
                    
                    // Returns the centre of the voxel in buffer space.
                    inline static Triple<double> voxel_centre(Index coord) {
                        return coord.to_point() + Triple<double>::Halves;
                    }
                    
                    // Returns the voxel coordinate that the point given in buffer space lies in.
                    static Index voxel_coord(Triple<double> point) {
                        return Index((int) point[X], (int) point[Y], (int) point[Z]);
                    }
                    
                protected:
                    
                    Triple<double> voxel_lengths;
                    Triple<double> corner_offsets;    // Spatial offsets of bottom-left-back corner of the image.
                    
                    Properties props;

                public:
                    
                    ~Buffer_tpl() {
                    }
                    
                    const Triple<double>& vox_lengths() const {
                        return voxel_lengths;
                    }
                    
                    double vox_length(size_t dim_index) const {
                        return voxel_lengths[dim_index];
                    }
                    
                    //! Offset of the lowest corner of the lowest voxel coordinate ([0,0,0]) of the image.
                    const Triple<double>& offsets() const {
                        return corner_offsets;
                    }
                    
                    //! Offset of the lowest corner of the lowest voxel coordinate ([0,0,0]) of the image.
                    double offset(size_t dim_index) const {
                        return corner_offsets[dim_index];
                    }
                    
                    std::ostream& to_stream(std::ostream& stream) const;

                    Properties& properties() {
                        return props;
                    }
                    
                    const Properties& properties() const {
                        return props;
                    }
                    
                    virtual void save(const std::string& location) const;

                    void reset(const Triple<size_t>& dims, const Triple<double>& voxel_sizes,
                               const Triple<double>& spatial_offset);

                    double rms(bool include_b0s) const;

                    virtual size_t num_encodings() const = 0;

                    virtual const Diffusion::Encoding& encoding(size_t index) const = 0;

                    double max_b0() const;

                protected:
                    
                    Buffer_tpl(bool enforce_bounds)
                            : Image::Buffer_tpl<T>(enforce_bounds) {
                    }
                    
                    Buffer_tpl(Triple<size_t> dimensions, Triple<double> voxel_lengths,
                               Triple<double> corner_offsets, bool enforce_bounds);

                    Buffer_tpl(const Buffer_tpl<T>& B)
                            : Image::Buffer_tpl<T>(B), voxel_lengths(B.voxel_lengths), corner_offsets(
                                      B.corner_offsets), props(B.props) {
                    }
                    
                    Buffer_tpl& operator=(const Buffer_tpl<T>& B) {
                        this->Image::Buffer_tpl<T>::operator=(B);
                        voxel_lengths = B.voxel_lengths;
                        corner_offsets = B.corner_offsets;
                        props = B.props;
                        return *this;
                    }
                    
            };
        
        }
    
    }

}

#include "bts/image/observed/voxel.h"

namespace BTS {
    
    namespace Image {
        
        namespace Observed {
            
            class Buffer: public Buffer_tpl<Voxel> {
                    
                protected:
                    
                    Diffusion::Encoding::Set encodings;

                public:
                    
                    Buffer()
                            : Buffer_tpl<Voxel>(true) {
                    }
                    
                    //TODO: Make default constructor, and open functions.
                    Buffer(const Triple<size_t>& dimensions, const Triple<double>& voxel_sizes,
                           const Triple<double>& spatial_offset,
                           const Diffusion::Encoding::Set& encodings)
                            : Buffer_tpl<Voxel>(dimensions, voxel_sizes, spatial_offset, true), encodings(
                                      encodings) {
                    }
                    
                    Buffer(const std::string& image_location,
                           const Diffusion::Encoding::Set& default_encodings =
                                   Diffusion::Encoding::Set())
                            : Buffer_tpl<Voxel>(true) {
                        load(image_location, default_encodings);
                    }
                    
                    Buffer(const Buffer& B)
                            : Buffer_tpl<Voxel>(B), encodings(B.encodings) {
                    }
                    
                    virtual ~Buffer() {
                    }
                    
                    Buffer& operator=(const Buffer& B) {
                        Buffer_tpl<Voxel>::operator=(B);
                        encodings = B.encodings;
                        return *this;
                    }
                    
                    void reset(const Triple<size_t>& dims, const Triple<double>& voxel_sizes,
                               const Triple<double>& spatial_offset,
                               const Diffusion::Encoding::Set& encodings)

                               {
                        Buffer_tpl<Voxel>::reset(dims, voxel_sizes, spatial_offset);
                        this->encodings = encodings;
                    }
                    
                    void load(const std::string& location,
                              const Diffusion::Encoding::Set& default_encodings =
                                      Diffusion::Encoding::Set());

                    Buffer* clone() const {
                        return new Buffer(*this);
                    }
                    
                    Voxel& operator()(Index coord) {
                        return Image::Buffer_tpl<Voxel>::operator()(coord);
                    }
                    
                    Voxel& operator()(int x, int y, int z) {
                        return Image::Buffer_tpl<Voxel>::operator()(x, y, z);
                    }
                    
                    const Voxel& operator()(Index coord) const {
                        return Image::Buffer_tpl<Voxel>::operator()(coord);
                    }
                    
                    const Voxel& operator()(int x, int y, int z) const {
                        return Image::Buffer_tpl<Voxel>::operator()(x, y, z);
                    }
                    
                    size_t num_encodings() const {
                        return encodings.num_encodings();
                    }
                    
                    const Diffusion::Encoding& encoding(size_t index) const {
                        return encodings[index];
                    }
                    
                    const Diffusion::Encoding::Set& get_encodings() const {
                        return encodings;
                    }
                    
                    size_t dim(size_t dim_index) const {
                        return Image::Buffer_tpl<Voxel>::dim(dim_index);
                    }
                    
                    double vox_length(size_t dim_index) const {
                        return Buffer_tpl<Voxel>::vox_length(dim_index);
                    }
                    
                    //! Offset of the lowest corner of the lowest voxel coordinate ([0,0,0]) of the image.
                    double offset(size_t dim_index) const {
                        return Buffer_tpl<Voxel>::offset(dim_index);
                    }
                    
                    const Triple<size_t>& dims() const {
                        return Image::Buffer_tpl<Voxel>::dims();
                    }
                    
                    const Triple<double>& vox_lengths() const {
                        return Buffer_tpl<Voxel>::vox_lengths();
                    }
                    
                    //! Offset of the lowest corner of the lowest voxel coordinate ([0,0,0]) of the image.
                    const Triple<double>& offsets() const {
                        return Buffer_tpl<Voxel>::offsets();
                    }
                    
                    Buffer& zero() {
                        this->Image::Buffer_tpl<Voxel>::zero();
                        return *this;
                    }
                    
                    Double::Buffer b0() const;

                    double average_b0() const;

                    double max_b0() const {
                        return Buffer_tpl<Voxel>::max_b0();
                    }
                    
                    Double::Buffer isotropic(bool include_b0s) const;

                    double rms(bool include_b0s) const {
                        return Buffer_tpl<Voxel>::rms(include_b0s);
                    }
                    
                    double mean_square(bool include_b0s) const;

                    virtual void remove_isotropic();

                    Buffer& operator+=(const Buffer& buff) {
                        this->Image::Buffer_tpl<Voxel>::operator+=(buff);
                        return *this;
                    }
                    
                    Buffer& operator-=(const Buffer& buff) {
                        this->Image::Buffer_tpl<Voxel>::operator-=(buff);
                        return *this;
                    }
                    
                    Buffer& operator*=(double M) {
                        this->Image::Buffer_tpl<Voxel>::operator*=(M);
                        return *this;
                    }
                    
                    Buffer& operator/=(double M) {
                        this->Image::Buffer_tpl<Voxel>::operator/=(M);
                        return *this;
                    }
                    
                    template<typename T> Buffer& operator+=(const Image::Buffer_tpl<T>& buff) {
                        this->Image::Buffer_tpl<Voxel>::operator+=(buff);
                        return *this;
                    }
                    
                    template<typename T> Buffer& operator-=(const Image::Buffer_tpl<T>& buff) {
                        this->Image::Buffer_tpl<Voxel>::operator-=(buff);
                        return *this;
                    }
                    
                    std::ostream& to_stream(std::ostream& stream) const {
                        return this->Buffer_tpl<Voxel>::to_stream(stream);
                    }
                    
                    Properties& properties() {
                        return this->Buffer_tpl<Voxel>::properties();
                    }
                    
                    const Properties& properties() const {
                        return this->Buffer_tpl<Voxel>::properties();
                    }
                    
                protected:
                    
                    virtual Voxel new_voxel(const Index& c) {
                        return Voxel(*this, c);
                    }
                    
            };
            
            std::ostream& operator<<(std::ostream& stream, const Buffer& buffer);
        
        }
    
    }

}

#undef LOOP

#endif
