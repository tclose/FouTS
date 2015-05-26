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

#ifndef __image_buffer_tpl_cpp_h__
#define __image_buffer_tpl_cpp_h__

#include "image/header.h"
#include "dataset/loop.h"

#define LOOP(op) \
for (int z = 0; z < dims[Z]; z++) { \
  for (int y = 0; y < dims[Y]; y++) { \
    for (int x = 0; x < dims[X]; x++) { \
      op \
    } \
  } \
}

namespace FTS {
    
    namespace Image {
        
        template<typename T> T& Buffer_tpl<T>::operator()(const Index& c) {
            
            if (enforce_bounds && !in_bounds(c))
                throw Exception(
                        "Coordinate " + str(c) + " is out of image bounds " + str(dimensions)
                        + " (set 'enforce_bounds=false' if this is intentional).");
            
            iterator it = voxels.find(c);
            
            if (it == voxels.end())
                return voxels[c] = new_voxel(c);
            else
                return it->second;
            
        }
        
        template<typename T> const T& Buffer_tpl<T>::operator()(const Index& c) const {
            
            if (enforce_bounds && !in_bounds(c))
                throw Exception(
                        "Coordinate " + str(c) + " is out of image bounds " + str(dimensions)
                        + " (set 'enforce_bounds=false' if this is intentional).");
            
            const_iterator it = voxels.find(c);
            
            if (it == voxels.end())
                throw Exception(
                        "Coordinate " + str(c)
                        + " has not been initialised and can't be since this buffer is const. Suggest using 'is_empty(const Coord&)' to check if it has been initialised first.");
            
            return it->second;
            
        }
        
        template<typename T> bool Buffer_tpl<T>::is_empty(const Index& c) const {
            
            return !voxels.count(c);
            
        }
        
        template<typename T> template<typename U> Buffer_tpl<T>& Buffer_tpl<T>::operator+=(
                const Buffer_tpl<U>& buff) {
            
            if (dims() != buff.dims())
                throw Exception(
                        "Buffer dimensions do not match, " + str(buff.dims()) + " and "
                        + str(dims()) + ".");
            
            for (typename Buffer_tpl<U>::const_iterator vox_it = buff.begin(); vox_it != buff.end();
                    ++vox_it) {
                if (!is_empty(vox_it->first))
                    this->operator()(vox_it->first) = vox_it->second;
                else
                    this->operator()(vox_it->first) += vox_it->second;
            }
            
            return *this;
            
        }
        
        template<typename T> template<typename U> Buffer_tpl<T>& Buffer_tpl<T>::operator-=(
                const Buffer_tpl<U>& buff) {
            
            if (dims() != buff.dims())
                throw Exception(
                        "Buffer dimensions do not match, " + str(buff.dims()) + " and "
                        + str(dims()) + ".");
            
            for (typename Buffer_tpl<U>::const_iterator vox_it = buff.begin(); vox_it != buff.end();
                    ++vox_it) {
                
                if (!is_empty(vox_it->first)) {
                    this->operator()(vox_it->first) = vox_it->second;
                    negate(this->operator()(vox_it->first));
                } else
                    this->operator()(vox_it->first) -= vox_it->second;
                
            }
            
            return *this;
            
        }
        
        template<typename T> Buffer_tpl<T>& Buffer_tpl<T>::operator+=(const Buffer_tpl<T>& buff) {
            
            if (dims() != buff.dims())
                throw Exception(
                        "Buffer dimensions do not match, " + str(buff.dims()) + " and "
                        + str(dims()) + ".");
            
            for (const_iterator vox_it = buff.begin(); vox_it != buff.end(); ++vox_it)
                this->operator()(vox_it->first) += vox_it->second;
            
            return *this;
            
        }
        
        template<typename T> Buffer_tpl<T>& Buffer_tpl<T>::operator-=(const Buffer_tpl<T>& buff) {
            
            if (dims() != buff.dims())
                throw Exception(
                        "Buffer dimensions do not match, " + str(buff.dims()) + " and "
                        + str(dims()) + ".");
            
            for (const_iterator vox_it = buff.begin(); vox_it != buff.end(); ++vox_it)
                this->operator()(vox_it->first) -= vox_it->second;
            
            return *this;
            
        }
        
        template<typename T> template<typename U> Buffer_tpl<T>& Buffer_tpl<T>::operator*=(
                const Buffer_tpl<U>& buff) {
            
            if (dims() != buff.dims())
                throw Exception(
                        "Buffer dimensions do not match, " + str(buff.dims()) + " and "
                        + str(dims()) + ".");
            
            for (typename Buffer_tpl<U>::const_iterator vox_it = buff.begin(); vox_it != buff.end();
                    ++vox_it) {
                
                if (!is_empty(vox_it->first))
                    this->operator()(vox_it->first) *= vox_it->second;
                
            }
            
            return *this;
            
        }
        
        template<typename T> template<typename U> Buffer_tpl<T>& Buffer_tpl<T>::operator/=(
                const Buffer_tpl<U>& buff) {
            
            if (dims() != buff.dims())
                throw Exception(
                        "Buffer dimensions do not match, " + str(buff.dims()) + " and "
                        + str(dims()) + ".");
            
            for (typename Buffer_tpl<U>::const_iterator vox_it = buff.begin(); vox_it != buff.end();
                    ++vox_it) {
                
                if (!is_empty(vox_it->first))
                    this->operator()(vox_it->first) /= vox_it->second;
                
            }
            
            return *this;
            
        }
        
        template<typename T> Buffer_tpl<T>& Buffer_tpl<T>::negate() {
            
            for (typename Buffer_tpl<T>::iterator vox_it = this->begin(); vox_it != this->end();
                    ++vox_it)
                vox_it->second.negate();
            
            return *this;
            
        }
        
        template<typename T> Buffer_tpl<T> Buffer_tpl<T>::operator+(
                const Buffer_tpl<T>& buff) const {
            
            Buffer_tpl<T> answer;
            
            if (this->num_not_empty_voxels() > buff.num_not_empty_voxels()) {
                answer = *this;
                answer += buff;
            } else {
                answer = buff;
                answer += *this;
            }
            
            return answer;
            
        }
        
        template<typename T> Buffer_tpl<T> Buffer_tpl<T>::operator-(
                const Buffer_tpl<T>& buff) const {
            
            Buffer_tpl<T> answer;
            
            if (this->num_not_empty_voxels() > buff.num_not_empty_voxels()) {
                answer = *this;
                answer -= buff;
            } else {
                answer = buff;
                answer -= *this;
            }
            
            return answer;
            
        }
        
        template<typename T> Buffer_tpl<T> Buffer_tpl<T>::operator*(
                const Buffer_tpl<T>& buff) const {
            
            Buffer_tpl<T> answer;
            
            if (this->num_not_empty_voxels() > buff.num_not_empty_voxels()) {
                answer = *this;
                answer *= buff;
            } else {
                answer = buff;
                answer *= *this;
            }
            
            return answer;
            
        }
        
        template<typename T> Buffer_tpl<T> Buffer_tpl<T>::operator/(
                const Buffer_tpl<T>& buff) const {
            
            Buffer_tpl<T> answer;
            
            if (this->num_not_empty_voxels() > buff.num_not_empty_voxels()) {
                answer = *this;
                answer /= buff;
            } else {
                answer = buff;
                answer /= *this;
            }
            
            return answer;
            
        }
        
        template<typename T> template<typename U> Buffer_tpl<T>& Buffer_tpl<T>::operator*=(
                const U& M) {
            
            for (iterator vox_it = begin(); vox_it != end(); ++vox_it)
                vox_it->second *= M;
            
            return *this;
            
        }
        
        template<typename T> template<typename U> Buffer_tpl<T>& Buffer_tpl<T>::operator/=(
                const U& M) {
            
            for (iterator vox_it = begin(); vox_it != end(); ++vox_it)
                vox_it->second /= M;
            
            return *this;
            
        }
        
        template<typename T> std::ostream& operator<<(std::ostream& stream,
                                                      const Image::Buffer_tpl<T>& B) {
            
            stream << "Dims: " << B.dims();
            
            if (B.enforce_bounds)
                stream << " (enforced)";
            
            stream << std::endl;
            
            for (typename Image::Buffer_tpl<T>::const_iterator vox_it = B.begin();
                    vox_it != B.end(); ++vox_it)
                stream << vox_it->first << ": " << vox_it->second;
            
            return stream;
            
        }
        
        template<typename T> std::set<Index> Buffer_tpl<T>::non_empty() const {
            
            std::set<Index> nempty;
            
            for (const_iterator vox_it = begin(); vox_it != end(); ++vox_it)
                nempty.insert(vox_it->first);
            
            return nempty;
            
        }
        
        template<typename T> std::set<Index> Buffer_tpl<T>::non_empty_or_inbounds() const {
            
            std::set<Index> nempty;
            
            for (const_iterator vox_it = begin(); vox_it != end(); ++vox_it)
                nempty.insert(vox_it->first);
            
            for (size_t x = 0; x < dim(X); ++x)
                for (size_t y = 0; y < dim(Y); ++y)
                    for (size_t z = 0; z < dim(Z); ++z)
                        nempty.insert(Index(x, y, z));
            
            return nempty;
            
        }
        
        template<typename T> std::set<Index> Buffer_tpl<T>::empty_inbounds() const {
            
            std::set<Index> empty;
            
            for (size_t x = 0; x < dim(X); ++x)
                for (size_t y = 0; y < dim(Y); ++y)
                    for (size_t z = 0; z < dim(Z); ++z)
                        if (is_empty(x, y, z))
                            empty.insert(Index(x, y, z));
            
            return empty;
            
        }
        
        template<typename T> void Buffer_tpl<T>::clear_and_enforce_bounds() {
            
            clear();
            
            enforce_bounds = true;
            
        }
    
        template<typename T> void Buffer_tpl<T>::load(const std::string& location) {

            MR::Image::Header header(location);

            if (header.ndim() != 3)
               throw Exception( "Image should contain 3 dimensions, found " + str(header.ndim()) + ".");

            //Get image dimensions and voxel lengths.
            Triple<size_t> dimensions(header.dim(X), header.dim(Y), header.dim(Z));
            Triple<double> voxel_lengths(header.vox(X), header.vox(Y), header.vox(Z));

            //Resize buffer to fit image.
            reset(dimensions, true);

            //Copy data from image to buffer.
            MR::Image::Voxel<double> voxel(header);
            MR::DataSet::Loop loop(0, 3);

            for (loop.start(voxel); loop.ok(); loop.next(voxel))
               operator()(voxel[X], voxel[Y], voxel[Z]) = voxel.value();

        }

    }

}
#undef LOOP

#endif /* __image_buffer_tpl_cpp_h__ */
