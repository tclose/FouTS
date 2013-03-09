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

#ifndef __bts_image_voxel_h__
#define __bts_image_voxel_h__

namespace BTS {
    namespace Image {
        
        template<typename T> class Voxel;
    
    }
}

#include <ostream>
#include <vector>

#include "bts/common.h"

#define LOOP(op) \
for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++) {\
  op \
}

namespace BTS {
    
    namespace Image {
        
        template<typename T> std::ostream& operator<<(std::ostream& stream, const Voxel<T>& voxel);
        
        template<typename T> class Voxel {
                
            protected:
                
                std::vector<T> intensities;

            public:
                
                Voxel() {
                }
                
                Voxel(size_t num_encodings)
                        : intensities(num_encodings) {
                }
                
                Voxel(size_t num_encodings, const T& default_value)
                        : intensities(num_encodings, default_value) {
                }
                
                Voxel(const Voxel& voxel)
                        : intensities(voxel.intensities) {
                }
                
                Voxel& operator=(const Voxel& voxel) {
                    intensities = voxel.intensities;
                    return *this;
                }
                
                template<typename U> Voxel& operator=(const U& vox) {
                    return operator=((const Voxel&) vox);
                }
                
                virtual ~Voxel() {
                }
                
                void resize(size_t new_size, const T& default_value = T()) {
                    intensities.resize(new_size, default_value);
                }
                
                T& operator[](size_t encode_i) {
                    return intensities[encode_i];
                }
                
                const T& operator[](size_t encode_i) const {
                    return intensities[encode_i];
                }
                
                Voxel operator-() const {
                    return Voxel(*this).negate();
                }
                
                Voxel& negate() {
                    LOOP(operator[](encode_i) = -operator[](encode_i);)
                    return *this;
                }
                
                Voxel& operator+=(const Voxel& v) {
                    LOOP(operator[](encode_i) += v[encode_i];)
                    return *this;
                }
                
                Voxel& operator+=(T f) {
                    LOOP(operator[](encode_i) += f;)
                    return *this;
                }
                
                Voxel& operator-=(const Voxel& v) {
                    LOOP(operator[](encode_i) -= v[encode_i];)
                    return *this;
                }
                
                Voxel& operator-=(T f) {
                    LOOP(operator[](encode_i) -= f;)
                    return *this;
                }
                
                Voxel& operator*=(T f) {
                    LOOP(operator[](encode_i) *= f;)
                    return *this;
                }
                
                Voxel& operator/=(T f) {
                    LOOP(operator[](encode_i) /= f;)
                    return *this;
                }
                
                Voxel operator+(const Voxel& v) const {
                    Voxel answer(*this);
                    answer += v;
                    return answer;
                }
                
                Voxel operator+(T f) const {
                    Voxel answer(*this);
                    answer += f;
                    return answer;
                }
                
                Voxel operator-(const Voxel& v) const {
                    Voxel answer(*this);
                    answer -= v;
                    return answer;
                }
                
                Voxel operator-(T f) const {
                    Voxel answer(*this);
                    answer -= f;
                    return answer;
                }
                
                Voxel operator*(T f) const {
                    Voxel answer(*this);
                    answer *= f;
                    return answer;
                }
                
                Voxel operator/(T f) const {
                    Voxel answer(*this);
                    answer /= f;
                    return answer;
                }
                
                Voxel& zero() {
                    LOOP(BTS::zero(operator[](encode_i));)
                    return *this;
                }
                
                size_t num_encodings() const {
                    return intensities.size();
                }
                
                size_t size() const {
                    return intensities.size();
                }
                
                friend std::ostream& operator<<<>(std::ostream& stream, const Voxel& voxel);
                
        };
        
        template<typename T> Voxel<T> operator+(T f, const Voxel<T>& v) {
            Voxel<T> answer(v);
            answer += f;
            return answer;
        }
        
        template<typename T> Voxel<T> operator*(T f, const Voxel<T>& v) {
            Voxel<T> answer(v);
            answer *= f;
            return answer;
        }
        
        template<typename T> std::ostream& operator<<(std::ostream& stream, const Voxel<T>& voxel) {
            
            for (size_t encode_i = 0; encode_i < voxel.num_encodings(); encode_i++) {
                stream << voxel[encode_i] << ", ";
            }
            stream << std::endl;
            return stream;
            
        }
    
    }

}

#undef LOOP

#endif
