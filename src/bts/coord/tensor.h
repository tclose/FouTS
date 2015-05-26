/*
 Copyright 2012 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 07/08/2010.

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

#ifndef __bts_coord_tensor_h__
#define __bts_coord_tensor_h__

#include "math/matrix.h"
#include "bts/coord.h"
#include "bts/fibre/base/tensor_writer.h"

namespace FTS {
    
    class Coord::Tensor: public MR::Math::Matrix<double> {
            
        public:
            
            typedef Fibre::Base::TensorWriter<Coord::Tensor> Writer;
            typedef Coord RowType;

        public:
            
            static Tensor create_diagonal(double value) {
                Tensor t;
                t.diagonal() = value;
                return t;
            }
            
        public:
            
            Tensor()
                    : MR::Math::Matrix<double>(3, 3) {
            }
            
            //Only included for template consistency in gradient tester
            Tensor(const Coord& c)
                    : MR::Math::Matrix<double>(3, 3) {
                invalidate();
            }
            
            Tensor(const MR::Math::Matrix<double>::View& view)
                    : MR::Math::Matrix<double>(view) {
                assert(rows()==3);
                assert(columns() == 3);
            }
            
            Coord operator[](size_t idx) {
                return Coord(row(idx));
            }
            
            const Coord operator[](size_t idx) const {
                return Coord(row(idx));
            }
            
            Tensor& invalidate() {
                for (size_t i = 0; i < 3; ++i)
                    for (size_t j = 0; j < 3; ++j)
                        operator()(i, j) = NAN;
                return *this;
            }
            
            Tensor& zero() {
                for (size_t i = 0; i < 3; ++i)
                    for (size_t j = 0; j < 3; ++j)
                        operator()(i, j) = NAN;
                return *this;
            }
            
            Tensor& operator*=(const Tensor& t) {
                MR::Math::Matrix<double>::operator*=(t);
                return *this;
            }
            
            Tensor& operator*=(double c) {
                MR::Math::Matrix<double>::operator*=(c);
                return *this;
            }
            
            Tensor& operator/=(const Tensor& t) {
                MR::Math::Matrix<double>::operator/=(t);
                return *this;
            }
            
            Tensor& operator/=(double c) {
                MR::Math::Matrix<double>::operator/=(c);
                return *this;
            }
            
            Tensor operator*(const Tensor& t) const {
                Tensor answer(*this);
                answer *= t;
                return answer;
            }
            
            Tensor operator*(double c) const {
                Tensor answer(*this);
                answer *= c;
                return answer;
            }
            
            Tensor operator/(const Tensor& t) const {
                Tensor answer(*this);
                answer /= t;
                return answer;
            }
            
            Tensor operator/(double c) const {
                Tensor answer(*this);
                answer /= c;
                return answer;
            }
            
            Tensor& operator+=(const Tensor& t) {
                MR::Math::Matrix<double>::operator+=(t);
                return *this;
            }
            
            Tensor& operator-=(const Tensor& t) {
                MR::Math::Matrix<double>::operator-=(t);
                return *this;
            }
            
            Tensor operator+(const Tensor& t) const {
                Tensor answer(*this);
                answer += t;
                return answer;
            }
            
            Tensor operator-(const Tensor& t) const {
                Tensor answer(*this);
                answer -= t;
                return answer;
            }
            
            double vsize() const {
                return 3;
            }
            
    };
    
    inline Coord::Tensor operator*(double c, const Coord::Tensor& t) {
        return t * c;
    }

}

#endif /*  #ifndef __bts_fibre_coord_tensor_h__*/

