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

#include "bts/coord.h"
#include "bts/coord/tensor.h"

namespace BTS {
    
    const Coord Coord::Invalid(NAN, NAN, NAN);
    const Coord Coord::Zeros(0.0, 0.0, 0.0);
    const Coord Coord::Ones(1.0, 1.0, 1.0);
    const Coord Coord::Halves(0.5, 0.5, 0.5);
    const std::string Coord::FILE_EXTENSION = "crd";
    
    std::ostream& operator<<(std::ostream& stream, const Coord& c) {
        stream << "[ " << c[X] << ", " << c[Y] << ", " << c[Z] << " ]";
        return (stream);
    }
    
    Coord& Coord::zero() {
        
        operator[](X) = 0.0;
        operator[](Y) = 0.0;
        operator[](Z) = 0.0;
        
        return *this;
    }
    
    Coord& Coord::invalidate() {
        
        operator[](X) = NAN;
        operator[](Y) = NAN;
        operator[](Z) = NAN;
        
        return *this;
    }
    
    bool Coord::operator==(const Coord& A) const {
        return (operator[](X) == A[X]) && (operator[](Y) == A[Y]) && (operator[](Z) == A[Z]);
    }
    ;
    
    bool Coord::operator!=(const Coord& A) const {
        return (operator[](X) != A[X]) || (operator[](Y) != A[Y]) || (operator[](Z) != A[Z]);
    }
    ;
    
    Coord operator*(double M, const Coord& c) {
        return Coord(c[X] * M, c[Y] * M, c[Z] * M);
    }
    
    double dist2(const Coord& a, const Coord& b) {
        return ((a - b).norm2());
    }
    
    double dist(const Coord& a, const Coord& b) {
        return ((a - b).norm());
    }
    
    double Coord::dot(const Coord& t, Coord& gradient) const {
        
        gradient = *this;
        
        return dot(t);
        
    }
    
    double Coord::angle(const Coord& t) const {
        
        return MR::Math::acos(dot(t) / (norm() * t.norm()));
        
    }
    
    double Coord::norm2(Coord& gradient) const {
        
        gradient = *this * 2.0;
        
        return norm2();
        
    }
    
    double Coord::norm(Coord& gradient) const {
        
        double nrm = norm();
        
        gradient = *this / nrm;
        
        return nrm;
        
    }
    
    double Coord::norm2(Coord& gradient, Coord::Tensor& hessian) const {
        
        gradient = *this * 2.0;
        
        hessian.diagonal() = 2.0;
        
        return norm2();
        
    }
    
    double Coord::norm(Coord& gradient, Coord::Tensor& hessian) const {
        
        double nrm = norm();
        
        gradient = *this / nrm;
        
        hessian = this->outer(*this);
        hessian *= -MR::Math::pow3(1 / nrm);
        hessian.diagonal() += 1 / nrm;
        
        return nrm;
        
    }
    
    Coord::Tensor Coord::outer(const Coord& c) const {
        
        Coord::Tensor tensor;
        
        tensor.row(X) = c * operator[](X);
        tensor.row(Y) = c * operator[](Y);
        tensor.row(Z) = c * operator[](Z);
        
        return tensor;
        
    }

}
