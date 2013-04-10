/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

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

#include "bts/triple.h"
#include "bts/coord.h"

namespace BTS {
    
    template<> const Triple<double> Triple<double>::Invalid(NAN, NAN, NAN);
    template<> const Triple<double> Triple<double>::Zeros(0.0, 0.0, 0.0);
    template<> const Triple<double> Triple<double>::Ones(1.0, 1.0, 1.0);
    template<> const Triple<double> Triple<double>::Halves(0.5, 0.5, 0.5);
    
    template<> const Triple<int> Triple<int>::Invalid(-999999999, -999999999, -999999999);
    template<> const Triple<int> Triple<int>::Zeros(0, 0, 0);
    template<> const Triple<int> Triple<int>::Ones(1, 1, 1);
    
    template<> const Triple<size_t> Triple<size_t>::Invalid(999999999, 999999999, 999999999);
    template<> const Triple<size_t> Triple<size_t>::Zeros(0, 0, 0);
    template<> const Triple<size_t> Triple<size_t>::Ones(1, 1, 1);
    
    template<> Triple<double>::Triple(const Coord& coord) {
        set(coord[0], coord[1], coord[2]);
    }
    
    template<> Triple<float>::Triple(const Coord& coord) {
        set(coord[0], coord[1], coord[2]);
    }

}
