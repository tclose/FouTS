/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 14/07/2010.

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

#include "bts/image/index.h"

namespace FTS {
    
    namespace Image {
        
        std::ostream& operator<<(std::ostream& stream, const Index& coord) {
            
            stream << "[" << coord[X] << ", " << coord[Y] << ", " << coord[Z] << "]";
            
            return stream;
            
        }
        
        void Index::bound(int image_dims) {
            
            this->p[X] = min2(max2(this->p[X], 0), image_dims);
            this->p[Y] = min2(max2(this->p[Y], 0), image_dims);
            this->p[Z] = min2(max2(this->p[Z], 0), image_dims);
            
        }
        
        bool Index::operator<(const Index& A) const {
            
            bool is_less_than;
            if (p[X] < A[X])
                is_less_than = true;
            else if (p[X] == A[X]) {
                if (p[Y] < A[Y])
                    is_less_than = true;
                else if (p[Y] == A[Y])
                    is_less_than = (p[Z] < A[Z]);
                else
                    is_less_than = false;
            } else
                is_less_than = false;
            return is_less_than;
        }
        
        bool Index::operator>(const Index& A) const {
            
            bool is_greater_than;
            if (p[X] > A[X])
                is_greater_than = true;
            else if (p[X] == A[X]) {
                if (p[Y] > A[Y])
                    is_greater_than = true;
                else if (p[Y] == A[Y])
                    is_greater_than = (p[Z] > A[Z]);
                else
                    is_greater_than = false;
            } else
                is_greater_than = false;
            return is_greater_than;
        }
    
    }
}
