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

#ifndef __bts_image_index_h__
#define __bts_image_index_h__

namespace BTS {

  namespace Image {

    class Index;

  }

}

#include "bts/common.h"

#include "bts/triple.h"



namespace BTS {

  namespace Image {

    class Index : public Triple<int> {

      public:
      
        Index () {}
        
        Index (int x_value, int y_value, int z_value)   
          : Triple<int>(x_value, y_value, z_value) {}
        
        Index(const Index& c) 
          : Triple<int>(c) {}
          
        Index(const Triple<int>& t)
          : Triple<int>(t) {}

        Index(const Triple<size_t>& t)
          : Triple<int>(t[X], t[Y], t[Z]) {}
      
        inline bool   positive()                       
          { return (this->p[X] > 0.0) && (this->p[Y] > 0.0) && (this->p[Z] > 0.0); }
          
        inline bool   non_negative()                  
          { return (this->p[X] >= 0.0) && (this->p[Y] >= 0.0) && (this->p[Z] >= 0.0); }
          
        inline bool   bounded_by(const Triple<int>& A) const     
          { return (this->p[X] < A[X]) && (this->p[Y] < A[Y]) && (this->p[Z] < A[Z]); }  

        inline bool   bounded_by(const Triple<size_t>& A) const     
          { return (this->p[X] < (int)A[X]) && (this->p[Y] < (int)A[Y]) && (this->p[Z] < (int)A[Z]); }

        void          bound(int image_dims);
        template <typename T> void bound(const Triple<T> &image_dims);

        void          bound(int low_bound, int up_bound);
        template <typename T> void bound(const Triple<T> &low_bound, const Triple<T> &up_bound);
                
        /*! Comparison between coordinates is implemented so that they can be sorted in std::maps and works by
         * placing the most significance on X, then Y, and the Z channels.
         */
        bool  operator< (const Index& A) const;
        
        bool  operator> (const Index& A) const;

        bool  operator<= (const Index& A) const
          { return this->operator<(A) ? true : this->operator == (A); }
        
        bool  operator>= (const Index& A) const
          { return this->operator>(A) ? true : this->operator == (A); }

                
        Triple<double>         to_point()
          { return Triple<double>((double)this->p[X] , (double)this->p[Y], (double)this->p[Z]); }        

    };
    

    std::ostream&        operator<< (std::ostream& stream, const Index& coord);

    template <typename T> void Index::bound(const Triple<T> &image_dims) {

      this->p[X] = min2(max2(this->p[X],0),image_dims[X]);
      this->p[Y] = min2(max2(this->p[Y],0),image_dims[Y]);
      this->p[Z] = min2(max2(this->p[Z],0),image_dims[Z]);

    }

    inline void Index::bound(int low_bound, int up_bound) {
      this->p[X] = min2(max2(this->p[X],low_bound),up_bound);
      this->p[Y] = min2(max2(this->p[Y],low_bound),up_bound);
      this->p[Z] = min2(max2(this->p[Z],low_bound),up_bound);
    }

    template <typename T> void Index::bound(const Triple<T> &low_bound, const Triple<T> &up_bound) {
      this->p[X] = min2(max2(this->p[X],low_bound[X]),up_bound[X]);
      this->p[Y] = min2(max2(this->p[Y],low_bound[Y]),up_bound[Y]);
      this->p[Z] = min2(max2(this->p[Z],low_bound[Z]),up_bound[Z]);
    }

  }

}








//#undef COPY_VALUES

#endif


