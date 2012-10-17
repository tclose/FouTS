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


#ifndef __bts_image_buffer_h__
#define __bts_image_buffer_h__

namespace BTS {

  namespace Image {

    template <typename T> class Buffer_tpl;

  }

}

#include <set>

#include "bts/common.h"

#include "bts/image/index.h"

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

    template <typename T> std::ostream&    operator<<(std::ostream& stream, const Buffer_tpl<T>& B);

    template <typename T> class Buffer_tpl {

      public:

        typedef typename std::map<Index, T>::iterator iterator;
        typedef typename std::map<Index, T>::const_iterator const_iterator;

        
      //Protected member variables
      protected:

        Triple<size_t>        dimensions;
        std::map<Index, T>  voxels;
        bool                enforce_bounds;
        
      //Public member functions
      public:

        Buffer_tpl (bool enforce_bounds = true)
          : dimensions(0,0,0), enforce_bounds(enforce_bounds) {}


        Buffer_tpl (const Triple<size_t>& dimensions, bool enforce_bounds = true)
         : dimensions(dimensions), enforce_bounds(enforce_bounds) {}


        Buffer_tpl<T> (const Buffer_tpl<T>& B)
          : dimensions(B.dims()), voxels(B.voxels), enforce_bounds(B.enforce_bounds) {}


        ~Buffer_tpl() {}


        Buffer_tpl<T>&                        operator= (const Buffer_tpl<T>& B)
          { dimensions = B.dimensions; voxels = B.voxels; enforce_bounds = B.enforce_bounds; return *this; }

        
        inline T&                             operator() (int x, int y, int z)
          { return operator()(Index(x,y,z)); }
          
          
        inline const T&                       operator() (int x, int y, int z) const
          { return operator()(Index(x,y,z)); }


        inline T&                             operator() (const Index& c);


        inline const T&                       operator() (const Index& c) const;


        bool                                  is_empty(size_t x, size_t y, size_t z) const
          { return is_empty(Index(x,y,z)); }


        bool                                  is_empty(const Index& c) const;


        std::set<Index>                       non_empty() const;


        std::set<Index>                       non_empty_or_inbounds() const;


        std::set<Index>                       empty_inbounds() const;


        iterator                              begin()
          { return voxels.begin(); }

        const_iterator                        begin() const
          { return voxels.begin(); }

        iterator                              end()
          { return voxels.end(); }

        const_iterator                        end() const
          { return voxels.end(); }


        virtual Buffer_tpl&                   clear()
          { voxels.clear(); return *this; }


        //FIXME: Would be better if it was virtual but that throws a compile error for some reason.
        Buffer_tpl&                           zero()
          { for (iterator vox_it = begin(); vox_it != end(); ++vox_it) BTS::zero(vox_it->second); return *this; }


        Buffer_tpl&                           invalidate()
          { for (iterator vox_it = begin(); vox_it != end(); ++vox_it) BTS::invalidate(vox_it->second); return *this; }

        /*! Resizes the dimensions of the image.  Note that is does not change the underlying voxels which are sparsely
         *  implemented in a std::map.  It does set enforce_bounds to false as it does not guarantee that the existing voxels
         *  are contained withinness the new bounds. */
        void                                  resize(const Triple<size_t>& dimensions)
          { this->dimensions = dimensions; relax_bounds(); }


        void                                  reset(const Triple<size_t>& dimensions)
          { reset(dimensions, enforce_bounds); }

        void                                  reset(const Triple<size_t>& dimensions, bool enforce_bounds)
          { this->dimensions = dimensions; this->enforce_bounds = enforce_bounds; voxels.clear(); }


        //!Relaxes the constraint that voxels need to be withinness the given bounds of the image.
        void                                  relax_bounds()
          { enforce_bounds = false; }

        void                                  clear_and_enforce_bounds();


        bool                                  bounds_are_enforced() const
          { return enforce_bounds; }


        Buffer_tpl<T>&                        operator+= (const Buffer_tpl<T>& buff);


        Buffer_tpl<T>&                        operator-= (const Buffer_tpl<T>& buff);


        template <typename U> Buffer_tpl<T>&  operator+= (const Buffer_tpl<U>& buff);


        template <typename U> Buffer_tpl<T>&  operator-= (const Buffer_tpl<U>& buff);


        template <typename U> Buffer_tpl<T>&  operator*= (const Buffer_tpl<U>& buff);


        template <typename U> Buffer_tpl<T>&  operator/= (const Buffer_tpl<U>& buff);


        template <typename U> Buffer_tpl<T>&  operator*= (const U& M);


        template <typename U> Buffer_tpl<T>&  operator/= (const U& M);


        Buffer_tpl<T>&                        negate();


        Buffer_tpl<T>                         operator+ (const Buffer_tpl<T>& buff) const;


        Buffer_tpl<T>                         operator- (const Buffer_tpl<T>& buff) const;


        Buffer_tpl<T>                         operator* (const Buffer_tpl<T>& buff) const;


        Buffer_tpl<T>                         operator/ (const Buffer_tpl<T>& buff) const;


        Buffer_tpl<T>                         operator* (double M) const
          { Buffer_tpl<T> mult = *this; mult *= M; return mult; }       


        Buffer_tpl<T>                         operator/ (double M) const
          { Buffer_tpl<T> mult = *this; mult /= M; return mult; }    

        
        const Triple<size_t>&                 dims () const
          { return dimensions; }
        

        size_t                                dim (size_t dim_index) const
          { return dimensions[dim_index]; }


        size_t                                num_voxels_in_bounds() const
          { return dim(X) * dim(Y) * dim(Z); }


        size_t                                num_not_empty_voxels() const
          { return voxels.size(); }


        bool                                  in_bounds(Index coord) const
          { return coord.non_negative() && coord.bounded_by(dims()); }

      protected:

        virtual T                             new_voxel(const Index& c)
          { return T(); }


      friend std::ostream& operator<<<>(std::ostream& stream, const Buffer_tpl<T>& B);

    };


    template <typename T> Buffer_tpl<T>       operator* (double M, Buffer_tpl<T> mult)
      { mult *= M; return mult; }

  }

}

#undef LOOP

#include "bts/image/voxel.h"

namespace BTS {

  namespace Image {

    class Buffer : public Buffer_tpl< Voxel<double> > {

      protected:

        size_t num_encodings;

      public:

        Buffer (size_t num_encodings = 0, bool enforce_bounds = true)
          : Buffer_tpl< Voxel<double> >(enforce_bounds), num_encodings(num_encodings) {}


        Buffer (const Triple<size_t>& dimensions, size_t num_encodings, bool enforce_bounds = true)
          : Buffer_tpl< Voxel<double> >(dimensions, enforce_bounds), num_encodings(num_encodings) {}


        Buffer (const Buffer& B)
          : Buffer_tpl< Voxel<double> >(B), num_encodings(B.num_encodings) {}

        template <typename T> Buffer (const Buffer_tpl<T>& B)
          : Buffer_tpl< Voxel<double> >(B.dims(), B.bounds_are_enforced()), num_encodings(0) {

          int read_num_encodings = -1;

          for (typename Buffer_tpl<T>::const_iterator vox_it = B.begin(); vox_it != B.end(); ++vox_it) {

            this->operator()(vox_it->first) = vox_it->second;

            if (read_num_encodings == -1)
              read_num_encodings = vox_it->second.num_encodings();
            else if (read_num_encodings != (int)vox_it->second.num_encodings())
              throw Exception ("Number of encodings in copied voxels do not match (" + str(read_num_encodings) + " and " + str(vox_it->second.num_encodings()) + ").");

          }

          num_encodings = read_num_encodings;

        }


        virtual ~Buffer() {}


        Buffer&                            operator= (const Buffer& B)
          { this->Buffer_tpl< Voxel<double> >::operator=(B); return *this; }


      protected:

        Voxel<double>                      new_voxel(const Index& coord)
          { return Voxel<double>(num_encodings); }


    };


    namespace Float {

      typedef Image::Buffer_tpl<float> Buffer;

    };

    namespace Double {

      typedef Image::Buffer_tpl<double> Buffer;

    };

    namespace Int {

      typedef Image::Buffer_tpl<int> Buffer;

    };

    namespace UInt {

      typedef Image::Buffer_tpl<size_t> Buffer;

    };

    namespace Bool {

      typedef Image::Buffer_tpl<bool> Buffer;

    };

  }
}

#undef LOOP


#endif
