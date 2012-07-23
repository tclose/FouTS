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


#ifndef __bts_fibre_segment_h__
#define __bts_fibre_segment_h__

#include <assert.h>
#include <vector>
#include <string>

#include "bts/triple.h"

#include "bts/fibre/base/object.h"

namespace BTS {

  namespace Fibre {

    class Segment : public Base::Object {

      public:

        class Tensor;
        friend class Tensor;

        typedef Coord Element;
        typedef Fibre::Base::Reader<Segment> Reader;
        typedef Fibre::Base::Writer<Segment> Writer;

      protected:

        const static char*                    INTENSITY_PROP;
        const static std::vector<const char*> PROPS;

      public:
        const static char*          PROPS_LIST[];

        const static std::string    FILE_EXTENSION;
        const static Coord          FILE_SEPARATOR;

      //Public member functions.
      public:

        Segment()
          : Base::Object((size_t)2, (size_t)6 + PROPS.size(), PROPS) { Base::Object::set(NAN); }

        Segment(const MR::Math::Vector<double>::View& view, std::vector<const char*>* props)
          : Base::Object((size_t)2, view, props) { assert(view.size()==7); }

      // Used in derived classes
      protected:

        Segment(size_t size, const std::vector<const char*>& props)
          : Base::Object(size, size * 3 + props.size(), props) { set(NAN); }

        Segment(size_t size, const MR::Math::Vector<double>::View& view, std::vector<const char*>* props)
          : Base::Object(size, view, props) {}

      public:

        virtual ~Segment() {}

        Coord           position()
          { return Coord(MR::Math::Vector<double>::sub(0,3)); }

        Coord           position() const
          { return Coord(MR::Math::Vector<double>::sub(0,3)); }

        Coord           tangent()
          { return Coord(MR::Math::Vector<double>::sub(3,6)); }

        Coord           tangent() const
          { return Coord(MR::Math::Vector<double>::sub(3,6)); }

        double&         position(size_t idx)
          { return MR::Math::Vector<double>::operator[](idx); }

        const double&   position(size_t idx) const
          { return MR::Math::Vector<double>::operator[](idx); }

        double&         tangent(size_t idx)
          { return MR::Math::Vector<double>::operator[](idx + 3); }

        const double&   tangent(size_t idx) const
          { return MR::Math::Vector<double>::operator[](idx + 3); }


        double&         intensity()
         { return prop(INTENSITY_PROP);}

        const double&   intensity() const
         { return prop(INTENSITY_PROP); }


        BASE_GENERAL_FUNCTIONS(Segment);

        BASE_MULT_DIVIDE_FUNCTIONS(Segment);

        BASE_ADD_SUBTRACT_FUNCTIONS(Segment);

        double                    length() const { return tangent().norm(); }

      friend std::ostream& operator<< (std::ostream& stream , const Segment& seg);

    };


    inline Segment     operator* (double scalar, const Segment& section)
      { return section * scalar; }


    std::ostream&                   operator<< (std::ostream& stream , const Segment& seg);


  }

}

#endif
