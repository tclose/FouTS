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

#ifndef __bts_fibre_track_h__
#define __bts_fibre_track_h__


#include "bts/common.h"
#include "bts/fibre/base/object.h"
#include "bts/coord.h"

#include "bts/fibre/base/writer.h"
#include "bts/fibre/base/reader.h"

namespace BTS {

  namespace Fibre {

    class Track : public Base::Object {

      public:

        typedef Coord Element;

        class Set;

        typedef Base::Reader<Track>              Reader;
        typedef Fibre::Base::Writer<Track>       Writer;

        typedef Base::TextReader<Track>          TextReader;
        typedef Base::TextWriter<Track>          TextWriter;

        friend class Object;
        friend class Strand;
        friend class Tractlet;
        friend class Set;
        friend class Base::Set<Track>;
        friend class Base::Reader<Track>;
        friend class Base::Writer<Track>;

      public:

        const static Coord                      FILE_SEPARATOR;
        const static std::string                FILE_EXTENSION;
        const static std::string                RADIUS_PROP;
        const static std::string                BUNDLE_INDEX_EPROP;
        const static std::string                RADIUS_PROP_DEFAULT;

      protected:

        const static char*          PROPS_LIST[];

      public:

        const static size_t NUM_LENGTH_SECTIONS_DEFAULT = 100;

      public:

        Track (size_t num_points = 0, const std::vector<const char*> props = std::vector<const char*>())
          : Base::Object (num_points, num_points * 3 + props.size(), props) { set(NAN); }

        Track (const Track& t) 
          : Base::Object(t) {}

        ~Track() {}

        Track&        operator=(const Track& t)
          { Base::Object::operator=(t); return *this; }


        Track (const Strand& s, size_t num_points, bool include_endpoints = 0);

      protected:

        /*! Used by Base::Set<Track>::operator[]() to initiate a from a section of a larger state vector or
         * by Base::Tensor<Track>::row/column() for a slice of the tensor matrix
         *
         * @param num_points The number of points in the track
         * @param view The view onto the larger vector or matrix
         * @param props The properties stored in the set or tensor
         */
        Track(size_t num_points, const MR::Math::Vector<double>::View& view, std::vector<const char*>* props)
          : Base::Object(num_points, view, props) {}

      public:

        BASE_GENERAL_FUNCTIONS(Track);

        BASE_ADD_SUBTRACT_FUNCTIONS(Track);

        BASE_MULT_DIVIDE_FUNCTIONS(Track);

        /*! Used in Base::Set<T> functions that need to differentiate between a tractlets degree and its size, which for
         *  Strand's and Track's they are equivlaent.
         *
         * @return The size of the Track
         */
        size_t                            degree() const
          { return size(); };

        double                              acs() const
          { return has_prop(ACS_SQRT_PROP) ? prop(ACS_SQRT_PROP) : 1.0; }

        void                                add_acs(double acs = 1.0)
          { add_prop(ACS_SQRT_PROP, acs); }

        void                                set_acs(double acs)
          { assert(acs >= 0); prop(ACS_SQRT_PROP) = MR::Math::sqrt(acs); }

        void                                remove_acs()
          { remove_prop(ACS_SQRT_PROP); }

        Track                               base()
          { return Track(sze, sub(0,bsize()), &EMPTY_PROPS); }

        const Track                         base() const
          { return Track(sze, sub(0,bsize()), &EMPTY_PROPS); }

        Track&                              operator+=  (const Coord& p)
          { for (size_t point_i = 0; point_i < size(); point_i++) operator[](point_i) += p; return *this; }

        Track&                              operator-=  (const Coord& p)
          { for (size_t point_i = 0; point_i < size(); point_i++) operator[](point_i) -= p; return *this; }

        Track                               operator+  (const Coord& p)  const
          { Track answer(*this); answer += p; return answer; }

        Track                               operator-  (const Coord& p) const
          { Track answer(*this); answer -= p; return answer; }

        size_t                              num_points() const
          { return size(); }

        Strand                              to_strand(size_t degree) const;

        Strand                              to_strand(size_t degree, const std::vector<double>& timeline) const;

        Coord                               left_product(const MR::Math::Vector<double>::View& row_vector) const;


    };


   std::ostream& operator<< (std::ostream& stream, const Track& tck);

  }

}

#undef LOOP

#endif



