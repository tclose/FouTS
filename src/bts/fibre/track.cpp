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


#include "bts/fibre/track.h"
#include "bts/fibre/strand.h"
#include "math/matrix.h"
//#include "math/matrix_ops.h"
#include "math/least_squares.h"



namespace BTS {

  namespace Fibre {

    const std::string                 Track::RADIUS_PROP = "track_radius";
    const std::string                 Track::BUNDLE_INDEX_EPROP = "bundle_index";
    const std::string                 Track::RADIUS_PROP_DEFAULT = "0.02";
    
    const Coord                       Track::FILE_SEPARATOR = Triple<double> (NAN, NAN, NAN);
    const std::string                 Track::FILE_EXTENSION = "tck";

    const char*                       Track::PROPS_LIST[] = { Object::ACS_SQRT_PROP, PROPS_LIST_END };

//    const char*                       Track::INTRINSIC_PROPERTIES[Track::NUM_INTRINSIC_PROPERTIES] = {};



    Track::Track(const Strand& s, size_t num_points, bool include_endpoints)
      :  Base::Object (num_points, num_points * 3 + s.num_props(), select_props<Track>(*s.props)) {

      if (num_points == 0)
        throw Exception ("'num_length_sections' option is required for conversion to tracks.");

      const MR::Math::Matrix<double>& position_matrix  = Strand::position_matrix(num_points, s.degree(), include_endpoints);

      for (size_t point_i = 0; point_i < num_points; point_i++)
        operator[](point_i) = s.left_product(position_matrix.row(point_i));

      for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
        prop(prop_i) = s.prop(prop_key(prop_i));

    }


    Strand      Track::to_strand(size_t degree) const {

      return Strand(*this,degree);

    }


    Strand      Track::to_strand(size_t degree, const std::vector<double>& timeline) const {

      return Strand(*this,degree,timeline);

    }


    Coord       Track::left_product(const MR::Math::Vector<double>::View& row_vector) const {

      assert(row_vector.size() >= num_points());

      Coord product(0,0,0);

      for (size_t point_i = 0; point_i < num_points(); point_i++)
        product += operator[](point_i) * row_vector[point_i];

      return product;

    }


    std::ostream& operator<< (std::ostream& stream, const Track& tck) {
      stream << std::endl;
      for (size_t point_i = 0; point_i < tck.num_points(); point_i++) {
        stream << tck[point_i] << std::endl;
      }

      return (stream);
    }


  }

}
