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

#include <set>

#include "bts/common.h"
#include "bts/math/svd.h"

#include "bts/fibre/tractlet.h"
#include "bts/fibre/track/set.h"
#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/strand/set.h"



#include "bts/fibre/base/reader.cpp.h"
#include "bts/fibre/base/writer.cpp.h"

namespace BTS {

  namespace Fibre {

    const std::string                         Tractlet::Set::FILE_EXTENSION = "tst";
    const char*                               Tractlet::Set::LENGTH_EPSILON_PROP = "length_epsilon";
    const char*                               Tractlet::Set::WIDTH_EPSILON_PROP = "width_epsilon";

    const double                              Tractlet::Set::LENGTH_EPSILON_DEFAULT = 0.01;
    const double                              Tractlet::Set::WIDTH_EPSILON_DEFAULT = 0.01;
    const char*                               Tractlet::Set::PROPS_LIST[] = { Base::Set<Tractlet>::BASE_INTENSITY_PROP,
                                                                                Tractlet::Set::LENGTH_EPSILON_PROP,
                                                                                Tractlet::Set::WIDTH_EPSILON_PROP,
                                                                                PROPS_LIST_END };



    Tractlet::Set::Set(const Strand::Set& strands)
      : Fibre::Base::Set<Tractlet>( select_props<Set>(*strands.props),
                                    select_props<Tractlet>(*strands.elem_props),
                                    strands.get_extend_props()) {

      from_strands(strands);

    }


    //TODO: Look into preallocating for speed
    Tractlet::Set::Set(const Strand::Set& strands, double width)
      : Fibre::Base::Set<Tractlet>( select_props<Set>(*strands.props),
                                    select_props<Tractlet>(*strands.elem_props),
                                    strands.get_extend_props()) {

      from_strands(strands,width);

    }


    //TODO: Look into preallocating for speed
    Tractlet::Set::Set(const Strand::Set& strands, size_t num_tractlets, const std::string& cluster_save_location)
      : Fibre::Base::Set<Tractlet>( select_props<Set>(*strands.props),
                                    select_props<Tractlet>(*strands.elem_props),
                                    strands.get_extend_props()) {

      from_strands(strands,num_tractlets,cluster_save_location);

    }


    Tractlet::Set::Set (const Track::Set& strands, size_t degree)
      : Fibre::Base::Set<Tractlet>(select_props<Set>(*strands.props),
                                  select_props<Tractlet>(*strands.elem_props),
                                  strands.get_extend_props()) {

      from_strands(Strand::Set(strands, degree));

    }


    //TODO: Look into preallocating for speed
    Tractlet::Set::Set (const Track::Set& tcks, size_t degree, size_t num_tractlets, const std::string& cluster_save_location)
      : Fibre::Base::Set<Tractlet>(select_props<Set>(*tcks.props),
                                  select_props<Tractlet>(*tcks.elem_props),
                                  tcks.get_extend_props()) {

      from_strands(Strand::Set(tcks, degree), num_tractlets, cluster_save_location);

    }


    //TODO: Look into preallocating for speed
    Tractlet::Set::Set (const Track::Set& tcks, size_t degree, double width)
      : Fibre::Base::Set<Tractlet>( select_props<Set>(*tcks.props),
                                    select_props<Tractlet>(*tcks.elem_props),
                                    tcks.get_extend_props()) {

      from_strands(Strand::Set(tcks),width);


    }



    void        	                              Tractlet::Set::from_strands(const Strand::Set& strands) {

      if (!strands.has_extend_elem_prop(BUNDLE_INDEX_EPROP))
        throw Exception ("'" + BUNDLE_INDEX_EPROP + "' property was not found in strand set, which is required for conversion to tractlet set.");

      add_extend_elem_prop(BUNDLE_INDEX_EPROP,"-1");

      std::set<size_t> bundle_indices;

      //Get a list of bundle indices present in the current set.
      for (size_t strand_i = 0; strand_i < strands.size(); ++strand_i)
        bundle_indices.insert(to<double>(strands.get_extend_elem_prop(BUNDLE_INDEX_EPROP, strand_i)));

      //Loop through the strand set and convert each bundle into a separate tractlet.
      for (std::set<size_t>::iterator bundle_it = bundle_indices.begin(); bundle_it != bundle_indices.end(); ++bundle_it) {

        Strand::Set bundle_strands(std::vector<const char*>(), *strands.elem_props);

        //TODO: Maybe look into preallocating the size of the strand set to avoid repeated 'push_back' calls.
        for (size_t strand_i = 0; strand_i < strands.size(); ++strand_i)
          if (to<size_t>(strands.get_extend_elem_prop(BUNDLE_INDEX_EPROP, strand_i)) == *bundle_it)
            bundle_strands.push_back(strands[strand_i]);

        size_t num_strands = bundle_strands.size();

        bundle_strands = bundle_strands.shortest_span();

        Strand centroid;
        centroid = bundle_strands[0];
        centroid.zero();

        for (size_t strand_i = 0; strand_i < bundle_strands.size(); ++strand_i)
          centroid += bundle_strands[strand_i];

        centroid /= (double)num_strands;

        Tractlet tractlet(centroid.degree(), select_props<Tractlet>(centroid.prop_keys()));

        for (size_t prop_i = 0; prop_i < tractlet.num_props(); ++prop_i)
          tractlet.prop(prop_i) = centroid.prop(tractlet.prop_key(prop_i));

        tractlet[0] = centroid.base();


        if (num_strands == 1) {
          //Get width axes perpendicular to the main orientation and each other.
          Coord ax1_basewidth = centroid[1].cross(centroid[1].min_axis()).normalise();
          Coord ax2_basewidth = centroid[1].cross(ax1_basewidth).normalise();

          tractlet[1][0] = ax1_basewidth * Tractlet::REASONABLE_WIDTH;
          tractlet[2][0] = ax2_basewidth * Tractlet::REASONABLE_WIDTH;


        } else if (num_strands == 2) {

          Strand ax1 = bundle_strands[0] - centroid;

          Coord ax2_basewidth = centroid[1].cross(ax1[1]).normalise() * Tractlet::REASONABLE_WIDTH;

          tractlet[1] = ax1;
          tractlet[2][0] = ax2_basewidth;

        } else {

          //Set up the difference matrix which holds the difference between each strand and the centroid.

          size_t num_rows = num_strands;
          if (centroid.bsize() > num_strands)
            num_rows = centroid.bsize(); //If vector size of strands is less than number in bundle the difference matrix is padded out with zeros;

          MR::Math::Matrix<double> diff_matrix (num_rows, centroid.bsize());
          diff_matrix = 0;

          for (size_t strand_i = 0; strand_i < num_strands; ++strand_i) {

            diff_matrix.row(strand_i) = bundle_strands[strand_i].base();
            diff_matrix.row(strand_i) -= centroid.base();

          }

          //Perform singular value decomposition to get the two axes that describe the most variance.
          Math::USV usv = Math::svd(diff_matrix);

          tractlet[1] = usv.v.column(0);
          tractlet[2] = usv.v.column(1);

          //Calculate the RMS of projections onto the two axes.
          double rms_proj1 = 0.0;
          double rms_proj2 = 0.0;

          for (size_t strand_i = 0; strand_i < num_strands; ++strand_i) {

            rms_proj1 += MR::Math::pow2(dot(diff_matrix.row(strand_i), tractlet[1]));
            rms_proj2 += MR::Math::pow2(dot(diff_matrix.row(strand_i), tractlet[2]));

          }

          rms_proj1 /= (double)num_strands;
          rms_proj2 /= (double)num_strands;

          rms_proj1 = MR::Math::sqrt(rms_proj1);
          rms_proj2 = MR::Math::sqrt(rms_proj2);

          // Multiply the RMS projections by two to get the extent of the tractlet along each axes
          // (Assumes that the strand bundle is rougly circular w.r.t to the two axes)
          tractlet[1] *= rms_proj1 * 2.0;
          tractlet[2] *= rms_proj2 * 2.0;

        }

        tractlet.sanitize();

        //Keep the bundle index property (so it keeps the same colour in the plotting routines).
        std::map<std::string,std::string> properties_row;
        properties_row[BUNDLE_INDEX_EPROP] = str(*bundle_it);

        push_back(tractlet, properties_row);

      }

      for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
        prop(prop_i) = strands.prop(prop_key(prop_i));


      if (!strands.var_elem_degrees())
        freeze_elem_degree();


      //TODO: Look into automatically trying to fix degree of tractlet set after conversion.

    }
    

    void   	                                    Tractlet::Set::from_strands(const Strand::Set& set, size_t num_tractlets,
                                                                            const std::string& cluster_save_location) {

      if (num_tractlets) {
        Strand::Set set_copy(set);
        set_copy.cluster(num_tractlets);

        if (cluster_save_location.size())
          set_copy.save(cluster_save_location);

        from_strands(set_copy);
      } else
        from_strands(set);

    }


    void 	                                      Tractlet::Set::from_strands (const Strand::Set& set, double width) {

      for (size_t strand_i = 0; strand_i < set.size(); ++strand_i)
        push_back(Tractlet(set[strand_i], width));

      for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
        prop(prop_i) = set.prop(prop_key(prop_i));

      copy_extend_elem_props(set);

      if (has_extend_elem_prop(Fibre::Track::RADIUS_PROP))
        remove_extend_elem_prop(Fibre::Track::RADIUS_PROP);

      if (!set.var_elem_degrees())
        freeze_elem_degree();


    }


//
//    //FIXME: Will not work yet
//    Tractlet::Set::Set(const Strand::Set& set) {
//
//
//      if (!set.has_extend_elem_prop(BUNDLE_INDEX_EPROP))
//        throw Exception ("'" + BUNDLE_INDEX_EPROP + "' property was not found in strand set, which is required for conversion to tractlet set.");
//
//      add_extend_elem_prop(BUNDLE_INDEX_EPROP,"-1");
//
//      std::set<size_t> bundle_indices;
//
//      //Get a list of bundle indices present in the current set.
//      for (size_t strand_i = 0; strand_i < set.size(); ++strand_i)
//        bundle_indices.insert(to<double>(set.get_extend_elem_prop(BUNDLE_INDEX_EPROP, strand_i)));
//
//      //Loop through the strand set and convert each bundle into a separate tractlet.
//      for (std::set<size_t>::iterator bundle_it = bundle_indices.begin(); bundle_it != bundle_indices.end(); ++bundle_it) {
//
//        Fibre::Strand::Set bundle_strands;
//
//        for (size_t strand_i = 0; strand_i < set.size(); ++strand_i)
//          if (to<size_t>(set.get_extend_elem_prop(BUNDLE_INDEX_EPROP, strand_i)) == *bundle_it)
//            bundle_strands.push_back(operator[](strand_i));
//
//        size_t num_strands = bundle_strands.size();
//
//        bundle_strands = bundle_strands.shortest_span();
//
//        Fibre::Strand centroid;
//        centroid.zero();
//
//        for (size_t strand_i = 0; strand_i < bundle_strands.size(); ++strand_i)
//          centroid += bundle_strands[strand_i];
//
//        centroid /= (double)num_strands;
//
//        Fibre::Tractlet tractlet(centroid.degree());
//        tractlet.zero();
//        tractlet.set_acs(centroid.acs());
//        tractlet[0] = centroid;
//
//
//        if (num_strands == 1) {
//          //Get width axes perpendicular to the main orientation and each other.
//          Coord ax1_basewidth = centroid[1].cross(centroid[1].min_axis()).normalise();
//          Coord ax2_basewidth = centroid[1].cross(ax1_basewidth).normalise();
//
//          tractlet[1][0] = ax1_basewidth * Tractlet::REASONABLE_WIDTH;
//          tractlet[2][0] = ax2_basewidth * Tractlet::REASONABLE_WIDTH;
//
//
//        } else if (num_strands == 2) {
//
//          Strand ax1 = bundle_strands[0] - centroid;
//
//          Coord ax2_basewidth = centroid[1].cross(ax1[1]).normalise() * Tractlet::REASONABLE_WIDTH;
//
//          tractlet[1] = ax1;
//          tractlet[2][0] = ax2_basewidth;
//
//        } else {
//
//          //Set up the difference matrix which holds the difference between each strand and the centroid.
//
//          size_t vector_size = centroid.vsize();
//
//          size_t num_rows = num_strands;
//          if (vector_size > num_strands)
//            num_rows = vector_size; //If vector size of strands is less than number in bundle the difference matrix is padded out with zeros;
//
//          MR::Math::Matrix<double> diff_matrix (num_rows, vector_size);
//          diff_matrix = 0;
//
//          for (size_t strand_i = 0; strand_i < num_strands; ++strand_i) {
//
//            diff_matrix.row(strand_i) = bundle_strands[strand_i];
//            diff_matrix.row(strand_i) -= centroid;
//
//          }
//
//          //Perform singular value decomposition to get the two axes that describe the most variance.
//          Math::USV usv = Math::svd(diff_matrix);
//
//          MR::Math::Vector<double> axis1 = usv.v.column(0);
//          MR::Math::Vector<double> axis2 = usv.v.column(1);
//
//
//          //Calculate the RMS of projections onto the two axes.
//          double rms_proj1 = 0.0;
//          double rms_proj2 = 0.0;
//
//          for (size_t strand_i = 0; strand_i < num_strands; ++strand_i) {
//
//            rms_proj1 += MR::Math::pow2(dot(diff_matrix.row(strand_i), axis1));
//            rms_proj2 += MR::Math::pow2(dot(diff_matrix.row(strand_i), axis2));
//
//          }
//
//          rms_proj1 /= (double)num_strands;
//          rms_proj2 /= (double)num_strands;
//
//          rms_proj1 = MR::Math::sqrt(rms_proj1);
//          rms_proj2 = MR::Math::sqrt(rms_proj2);
//
//          // Multiply the RMS projections by two to get the extent of the tractlet along each axes
//          // (Assumes that the strand bundle is rougly circular w.r.t to the two axes)
//          tractlet[1] *= rms_proj1 * 2.0;
//          tractlet[2] *= rms_proj2 * 2.0;
//
//        }
//
//        tractlet.sanitize();
//
//        //Keep the bundle index property (so it keeps the same colour in the plotting routines).
//        std::map<std::string,std::string> properties_row;
//        properties_row[BUNDLE_INDEX_EPROP] = str(*bundle_it);
//
//        push_back(tractlet, properties_row);
//
//      }
//
//      copy_props(set);
//      copy_relevant_elem_props(set);
//
//
//    }



    Track::Set                                  Tractlet::Set::to_tracks(size_t num_length_sections, size_t num_width_sections) const {

      return Track::Set(*this, num_length_sections, num_width_sections);

    }


    Track::Set                                  Tractlet::Set::to_tracks(size_t num_length_sections, double strands_per_acs) const {

      return Track::Set(*this, num_length_sections, strands_per_acs);

    }


    Strand::Set                                 Tractlet::Set::to_strands(size_t num_samples) const {

      return Strand::Set(*this, num_samples);

    }


    Strand::Set                                 Tractlet::Set::to_strands(double strands_per_vol) const {

      return Strand::Set(*this, strands_per_vol);

    }


    std::vector<std::string>&                   Tractlet::Set::append_characteristic_keys(std::vector<std::string>& header) {

      header.push_back(Tractlet::ACS_EXT_PROP);
      header.push_back(Tractlet::LENGTH_EPSILON_COMPONENT_EXT_PROP);
      header.push_back(Tractlet::WIDTH_EPSILON_COMPONENT_EXT_PROP);
      header.push_back(Tractlet::AVG_DENSITY_EXT_PROP);
      return header;

    }


    void                                        Tractlet::Set::set_characteristics() {

      for (size_t tractlet_i = 0; tractlet_i < size(); ++tractlet_i) {
        set_extend_elem_prop(Tractlet::ACS_EXT_PROP,str(operator[](tractlet_i).acs()), tractlet_i);
        set_extend_elem_prop(Tractlet::LENGTH_EPSILON_COMPONENT_EXT_PROP,str(prop(LENGTH_EPSILON_PROP) * MR::Math::sqrt(operator[](tractlet_i)(0,1).norm())), tractlet_i);
        set_extend_elem_prop(Tractlet::WIDTH_EPSILON_COMPONENT_EXT_PROP,str(prop(WIDTH_EPSILON_PROP) * (operator[](tractlet_i)(1,0).norm() + operator[](tractlet_i)(2,0).norm())), tractlet_i);
        set_extend_elem_prop(Tractlet::AVG_DENSITY_EXT_PROP, str(operator[](tractlet_i).average_density()), tractlet_i);
      }

    }


    Tractlet::Set&                              Tractlet::Set::sanitize() {

      for (size_t tract_i = 0; tract_i < size(); ++tract_i)
        operator[](tract_i).sanitize();

      return *this;

    }


    void                                        Tractlet::Set::load(const std::string& location, size_t degree, size_t num_tractlets) {

      this->clear();

      if (File::has_or_txt_extension<Tractlet>(location)) {

        this->Base::Set<Tractlet>::load(location);

      } else if (File::has_or_txt_extension<Strand>(location) || File::has_or_txt_extension<Track>(location)) {

        Strand::Set strands (location, degree);

        if (num_tractlets)
          operator=(strands.to_tractlets(num_tractlets));
        else if (strands.has_extend_elem_prop(Track::BUNDLE_INDEX_EPROP))
          operator=(strands.to_tractlets());
        else
          throw Exception(" \"-width\" parameter is required for implicit conversion to tractlet set ('" + location + "') if '" + Track::BUNDLE_INDEX_EPROP + "' property is not present.");

      } else
         throw Exception("The extension of file \"" + location + "\" is not a recognised type (\"" + Tractlet::FILE_EXTENSION + "\", \"" + Strand::FILE_EXTENSION + "\" or \"" + Track::FILE_EXTENSION + "\").");

    }


    void                                        Tractlet::Set::save(const std::string& location, size_t num_length_sections, double strands_per_acs) const {

      if (File::has_or_txt_extension<Tractlet>(location)) {

        this->Base::Set<Tractlet>::save(location);

      } else if (File::has_or_txt_extension<Strand>(location)) {

        if (strands_per_acs <= 0.0)
          throw Exception(" \"-strands_per_acs\" parameter is required for implicit conversion to strand set ('" + location + "').");

        Strand::Set strands = this->to_strands(strands_per_acs);

        strands.save(location);

      } else if (File::has_or_txt_extension<Strand>(location)) {

        if (num_length_sections == 0)
          throw Exception(" \"-num_length_sections\" parameter is required for implicit conversion to track set ('" + location + "').");

        if (strands_per_acs <= 0.0)
          throw Exception(" \"-strands_per_acs\" parameter is required for implicit conversion to track set ('" + location + "').");

        Track::Set tcks = this->to_tracks(num_length_sections, strands_per_acs);

        tcks.save(location);

      } else
         throw Exception("The extension of file \"" + location + "\" is not a recognised type (\"" + Tractlet::FILE_EXTENSION + "\" or \"" + Strand::FILE_EXTENSION + "\").");

    }



    void                                        Tractlet::Set::project_to_sphere(double sphere_radius) {

      for (size_t tract_i = 0; tract_i < size(); ++tract_i)
        operator[](tract_i).project_to_sphere(sphere_radius);

    }



    Strand::Set                                 Tractlet::Set::backbones() const {

      Strand::Set strands(size(), degree());
      for (size_t tractlet_i = 0; tractlet_i < size(); tractlet_i++)
        strands[tractlet_i] = operator[](tractlet_i).backbone();

      strands.copy_props(*this);
      strands.copy_relevant_elem_props(*this);

      return strands;
    }


    Tractlet::Set                               Tractlet::Set::permute(const std::vector<size_t>& indices) const {

      if (var_elem_degrees())
        throw Exception ("'permute' is currently not implemented for variable element size sets.");

      Set permuted_set (size(), degree(), *props, *elem_props, get_extend_props());

      this->Fibre::Base::Set<Tractlet>::permute(permuted_set, indices);

      return permuted_set;

    }


    size_t                                      Tractlet::Set::max_degree() const {

      size_t max_degree = 0;

      for (size_t tract_i = 0; tract_i < size(); ++tract_i)
        if (operator[](tract_i).degree() > max_degree)
          max_degree = operator[](tract_i).degree();

      return max_degree;

    }



    std::string                                 Tractlet::Set::matlab_str() const {

      std::stringstream ss;

      ss << "{ 'b_intens', " << base_intensity() << ", ";

      for (size_t tractlet_i = 0; tractlet_i < size(); tractlet_i++) {

        ss << "{ 'acs', " << operator[](tractlet_i).acs() << ", ";

        for (size_t ax_i = 0; ax_i < 3; ax_i++) {

          ss << "[";

          for (size_t degree_i = 0; degree_i < operator[](tractlet_i)[ax_i].size(); degree_i++) {

            for (size_t dim_i = 0; dim_i < 3; dim_i++)
              ss << operator[](tractlet_i)[ax_i][degree_i][dim_i] << " ";

            if (degree_i != 2)
              ss << "; ";
          }

          ss << "]";

          if (tractlet_i != 2)
            ss << ", ";

        }

        ss << "}";

        if (tractlet_i != (size() -1))
          ss << ", ";
      }

      ss << "}";

      return ss.str();

    }

    Tractlet::Set                               Tractlet::Set::peel(double distance, size_t num_length_sections, size_t num_width_sections, bool to_cube) const {

      Tractlet::Set include;

      include.set_extend_prop("peel_distance", str(distance));

      if (to_cube)
        include.set_extend_prop("peel_shape", "cube");
      else
        include.set_extend_prop("peel_shape", "sphere");

      include.add_extend_elem_props(*this);


      for (size_t tract_i = 0; tract_i < size(); ++tract_i) {

        Track::Set tcks = operator[](tract_i).to_strands(num_width_sections).to_tracks(num_length_sections);

        bool in_roi = false;

        for (size_t tck_i = 0; tck_i < tcks.size(); ++tck_i) {
          for (size_t point_i = 0; point_i < tcks[tck_i].num_points(); point_i++) {

            Track tck = tcks[tck_i];

            if (to_cube) {
              if ((abs(tck[point_i][X]) <= distance) && (abs(tck[point_i][Y]) <= distance) && (abs(tck[point_i][Z]) <= distance)) {
                in_roi = true;
                break;
              }
            } else {
              if (tck[point_i].norm() <= distance) {
                in_roi = true;
                break;
              }
            }
          }

          if (in_roi)
            break;
        }

        if (in_roi)
          include.push_back(operator[](tract_i), get_extend_elem_prop_row(tract_i));

      }

      include.copy_props(*this);
      include.copy_relevant_elem_props(*this);

      return include;
    }


    void                                        Tractlet::Set::set_width_epsilon(double width_epsilon) {
      if (!has_prop(WIDTH_EPSILON_PROP))
        add_prop(WIDTH_EPSILON_PROP, 0.0);
      std::vector<double> old_acss;
      for (size_t tractlet_i = 0; tractlet_i < size(); ++tractlet_i)
        old_acss.push_back(operator[](tractlet_i).acs());
      prop(WIDTH_EPSILON_PROP) = width_epsilon;
      for (size_t tractlet_i = 0; tractlet_i < size(); ++tractlet_i)
        operator[](tractlet_i).set_acs(old_acss[tractlet_i]);
    }

    void                                        Tractlet::Set::set_length_epsilon(double length_epsilon) {
      if (!has_prop(LENGTH_EPSILON_PROP))
        add_prop(LENGTH_EPSILON_PROP, 0.0);
      std::vector<double> old_acss;
      for (size_t tractlet_i = 0; tractlet_i < size(); ++tractlet_i)
        old_acss.push_back(operator[](tractlet_i).acs());
      prop(LENGTH_EPSILON_PROP) = length_epsilon;
      for (size_t tractlet_i = 0; tractlet_i < size(); ++tractlet_i)
        operator[](tractlet_i).set_acs(old_acss[tractlet_i]);
    }


  }

}
