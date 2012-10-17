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


#include "math/matrix.h"

#include "bts/triple.h"

#include "bts/common.h"
#include "bts/math/common.h"

#include "bts/math/svd.h"
#include "bts/fibre/tractlet.h"
#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"
#include "bts/fibre/tractlet/section.h"

#include "phantom/subdiv/subdiv.h"
#include "phantom/resample/resample.h"

#include "phantom/interface.h"

//#include "bts/image/expected/buffer.h"

namespace BTS {

  namespace Fibre {


    const Coord         Tractlet::FILE_SEPARATOR = Triple<double> (-INFINITY, -INFINITY, -INFINITY);
    const std::string  Tractlet::FILE_EXTENSION = "tct";
    
    const char*         Tractlet::PROPS_LIST[] = { Object::ALPHA_PROP, Tractlet::LENGTH_ACS_PROP, Tractlet::WIDTH_ACS_PROP, PROPS_LIST_END };

    const char*         Tractlet::LENGTH_ACS_PROP = "length_acs";
    const char*         Tractlet::WIDTH_ACS_PROP = "width_acs";

    const double        Tractlet::STRANDS_PER_AREA_DEFAULT = 1000;

    const double        Tractlet::REASONABLE_WIDTH = 0.1;

    const double        Tractlet::LENGTH_EPSILON_DEFAULT = 0.01;
    const double        Tractlet::WIDTH_EPSILON_DEFAULT = 0.01;

//    Tractlet::Tractlet (std::vector<Strand> axes, double base_width, double acs)
//
//      : acs(acs), base_width(base_width), axes(axes) {
//
//      if (axes.size() != 3)
//        throw Exception ("axes must be of length 3 (one for each spatial dimension).");
//
//      if (axes[1].degree() != axes[0].degree() || axes[2].degree() != axes[0].degree())
//        throw Exception ("axes degrees do not match {" + str(axes[0].degree()) + ", " + str(axes[1].degree()) + ", " + str(axes[2].degree()) +"}");
//
//      this->axes[1][0].normalise();
//      this->axes[2][0].normalise();
//
//      init();
//
//    }



    Tractlet::Tractlet (Track t, size_t degree, double width)
      : Base::Object((size_t)3, 3 * degree, select_props<Tractlet>(*t.props)), dgree(degree) {

      from_track(t, degree, width);

      for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
        prop(prop_i) = t.prop(t.prop_key(prop_i));

    }


    Tractlet::Tractlet(const Strand& s, double width)
      : Base::Object((size_t)3, 3 * s.degree(), select_props<Tractlet>(*s.props)), dgree(s.degree()) {

      from_track(Track(s, Track::NUM_LENGTH_SECTIONS_DEFAULT), s.degree(), width);

      for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
        prop(prop_i) = s.prop(s.prop_key(prop_i));

    }


    void                          	Tractlet::from_track(Track primary_axis, size_t degree, double width) {

      //acs(acs), base_width(1.0), axes(3) {

      //Requires a bit of stuffing around changing formats to use subdiv function from NFG suite.
      Fibre::Track::Set tcks(1, primary_axis.num_points()), subdivided;
      tcks[0] = primary_axis.base();

      std::vector<Triple<double> > pre_points;
      std::vector<Triple<double> > post_points;

      ::Strand_collection c, subdivided_c;//;

      ::generate_pre_points(tcks, pre_points);

      ::generate_post_points(tcks, post_points);

      //Assign the radius of the primary axis to be that of the desired width.
      tcks.add_extend_elem_prop(BTS::Fibre::Track::RADIUS_PROP, str(width));

      ::convert_mr_to_nfg(&c, tcks, pre_points, post_points);

      //Subdivide the strand so that 7 subdivided strands are in place of the original axis.
      ::subdivide_collection(&subdivided_c, &c, width/2.0001);

      ::convert_nfg_to_mr(subdivided, pre_points, post_points, &subdivided_c);

      //Convert subdivided back into strands.
      Strand::Set subdiv_strands = subdivided.to_strands(degree);

      Strand ax1_max = (subdiv_strands[0] + subdiv_strands[1])/2.0;
      Strand ax2_max = subdiv_strands[2];

      operator[](0) = primary_axis.base().to_strand(degree);

      operator[](1) = ax1_max - operator[](0);
      operator[](2) = ax2_max - operator[](0);

      for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
        prop(prop_i) = primary_axis.prop(prop_key(prop_i));

    }


    std::vector<double>             Tractlet::cross_sectional_areas(size_t num_points) const {

      std::vector<double> areas(num_points);

      // Get tangents to the primary path of the tractlet at 'num_points' along its path.
      Fibre::Track primary_tang = operator[](0).to_tangents(num_points);

      // Get the extent of the auxiliary axis at 'num_points' along its path.
      Fibre::Track aux_disp_1 = operator[](1).to_track(num_points);
      Fibre::Track aux_disp_2 = operator[](2).to_track(num_points);

      for (size_t point_i = 0; point_i < num_points; ++point_i) {

        // Get the normal to the oval defined by the two auxiliary axes at tau=point_i/num_points along the tractlet,
        // scaled by the oval's area.
        Coord area_normal = (M_PI / 4.0) * aux_disp_1[point_i].cross(aux_disp_2[point_i]);

        // Normalise the tangent of the primary axis at the current point (!note that it is done in place as it is not used again).
        const Coord norm_tangent = primary_tang[point_i].normalise();

        // Get the component of the ovals area that is perpendicular to the tractlet's primary
        areas[point_i] = MR::Math::abs(area_normal.dot(norm_tangent));

      }

      return areas;

    }





    //Ensure that tractlets do not pass through themselves
    Tractlet&                     	Tractlet::sanitize() {

      Coord primary = operator[](0)[1];
      Coord norm_primary = primary;

      if (norm_primary.norm() != 0) {

        norm_primary.normalise();


        std::vector<Triple<double> > perps(2);

        //Get the component of the first axis that is perpendicular to the primary orientation
        perps[0] = operator()(0,0) - norm_primary * operator()(0,0).dot(norm_primary);

        //Get the component of the first axis that is perpendicular to both the primary orientation and first axis
        Triple<double> perp2_orient = perps[0].cross(primary).normalise();

        perps[1] = perp2_orient * operator[](2)[0].dot(perp2_orient);

        for (size_t ax_i = 1; ax_i < 3; ++ax_i) {

          Coord base = operator[](ax_i)[0];
          Coord perp = perps[ax_i-1];

          if (perp.norm() < REASONABLE_WIDTH) {
            base -= perp;
            perp *= REASONABLE_WIDTH / perp.norm();
            base += perp;
          }

          for (size_t degree_i = 1; degree_i < degree(); ++degree_i)
            for (size_t dim_i = 0; dim_i < 3; ++dim_i) {

              double& higher = operator[](ax_i)[degree_i][dim_i];
              double ref = abs(perp[dim_i] / M_SQRT2);

              //If the higher degree is larger than the base degree reset it to be less than the base degree,
              //while preserving end points.
              if (abs(higher) > 1.8 * ref)
                higher = Math::sign(higher) * ref * 0.9;
              else if (abs(higher) > abs(perp[dim_i]))
                higher = Math::sign(higher) * (abs(higher) - abs(perp[dim_i]));

            }

        }

      }

      return *this;

    }


    void                          	Tractlet::add_section_gradient(const Tractlet& tractlet, const Section& section, const Strand::BasicSection& gradient) {


      Strand pos_gradient = Strand::outer_product(section.position_coeffs, gradient.position());

      //Add position gradient.
      operator[](0) += pos_gradient;
      operator[](1) += pos_gradient * section.ax1_fraction;
      operator[](2) += pos_gradient * section.ax2_fraction;

      //Add tangent gradient.

      Strand tang_gradient = Strand::outer_product(section.tangent_coeffs,  gradient.tangent() * section.length_fraction);

      operator[](0) += tang_gradient;
      operator[](1) += tang_gradient * section.ax1_fraction;
      operator[](2) += tang_gradient * section.ax2_fraction;


      if (has_var_acs()) {

        if (true) {

          throw Exception ("tied_width gradient needs to be adjusted after accounting for tractlet sheer.");

          alpha() += gradient.intensity() * (tractlet[1][0].norm() * tractlet[2][0].norm() - tractlet[1][0].dot(tractlet[2][0]));
          operator()(1,0) += gradient.intensity() * tractlet.acs() * (tractlet[1][0] * (tractlet[2][0].norm() / tractlet[1][0].norm()) - tractlet[2][0]);
          operator()(2,0) += gradient.intensity() * tractlet.acs() * (tractlet[2][0] * (tractlet[1][0].norm() / tractlet[2][0].norm()) - tractlet[1][0]);

        } else
          alpha() += gradient.intensity();
      }


    }


    size_t                          Tractlet::num_width_strands(size_t num_width_sections) {

      double width_fraction = 1.0 / (double)num_width_sections;

      size_t num_strands = 0;

      for (double ax1_frac = (-1.0 + width_fraction); ax1_frac < 1.0; ax1_frac += 2.0 * width_fraction)
        for (double ax2_frac = (-1.0 + width_fraction); ax2_frac < 1.0; ax2_frac += 2.0 * width_fraction)
          if (MR::Math::pow2(ax1_frac) + MR::Math::pow2(ax2_frac) <= (1.0 - width_fraction / 2.0))
            ++num_strands;

      return num_strands;

    }


    std::vector<Tractlet::Section>& Tractlet::sections(std::vector<Tractlet::Section>& sections,
                                                          size_t num_length_sections,
                                                          size_t num_width_sections,
                                                          const Triple<double>& vox_lengths,
                                                          const Triple<double>& offsets,
                                                          size_t num_encodings) const {


      sections.resize(total_num_sections(num_length_sections, num_width_sections), Tractlet::Section(num_encodings));

      double length_fraction = 1.0 / (double)num_length_sections;
      double width_fraction = 1.0 / (double)num_width_sections;

      const MR::Math::Matrix<double>& position_matrix  = Strand::position_matrix(num_length_sections, this->degree());
      const MR::Math::Matrix<double>& tangent_matrix   = Strand::tangent_matrix(num_length_sections, this->degree());

      size_t section_count = 0;
      //Loop over length of tractlet.
      for (size_t section_i = 0; section_i < num_length_sections; section_i++) {

        // Loop across cross-section of tractlet. Loop each perturbation axis from a 'width_fraction' from -1 to a width_fraction
        // from +1, in intervals of 2*width_fraction's. It is performed this way because each section is to represent the
        // tractlet +- a width fraction in both directions about its centre ('position').
        for (double ax1_frac = (-1.0 + width_fraction); ax1_frac < 1.0; ax1_frac += 2.0 * width_fraction) {

          for (double ax2_frac = (-1.0 + width_fraction); ax2_frac < 1.0; ax2_frac += 2.0 * width_fraction) {

            // Include section if the "coordinate" consisting of the combined fractions along each perturbation axes lies
            // withinness the unit circle minus a half a width_fraction. This will give the resulting cross-sections an
            // elliptical shape (or circular if axes are equal).
            if (MR::Math::pow2(ax1_frac) + MR::Math::pow2(ax2_frac) <= (1.0 - width_fraction / 2.0)) {

              sections[section_count].set(*this, position_matrix.row(section_i), tangent_matrix.row(section_i), ax1_frac, ax2_frac, length_fraction, width_fraction);
              sections[section_count].normalize(vox_lengths, offsets);

              ++section_count;

            }

          }

        }

      }

      return sections;

    }


    Strand::Set                     Tractlet::to_strands(size_t num_width_sections) const {

      Strand::Set strands(std::vector<const char*>(), this->prop_keys());

      double width_fraction = 1.0 / (double)num_width_sections;

      // Loop across cross-section of tractlet. Loop each perturbation axis from a 'width_fraction' from -1 to a width_fraction
      // from +1, in intervals of 2*width_fraction's. It is performed this way because each strand is to represent the
      // tractlet +- a width_fraction radius about it.
      for (double ax1_frac = (-1.0 + width_fraction); ax1_frac < 1.0; ax1_frac += 2.0 * width_fraction) {

        for (double ax2_frac = (-1.0 + width_fraction); ax2_frac < 1.0; ax2_frac += 2.0 * width_fraction) {

          //Include section if coord consisting of the combined fractions along each perturbation axes lies withinness the
          // unit circle minus a half a width_fraction. This will give the resulting cross-sections an elliptical
          // shape (or circular if axes are equal).
          if (MR::Math::pow2(ax1_frac) + MR::Math::pow2(ax2_frac) <= (1.0 - width_fraction /2.0)) {
            Strand s(degree(), this->prop_keys());
            s.base() = operator[](0) + operator[](1) * ax1_frac + operator[](2) * ax2_frac;
            for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
              s.prop(prop_i) = prop(prop_i);
            strands.push_back(s);
          }

        }

      }

      return strands;

    }


    Strand::Set                    	Tractlet::to_strands(double strands_per_acs) const {

      double num_strands = strands_per_acs * acs();

      size_t num_width_sections = (size_t)MR::Math::round(MR::Math::sqrt(num_strands * 4.0 / M_PI));

      return to_strands(num_width_sections);

    }



    Tractlet                       	Tractlet::flip() const {

      Tractlet t(*this);

      for (size_t ax_i = 0; ax_i < 3; ax_i++)
        t[ax_i] = operator[](ax_i).flip();

      return t;

    }


    Tractlet                       	Tractlet::switch_axes() const {

      Tractlet t(*this);

      t[1] = operator[](2);
      t[2] = operator[](1);

      return t;

    }


    Tractlet                       	Tractlet::invert_axis(size_t ax_i) const {

      assert(ax_i == 1 || ax_i == 2);

      Tractlet t(*this);

      t[ax_i] = operator[](ax_i) * -1;

      return t;

    }


    double                          Tractlet::distance(const Tractlet& reference, double strands_per_acs) const {

      double dist;

      if (strands_per_acs <= 0.0) {

        bool flipped;
        bool switched;
        bool invert1;
        bool invert2;

        dist = distance(reference,flipped,switched,invert1, invert2);

      } else {

        Strand::Set this_strands = to_strands(strands_per_acs);
        Strand::Set reference_strands = reference.to_strands(strands_per_acs);

        if (this_strands.size() > reference_strands.size())
          dist = reference_strands.distance(this_strands);
        else
          dist = this_strands.distance(reference_strands);

      }

      return dist;

    }


    double                        	Tractlet::distance(const Tractlet& reference, bool& flipped, bool& switched, bool& invert1, bool& invert2) const {

      double min_dist = INFINITY;

      //There are four ambiguities in the Tractlet model:
      //  - Flipping the direction of the sections
      //  - Switching the two axes
      //  - Inverting the first axis
      //  - Inverting the second axis
      // Each combination is checked to see which is the minimal distance.

      flipped = false;
      switched = false;
      invert1 = false;
      invert2 = false;

      for (size_t invt1 = 0; invt1 <= 1; ++invt1)
        for (size_t invt2 = 0; invt2 <= 1; ++invt2)
          for (size_t flip = 0; flip <= 1; ++flip)
            for (size_t swtch = 0; swtch <= 1; ++swtch) {

              Tractlet modified (reference);

              if (invt1)
                modified = modified.invert_axis(1);

              if (invt2)
                modified = modified.invert_axis(2);

              if (flip)
                modified = modified.flip();

              if (swtch)
                modified = modified.switch_axes();

              double dist = (*this - modified).norm2();

              if (dist < min_dist) {
                min_dist = dist;
                flipped = flip;
                switched = swtch;
                invert1 = invt1;
                invert2 = invt2;
              }

            }

      return min_dist;

    }


    Tractlet                       	Tractlet::smallest_distance_set(const Tractlet& reference) const {

      bool flipped;
      bool switched;
      bool invert1;
      bool invert2;

      distance(reference, flipped, switched, invert1, invert2);

      Tractlet t (*this);

      if (invert1)
        t = t.invert_axis(1);

      if (invert2)
        t = t.invert_axis(2);

      if (flipped)
        t = t.flip();

      if (switched)
        t = t.switch_axes();

      return t;

    }


    double                          Tractlet::rotation() const {


      std::pair< Triple<double>, Triple<double> > ax1_endpoints = Strand(operator[](1)).endpoints();
      std::pair< Triple<double>, Triple<double> > ax2_endpoints = Strand(operator[](2)).endpoints();

      Triple<double> start = (ax1_endpoints.first.normalise() + ax2_endpoints.first.normalise()) / 2.0;
      Triple<double> finish = (ax1_endpoints.second.normalise() + ax2_endpoints.second.normalise()) / 2.0;

      Triple<double> primary_orient = operator[](0)[1];
      primary_orient.normalise();

      start = start - start.dot(primary_orient) * primary_orient;
      finish = finish - finish.dot(primary_orient) * primary_orient;

      start.normalise();
      finish.normalise();

      Triple<double> sign_vec =  start.cross(finish).normalise() - primary_orient;

      double sign;
      if (sign_vec.norm() > 1.0)
        sign = 1.0;
      else
        sign = -1.0;

      return MR::Math::acos(start.dot(finish)/(finish.norm() * start.norm())) * sign;


      return start.angle(finish);

    }


    //!Resizes each of the 3 axes to the new degree value.
    void                            Tractlet::redegree(size_t new_degree, double default_value) {

      Tractlet new_tractlet (new_degree, *this->props);
      new_tractlet.zero();

      for (size_t ax_i = 0; ax_i < degree(); ++ax_i)
        for (size_t degree_i = 0; degree_i < degree(); ++degree_i)
          new_tractlet(ax_i,degree_i) = operator()(ax_i,degree_i);

      *this = new_tractlet;

    }


    Tractlet                        operator+ (double c, Tractlet t) { return t + c; }


    Tractlet                        operator* (double c, Tractlet t) { return t * c; }


    std::ostream&                  	operator<< (std::ostream& stream, const Tractlet& tractlet) {

      if (tractlet.has_var_acs()) {
        stream << "ACS: " << tractlet.acs() << std::endl;
        // Output 10 densities from t=0 to t=1
        stream << "densities:" << std::endl;
        std::vector<double> areas = tractlet.cross_sectional_areas(10);
        for (size_t i = 0; i < 10; ++i)
          std::cout << tractlet.acs() / areas[i] << " ";
        std::cout << std::endl;
        stream << "sqrt(ACS) for density=1:" << std::endl;
        for (size_t i = 0; i < 10; ++i)
          std::cout << "(" << MR::Math::sqrt(areas[i]) << ") ";
        std::cout << std::endl;
      }

      stream << std::endl;

      for (size_t axis_i = 0; axis_i < 3; axis_i++)
        stream << "axis " + str(axis_i) + ": " << std::endl << tractlet[axis_i] << std::endl;

      return stream;

    }


    //!Only implemented for resizing empty tractlets and for use in general template functions.
    void                            Tractlet::resize(size_t num_axes, const Strand& default_value) {

      size_t old_size = size();

      resize(num_axes);

      for (size_t ax_i = old_size; ax_i < num_axes; ++ax_i)
        operator[](ax_i) = default_value;

    }


    Tractlet&                       Tractlet::push_back(const Strand& strand) {

      if (!dgree)
        dgree = strand.degree();
      else if (dgree != strand.degree())
        throw Exception ("Attempting to push a strand of degree " + str(strand.degree()) + " on a tractlet of degree " + str(dgree) + ".");

      assert(bsize() % dgree == 0);

      size_t strand_i = bsize() / (dgree * 3);

      if (strand_i > 2)
        throw Exception ("Attempting to push a fourth strand onto an already created tractlet.");

      resize(strand_i+1);
      MR::Math::Vector<double>::sub(strand_i * dgree * 3, (strand_i+1) * dgree * 3) = strand;

      return *this;

    }


  }

}
