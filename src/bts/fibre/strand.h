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


#ifndef __bts_fibre_strand_h__
#define __bts_fibre_strand_h__

#include "bts/fibre/base/object.h"
#include "bts/coord.h"


#include "bts/fibre/base/writer.h"
#include "bts/fibre/base/reader.h"
#include "bts/fibre/base/set.h"

namespace BTS {

	namespace Fibre {

	  class Strand : public Base::Object {

      //Public nested classes and typedefs
      public:

        typedef Coord Element;

        class BasicSection;
        class Section;
        class Set;
        class Tensor;
        class Walker;

        typedef Base::Writer<Strand> Writer;
        typedef Base::Reader<Strand> Reader;

        typedef Base::TextWriter<Strand> TextWriter;
        typedef Base::TextReader<Strand> TextReader;

        friend class Track;
        friend class Tractlet;
        friend class Set;
        friend class Base::Set<Strand>;
        friend class Tensor;
        friend class Base::Tensor<Strand>;

      //Public static constants
      public:

        const static std::string    FILE_EXTENSION;
        const static Coord          FILE_SEPARATOR;

        const static int            INCLUDE_ENDPOINTS = 1;
        const static int            DONT_INCLUDE_ENDPOINTS = 0;

        const static size_t         DEFAULT_DEGREE;

        const static char*          PROPS_LIST[];


      //Public static variables
      public:

        static std::map<size_t,std::map<size_t,MR::Math::Matrix<double> > > position_matrices;
        static std::map<size_t,std::map<size_t,MR::Math::Matrix<double> > > tangent_matrices;
        static std::map<size_t,std::map<size_t,MR::Math::Matrix<double> > > inverse_position_matrices;

        static std::map<size_t,std::map<size_t,MR::Math::Matrix<double> > > position_matrices_w_endpoints;
        static std::map<size_t,std::map<size_t,MR::Math::Matrix<double> > > tangent_matrices_w_endpoints;
        static std::map<size_t,std::map<size_t,MR::Math::Matrix<double> > > inverse_position_matrices_w_endpoints;

      //Public static methods.
      public:

        // Generate conversion matrix to convert Strand descriptors to path points.
        static MR::Math::Matrix<double>    create_position_matrix(size_t degree, size_t num_sections, bool include_endpoints = INCLUDE_ENDPOINTS);
        // Generate conversion matrix to convert Strand descriptors to path points from arbitrary 'time' points.
        static MR::Math::Matrix<double>    create_position_matrix(size_t degree, const std::vector<double>& timeline);
        // Generate conversion matrix to convert Strand descriptors to the tangents of the strand at the path points.
        static MR::Math::Matrix<double>    create_tangent_matrix(size_t degree, size_t num_sections, bool include_endpoints = INCLUDE_ENDPOINTS);
        static MR::Math::Matrix<double>    create_inverse_position_matrix(size_t num_sections, size_t degree, bool include_endpoints = INCLUDE_ENDPOINTS);
        static MR::Math::Matrix<double>    create_inverse_position_matrix(const std::vector<double>& timeline, size_t degree);

        static const MR::Math::Matrix<double>& position_matrix (size_t num_sections, size_t degree, bool include_endpoints = DONT_INCLUDE_ENDPOINTS);
        static const MR::Math::Matrix<double>& tangent_matrix (size_t num_sections, size_t degree);
        static const MR::Math::Matrix<double>& inverse_position_matrix (size_t num_sections, size_t degree);

        static const MR::Math::Matrix<double>& position_matrix_w_endpoints (size_t num_sections, size_t degree);
        static const MR::Math::Matrix<double>& tangent_matrix_w_endpoints (size_t num_sections, size_t degree);
        static const MR::Math::Matrix<double>& inverse_position_matrix_w_endpoints (size_t num_sections, size_t degree);

        static Strand                     outer_product(const MR::Math::Vector<double>& column_vector, const Coord& coord);

      protected:

        const Base::Set<Strand>* parent;

      //Public methods
      public:

        Strand (size_t degree = 0, const std::vector<const char*>& props = std::vector<const char*>())
          : Base::Object(degree, degree * 3 + props.size(), props), parent(0) {}


        Strand (int degree, double default_value, const std::vector<const char*>& props = std::vector<const char*>())
          : Base::Object(degree, degree * 3 + props.size(), props), parent(0) { set(default_value); }


        Strand (const Strand& s)
          : Base::Object(s), parent(0) {}


        ~Strand () {}


        Strand                            base()
          { return Strand(sze, sub(0,bsize()), &EMPTY_PROPS, parent); }

        const Strand                      base() const
          { return Strand(sze, sub(0,bsize()), &EMPTY_PROPS, parent); }

        Strand&                           operator=(const Strand& strand)
          { Base::Object::operator=(strand); return *this; }

        //        Strand (int degree, const Coord& default_value, const std::vector<const char*>& props = std::vector<const char*>())
        //          : Base::Object(degree, (size_t)3, props) { set(default_value); }


        Strand(const Track& t, size_t degree);


        Strand(const Track& t, size_t degree, const std::vector<double>& timeline);


      protected:

        /*! Used by Base::Set<Strand>::operator[]() to initiate a from a section of a larger state vector or
         * by Base::Tensor<Strand>::row/column() for a slice of the tensor matrix
         *
         * @param degree Degree of the strand
         * @param view The view onto the larger vector or matrix
         * @param props The properties stored in the set or tensor
         */
        Strand(size_t degree, const MR::Math::Vector<double>::View& view, std::vector<const char*>* props,
                                                                                      const Base::Set<Strand>* parent)
          : Base::Object(degree, view, props), parent(parent) {}


      public:


        size_t                     	      degree() const
          { return size(); };


        double                        	  acs() const
          { return has_prop(ALPHA_PROP) ? MR::Math::pow2(prop(ALPHA_PROP)) : 1.0; /* The ACS is stored in the state vector via its squareroot to prevent it from becoming < 0 */}


        void                          	  add_acs(double acs = 1.0)
          { add_prop(ALPHA_PROP, acs); }


        void                          	  set_acs(double acs)
          { assert(acs >= 0); prop(ALPHA_PROP) = MR::Math::sqrt(acs); }


        void                          	  remove_acs()
          { remove_prop(ALPHA_PROP); }


        bool                              has_var_acs() const
          { return has_prop(ALPHA_PROP); }


        double&                           var_acs()
          { assert(has_var_acs()); return prop(ALPHA_PROP); }


        BASE_GENERAL_FUNCTIONS(Strand);


        BASE_MULT_DIVIDE_FUNCTIONS(Strand);


        //Conversion methods.
        Track                             to_track(size_t num_points, bool include_endpoints = DONT_INCLUDE_ENDPOINTS) const;

        Track                             to_tangents(size_t num_points, bool include_endpoints = DONT_INCLUDE_ENDPOINTS) const;


        //Projects extends or trims the strand so that its endpoints lie on the sphere.
        void                              project_to_sphere(double sphere_radius);


        //! Flips the strand so that runs in the opposite direction but the path taken is the same.
        Strand                            flip() const;


        Strand&                           operator+=(const Strand& s);


        Strand&                           operator-=(const Strand& s);


        Strand                            operator+(const Strand& s) const
          { Strand answer(*this); answer += s; return answer; }


        Strand                            operator-(const Strand& s) const
          { Strand answer(*this); answer -= s; return answer; }


        Tractlet                          to_tractlet(double width) const;


        /*!Split strand in to many 'Strand::Sections'. Note num_width_sections is ignored and is only included to get a template
         * function to work.*/
        std::vector<Section>&             sections(std::vector<Section>& sections, size_t num_length_sections,
                                                    size_t num_width_sections = 0,
                                                    const Triple<double>& vox_lengths = Triple<double>::Ones,
                                                    const Triple<double>& offsets = Triple<double>::Zeros,
                                                    size_t num_encodings = 0,
                                                    double width_epsilon = 0.0,
                                                    double length_epsilon = 0.0) const;

        //Calculates the contribution of the section gradient back to the original strand's gradient.
        void                              add_section_gradient( const Strand& strand,
                                                                  const Section& section,
                                                                  const BasicSection& gradient);

        // Get the endpoints of the strand.
        std::pair<Coord, Coord>           endpoints() const;


        /*! Get the gradient of the endpoints of the strand. Note that the gradient is the same for X, Y, & Z and each
         * dimension is only related to its corresponding fourier coefficient so a vector is sufficient to capture the
         * relationships between every Fourier coefficient and the each endpoint. */
        std::pair< MR::Math::Vector<double>, MR::Math::Vector<double> >           endpoint_coeffs() const;

        // Gets the point at the top of the 'arc', which is just the sum of even vectors.
        Coord                             midpoint() const;

        /*! Get the gradient of the midpoint of the strand. Note that the gradient is the same for X, Y, & Z and each
         * dimension is only related to its corresponding fourier coefficient so a vector is sufficient to capture the
         * relationships between every Fourier coefficient and the midpoint.*/
        MR::Math::Vector<double>          midpoint_gradient() const;


        //FIXME: Need to sort out the names of these functions.
        double                            distance(const Strand& reference, bool& flipped) const;


        double                            distance(const Strand& reference, double dummy = 0) const
          { bool flipped; return distance(reference, flipped); }


        Strand                            smallest_distance_set(const Strand& reference) const;


        // Join two strands together into one longer strand of degree 'new_degree'.  If new_degree = 0 or is not specified it defaults to
        // the maximum degree of the strands + 1.
        static Strand                     merge(std::vector<const Strand*>& strands, size_t num_sections, size_t new_degree = 0);


        Strand                            merge(const Strand& strand, size_t num_sections, size_t new_degree = 0) const;


        //! Return a MR::Math::Matrix<double> window to the underlying state vector interpreted as a matrix (without properties).
        MR::Math::Matrix<double>          bmatrix()
          { return Base::Object::bmatrix(3); }


        std::vector<double>               cross_sectional_areas(size_t num_points)
          { throw Exception ("Should not be used, just added to debug template function with tractlets."); }

      protected:

        void                          set_parent(const Base::Set<Strand>* prent)
          { parent = prent; }
    };


    std::ostream&         operator<< (std::ostream& stream, const Strand& strand);

	}

}

#undef LOOP

#endif /* __bts_fibre_strand_h__ */
