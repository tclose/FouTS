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

#ifndef __bts_fibre_tractlet_h__
#define __bts_fibre_tractlet_h__

#include <assert.h>

#include "math/matrix.h"

#include "bts/common.h"
#include "bts/file.h"

#include "bts/coord.h"
#include "bts/fibre/strand.h"
#include "bts/fibre/base/reader.h"
#include "bts/fibre/base/writer.h"
#include "bts/fibre/base/set.h"

namespace BTS {
    
    namespace Fibre {
        
        class Tractlet: public Base::Object {
                
                // Nested classes and typedefs.
            public:
                
                typedef Strand Element;

                class Set;
                class Section;
                class Tensor;
                class Walker;

                typedef Base::Reader<Tractlet> Reader;
                typedef Base::Writer<Tractlet> Writer;

                typedef Base::TextReader<Tractlet> TextReader;
                typedef Base::TextWriter<Tractlet> TextWriter;

                friend class Object;
                friend class Track;
                friend class Strand;
                friend class Set;
                friend class Base::Set<Tractlet>;
                friend class Tensor;
                friend class Base::Tensor<Tractlet>;
                friend class Base::Reader<Tractlet>;
                friend class Base::Writer<Tractlet>;
                friend class Base::Reader<Tractlet::Set>;
                friend class Base::Writer<Tractlet::Set>;

                //Public const static members,
            public:
                
                const static Coord FILE_SEPARATOR;
                const static std::string FILE_EXTENSION;

                const static char* ACS_EXT_PROP;
                const static char* WIDTH_EPSILON_COMPONENT_EXT_PROP;
                const static char* LENGTH_EPSILON_COMPONENT_EXT_PROP;
                const static char* AVG_DENSITY_EXT_PROP;

                const static double STRANDS_PER_AREA_DEFAULT;

                const static double REASONABLE_WIDTH;

                const static char* PROPS_LIST[];

                //Public static functions
            public:
                
                static size_t num_width_strands(size_t num_width_sections);

                static size_t total_num_sections(size_t num_length_sections,
                                                 size_t num_width_sections) {
                    return num_width_strands(num_width_sections) * num_length_sections;
                }
                
                //Protected member variables.
            protected:
                
                size_t dgree;
                const Base::Set<Tractlet>* parent;

                //Public member functions.
            public:
                
                Tractlet(size_t degree,
                         const std::vector<const char*>& props = std::vector<const char*>())
                        : Base::Object(3, degree * 9 + props.size(), props), dgree(degree), parent(
                                  0) {
                }
                
                //TODO: Try to make protected by making Base::Reader<Tractlet> a friend
            public:
                
                Tractlet(const std::vector<const char*>& props = std::vector<const char*>())
                        : Base::Object(3, 0, props), dgree(0), parent(0) {
                }
                
            public:
                
                Tractlet(const Tractlet& t)
                        : Base::Object(t), dgree(t.dgree), parent(0) {
                    set_acs(t.acs());
                }
                
                ~Tractlet() {
                }
                
                Tractlet& operator=(const Tractlet& t) {
                    Base::Object::operator=(t);
                    dgree = t.dgree;
                    if (has_var_acs() || t.has_var_acs())
                        set_acs(t.acs());
                    return *this;
                }
                
                Tractlet(const Strand& s, double width);

                Tractlet(Track primary_axis, size_t degree, double width);

            protected:
                
                /*! Used by Base::Set<Tractlet>::operator[]() to initiate a from a section of a larger state vector or
                 * by Base::Tensor<Tractlet>::row/column() for a slice of the tensor matrix.
                 *
                 * @param degree Degree of the tractlet
                 * @param view The view onto the larger vector or matrix
                 * @param props The properties stored in the set or tensor
                 */
                Tractlet(size_t degree, const MR::Math::Vector<double>::View& view,
                         std::vector<const char*>* props, const Base::Set<Tractlet>* parent)
                        : Base::Object((size_t) 3, view, props), dgree(degree), parent(parent) {
                }
                
                void from_track(Track primary_axis, size_t degree, double width);

            public:
                
                Strand operator[](size_t index) {
                    assert(index < 3);
                    return Strand(dgree, sub(index * degree() * 3, (index + 1) * degree() * 3),
                            &EMPTY_PROPS, 0);
                }
                
                Strand operator[](size_t index) const {
                    assert(index < 3);
                    return Strand(dgree, sub(index * degree() * 3, (index + 1) * degree() * 3),
                            &EMPTY_PROPS, 0);
                }
                
                Coord operator()(size_t axis_i, size_t degree_i) {
                    assert(axis_i < 3);
                    assert(bsize() % dgree == 0);
                    assert(bsize() / 9 == dgree);
                    return Coord(
                            sub((axis_i * degree() + degree_i) * 3,
                                    (axis_i * degree() + degree_i) * 3 + 3));
                }
                
                const Coord operator()(size_t axis_i, size_t degree_i) const {
                    assert(axis_i < 3);
                    assert(bsize() % dgree == 0);
                    assert(bsize() / 9 == dgree);
                    return Coord(
                            sub((axis_i * degree() + degree_i) * 3,
                                    (axis_i * degree() + degree_i) * 3 + 3));
                }
                
                //!Resizes each of the 3 axes to the new degree value.
                void redegree(size_t new_degree, double default_value = NAN);

                double acs() const;

                void set_acs(double acs);

                void remove_acs() {
                    remove_prop(ALPHA_PROP);
                }
                
                bool has_var_acs() const {
                    return has_prop(ALPHA_PROP);
                }
                
                double average_area(size_t num_points = 100);

                double average_density(size_t num_points = 100) {
                    return acs() / average_area(num_points);
                }
                
                void normalise_density(size_t num_points = 100) {
                    set_acs(average_area(num_points));
                }
                
                double& alpha() {
                    assert(has_var_acs());
                    return prop(ALPHA_PROP);
                }
                
                double alpha() const {
                    assert(has_var_acs());
                    return prop(ALPHA_PROP);
                }
                
                //Ensure that tractlets do not pass through themselves
                Tractlet& sanitize();

                std::vector<double> cross_sectional_areas(size_t num_points) const;

                Tractlet base() {
                    return Tractlet(dgree, sub(0, bsize()), &EMPTY_PROPS, 0);
                }
                
                const Tractlet base() const {
                    return Tractlet(dgree, sub(0, bsize()), &EMPTY_PROPS, 0);
                }
                
                Strand backbone() const {
                    return Strand(dgree, sub(0, dgree * 3), &EMPTY_PROPS, 0);
                }
                
                std::vector<Section>& sections(
                        std::vector<Section>& sections, size_t num_length_sections,
                        size_t num_width_sections, const Triple<double>& vox_lengths =
                                Triple<double>::Ones,
                        const Triple<double>& offsets = Triple<double>::Zeros,
                        size_t num_encodings = 0) const;

                /*!Section positions are normalized to the image so that the bottom-left-back (the minimal coordinate point)
                 corner is mapped to [0,0,0] and the voxel sizes are normalized to [1,1,1] so the inverse voxel lengths are required to calculate the gradient*/
                void add_section_gradient(const Tractlet& tractlet, const Section& section,
                                          const Strand::BasicSection& gradient);

                BASE_GENERAL_FUNCTIONS(Tractlet)
                ;

                BASE_MULT_DIVIDE_FUNCTIONS(Tractlet)
                ;

                BASE_ADD_SUBTRACT_FUNCTIONS(Tractlet)
                ;

                /*! Returns the degree of the tractlet (as determined by the size of the base vector divided by 3)
                 *
                 * @return The degree of the tractlet
                 */
                size_t degree() const {
                    return dgree;
                }
                
                /*! Projects extends or trims the primary axis so that its endpoints lie on the sphere.
                 *
                 * @return The degree of the tractlet
                 */
                void project_to_sphere(double sphere_radius) {
                    operator[](0).project_to_sphere(sphere_radius);
                }
                
                Strand::Set to_strands(size_t num_samples) const;

                Strand::Set to_strands(double strands_per_acs) const;

                double distance(const Tractlet& reference, bool& flipped, bool& switched,
                                bool& invert1, bool& invert2) const;

                /*!If 'num_strands' argument is non-zero, instead of measuring the distance between parameters,
                 the distance between generated strands is calculated instead*/
                double distance(const Tractlet& reference, double strands_per_acs = 0.0) const;

                Tractlet smallest_distance_set(const Tractlet& reference) const;

                //Doesn't actually effect the shape of the tract only the direction it travels.
                Tractlet flip() const;

                //Doesn't actually effect the shape of the tract only the order of its secondary axes.
                Tractlet switch_axes() const;

                //Doesn't actually effect the shape of the tract only the sign of its secondary axes.
                Tractlet invert_axis(size_t axis_i) const;

                double rotation() const;

                //! Only included to get template function to work (a strand will never have element properties.
                Tractlet& push_back(const Strand& strand, bool without_properties) {
                    return push_back(strand);
                }
                
                double path_length(size_t num_points = 1000) const;

            protected:
                
                void set_parent(const Base::Set<Tractlet>* prent) {
                    parent = prent;
                }
                
                void clear() {
                    Object::clear();
                    dgree = 0;
                }
                
                void clear(const std::vector<const char*>& properties) {
                    Object::clear(properties);
                    dgree = 0;
                }
                
                /*! Only used in reader as tractlets should always only have 3 strands. Is only named push_back to fit in
                 *  with other object classes. NB: This function doesn't actually resize the tractlet class (because, as
                 *  mentioned above, it always has 3 strands) but checks to see which of its strands are still invalidated (set
                 *  to NAN)
                 *
                 * @param strand The strand to be "pushed" onto the back of the strand.
                 * @return
                 */

                Tractlet& push_back(const Strand& strand);

                //TODO: work out how to remove this function
                //!Only implemented for resizing empty tractlets and for use in general template functions.
                void resize(size_t num_axes, const Strand& default_value);

                //TODO: work out how to remove this function
                //!Only implemented for resizing empty tractlets and for use in general template functions.
                void resize(size_t num_axes, double fill_value = NAN) {
                    Base::Object::resize(num_axes, dgree * 3, fill_value);
                }
                
                friend std::ostream& operator<<(std::ostream& stream, const Tractlet& tractlet);
                
        };
        
        Tractlet operator+(double c, Tractlet t);
        Tractlet operator*(double c, Tractlet t);
    
    }
}

#endif
