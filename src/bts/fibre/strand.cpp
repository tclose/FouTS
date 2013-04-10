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

#include "math/matrix.h"
#include "math/least_squares.h"

#include "bts/math/odr.h"
#include "bts/fibre/track.h"
#include "bts/fibre/tractlet.h"
#include "bts/fibre/strand/section.h"

namespace BTS {
    
    namespace Fibre {
        
        const Coord Strand::FILE_SEPARATOR = Triple<double>(NAN, NAN, NAN);
        
        const char* Strand::PROPS_LIST[] = { Object::ALPHA_PROP, PROPS_LIST_END };
        
        const std::string Strand::FILE_EXTENSION = "str";
        const size_t Strand::DEFAULT_DEGREE = 3;
        
        std::map<size_t, std::map<size_t, MR::Math::Matrix<double> > > Strand::position_matrices;
        std::map<size_t, std::map<size_t, MR::Math::Matrix<double> > > Strand::position_matrices_w_endpoints;
        std::map<size_t, std::map<size_t, MR::Math::Matrix<double> > > Strand::tangent_matrices;
        std::map<size_t, std::map<size_t, MR::Math::Matrix<double> > > Strand::inverse_position_matrices;
        std::map<size_t, std::map<size_t, MR::Math::Matrix<double> > > Strand::tangent_matrices_w_endpoints;
        std::map<size_t, std::map<size_t, MR::Math::Matrix<double> > > Strand::inverse_position_matrices_w_endpoints;
        
        const double MERGE_FUDGE_FACTOR = 1.2;
        
        Strand::Strand(const Track& t, size_t degree)
                : Base::Object(degree, (size_t) 3, select_props<Strand>(*t.props)) {
            
            size_t effective_degree = degree;
            
            if (t.num_points() < degree)
                effective_degree = t.num_points();
            
            const MR::Math::Matrix<double>& inverse_position_matrix =
                    Strand::inverse_position_matrix(t.num_points(), effective_degree);
            
            for (size_t degree_i = 0; degree_i < effective_degree; ++degree_i)
                operator[](degree_i) = t.left_product(inverse_position_matrix.row(degree_i));
            
            for (size_t degree_i = effective_degree; degree_i < degree; ++degree_i)
                operator[](degree_i).zero();
            
            for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
                prop(prop_i) = t.prop(prop_key(prop_i));
            
        }
        
        Strand::Strand(const Track& t, size_t degree, const std::vector<double>& timeline)
                : Base::Object(degree, (size_t) 3, select_props<Strand>(*t.props)) {
            
            size_t effective_degree = degree;
            
            if (t.num_points() < degree)
                effective_degree = t.num_points();
            
            MR::Math::Matrix<double> inverse_position_matrix =
                    Strand::create_inverse_position_matrix(timeline, effective_degree);
            
            for (size_t degree_i = 0; degree_i < effective_degree; degree_i++)
                operator[](degree_i) = t.left_product(inverse_position_matrix.row(degree_i));
            
            for (size_t degree_i = effective_degree; degree_i < degree; ++degree_i)
                operator[](degree_i).zero();
            
            for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
                prop(prop_i) = t.prop(prop_key(prop_i));
            
        }
        
        //Projects extends or trims the strand so that its endpoints lie on the sphere.
        void Strand::project_to_sphere(double sphere_radius) {
            
            Track tck = this->to_track(this->degree());
            
            tck[0].normalise();
            tck[0] *= sphere_radius;
            tck[tck.size() - 1].normalise();
            tck[tck.size() - 1] *= sphere_radius;
            
            *this = tck.to_strand(this->degree());
            
        }
        
        Track Strand::to_track(size_t num_points, bool include_endpoints) const {
            return Track(*this, num_points, include_endpoints);
        }
        
        Track Strand::to_tangents(size_t num_points, bool include_endpoints) const {
            
            if (num_points == 0)
                throw Exception("'num_points' option is required for conversion to tracks.");
            
            Track tck(num_points);
            const MR::Math::Matrix<double>& position_matrix = Strand::tangent_matrix(num_points,
                    degree());
            
            for (size_t point_i = 0; point_i < num_points; point_i++)
                tck[point_i] = left_product(position_matrix.row(point_i));
            
            return tck;
            
        }
        
        //! Flips the strand so that runs in the opposite direction but the path taken is the same.
        Strand Strand::flip() const {
            
            Strand s(*this);
            
            for (size_t degree_i = 1; degree_i < degree(); degree_i += 2)
                s[degree_i] = -operator[](degree_i);
            
            return s;
            
        }
        
        Strand& Strand::operator+=(const Strand& s) {
            
            assert(props_match(s));
            
            if (size() == s.size())
                
                MR::Math::Vector<double>::operator+=(s);
            
            else if (size() <= s.size()) {
                
                resize(s.degree(), 0.0);
                
                MR::Math::Vector<double>::operator+=(s);
                
            } else {
                
                for (size_t degree_i = 0; degree_i < s.size(); ++degree_i)
                    operator[](degree_i) += s[degree_i];
                
                for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
                    prop(prop_i) += s.prop(prop_i);
                
            }
            
            return *this;
            
        }
        
        Strand& Strand::operator-=(const Strand& s) {
            
            assert(props_match(s));
            
            if (size() == s.size())
                
                MR::Math::Vector<double>::operator-=(s);
            
            else if (size() <= s.size()) {
                
                resize(s.degree(), 0.0);
                
                MR::Math::Vector<double>::operator-=(s);
                
            } else {
                
                for (size_t degree_i = 0; degree_i < s.size(); ++degree_i)
                    operator[](degree_i) -= s[degree_i];
                
                for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
                    prop(prop_i) -= s.prop(prop_i);
                
            }
            
            return *this;
            
        }
        
        Strand Strand::outer_product(const MR::Math::Vector<double>& column_vector,
                                     const Coord& row) {
            
            Strand strand(column_vector.size());
            
            for (size_t row_i = 0; row_i < column_vector.size(); row_i++)
                strand[row_i] = column_vector[row_i] * row;
            
            return strand;
            
        }
        
        Strand Strand::merge(const Strand& strand, size_t num_sections, size_t new_degree) const {
            
            std::vector<const Strand*> to_merge;
            
            to_merge.push_back(this);
            
            to_merge.push_back(&strand);
            
            return merge(to_merge, num_sections, new_degree);
            
        }
        
        //!Generates, stores and returns position conversion matrices when required.
        const MR::Math::Matrix<double>& Strand::position_matrix(size_t num_sections, size_t degree,
                                                                bool include_endpoints) {
            
            assert(num_sections != 0);
            assert(degree != 0);
            
            if (include_endpoints) {
                
                MR::Math::Matrix<double>& position_matrix =
                        Strand::position_matrices_w_endpoints[num_sections][degree];
                
                if (position_matrix.rows())
                    return position_matrix;
                else
                    return position_matrix = Strand::create_position_matrix(degree, num_sections,
                            INCLUDE_ENDPOINTS);
                //TODO: Change dont_include_endpoints to enum.
                
            } else {
                
                MR::Math::Matrix<double>& position_matrix =
                        Strand::position_matrices[num_sections][degree];
                
                if (position_matrix.rows())
                    return position_matrix;
                else
                    return position_matrix = Strand::create_position_matrix(degree, num_sections,
                            DONT_INCLUDE_ENDPOINTS);
                
            }
            
        }
        
        //!Generates, stores and returns position conversion matrices when required.
        const MR::Math::Matrix<double>& Strand::tangent_matrix(size_t num_sections, size_t degree) {
            
            assert(num_sections != 0);
            assert(degree != 0);
            
            MR::Math::Matrix<double>& tangent_matrix =
                    Strand::tangent_matrices[num_sections][degree];
            
            if (tangent_matrix.rows())
                return tangent_matrix;
            else
                return tangent_matrix = Strand::create_tangent_matrix(degree, num_sections,
                        DONT_INCLUDE_ENDPOINTS);
            
        }
        
        //!Generates, stores and returns position conversion matrices when required.
        const MR::Math::Matrix<double>& Strand::inverse_position_matrix(size_t num_sections,
                                                                        size_t degree) {
            
            assert(num_sections != 0);
            assert(degree != 0);
            
            MR::Math::Matrix<double>& position_matrix =
                    Strand::inverse_position_matrices[degree][num_sections];
            
            if (position_matrix.rows())
                return position_matrix;
            else
                return position_matrix = Strand::create_inverse_position_matrix(num_sections,
                        degree, DONT_INCLUDE_ENDPOINTS);
            
        }
        
        //!Generates, stores and returns position conversion matrices when required.
        const MR::Math::Matrix<double>& Strand::position_matrix_w_endpoints(size_t num_sections,
                                                                            size_t degree) {
            
            assert(num_sections != 0);
            assert(degree != 0);
            
            MR::Math::Matrix<double>& position_matrix =
                    Strand::position_matrices_w_endpoints[num_sections][degree];
            
            if (position_matrix.rows())
                return position_matrix;
            else
                return position_matrix = Strand::create_position_matrix(degree, num_sections,
                        INCLUDE_ENDPOINTS);
            
        }
        
        //Generates, stores and returns position conversion matrices_w_endpoints when required.
        const MR::Math::Matrix<double>& Strand::tangent_matrix_w_endpoints(size_t num_sections,
                                                                           size_t degree) {
            
            assert(num_sections != 0);
            assert(degree != 0);
            
            MR::Math::Matrix<double>& tangent_matrix =
                    Strand::tangent_matrices_w_endpoints[num_sections][degree];
            
            if (tangent_matrix.rows())
                return tangent_matrix;
            else
                return tangent_matrix = Strand::create_tangent_matrix(degree, num_sections,
                        INCLUDE_ENDPOINTS);
        }
        
        //Generates, stores and returns position conversion matrices_w_endpoints when required.
        const MR::Math::Matrix<double>& Strand::inverse_position_matrix_w_endpoints(
                size_t num_sections, size_t degree) {
            
            assert(num_sections != 0);
            assert(degree != 0);
            
            MR::Math::Matrix<double>& position_matrix =
                    Strand::inverse_position_matrices_w_endpoints[degree][num_sections];
            
            if (position_matrix.rows())
                return position_matrix;
            else
                return position_matrix = Strand::create_inverse_position_matrix(num_sections,
                        degree, INCLUDE_ENDPOINTS);
            
        }
        
        void Strand::add_section_gradient(const Strand& strand, const Section& section,
                                          const BasicSection& gradient) {
            
            this->operator+=(outer_product(section.position_coeffs, gradient.position()));
            this->operator+=(
                    outer_product(section.tangent_coeffs,
                            gradient.tangent() * section.length_fraction));
            
            if (has_var_acs())
                var_acs() += gradient.intensity();
            
        }
        
        std::vector<Strand::Section>& Strand::sections(std::vector<Strand::Section>& sections,
                                                       size_t num_length_sections,
                                                       size_t num_width_sections,
                                                       const Triple<double>& vox_lengths,
                                                       const Triple<double>& offsets,
                                                       size_t num_encodings) const {
            
            sections.resize(num_length_sections, Strand::Section(num_encodings));
            
            double length_fraction = 1.0 / num_length_sections;
            
            const MR::Math::Matrix<double>& position_matrix = Strand::position_matrix(
                    num_length_sections, this->degree());
            const MR::Math::Matrix<double>& tangent_matrix = Strand::tangent_matrix(
                    num_length_sections, this->degree());
            
            for (size_t section_i = 0; section_i < num_length_sections; section_i++) {
                sections[section_i].set(*this, position_matrix.row(section_i),
                        tangent_matrix.row(section_i), length_fraction);
                sections[section_i].normalize(vox_lengths, offsets);
            }
            
            return sections;
            
        }
        
        // Used to convert Strand descriptors to a path of points
        MR::Math::Matrix<double> Strand::create_position_matrix(size_t degree, size_t num_sections,
                                                                bool include_endpoints) {
            
            if (!degree)
                throw Exception("The degree (" + str(degree) + ") cannot be 0.");
            
            if (!num_sections)
                throw Exception(
                        "The size of the number of length samples (" + str(num_sections)
                        + ") cannot be 0.");
            
            MR::Math::Matrix<double> Psi(num_sections, degree);
            
            double t_incr;
            
            if (include_endpoints)
                t_incr = 1.0 / (double) (num_sections - 1);
            else
                t_incr = 1.0 / (double) (num_sections + 1);
            
            for (size_t row_i = 0; row_i < num_sections; row_i++) {
                Psi(row_i, 0) = 1.0;
            }
            
            for (size_t degree_i = 1; degree_i < degree; degree_i++) {
                
                double t;
                if (include_endpoints)
                    t = 0.0;
                else
                    t = t_incr;
                
                for (size_t row_i = 0; row_i < num_sections; row_i++) {
                    
                    Psi(row_i, degree_i) = M_SQRT2 * cos(degree_i * M_PI * t);
                    t += t_incr;
                }
                
            }
            
            return Psi;
            
        }
        
        // Generate conversion matrix to convert Strand descriptors to path points from arbitrary 'time' points.
        MR::Math::Matrix<double> Strand::create_position_matrix(
                size_t degree, const std::vector<double>& timeline) {
            
            if (degree == 0 || timeline.size() == 0)
                throw Exception(
                        "The degree (" + str(degree) + ") or the size of the timeline vector ("
                        + str(timeline.size()) + ") cannot be 0.");
            
            size_t num_sections = timeline.size();
            
            MR::Math::Matrix<double> Psi(num_sections, degree);
            
            for (size_t row_i = 0; row_i < num_sections; row_i++) {
                Psi(row_i, 0) = 1.0;
            }
            
            for (size_t degree_i = 1; degree_i < degree; degree_i++) {
                
                for (size_t row_i = 0; row_i < num_sections; row_i++) {
                    Psi(row_i, degree_i) = M_SQRT2 * cos(degree_i * M_PI * timeline[row_i]);
                }
                
            }
            
            return Psi;
            
        }
        
        //Used to convert a path of points to Strand descriptors
        MR::Math::Matrix<double> Strand::create_inverse_position_matrix(size_t num_sections,
                                                                        size_t degree,
                                                                        bool include_endpoints) {
            
            MR::Math::Matrix<double> Psi = Strand::create_position_matrix(degree, num_sections,
                    include_endpoints);
            
            MR::Math::Matrix<double> inv_Psi(degree, num_sections);
            
            MR::Math::pinv(inv_Psi, Psi);
            
            return inv_Psi;
            
        }
        
        MR::Math::Matrix<double> Strand::create_inverse_position_matrix(
                const std::vector<double>& timeline, size_t degree) {
            
            MR::Math::Matrix<double> Psi = Strand::create_position_matrix(degree, timeline);
            
            MR::Math::Matrix<double> inv_Psi(degree, timeline.size());
            
            MR::Math::pinv(inv_Psi, Psi);
            
            return inv_Psi;
            
        }
        
        // Used to get the parameterization gradient at the points along the strand.
        MR::Math::Matrix<double> Strand::create_tangent_matrix(size_t degree, size_t num_sections,
                                                               bool include_endpoints) {
            
            MR::Math::Matrix<double> d_Psi(num_sections, degree);
            
            double t_incr;
            
            if (include_endpoints)
                t_incr = 1.0 / (double) (num_sections - 1);
            else
                t_incr = 1.0 / (double) (num_sections + 1);
            
            for (size_t row_i = 0; row_i < num_sections; row_i++) {
                d_Psi(row_i, 0) = 0.0;
            }
            
            for (size_t degree_i = 1; degree_i < degree; degree_i++) {
                
                double t;
                if (include_endpoints)
                    t = 0.0;
                else
                    t = t_incr;
                
                for (size_t row_i = 0; row_i < num_sections; row_i++) {
                    
                    d_Psi(row_i, degree_i) = -M_SQRT2 * M_PI * degree_i * sin(degree_i * M_PI * t);
                    t += t_incr;
                }
                
            }
            
            return d_Psi;
            
        }
        
        // Gets the endpoints of the strand.
        std::pair<Coord, Coord> Strand::endpoints() const {
            
            std::pair<Coord, Coord> endpoints;
            
            endpoints.first = operator[](0);
            endpoints.second = operator[](0);
            
            for (size_t odd_i = 1; odd_i < degree(); odd_i += 2) {
                endpoints.first += operator[](odd_i) * M_SQRT2;
                endpoints.second -= operator[](odd_i) * M_SQRT2;
            }
            
            for (size_t even_i = 2; even_i < degree(); even_i += 2) {
                endpoints.first += operator[](even_i) * M_SQRT2;
                endpoints.second += operator[](even_i) * M_SQRT2;
            }
            
            return endpoints;
            
        }
        
        // Get the gradient of the endpoints of the strand. Note that the gradient is the same for X, Y, & Z and each
        // dimension is only related to its corresponding fourier coefficient so a vector is sufficient to capture the
        // relationships between every Fourier coefficient and the each endpoint.
        std::pair<MR::Math::Vector<double>, MR::Math::Vector<double> > Strand::endpoint_coeffs() const {
            
            std::pair<MR::Math::Vector<double>, MR::Math::Vector<double> > endpoints_gradient(
                    MR::Math::Vector<double>(this->degree()),
                    MR::Math::Vector<double>(this->degree()));
            
            endpoints_gradient.first[0] = 1.0;
            endpoints_gradient.second[0] = 1.0;
            
            for (size_t odd_i = 1; odd_i < degree(); odd_i += 2) {
                endpoints_gradient.first[odd_i] = M_SQRT2;
                endpoints_gradient.second[odd_i] = -M_SQRT2;
            }
            
            for (size_t even_i = 2; even_i < degree(); even_i += 2) {
                endpoints_gradient.first[even_i] = M_SQRT2;
                endpoints_gradient.second[even_i] = M_SQRT2;
            }
            
            return endpoints_gradient;
            
        }
        
        // Gets the midpoint of the strand, which is the sum of even vectors.
        Coord Strand::midpoint() const {
            
            Coord midpoint = operator[](0);
            
            for (size_t degree_i = 2; degree_i < degree(); degree_i += 4)
                midpoint -= (double) M_SQRT2 * operator[](degree_i);
            
            for (size_t degree_i = 4; degree_i < degree(); degree_i += 4)
                midpoint += (double) M_SQRT2 * operator[](degree_i);
            
            return midpoint;
            
        }
        
        // Get the gradient of the midpoint of the strand. Note that the gradient is the same for X, Y, & Z and each
        // dimension is only related to its corresponding Fourier coefficient so a vector is sufficient to capture the
        // relationships between every Fourier coefficient and the midpoint.
        MR::Math::Vector<double> Strand::midpoint_gradient() const {
            
            MR::Math::Vector<double> gradient(this->degree());
            
            gradient.zero();
            
            gradient[0] = 1.0;
            
            for (size_t degree_i = 2; degree_i < this->degree(); degree_i += 4)
                gradient[degree_i] = -M_SQRT2;
            
            for (size_t degree_i = 4; degree_i < this->degree(); degree_i += 4)
                gradient[degree_i] = M_SQRT2;
            
            return gradient;
        }
        
        //Gets the squared distance between this strand and the reference strand.  Due to the Strand descriptor
        //parameterization being orthonormal, by Parseval's theorem it corresponds to the distance between identically-spaced
        //trains of path points along each strand.
        double Strand::distance(const Strand& reference, bool& flipped) const {
            
            double dist = MR::Math::sqrt((*this - reference).norm2());
            double flipped_dist = MR::Math::sqrt((*this - reference.flip()).norm2());
            
            if (flipped_dist < dist) {
                dist = flipped_dist;
                flipped = true;
            } else
                flipped = false;
            
            return dist;
            
        }
        
        Strand Strand::smallest_distance_set(const Strand& reference) const {
            
            bool flipped;
            
            distance(reference, flipped);
            
            Strand s(*this);
            
            if (flipped)
                s = s.flip();
            
            return s;
            
        }
        
        // Join two strands together into one longer strand of degree 'new_degree'.  If new_degree = 0 or is not specified it defaults to
        // the maximum degree of the strands + 1.
        Strand Strand::merge(std::vector<const Strand*>& strands, size_t num_points,
                             size_t new_degree) {
            
            std::vector<Coord> landmarks;
            landmarks.reserve(strands.size() * 3);
            
            size_t max_degree = 0;
            
            double avg_acs = 0;
            
            for (size_t strand_i = 0; strand_i < strands.size(); strand_i++) {
                
                const Strand& strand = *strands[strand_i];
                
                std::pair<Coord, Coord> e_points = strand.endpoints();
                
                Track tcks = strand.to_track(num_points, INCLUDE_ENDPOINTS);
                
                for (size_t tck_i = 0; tck_i < tcks.size(); ++tck_i)
                    landmarks.push_back(tcks[tck_i]);
                
                if (strand.degree() > max_degree)
                    max_degree = strand.degree();
                
                avg_acs += strand.acs();
            }
            
            avg_acs /= (double) strands.size();
            
            if (new_degree == 0)
                new_degree = max_degree + strands.size() - 1;
            
            // Get the orthogonal-distance regressor line between the path landmarks.
            Math::Line odr_line = Math::odr(landmarks);
            
            Track complete_tck(num_points * strands.size());
            
            std::vector<double> timeline;
            timeline.reserve(sizeof(double) * num_points * strands.size());
            
            double min_cos_time = INFINITY;
            double max_cos_time = -INFINITY;
            
            for (size_t strand_i = 0; strand_i < strands.size(); strand_i++) {
                
                Track tck = strands[strand_i]->to_track(num_points, INCLUDE_ENDPOINTS);
                
                for (size_t point_i = 0; point_i < num_points; point_i++) {
                    
                    Coord point = tck[point_i];
                    Coord disp = point - odr_line.centre;
                    double cos_time = odr_line.direction.dot(disp);
                    
                    complete_tck[strand_i * num_points + point_i] = point;
                    timeline.push_back(cos_time);
                    
                    if (cos_time > max_cos_time)
                        max_cos_time = cos_time;
                    
                    if (cos_time < min_cos_time)
                        min_cos_time = cos_time;
                    
                }
                
            }
            
            double cos_time_range = max_cos_time - min_cos_time;
            
            for (size_t time_i = 0; time_i < timeline.size(); time_i++)
                timeline[time_i] = (MR::Math::cos(
                                            (timeline[time_i] - min_cos_time) * (M_PI
                                                    / cos_time_range))
                                    + 1)
                                   * 0.5;
            
            Strand new_strand = complete_tck.to_strand(new_degree, timeline);
            
            if (strands[0]->has_var_acs())
                new_strand.set_acs(avg_acs);
            
            new_strand[1] *= MERGE_FUDGE_FACTOR;
            
            return new_strand;
            
        }
        
        std::ostream& operator<<(std::ostream& stream, const Strand& strand) {
            
            stream << std::endl;
            
            if (strand.has_var_acs())
                stream << "acs: " << strand.acs() << std::endl;
            
            for (size_t degree_i = 0; degree_i < strand.degree(); degree_i++) {
                stream << strand[degree_i] << std::endl;
            }
            return (stream);
            
        }
    
    }

}
