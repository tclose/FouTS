/*
 Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

 Written by Thomas G Close, 5/05/09.

 This file is part of Fourier Tract Sampling (FouTS).

 FouTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FouTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FTS.  If not, see <http://www.gnu.org/licenses/>.

 */

#include "bts/fibre/track.h"
#include "bts/fibre/strand.h"
#include "math/matrix.h"
//#include "math/matrix_ops.h"
#include "math/least_squares.h"

namespace FTS {
    
    namespace Fibre {
        
        const std::string Track::RADIUS_PROP = "track_radius";
        const std::string Track::BUNDLE_INDEX_EPROP = "bundle_index";
        const std::string Track::RADIUS_PROP_DEFAULT = "0.02";
        
        const Coord Track::FILE_SEPARATOR = Triple<double>(NAN, NAN, NAN);
        const std::string Track::FILE_EXTENSION = "tck";
        
        const char* Track::PROPS_LIST[] = { Object::ALPHA_PROP, PROPS_LIST_END };
        
//    const char*                       Track::INTRINSIC_PROPERTIES[Track::NUM_INTRINSIC_PROPERTIES] = {};
        
        Track::Track(const Strand& s, size_t num_points, bool include_endpoints)
                : Base::Object(num_points, num_points * 3 + s.num_props(),
                        select_props<Track>(*s.props)) {
            
            if (num_points == 0)
                throw Exception(
                        "'num_length_sections' option is required for conversion to tracks.");
            
            const MR::Math::Matrix<double>& position_matrix = Strand::position_matrix(num_points,
                    s.degree(), include_endpoints);
            
            for (size_t point_i = 0; point_i < num_points; point_i++)
                operator[](point_i) = s.left_product(position_matrix.row(point_i));
            
            for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
                prop(prop_i) = s.prop(prop_key(prop_i));
            
        }
        
        Strand Track::to_strand(size_t degree) const {
            
            return Strand(*this, degree);
            
        }
        
        Strand Track::to_strand(size_t degree, const std::vector<double>& timeline) const {
            
            return Strand(*this, degree, timeline);
            
        }
        
        Coord Track::left_product(const MR::Math::Vector<double>::View& row_vector) const {
            
            assert(row_vector.size() >= num_points());
            
            Coord product(0, 0, 0);
            
            for (size_t point_i = 0; point_i < num_points(); point_i++)
                product += operator[](point_i) * row_vector[point_i];
            
            return product;
            
        }
        
        std::ostream& operator<<(std::ostream& stream, const Track& tck) {
            stream << std::endl;
            for (size_t point_i = 0; point_i < tck.num_points(); point_i++) {
                stream << tck[point_i] << std::endl;
            }
            
            return (stream);
        }


        Track Track::resample(double step_size) {

            if (!size())
                throw Exception("Track has no points");

            Fibre::Track new_tck;
            // Add the first point of the track to the new track
            // The distance along the previous interval to start from
            // (a step_size away) is saved in in this variable.
            double t = 0.0;
            for (size_t point_i = 0; point_i + 1 < size(); ++point_i) {

                // Get the orientation and length of the current interval.
                const Triple<double> start = operator[](point_i);
                const Triple<double> end = operator[](point_i + 1);
                Triple<double> disp = end - start;
                double length = disp.norm();

                // Add the resampled points along the current inverval
                double step_incr = step_size / length;
                for (; t < 1.0; t += step_incr)
                    new_tck.push_back(start + disp * t);

                // Find the starting "time" on the next interval
                if (point_i + 2 < size()) {
                    t = 1.0;
                    while (t >= 1.0) {
                        // Get the first step on the current interval of the original track.
                        const Triple<double>& new_point = new_tck[new_tck.size() - 1];
                        const Triple<double>& next = operator[](point_i + 2);
                        // Solve for the "time" on the line passing through the next interval that is a 'step_size' away from the new_point
                        double a = next.norm2() - 2.0 * end.dot(next) + end.norm2();
                        double b = 2.0 * end.dot(next) - 2.0 * next.dot(new_point) - 2.0 * end.norm2() + 2.0 * end.dot(new_point);
                        double c = end.norm2() - 2.0 * end.dot(new_point) + new_point.norm2() - MR::Math::pow2(step_size);
                        // Get the times along the line joining the interval
                        double t1 = (-b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
                        double t2 = (-b - sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
                        // Find if the points actually lie on the interval
                        if (t1 >= 0.0 && t1 < 1.0)
                            t = t1;
                        else if (t2 >= 0.0 && t2 < 1.0)
                            t = t2;
                        else
                            ++point_i; // Skip to the next point if the next step is outside the next interval
                    }
                }

            }
            return new_tck;
        }
    }

}
