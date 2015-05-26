/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

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

#include "bts/prob/prior_component/hook.h"

#include "bts/prob/inline_functions.h"

namespace FTS {
    
    namespace Prob {
        
        namespace PriorComponent {
            
            const double Hook::SCALE_DEFAULT = 100000.0;
            const std::string Hook::NAME = "hook";
            const size_t Hook::NUM_POINTS_DEFAULT = 100;
            const size_t Hook::NUM_WIDTH_SECTIONS_DEFAULT = 15;
            
            double Hook::log_prob(const Fibre::Strand strand, Fibre::Strand gradient) {
                
                double lprob = 0.0;
                gradient.invalidate();
                
                std::pair<Coord, Coord> endpoints = strand.endpoints();
                Coord end_to_end = endpoints.second - endpoints.first;
                double end_to_end_length2 = end_to_end.norm2();
                end_to_end.normalise();
                
                Fibre::Track track(strand, num_points);
                
                for (size_t point_i = 0; point_i < num_points; ++point_i) {
                    
                    Coord start_diff = track[point_i] - endpoints.first;
                    double start_dot = end_to_end.dot(start_diff);
                    
                    //If the point on the track is before the start point of the strand w.r.t. the ultimate endpoint of the strand
                    if (start_dot < 0.0) {
                        double start_dist2 = start_diff.norm2();
                        lprob -= 0.5 * MR::Math::pow2(start_dot)
                                / (start_dist2 * end_to_end_length2);
                    } else {
                        
                        Coord end_diff = track[point_i] - endpoints.second;
                        double end_dot = end_to_end.dot(end_diff);
                        
                        if (end_dot > 0.0) {
                            double end_dist2 = end_diff.norm2();
                            lprob -= 0.5 * MR::Math::pow2(end_dot)
                                    / (end_dist2 * end_to_end_length2);
                            
                        }
                        
                    }
                    
                }
                
                lprob *= scale / (double) num_points;
                return lprob;
            }
            
            double Hook::log_prob(const Fibre::Tractlet tractlet, Fibre::Tractlet gradient) {
                
                double lprob = 0.0;
                gradient.invalidate();
                
                Fibre::Strand::Set strands = tractlet.to_strands(num_width_sections);
                
                for (size_t strand_i = 0; strand_i < strands.size(); ++strand_i) {
                    Fibre::Strand strand_gradient(strands[strand_i]);
                    lprob += log_prob(strands[strand_i], strand_gradient);
                }
                
                lprob /= (double) strands.size();
                return lprob;
            }
            
            double Hook::log_prob(const Fibre::Strand& strand, Fibre::Strand& gradient,
                                  Fibre::Strand::Tensor& hessian) {
                
                throw Exception("Not implemented yet");
                
            }
        
        }
    
    }

}
