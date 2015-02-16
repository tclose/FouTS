/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

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

#include <bts/prob/prior_component/in_image.h>
#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/tractlet/section.h"
#include "bts/prob/inline_functions.h"
#include "bts/image/inline_functions.h"

namespace BTS {
    
    namespace Prob {
        
        namespace PriorComponent {
            
            const double InImage::SCALE_DEFAULT = 1e6;
            const size_t InImage::POWER_DEFAULT = 10;
            const std::string InImage::NAME = "inside_image";
            const size_t InImage::NUM_LENGTH_SECTIONS_DEFAULT = 100;
            const size_t InImage::NUM_WIDTH_SECTIONS_DEFAULT = 7;
            const double InImage::BORDER_DEFAULT = 0.5;
            
            double InImage::log_prob(const Fibre::Tractlet tractlet, Fibre::Tractlet gradient) {
                
                double log_prob = 0.0;
                gradient.invalidate();

                size_t num_sections = Fibre::Tractlet::total_num_sections(num_length_sections, num_width_sections);
                std::vector<Fibre::Tractlet::Section> sections = std::vector<Fibre::Tractlet::Section>(num_sections);
                tractlet.sections(sections, num_length_sections, num_width_sections,
                                  extent, offset);

                for (size_t section_i = 0; section_i < sections.size(); ++section_i) {
                    const Triple<double>& pos = sections[section_i].position();
                    const Triple<double>& x = (pos - extent / 2.0) / extent;
                    for (size_t dim_i = 0; dim_i < 3; ++dim_i)
                        log_prob -= scale * MR::Math::pow(x[dim_i], (double)power * 2.0);
                }

                return log_prob / (double) num_sections;
                
            }


            Triple<double> InImage::get_offset(const Image::Observed::Buffer& obs_image,
                                               double border) {

                return obs_image.offsets() - obs_image.vox_lengths() * border;

            }


            Triple<double> InImage::get_extent(const Image::Observed::Buffer& obs_image,
                                               double border) {
                return obs_image.vox_lengths() * (obs_image.dims() + 2.0 * border);
            }
        
        }
    
    }

}
