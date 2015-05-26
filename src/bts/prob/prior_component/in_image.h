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

#ifndef __bts_prob_prior_in_image__
#define __bts_prob_prior_in_image__

#include "bts/fibre/tractlet.h"
#include "bts/fibre/tractlet/tensor.h"
#include "bts/image/observed/buffer.h"

namespace FTS {
    
    namespace Prob {
        
        namespace PriorComponent {
            
            class InImage {
                    
                public:
                    
                    static Triple<double> get_offset(const Image::Observed::Buffer& obs_image,
                                              double border);
                    static Triple<double> get_extent(const Image::Observed::Buffer& obs_image,
                                                                  double border);

                public:

                    const static double SCALE_DEFAULT;
                    const static size_t POWER_DEFAULT;
                    const static std::string NAME;
                    const static size_t NUM_LENGTH_SECTIONS_DEFAULT;
                    const static size_t NUM_WIDTH_SECTIONS_DEFAULT;
                    const static double BORDER_DEFAULT;


                protected:
                    
                    double scale;
                    double double_power;
                    Triple<double> offset;
                    Triple<double> extent;
                    size_t num_length_sections;
                    size_t num_width_sections;

                public:
                    
                    InImage()
                     : scale(0.0), double_power(0.0), num_length_sections(0),
                       num_width_sections(0) {
                    }
                    
                    InImage(double scale, size_t power,
                            const Triple<double>& offset, const Triple<double>& extent,
                            size_t num_length_sections = NUM_LENGTH_SECTIONS_DEFAULT,
                            size_t num_width_sections = NUM_WIDTH_SECTIONS_DEFAULT)
                       : scale(scale), double_power((double)(power * 2)),
                         offset(offset), extent(extent),
                         num_length_sections(num_length_sections),
                         num_width_sections(num_width_sections) {

                    }
                    
                    InImage(const InImage& ii)
                       : scale(ii.scale), double_power(ii.double_power),
                         offset(ii.offset), extent(ii.extent),
                         num_length_sections(ii.num_length_sections),
                         num_width_sections(ii.num_width_sections) {
                    }
                    
                    InImage& operator=(const InImage& ii) {
                        scale = ii.scale;
                        double_power = ii.double_power;
                        offset = ii.offset;
                        extent = ii.extent;
                        num_length_sections = ii.num_length_sections;
                        num_width_sections = ii.num_width_sections;
                        return *this;
                    }
                    
                    InImage* clone() const {
                        return new InImage(*this);
                    }
                    
                    double log_prob(const Fibre::Tractlet, Fibre::Tractlet gradient);

                    double log_prob(const Fibre::Tractlet tractlet, Fibre::Tractlet gradient,
                                    Fibre::Tractlet::Tensor hessian) {
                        throw Exception("Not implemented");
                    }
                    
                    double log_prob_and_fisher(const Fibre::Tractlet tractlet,
                                               Fibre::Tractlet gradient,
                                               Fibre::Tractlet::Tensor hessian) {
                        throw Exception("Not implemented");
                    }
                    
                    double log_prob_and_fisher(
                            const Fibre::Tractlet tractlet, Fibre::Tractlet gradient,
                            Fibre::Tractlet::Tensor fisher,
                            std::vector<Fibre::Tractlet::Tensor>& fisher_gradient) {
                        throw Exception("Not implemented");
                    }
                    
                    const std::string& get_name() {
                        return NAME;
                    }
                    
            };
        
        }
    
    }

}

#endif
