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

#ifndef __bts_prob_prior_tractlet_density__
#define __bts_prob_prior_tractlet_density__

#include "bts/fibre/tractlet.h"
#include "bts/fibre/tractlet/tensor.h"
#include "bts/image/observed/buffer.h"

namespace BTS {
    
    namespace Prob {
        
        namespace PriorComponent {
            
            class InImage {
                    
                public:
                    
                    const static double SCALE_DEFAULT;
                    const static double POWER_DEFAULT;
                    const static std::string NAME;
                    const static size_t NUM_LENGTH_SECTIONS_DEFAULT;
                    const static size_t NUM_WIDTH_SECTIONS_DEFAULT;
                    const static double BORDER_DEFAULT;

                protected:
                    
                    double scale;
                    double power;
                    Image::Observed::Buffer obs_image;
                    size_t num_length_sections;
                    size_t num_width_sections;
                    Triple<double> centre;
                    Triple<double> lengths;

                public:
                    
                    InImage()
                     : scale(0.0), power(0.0), num_length_sections(0),
                       num_width_sections(0) {
                    }
                    
                    InImage(double scale, double power,
                            Image::Observed::Buffer obs_image,
                            size_t num_length_sections = NUM_LENGTH_SECTIONS_DEFAULT,
                            size_t num_width_sections = NUM_WIDTH_SECTIONS_DEFAULT,
                            double border = BORDER_DEFAULT)
                       : scale(scale), power(power), obs_image(obs_image)
                         num_length_sections(num_length_sections),
                         num_width_sections(num_width_sections) {
                        for (size_t dim_i = 0; dim_i < 3; ++dim_i) {
                            centre[dim_i] = (double)obs_image.dim(dim_i) / 2.0;
                            lengths[dim_i] = (double)obs_image.dim(dim_i) / 2.0 + border;
                        }

                    }
                    
                    InImage(const InImage& ii)
                       : scale(ii.scale), power(ii.power), obs_image(ii.obs_image),
                         num_length_sections(ii.num_length_sections),
                         num_width_sections(ii.num_width_sections),
                         centre(centre), lengths(lengths)) {
                    }
                    
                    InImage& operator=(const InImage& ii) {
                        scale = ii.scale;
                        power = ii.power;
                        obs_image = ii.obs_image;
                        num_length_sections = ii.num_length_sections;
                        num_width_sections = ii.num_width_sections;
                        centre = ii.centre;
                        lengths = ii.lengths;
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
