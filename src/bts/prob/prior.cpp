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

#include "bts/prob/prior.h"
#include "bts/prob/prior_component/frequency.h"
#include "bts/prob/prior_component/hook.h"
#include "bts/prob/prior_component/density.h"
#include "bts/prob/prior_component/in_image.h"
#include "bts/image/observed/buffer.h"

#include "bts/image/inline_functions.h"
#include "bts/prob/inline_functions.h"

namespace FTS {
    
    namespace Prob {
        
        const char* Prior::TYPE_DEFAULT = "default";
        const double Prior::SCALE_DEFAULT = 1;
        
        Prior::Prior(double scale, double freq_scale, double freq_aux_scale, double hook_scale,
                     double hook_num_points, double hook_num_width_sections,
                     double density_high_scale, double density_low_scale, double density_num_points,
                     double acs_scale, double acs_mean, double length_scale, double length_mean,
                     double in_image_scale, double in_image_power, const Triple<double>& in_image_offset,
                     const Triple<double>& in_image_extent, size_t in_image_num_length_sections,
                     size_t in_image_num_width_sections)
                : scale(scale), frequency(freq_scale, freq_aux_scale), hook(hook_scale,
                  hook_num_points, hook_num_width_sections), density(density_high_scale,
                  density_low_scale, density_num_points), acs(acs_scale, acs_mean), length(
                  length_scale, length_mean),
                  in_image(in_image_scale, in_image_power, in_image_offset,
                           in_image_extent, in_image_num_length_sections,
                           in_image_num_width_sections) { }
        
        std::map<std::string, double> Prior::get_component_values(const Fibre::Strand fibres) {
            Fibre::Strand gradient;
            gradient = fibres;
            std::map<std::string, double> component_map;
            component_map[PriorComponent::Frequency::NAME] = frequency.log_prob(fibres, gradient);
            component_map[PriorComponent::Hook::NAME] = hook.log_prob(fibres, gradient);
            component_map[PriorComponent::Length::NAME] = length.log_prob(fibres, gradient);
            component_map[PriorComponent::Density::NAME] = 0.0;
            component_map[PriorComponent::ACS::NAME] = 0.0;
            return component_map;
        }
        
        std::map<std::string, double> Prior::get_component_values(const Fibre::Tractlet fibres) {
            Fibre::Tractlet gradient;
            gradient = fibres;
            std::map<std::string, double> component_map;
            component_map[PriorComponent::Frequency::NAME] = frequency.log_prob(fibres, gradient);
            component_map[PriorComponent::Hook::NAME] = hook.log_prob(fibres, gradient);
            component_map[PriorComponent::Length::NAME] = length.log_prob(fibres[0], gradient[0]);
            component_map[PriorComponent::InImage::NAME] = in_image.log_prob(fibres, gradient);
            component_map[PriorComponent::Density::NAME] = density.log_prob(fibres, gradient);
            component_map[PriorComponent::ACS::NAME] = acs.log_prob(fibres);
            return component_map;
        }
        
        double Prior::log_prob(const Fibre::Strand strand, Fibre::Strand gradient) {
            
            double lprob = 0.0;
            
            lprob += frequency.log_prob(strand, gradient);
            lprob += hook.log_prob(strand, gradient);
            lprob += length.log_prob(strand, gradient);
            
            return lprob;
            
        }
        
        double Prior::log_prob(const Fibre::Tractlet tractlet, Fibre::Tractlet gradient) {
            
            double lprob = 0.0;
            lprob += frequency.log_prob(tractlet, gradient);
            lprob += hook.log_prob(tractlet, gradient);
            lprob += length.log_prob(tractlet[0], gradient[0]);
            lprob += density.log_prob(tractlet, gradient);
            lprob += acs.log_prob(tractlet, gradient);
            lprob += in_image.log_prob(tractlet, gradient);
            
            return lprob;
            
        }
    
    }
}
