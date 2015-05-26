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

#include "math/math.h"
#include "bts/prob/prior_component/acs.h"

namespace FTS {
    
    namespace Prob {
        
        namespace PriorComponent {
            
            const double ACS::SCALE_DEFAULT = 0.0;    //1000
            const double ACS::MEAN_DEFAULT = 0.05;
            const std::string ACS::NAME = "acs";
            
            double ACS::log_prob(const Fibre::Tractlet tractlet) {
                
                return -scale * MR::Math::pow2(tractlet.acs() - mean);
                
            }
            
            double ACS::log_prob(const Fibre::Tractlet tractlet, Fibre::Tractlet gradient) {
                
                // If gradient hasn't been initialised, initialise it to the size of the tractlet, otherwise check its degree.
                if (!gradient.degree()) {
                    gradient = tractlet;
                    gradient.zero();
                } else if (gradient.degree() != tractlet.degree())
                    throw Exception(
                            "Tractlet degree (" + str(tractlet.degree())
                            + ") does not match that of the supplied gradient ("
                            + str(gradient.degree()) + ").");
                
                double lprob = -scale * MR::Math::pow2(tractlet.acs() - mean);
                
                if (tractlet.has_var_acs())
                    gradient.alpha() -= 4.0 * scale * tractlet.alpha() * (tractlet.acs() - mean);
                
                return lprob;
            }
        
        }
    
    }

}
