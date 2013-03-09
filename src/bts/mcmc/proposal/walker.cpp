/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Jul 28, 2010.

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

#include "bts/common.h"
#include "bts/mcmc/proposal/walker.h"

namespace BTS {
    
    namespace MCMC {
        
        namespace Proposal {
            
            const std::string Walker::TYPE_DEFAULT = "standard";
            const double Walker::STEP_SCALE_DEFAULT = 1.0;
            const std::string Walker::STEP_LOCATION_DEFAULT =
                    "/home/tclose/data/mcmc/parameters/walker.sta";
            
            void Walker::set(Distribution* const proposal_distribution,
                             MR::Math::Vector<double> relative_step_sizes) {
                
                if (prop_distr)
                    delete prop_distr;
                
                this->prop_distr = proposal_distribution->clone();
                
                step_sizes = relative_step_sizes;
                
            }
            
            void Walker::step(const MR::Math::Vector<double>& current,
                              MR::Math::Vector<double>& proposed, double scalar) const {
                
                assert(current.size() == proposed.size());
                
                for (size_t elem_i = 0; elem_i < current.size(); elem_i++)
                    proposed[elem_i] = prop_distr->sample(current[elem_i],
                            step_sizes[elem_i] * scalar);
                
            }
        
        }
    
    }

}
