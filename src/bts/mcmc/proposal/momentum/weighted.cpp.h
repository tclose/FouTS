/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 22/11/2010.

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

#ifndef __bts_mcmc_proposal_momentum_weighted_cpp_h__
#define __bts_mcmc_proposal_momentum_weighted_cpp_h__

namespace FTS {
    
    namespace MCMC {
        
        namespace Proposal {
            
            template<typename T> Momentum::Weighted Momentum::Weighted::factory(
                    const T& state, double step_scale, const std::string& step_location,
                    MCMC::Proposal::Distribution* const proposal_distribution) {
                
                Momentum::Weighted momentum;
                
                T step;
                
                if (step_location.size()) {
                    
                    T loaded_step;
                    
                    loaded_step.load(step_location);
                    
                    if (loaded_step.size() == 1)
                        step.resize(state.size(), loaded_step[0]);
                    else if (loaded_step.size() == state.size())
                        step = loaded_step;
                    else
                        throw Exception(
                                "Size of loaded step (" + str(loaded_step.size()) + ") from file '"
                                + step_location
                                + "' is not singular and does not match that of state ("
                                + str(state.size()) + ").");
                    
                    step *= step_scale;
                    
                } else
                    step.resize(state.size(), step_scale);
                
//        step.base_intensity = 0.0;
//
//        //TODO: Remove from production version
//        for (size_t i = 0; i < step.size(); ++i)
//          step[i].acs = 0.0;
                
                MR::Math::Vector<double>& step_vector = step;
                
                return Momentum::Weighted(proposal_distribution, step_vector);
                
            }
        
        }
    
    }

}

#endif /* __bts_mcmc_proposal_momentum_weighted_cpp_h__ */
