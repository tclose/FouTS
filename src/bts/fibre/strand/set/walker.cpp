/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Aug 2, 2010.

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

#include "bts/fibre/strand/set/walker.h"

namespace BTS {
    
    namespace Fibre {
        
        const char* Strand::Set::Walker::TYPE_DEFAULT = "standard";
        const char* Strand::Set::Walker::STEP_LOCATION_DEFAULT =
                "/home/tclose/data/mcmc/metropolis/params/walker/scaling.str.txt";
        const double Strand::Set::Walker::STEP_SCALE_DEFAULT = 0.01;
        
        Strand::Set::Walker* Strand::Set::Walker::factory(
                const Strand::Set& state, const std::string& type, double step_scale,
                const std::string& step_location, MCMC::Proposal::Distribution* distribution,
                double base_intensity_scale) {
            
            Walker* walker;
            
            if (type == "standard") {
                Fibre::Strand::Set step;
                
                if (step_location.size()) {
                    step.load(step_location);
                    
                    if (!step.props_match(state))
                        throw Exception(
                                "Properties do not match between loaded properties and state properties.");
                    
                    if (step.size() == 1) {
                        
                        Fibre::Strand step_template = step[0];
                        
                        step.resize(state.size());
                        
                        for (size_t tract_i = 0; tract_i < state.size(); ++tract_i) {
                            
                            if (state[tract_i].size() > step_template.size())
                                throw Exception(
                                        "Strand (" + str(tract_i)
                                        + ") is larger than step template.");
                            
                            step[tract_i] = step_template;
                            step[tract_i].resize(state[tract_i].size());
                            
                        }
                        
                    } else if (step.size() != state.size())
                        throw Exception(
                                "Loaded step sizes' size (" + str(step.size())
                                + ") does not match state (" + str(state.size()) + ").");
                    
                    step *= step_scale;
                    
                } else {
                    step = state;
                    step.set(step_scale);
                }
                
                walker = new Walker(distribution, step);
                
            } else
                
                throw Exception("Unrecognised Fibre::Strand::Set::Walker type '" + type + "'.");
            
            return walker;
            
        }
        
        void Strand::Set::Walker::step(const Strand::Set& current, Strand::Set& proposed,
                                       double scalar) {
            
            MCMC::Proposal::Walker::step(current, proposed, scalar);
            
        }
    
    }

}
