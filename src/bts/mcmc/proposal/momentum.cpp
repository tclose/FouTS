/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Jul 20, 2010.

 This file is part of MRtrix.

 MRtrix is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version size() of the License, or
 (at your option) any later version.

 MRtrix is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef __bts_proposal_momentum_cpp_h__
#define __bts_proposal_momentum_cpp_h__

#include "bts/mcmc/proposal/momentum.h"
#include "math/vector.h"

namespace BTS {
    
    namespace MCMC {
        
        namespace Proposal {
            
            const double Momentum::STEP_SCALE_DEFAULT = 0.05;
            const std::string Momentum::STEP_LOCATION_DEFAULT =
                    "/home/tclose/data/mcmc/hamiltonian/params/momentum/scaling.sta";
            
            void Momentum::set(Distribution* const proposal_distribution,
                               MR::Math::Vector<double> relative_step) {
                
                momen = relative_step;
                momen *= NAN;
                step = relative_step;
                
                if (prop_distr)
                    delete prop_distr;
                
                this->prop_distr = proposal_distribution->clone();
                
            }
            
            void Momentum::randomize() {
                
                for (size_t elem_i = 0; elem_i < size(); elem_i++)
                    momen[elem_i] = prop_distr->sample(0.0, 1.0);
                
            }
            
            void Momentum::half_update_momentum(const MR::Math::Vector<double>& gradient,
                                                double time_direction) {
                
                assert((time_direction == -1.0) || (time_direction == 1.0));
                
                for (size_t elem_i = 0; elem_i < size(); elem_i++)
                    momen[elem_i] += time_direction * gradient[elem_i] * step[elem_i] / 2.0;
                
            }
            
            void Momentum::update_state(MR::Math::Vector<double>& state,
                                        double time_direction) const {
                
                assert((time_direction == -1.0) || (time_direction == 1.0));
                
                for (size_t elem_i = 0; elem_i < size(); elem_i++)
                    state[elem_i] += time_direction * momen[elem_i] * step[elem_i];
                
            }
            
            double Momentum::predicted_change(const MR::Math::Vector<double>& gradient,
                                              double time_direction) const {
                
                assert((time_direction == -1.0) || (time_direction == 1.0));
                
                double pred_change = 0.0;
                
                for (size_t elem_i = 0; elem_i < size(); elem_i++)
                    pred_change += time_direction * momen[elem_i] * gradient[elem_i] * step[elem_i];
                
                return pred_change;
                
            }
            
            std::ostream& operator<<(std::ostream& stream, const Momentum& momen) {
                
                stream << "Momentum: " << momen.momen << std::endl;
                stream << "Step sizes: " << momen.step << std::endl;
                
                return stream;
                
            }
        
        }
    
    }

}

#endif /* __bts_proposal_momentum_cpp_h__ */
