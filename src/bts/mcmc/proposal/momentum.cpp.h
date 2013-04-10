/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 15/09/2010.

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

#ifndef __bts_mcmc_proposal_momentum_cpp_h__
#define __bts_mcmc_proposal_momentum_cpp_h__

namespace BTS {
    
    namespace MCMC {
        
        namespace Proposal {
            
            template<typename T> Momentum Momentum::factory(
                    const T& state, double step_scale, const std::string& step_location,
                    const std::string& weight_location,
                    MCMC::Proposal::Distribution* const proposal_distribution) {
                
                Momentum momentum;
                
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
                
                MR::Math::Vector<double>& step_vector = step;
                
                //TODO: Reinstate this and create fixed weighted class.
//        if (weight_location.size()) {
//
//          size_t state_vector_size = state.vsize();
//          MR::Math::Matrix<double> weights(state_vector_size, state_vector_size);
//          weights = 0.0;
//
//          MR::Math::Matrix<double> loaded_weights (weight_location);
//
//          size_t loaded_weights_size = loaded_weights.rows();
//
//          if (loaded_weights_size != loaded_weights.columns())
//            throw Exception ("Loaded matrix is not square (" + str(loaded_weights_size) + " x " + str(loaded_weights.columns()) + ").");
//
//          if (loaded_weights_size == state_vector_size)
//            weights = loaded_weights;
//
//          else if (!(state_vector_size % loaded_weights_size)) {
//
//            for (size_t block_i = 0; block_i < state_vector_size; block_i += loaded_weights_size)
//              weights.sub(block_i, block_i+loaded_weights_size, block_i, block_i+loaded_weights_size) = loaded_weights;
//
//          } else
//            throw Exception ("Loaded weights size (" + str(loaded_weights_size) + ") does not divide evenly into state vector size (" + str(state_vector_size) + ").");
//
//
//          momentum = Momentum::Weighted(proposal_distribution, step_vector, weights);
//        } else
                momentum = Momentum(proposal_distribution, step_vector);
                
                VAR(momentum);
                
                return momentum;
                
            }
        
        }
    
    }

}

#endif /* __bts_mcmc_proposal_momentum_cpp_h__ */
