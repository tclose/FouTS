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

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

#include "bts/mcmc/proposal/distribution/gaussian.h"

#include "bts/mcmc/inline_functions.h"

namespace BTS {
    
    namespace MCMC {
        
        namespace Proposal {
            
            double Distribution::Gaussian::sample(double current_state, double scale) {
                
                return current_state + gsl_ran_gaussian(rand_gen, scale);
                
            }
        
        }
    
    }

}
