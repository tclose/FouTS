/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Sep 28, 2010.

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

#ifndef __bts_prob_posterior_h__
#define __bts_prob_posterior_h__

namespace FTS {
    
    namespace Prob {
        
        template<typename State> class Posterior {
                
                //Public static variables, nested classes and typedefs
            public:
                
                //Protected member variables
            protected:
                
                Prior* prior;
                Likelihood* likelihood;

                //Public member functions
            public:
                
                Posterior(const State& state, Prior* prior, Likelihood* likelihood)
                        : prior(prior->clone()), likelihood(likelihood->clone()) {
                }
                
                //Protected member functions
            protected:
                
        };
    
    }

}

#endif /* __bts_prob_posterior_h__ */
