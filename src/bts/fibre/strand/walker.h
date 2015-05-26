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

#ifndef __bts_fibre_strand_walker_h__
#define __bts_fibre_strand_walker_h__

#include "bts/fibre/strand.h"
#include "bts/mcmc/proposal/walker.h"
#include "bts/mcmc/proposal/distribution.h"

namespace FTS {
    
    namespace Fibre {
        
        class Strand::Walker: public MCMC::Proposal::Walker {
                
                //Public member functions
            public:
                
                Walker(MCMC::Proposal::Distribution* const prop_distr)
                        : MCMC::Proposal::Walker(prop_distr) {
                }
                
                Walker(MCMC::Proposal::Distribution* const proposal_distribution,
                       const Strand& relative_step_sizes)
                        
                        : MCMC::Proposal::Walker(proposal_distribution, relative_step_sizes) {
                }
                
                virtual ~Walker() {
                }
                
        };
    
    }

}

#endif /* __bts_fibre_strand_walker_h__ */
