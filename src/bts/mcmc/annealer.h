/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 24/01/2011.

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

#ifndef __bts_mcmc_annealer_h__
#define __bts_mcmc_annealer_h__

#include "bts/common.h"

namespace BTS {
    
    namespace MCMC {
        
        class Annealer {
                
                //Public static variables, nested classes and typedefs
            public:
                
                const static double START_FRACTION_DEFAULT;

                //Protected member variables
            protected:
                
                double t;
                double t_inc;

                //Public member functions
            public:
                
                Annealer(size_t num_steps, double start_fraction);

                double factor() const;

                void increment();
                
        };
    
    }

}

#endif /* __bts_mcmc_annealer_h__ */
