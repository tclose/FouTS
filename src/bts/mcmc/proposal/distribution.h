/*
 Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

 Written by Thomas G Close, 5/05/09.

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

#ifndef __bts_mcmc_proposal_distribution_h__
#define __bts_mcmc_proposal_distribution_h__

#define PROPOSAL_DISTRIBUTION_PARAMETERS \
  Option ("prop_distr_type", "The 'Proposal::Distribution' object used to calculate the prior probability.") \
   + Argument ("prop_distr_type", "").type_text(MCMC::Proposal::Distribution::TYPE_DEFAULT) \
\


//Loads the 'proposal' parameters into variables
#define SET_PROPOSAL_DISTRIBUTION_PARAMETERS \
  std::string prop_distr_type         = MCMC::Proposal::Distribution::TYPE_DEFAULT; \
\
  Options prop_opt = get_options("prop_distr_type"); \
  if (prop_opt.size()) \
    prop_distr_type = prop_opt[0][0].c_str(); \
\


//Adds the 'proposal distribution' parameters to the properties to be saved with the data.
#define ADD_PROPOSAL_DISTRIBUTION_PROPERTIES(properties) \
  properties["prop_distr_type"]               = prop_distr_type; \


#include <vector>
#include <string>

namespace BTS {
    
    namespace MCMC {
        
        namespace Proposal {
            
            class Distribution;
        
        }
    
    }

}

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

namespace BTS {
    
    namespace MCMC {
        
        namespace Proposal {
            
            class Distribution {
                    
                public:
                    
                    class Gaussian;
                    class HeavyTailed;

                public:
                    
                    const static double STEP_SIZE_DEFAULT;
                    const static double SHEER_SCALE_DEFAULT;
                    const static double SPREAD_SCALE_DEFAULT;
                    const static double BASEWIDTH_SCALE_DEFAULT;
                    const static double ENDPOINT_SCALE_DEFAULT;
                    const static double MIDPOINT_SCALE_DEFAULT;
                    const static double INTENSITY_STEP_SIZE_DEFAULT;
                    const static double BASE_INTENSITY_STEP_SIZE_DEFAULT;

//          const static std::vector<double> DEGREE_DEFAULTS;
//          const static std::vector<double> SECONDARY_DEGREE_DEFAULTS;
                    
//Public static constants
                public:
                    
                    const static char* TYPE_DEFAULT;

                protected:
                    
                    gsl_rng* rand_gen;

                    //Public static methods
                public:
                    
                    static Distribution* factory(const std::string& type, gsl_rng* rand_gen);

                public:
                    
                    virtual double sample(double current_state, double scale) = 0;

                    virtual Distribution* clone() const = 0;

                    virtual ~Distribution() {
                    }
                    
                protected:
                    
                    Distribution(gsl_rng* rand_gen)
                            : rand_gen(rand_gen) {
                    }
                    
                    Distribution(const Distribution& g)
                            : rand_gen(g.rand_gen) {
                    }
                    
                    Distribution& operator=(const Distribution& g) {
                        rand_gen = g.rand_gen;
                        return *this;
                    }
                    
            };
        
        }
    
    }

}

#endif
