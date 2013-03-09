/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 26/07/2010.

 This file is part of MRtrix.

 MRtrix is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your walk_option) any later version.

 MRtrix is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef __bts_mcmc_proposal_walker_h__
#define __bts_mcmc_proposal_walker_h__

namespace BTS {
    
    namespace MCMC {
        
        namespace Proposal {
            
            class Walker;
        
        }
    
    }

}

#define PROPOSAL_WALKER_PARAMETERS \
  Option ("walk_type", "The 'Proposal Walker' object used to propose new fibres in MCMC sampling.") \
   + Argument ("walk_type", "").type_text(), \
\
  Option ("walk_step_scale", "Overall scaling, every scaling is multiplied by this value.") \
   + Argument ("walk_step_scale", "").type_float (0.0, 1e5, NAN), \
\
  Option ("walk_step_location", "The location of the relative scalings for each of the parameters.") \
   + Argument ("walk_step_location", "").type_text(), \
\
  Option ("walk_intens_scale", "Intensity step size scaling.") \
   + Argument ("walk_intens_scale", "").type_float (0.0, 1e5, NAN)

//Loads the 'proposal' parameters into variables
#define SET_PROPOSAL_WALKER_PARAMETERS(state_location) \
  std::string walk_type; \
  double walk_step_scale; \
  std::string walk_step_location; \
  double walk_base_intens_scale; \
  \
  if (File::has_extension<Fibre::Strand>(state_location)) { \
   walk_type                     = Fibre::Strand::Set::Walker::TYPE_DEFAULT; \
   walk_step_scale               = Fibre::Strand::Set::Walker::STEP_SCALE_DEFAULT; \
   walk_base_intens_scale        = 0.01; \
   walk_step_location            = Fibre::Strand::Set::Walker::STEP_LOCATION_DEFAULT; \
  \
  } else if (File::has_extension<Fibre::Tractlet>(state_location)) { \
   walk_type                     = Fibre::Tractlet::Set::Walker::TYPE_DEFAULT; \
   walk_step_scale               = Fibre::Tractlet::Set::Walker::STEP_SCALE_DEFAULT; \
   walk_base_intens_scale        = 0.01; \
   walk_step_location            = Fibre::Tractlet::Set::Walker::STEP_LOCATION_DEFAULT; \
  \
  } else { \
   walk_type                     = MCMC::Proposal::Walker::TYPE_DEFAULT; \
   walk_step_scale               = MCMC::Proposal::Walker::STEP_SCALE_DEFAULT; \
   walk_step_location            = MCMC::Proposal::Walker::STEP_LOCATION_DEFAULT; \
  } \
  \
  Options walk_opt = get_options("walk_type"); \
  if (walk_opt.size()) \
   walk_type = walk_opt[0][0].c_str(); \
  \
  walk_opt = get_options("walk_step_scale"); \
  if (walk_opt.size()) \
   walk_step_scale = walk_opt[0][0]; \
  \
  walk_opt = get_options("walk_step_location"); \
  if (walk_opt.size()) \
   walk_step_location = walk_opt[0][0].c_str(); \
  \
  walk_opt = get_options("walk_base_intens_scale"); \
  if (walk_opt.size()) \
   walk_base_intens_scale = walk_opt[0][0];

//Adds the 'proposal' parameters to the properties to be saved with the data.
#define ADD_PROPOSAL_WALKER_PROPERTIES(properties) \
  properties["walk_type"]                       = walk_type; \
  properties["walk_step_scale"]                 = str(walk_step_scale); \
  properties["walk_base_intens_scale"]          = str(walk_base_intens_scale); \
\
  if (walk_type == "standard") { \
    properties["walk_step"]                     = Fibre::Base::Object::load_matlab_str(walk_step_location, walk_step_scale); \
    properties["walk_step_location"]            = walk_step_location; \
\
  }

#include "bts/mcmc/proposal/distribution.h"
#include "bts/mcmc/state.h"

namespace BTS {
    
    namespace MCMC {
        
        namespace Proposal {
            
            class Walker {
                    
                public:
                    
                    const static std::string TYPE_DEFAULT;
                    const static double STEP_SCALE_DEFAULT;
                    const static std::string STEP_LOCATION_DEFAULT;

                    //Member variables
                protected:
                    
                    MR::Math::Vector<double> step_sizes;

                    Distribution* prop_distr;

                    //Member functions
                public:
                    
                    Walker(Distribution* const prop_distr = 0)
                            : prop_distr(prop_distr) {
                    }
                    
                    Walker(Distribution* const proposal_distribution,
                           const MR::Math::Vector<double>& relative_step_sizes)
                            : prop_distr(0)

                    {
                        set(proposal_distribution, relative_step_sizes);
                    }
                    
                    ~Walker() {
                        if (prop_distr)
                            delete prop_distr;
                    }
                    
                    Walker(const Walker& m)
                            : step_sizes(m.step_sizes) {
                        if (m.prop_distr)
                            prop_distr = m.prop_distr->clone();
                    }
                    
                    Walker& operator=(const Walker& m) {
                        step_sizes = m.step_sizes;
                        if (m.prop_distr)
                            prop_distr = m.prop_distr->clone();
                        return *this;
                    }
                    
                    void step(const MR::Math::Vector<double>& current,
                              MR::Math::Vector<double>& proposed, double scalar = 1.0) const;

                    void set(Distribution* const proposal_distribution,
                             MR::Math::Vector<double> relative_step_sizes);
                    
            };
        
        }
    
    }

}

#endif /* __bts_mcmc_proposal_walker_h__ */
