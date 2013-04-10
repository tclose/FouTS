/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Jun 28, 2010.

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

#ifndef __bts_mcmc_proposal_momentum_h__
#define __bts_mcmc_proposal_momentum_h__

namespace BTS {
    
    namespace MCMC {
        
        namespace Proposal {
            
            class Momentum;
        
        }
    
    }

}

#include "math/vector.h"
#include "bts/common.h"

#define PROPOSAL_MOMENTUM_PARAMETERS \
  Option ("momen_step_scale", "The overall scale of the momen step sizes.") \
   + Argument ("momen_step_scale", "").type_float(SMALL_FLOAT, MCMC::Proposal::Momentum::STEP_SCALE_DEFAULT, LARGE_FLOAT), \
\
  Option ("momen_step_location", "The relative scales of the momen step sizes.") \
   + Argument ("momen_step_location", "").type_text(), \
  \
    Option ("momen_weight_location", "The weight matrix used for the update of the momen.") \
     + Argument ("momen_weight_location", "").type_text() \

//Loads the 'proposal' parameters into variables
#define SET_PROPOSAL_MOMENTUM_PARAMETERS(state_location) \
  double momen_step_scale; \
  std::string momen_step_location; \
  std::string momen_weight_location; \
  if (File::has_extension<Fibre::Strand>(state_location)) { \
    momen_step_scale    = Fibre::Strand::Set::Momentum::STEP_SCALE_DEFAULT; \
    momen_step_location = Fibre::Strand::Set::Momentum::STEP_LOCATION_DEFAULT; \
  } else if (File::has_extension<Fibre::Tractlet>(state_location)) { \
    momen_step_scale    = Fibre::Tractlet::Set::Momentum::STEP_SCALE_DEFAULT; \
    momen_step_location = Fibre::Tractlet::Set::Momentum::STEP_LOCATION_DEFAULT; \
  } else { \
    momen_step_scale    = MCMC::Proposal::Momentum::STEP_SCALE_DEFAULT; \
    momen_step_location = MCMC::Proposal::Momentum::STEP_LOCATION_DEFAULT; \
  } \
\
  opt = get_options("momen_step_scale"); \
  if (opt.size()) \
    momen_step_scale = opt[0][0]; \
\
  opt = get_options("momen_step_location"); \
  if (opt.size()) \
    momen_step_location = opt[0][0].c_str(); \
\
  opt = get_options("momen_weight_location"); \
  if (opt.size()) \
    momen_weight_location = opt[0][0].c_str(); \

//Adds the 'prior' parameters to the momerties to be saved with the data.
#define ADD_PROPOSAL_MOMENTUM_PROPERTIES(properties) \
\
  properties["momen_step_scale"]      = str(momen_step_scale); \
  properties["momen_step_location"]   = momen_step_location; \
  properties["momen_step"]            = Fibre::Base::Object::load_matlab_str(momen_step_location, momen_step_scale); \
  if (momen_weight_location.size()) { \
    properties["momen_weight_location"] = momen_weight_location; \
    properties["momen_weight"]          = Fibre::Base::Object::load_matlab_str(momen_weight_location,1.0); \
  } \

#include "bts/mcmc/proposal/distribution.h"
#include "bts/mcmc/state.h"

namespace BTS {
    
    namespace MCMC {
        
        namespace Proposal {
            
            class Momentum {
                    
                    //Nested Classes and Typedefs
                public:
                    
                    class Weighted;

                    //Public static constants
                public:
                    
                    static const double STEP_SCALE_DEFAULT;
                    static const std::string STEP_LOCATION_DEFAULT;

                    //Public static functions.
                public:
                    
                    template<typename T> static Momentum factory(
                            const T& state, double step_scale, const std::string& step_location,
                            const std::string& weight_location,
                            MCMC::Proposal::Distribution* const proposal_distribution);

                    //Member variables
                protected:
                    
                    MR::Math::Vector<double> momen;
                    MR::Math::Vector<double> step;

                    Distribution* prop_distr;

                    //Member functions
                public:
                    
                    Momentum(Distribution* const prop_distr = 0)
                            : prop_distr(prop_distr) {
                    }
                    
                    Momentum(Distribution* const proposal_distribution,
                             const MR::Math::Vector<double>& relative_step)
                            : prop_distr(0) {
                        set(proposal_distribution, relative_step);
                    }
                    
                    ~Momentum() {
                        if (prop_distr)
                            delete prop_distr;
                    }
                    
                    Momentum(const Momentum& m)
                            : momen(m.momen), step(m.step) {
                        if (m.prop_distr)
                            delete prop_distr;
                        prop_distr = m.prop_distr->clone();
                    }
                    
                    Momentum& operator=(const Momentum& m) {
                        momen = m.momen, step = m.step;
                        if (m.prop_distr)
                            delete prop_distr;
                        prop_distr = m.prop_distr->clone();
                        return *this;
                    }
                    
                    void set(Distribution* const proposal_distribution,
                             MR::Math::Vector<double> relative_step);

                    size_t size() const {
                        return momen.size();
                    }
                    
                    double log_kinetic_energy() const {
                        return norm2(momen) / 2.0;
                    }
                    
                    const MR::Math::Vector<double>& momentum() const {
                        return momen;
                    }
                    
                    const MR::Math::Vector<double>& step_sizes() const {
                        return step;
                    }
                    
                    void randomize();

                    virtual void half_update_momentum(const MR::Math::Vector<double>& gradient,
                                                      double time_direction = 1.0);

                    virtual void update_state(MR::Math::Vector<double>& state,
                                              double time_direction = 1.0) const;

                    template<typename T> void half_update_momentum(const T& gradient,
                                                                   double time_direction = 1.0) {
                        half_update_momentum(gradient, time_direction);
                    }
                    
                    template<typename T> void update_state(T& state,
                                                           double time_direction = 1.0) const {
                        MR::Math::Vector<double>& state_vector = state;
                        update_state(state_vector, time_direction);
                    }
                    
                    void half_update_momentum(const MCMC::State& gradient, double time_direction =
                            1.0) {
                        half_update_momentum((const MR::Math::Vector<double>&) gradient,
                                time_direction);
                    }
                    
                    void update_state(MCMC::State& state, double time_direction = 1.0) const {
                        update_state((MR::Math::Vector<double>&) state, time_direction);
                    }
                    
                    //For debugging purposes.
                    double predicted_change(const MR::Math::Vector<double>& gradient,
                                            double time_direction = 1.0) const;

                    template<typename T> double predicted_change(
                            const T& gradient, double time_direction = 1.0) const {
                        return predicted_change(gradient, time_direction);
                    }
                    
                    friend std::ostream& operator<<(std::ostream& stream, const Momentum& momen);
                    
            };
            
            std::ostream& operator<<(std::ostream& stream, const Momentum& momen);
        
        }
    }
}

#endif /* __bts_proposal_momentum_h__ */
