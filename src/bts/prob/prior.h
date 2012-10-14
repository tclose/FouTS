/*
    Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

    Written by Thomas G Close, 5/05/09.

    This file is part of Bayesian Tractlet Sampling (BTS).

    BTS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your prior_option) any later version.

    BTS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BTS.  If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef __bts_prob_prior_h__
#define __bts_prob_prior_h__


//Defines the parameters that are required to initialise the prior objects.
#define PRIOR_PARAMETERS \
  Option ("prior_scale",  "Scales the weight of the whole prior distribution") \
   + Argument ("prior_scale", "").type_float(0, Prob::Prior::SCALE_DEFAULT,LARGE_FLOAT), \
\
  Option ("prior_hook", "The type of prior placed on the sheer of the secondary axes") \
   + Argument ("prior_hook_scale", "").optional().type_float(0, Prob::PriorComponent::Hook::SCALE_DEFAULT,LARGE_FLOAT) \
   + Argument ("prior_hook_num_points", "").optional().type_integer(1, Prob::PriorComponent::Hook::NUM_POINTS_DEFAULT,LARGE_INT), \
\
  Option ("prior_freq", "The prior placed on the tractlet frequency") \
   + Argument ("prior_freq_scale", "").optional().type_float(0, Prob::PriorComponent::Frequency::SCALE_DEFAULT,LARGE_FLOAT) \
   + Argument ("prior_freq_aux_scale", "").optional().type_float(0, Prob::PriorComponent::Frequency::AUX_SCALE_DEFAULT,LARGE_FLOAT), \
\
  Option ("prior_density", "The type of prior placed on the sheer of the secondary axes") \
   + Argument ("prior_density_high_scale", "").optional().type_float(0, Prob::PriorComponent::Density::HIGH_SCALE_DEFAULT,LARGE_FLOAT) \
   + Argument ("prior_density_low_scale", "").optional().type_float(0, Prob::PriorComponent::Density::LOW_SCALE_DEFAULT,LARGE_FLOAT) \
   + Argument ("prior_density_num_points", "").optional().type_integer(1, Prob::PriorComponent::Density::NUM_POINTS_DEFAULT,LARGE_INT), \
\
  Option ("prior_acs", "The prior placed on the tractlet acs") \
   + Argument ("prior_acs_scale", "").optional().type_float(0, Prob::PriorComponent::ACS::SCALE_DEFAULT,LARGE_FLOAT) \
   + Argument ("prior_acs_mean", "").optional().type_float(0, Prob::PriorComponent::ACS::MEAN_DEFAULT,LARGE_FLOAT), \
\
  Option ("prior_length", "The prior placed on the tractlet length") \
   + Argument ("prior_length_scale", "").optional().type_float(0, Prob::PriorComponent::Length::SCALE_DEFAULT,LARGE_FLOAT) \
   + Argument ("prior_length_mean", "").optional().type_float(0, Prob::PriorComponent::Length::MEAN_DEFAULT,LARGE_FLOAT), \
\
  Option ("prior_skinny", "The prior placed on the tractlet skinny") \
   + Argument ("prior_skinny_scale", "").optional().type_float(0, Prob::PriorComponent::Skinny::SCALE_DEFAULT,LARGE_FLOAT) \
   + Argument ("prior_skinny_power", "").optional().type_float(0, Prob::PriorComponent::Skinny::POWER_DEFAULT,LARGE_FLOAT) \
\

//Loads the 'prior' parameters into variables
#define SET_PRIOR_PARAMETERS \
  double prior_scale = Prob::Prior::SCALE_DEFAULT; \
 \
  Options prior_opt = get_options("prior_scale"); \
  if (prior_opt.size()) \
    prior_scale = prior_opt[0][0]; \
\
  double prior_hook_scale  = Prob::PriorComponent::Hook::SCALE_DEFAULT; \
  size_t prior_hook_num_points = Prob::PriorComponent::Hook::NUM_POINTS_DEFAULT; \
 \
  prior_opt = get_options("prior_hook"); \
  if (prior_opt.size()) { \
    if (prior_opt[0].size() >= 1) \
      prior_hook_scale = prior_opt[0][0]; \
    if (prior_opt[0].size() >= 2) \
      prior_hook_num_points = prior_opt[0][1]; \
  } \
 \
  double prior_freq_scale  = Prob::PriorComponent::Frequency::SCALE_DEFAULT; \
  double prior_freq_aux_scale  = Prob::PriorComponent::Frequency::AUX_SCALE_DEFAULT; \
 \
  prior_opt = get_options("prior_freq"); \
  if (prior_opt.size()) { \
    if (prior_opt[0].size() >= 1) \
      prior_freq_scale = prior_opt[0][0]; \
    if (prior_opt[0].size() >= 2) \
      prior_freq_aux_scale = prior_opt[0][1]; \
  } \
 \
   double prior_density_high_scale  = Prob::PriorComponent::Density::HIGH_SCALE_DEFAULT; \
   double prior_density_low_scale  = Prob::PriorComponent::Density::LOW_SCALE_DEFAULT; \
   size_t prior_density_num_points = Prob::PriorComponent::Density::NUM_POINTS_DEFAULT; \
\
  prior_opt = get_options("prior_density"); \
  if (prior_opt.size()) { \
    if (prior_opt[0].size() >= 1) \
      prior_density_high_scale = prior_opt[0][0]; \
    if (prior_opt[0].size() >= 1) \
      prior_density_low_scale = prior_opt[0][1]; \
    if (prior_opt[0].size() >= 2) \
      prior_density_num_points = prior_opt[0][2]; \
  } \
  \
  double prior_acs_scale  = Prob::PriorComponent::ACS::SCALE_DEFAULT; \
  double prior_acs_mean  = Prob::PriorComponent::ACS::MEAN_DEFAULT; \
 \
  prior_opt = get_options("prior_acs"); \
  if (prior_opt.size()) { \
    if (prior_opt[0].size() >= 1) \
      prior_acs_scale = prior_opt[0][0]; \
    if (prior_opt[0].size() >= 2) \
      prior_acs_mean = prior_opt[0][1]; \
  } \
   \
  double prior_length_scale  = Prob::PriorComponent::Length::SCALE_DEFAULT; \
  double prior_length_mean  = Prob::PriorComponent::Length::MEAN_DEFAULT; \
 \
  prior_opt = get_options("prior_length"); \
  if (prior_opt.size()) { \
    if (prior_opt[0].size() >= 1) \
      prior_length_scale = prior_opt[0][0]; \
    if (prior_opt[0].size() >= 2) \
      prior_length_mean = prior_opt[0][1]; \
  } \
 \
  double prior_skinny_scale  = Prob::PriorComponent::Skinny::SCALE_DEFAULT; \
  double prior_skinny_power  = Prob::PriorComponent::Skinny::POWER_DEFAULT; \
 \
  prior_opt = get_options("prior_skinny"); \
  if (prior_opt.size()) { \
    if (prior_opt[0].size() >= 1) \
      prior_skinny_scale = prior_opt[0][0]; \
    if (prior_opt[0].size() >= 2) \
      prior_skinny_power = prior_opt[0][1]; \
  } \
 \


//Adds the 'prior' parameters to the properties to be saved with the data.
#define ADD_PRIOR_PROPERTIES(properties) \
  if (prior_scale) { \
    properties["prior_scale"]                         = str(prior_scale); \
\
    if (prior_freq_scale) { \
      properties["prior_freq_scale"]                   = str(prior_freq_scale); \
      properties["prior_freq_aux_scale"]               = str(prior_freq_aux_scale); \
    } \
\
    if (prior_hook_scale) { \
      properties["prior_hook_scale"]                  = str(prior_hook_scale); \
      properties["prior_hook_num_points"]             = str(prior_hook_num_points); \
    } \
\
    if (prior_density_high_scale) { \
      properties["prior_density_high_scale"]          = str(prior_density_high_scale); \
      properties["prior_density_low_scale"]           = str(prior_density_low_scale); \
      properties["prior_density_num_points"]          = str(prior_density_num_points); \
    } \
\
    if (prior_acs_scale) { \
      properties["prior_acs_scale"]                   = str(prior_acs_scale); \
      properties["prior_acs_mean"]                    = str(prior_acs_mean); \
    } \
  \
    if (prior_length_scale) { \
      properties["prior_length_scale"]                = str(prior_length_scale); \
      properties["prior_length_mean"]                 = str(prior_length_mean); \
    } \
  \
    if (prior_skinny_scale) { \
      properties["prior_skinny_scale"]                = str(prior_skinny_scale); \
      properties["prior_skinny_power"]                = str(prior_skinny_power); \
    } \
  }

#include "bts/prob/prior_component/frequency.h"
#include "bts/prob/prior_component/hook.h"
#include "bts/prob/prior_component/density.h"
#include "bts/prob/prior_component/acs.h"
#include "bts/prob/prior_component/length.h"
#include "bts/prob/prior_component/skinny.h"

#include "bts/common.h"  

#include "bts/prob/prob.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"
#include "bts/mcmc/state.h"


namespace BTS {

  namespace Prob {

    class Prior {

      public:

        const static char*         TYPE_DEFAULT;
        const static double        SCALE_DEFAULT;

      protected:

        double scale;
        PriorComponent::Frequency frequency;
        PriorComponent::Hook hook;
        PriorComponent::Density density;
        PriorComponent::ACS acs;
        PriorComponent::Length length;
        PriorComponent::Skinny skinny;

      public:

        Prior(double scale,
              double freq_scale,
              double freq_aux_scale,
              double hook_scale,
              double hook_num_points,
              double density_high_scale,
              double density_low_scale,
              double density_num_points,
              double acs_scale,
              double acs_mean,
              double length_scale,
              double length_mean,
              double skinny_scale,
              size_t skinny_power);

        Prior(const Prior& p)
         : scale(p.scale), frequency(p.frequency), hook(p.hook), density(p.density), acs(p.acs), length(p.length),
                                                                                           skinny(p.skinny) {}

        Prior& operator=(const Prior& p)
          { scale =p.scale; frequency = p.frequency; hook = p.hook; density = p.density; acs = p.acs; length = p.length;
                                                                              skinny = p.skinny; return *this; }

        virtual ~Prior() {}

        std::vector<std::string>        list_components() const {
          std::vector<std::string> components;
          components.push_back(PriorComponent::Frequency::NAME);
          components.push_back(PriorComponent::Hook::NAME);
          components.push_back(PriorComponent::Density::NAME);
          components.push_back(PriorComponent::ACS::NAME);
          components.push_back(PriorComponent::Length::NAME);
          components.push_back(PriorComponent::Skinny::NAME);
          return components;
        }

        std::map<std::string, std::string> get_component_values(const Fibre::Strand strand);

        std::map<std::string, std::string> get_component_values(const Fibre::Tractlet tractlet);

        template <typename T> std::map<std::string, std::string> get_component_values(const T fibres) {
          std::map<std::string,std::string> overall_map, elem_map;
          std::vector<string> components = list_components();
          for (std::vector<std::string>::iterator comp_it = components.begin(); comp_it != components.end(); ++comp_it)
            overall_map[*comp_it] = 0.0;
          for (size_t fibre_i = 0; fibre_i < fibres.size(); ++fibre_i) {
            elem_map = get_component_values(fibres[fibre_i]);
            for (std::map<std::string,std::string>::iterator comp_it = elem_map.begin(); comp_it != elem_map.end(); ++comp_it)
              overall_map[comp_it->first] += comp_it->second;

          }
          return overall_map;
        }

        double                          log_prob(const Fibre::Strand strand, Fibre::Strand gradient);

        double                          log_prob(const Fibre::Tractlet tractlet, Fibre::Tractlet gradient);

        template <typename T> double    log_prob(const T fibres, T gradient) {

          double lprob = 0.0;

          if (scale)
            for (size_t fibre_i = 0; fibre_i < fibres.size(); ++fibre_i)
              lprob += log_prob(fibres[fibre_i], gradient[fibre_i]);

          return lprob;
        }

        template <typename T> double    log_prob(const T& fibre) {

          T dummy_gradient;
          dummy_gradient = fibre;

          return log_prob(fibre, dummy_gradient);

        }

        template <typename T> double    log_prob(const T strand, T gradient, typename T::Tensor hessian)
          { throw Exception ("Not implemented yet."); }


        template <typename T> double    log_prob_and_fisher(const T strand, T gradient, typename T::Tensor fisher)
          { throw Exception ("Not implemented yet."); }


        template <typename T> double    log_prob_and_fisher(const T strand,
                                                            T gradient, typename T::Tensor fisher,
                                                            std::vector<typename T::Tensor> fisher_gradient)
          { throw Exception ("Not implemented yet."); }

    };


  }

}




#endif
