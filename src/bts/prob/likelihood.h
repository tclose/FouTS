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

#ifndef __bts_prob_likelihood_h__
#define __bts_prob_likelihood_h__

namespace BTS {
    
    namespace Prob {
        
        class Likelihood;
    
    }

}

#define LIKELIHOOD_PARAMETERS \
  Option ("like_snr", "The assumed Signal-to-noise ratio (SNR) in the likelihood model.") \
   + Argument ("like_snr", "").type_float(0, Prob::Likelihood::SNR_DEFAULT,LARGE_FLOAT), \
\
  Option ("like_type", "The 'Likelihood' object used to calculate the likelihood probability.") \
   + Argument ("like_type", "").type_text(Prob::Likelihood::TYPE_DEFAULT), \
\
  Option ("like_b0_include", "Include b values in comparison between observed and expected images.") \
   + Argument ("like_b0_include", "").type_text(Prob::Likelihood::B0_INCLUDE_DEFAULT), \
\
  Option ("like_outside_scale", "The scale of the outside signal contribution to the probability.") \
   + Argument ("like_outside_scale", "").type_float(0.0, Prob::Likelihood::OUTSIDE_SCALE_DEFAULT,LARGE_FLOAT), \
\
  Option ("like_ref_b0", "Speicify whether to use the 'max' or 'average' b0 signal in the observed image for the calculation of the noise variance from the SNR.") \
   + Argument ("like_ref_b0", "").type_text(Prob::Likelihood::REF_B0_DEFAULT), \
\
  Option ("like_ref_signal", "Manually speicify b0 signal for the calculation of the noise variance from the SNR.") \
   + Argument ("like_ref_signal", "").type_float(SMALL_FLOAT, NAN, LARGE_FLOAT), \
\
 Option ("like_noise_map", "An image containing an explicit map of the estimate noise at each voxel.") \
  + Argument ("like_noise", "").type_image_in() \


\
//Loads the 'prior' parameters into variables
#define SET_LIKELIHOOD_PARAMETERS \
  double like_snr             = Prob::Likelihood::SNR_DEFAULT; \
  std::string like_type       = Prob::Likelihood::TYPE_DEFAULT; \
  std::string  like_b0_include = Prob::Likelihood::B0_INCLUDE_DEFAULT; \
  double like_outside_scale   = Prob::Likelihood::OUTSIDE_SCALE_DEFAULT; \
  std::string like_ref_b0     = Prob::Likelihood::REF_B0_DEFAULT; \
  double like_ref_signal      = NAN; \
  std::string like_noise_map_name = ""; \
  Image::Double::Buffer like_noise_map; \
\
  Options like_opt = get_options("like_snr"); \
  if (like_opt.size()) \
    like_snr = like_opt[0][0]; \
\
  like_opt = get_options("like_type"); \
  if (like_opt.size()) \
    like_type = like_opt[0][0].c_str(); \
\
  like_opt = get_options("like_b0_include"); \
  if (like_opt.size()) \
    like_b0_include = like_opt[0][0].c_str(); \
\
  like_opt = get_options("like_outside_scale"); \
  if (like_opt.size()) \
    like_outside_scale = like_opt[0][0]; \
\
  like_opt = get_options("like_ref_b0"); \
  if (like_opt.size()) \
    like_ref_b0 = like_opt[0][0].c_str(); \
\
  like_opt = get_options("like_ref_signal"); \
  if (like_opt.size()) \
    like_ref_signal = like_opt[0][0]; \
\
  like_opt = get_options("like_noise_map"); \
  if (like_opt.size()) { \
    like_noise_map_name = like_opt[0][0].c_str(); \
    like_noise_map.load(like_noise_map_name); \
  }

//Adds the 'prior' parameters to the properties to be saved with the data.
#define ADD_LIKELIHOOD_PROPERTIES(properties) \
  properties["like_type"]      = like_type; \
  if (like_type != "uniform") {\
    properties["like_snr"]            = str(like_snr); \
    properties["like_b0_include"]     = like_b0_include; \
    properties["like_outside_scale"]  = str(like_outside_scale); \
    if (!isnan(like_ref_signal)) \
      properties["like_ref_signal"]   = str(like_ref_signal); \
    else \
      properties["like_ref_b0"]       = like_ref_b0; \
    if (!like_noise_map) \
        properties["like_noise_map"] = like_noise_map_name; \
  } \

#include "bts/image/expected/buffer.h"
#include "bts/image/expected/trilinear/buffer.h"
#include "bts/image/observed/buffer.h"

#include "bts/prob/prob.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"
#include "bts/mcmc/state.h"

namespace BTS {
    
    namespace Prob {
        
        class Likelihood {
                
                //Public nested classes.
            public:
                
                class Gaussian;
                class OneSidedGaussian;
                class Rician;

                //Public static constants
            public:
                
                const static char* TYPE_DEFAULT;
                const static double SNR_DEFAULT;
                const static char* B0_INCLUDE_DEFAULT;
                const static double OUTSIDE_SCALE_DEFAULT;
                const static char* REF_B0_DEFAULT;

                const static Image::Observed::Buffer dummy_observed_buffer;
                const static Image::Expected::Trilinear::Buffer dummy_expected_buffer;

            public:
                
                static Likelihood* factory(const std::string& type,
                                           const Image::Observed::Buffer& obs_image,
                                           Image::Expected::Buffer* exp_image, double assumed_snr,
                                           const std::string& b0_include, double outside_scale,
                                           const std::string& ref_b0, double ref_signal,
                                           const Image::Double::Buffer& noise_map);

            protected:
                
                double sigma2;

                Image::Observed::Buffer obs_image;
                Image::Expected::Buffer* exp_image;

                Image::Container::Buffer<Fibre::Strand>::Set recycled_strand_gradients;
                Image::Container::Buffer<Fibre::Tractlet>::Set recycled_tractlet_gradients;
                Image::Container::Buffer<Fibre::Strand::Tensor>::Set recycled_strand_hessians;
                Image::Container::Buffer<Fibre::Tractlet::Tensor>::Set recycled_tractlet_hessians;

                size_t difference_mode;

                std::string b0_include;

            public:
                
                Likelihood(const Image::Observed::Buffer& observed_image,
                           Image::Expected::Buffer* const expected_image, double assumed_snr,
                           const std::string& b0_include, double outside_scale,
                           const std::string& ref_b0, double ref_signal);

                Likelihood(const Likelihood& l);

                virtual ~Likelihood() {
                    delete this->exp_image;
                }
                
                Likelihood& operator=(const Likelihood& l);

                virtual double log_prob(Image::Expected::Buffer& image);

                double log_prob(const Fibre::Strand& strand) {
                    Fibre::Strand::Set set;
                    set.push_back(strand);
                    return log_prob(set);
                }
                
                double log_prob(const Fibre::Tractlet& tractlet) {
                    Fibre::Tractlet::Set set;
                    set.push_back(tractlet);
                    return log_prob(set);
                }
                
                double log_prob(const Fibre::Strand& strand, Fibre::Strand& gradient) {
                    Fibre::Strand::Set set, set_gradient;
                    set.push_back(strand);
                    set_gradient = set;
                    double lprob = log_prob(set, set_gradient);
                    gradient = set_gradient[0];
                    return lprob;
                }
                
                double log_prob(const Fibre::Tractlet& tractlet, Fibre::Tractlet& gradient) {
                    Fibre::Tractlet::Set set, set_gradient;
                    set.push_back(tractlet);
                    set_gradient = set;
                    double lprob = log_prob(set, set_gradient);
                    gradient = set_gradient[0];
                    return lprob;
                }
                
                virtual double log_prob(const Fibre::Strand::Set& strands) {
                    return log_prob_tpl<Fibre::Strand>(strands);
                }
                
                virtual double log_prob(const Fibre::Tractlet::Set& tractlets) {
                    return log_prob_tpl<Fibre::Tractlet>(tractlets);
                }
                
                virtual double log_prob(const Fibre::Strand::Set& strands,
                                        Fibre::Strand::Set& gradient) {
                    return log_prob_tpl<Fibre::Strand>(strands, gradient);
                }
                
                virtual double log_prob(const Fibre::Tractlet::Set& tractlets,
                                        Fibre::Tractlet::Set& gradient) {
                    return log_prob_tpl<Fibre::Tractlet>(tractlets, gradient);
                }
                
                virtual double log_prob(const Fibre::Strand::Set& strands,
                                        Fibre::Strand::Set& gradient,
                                        Fibre::Strand::Set::Tensor& hessian) {
                    return log_prob_tpl<Fibre::Strand>(strands, gradient, hessian);
                }
                
                virtual double log_prob(const Fibre::Tractlet::Set& tractlets,
                                        Fibre::Tractlet::Set& gradient,
                                        Fibre::Tractlet::Set::Tensor& hessian) {
                    return log_prob_tpl<Fibre::Tractlet>(tractlets, gradient, hessian);
                }
                
                virtual double log_prob_and_fisher(const Fibre::Strand::Set& strands,
                                                   Fibre::Strand::Set& gradient,
                                                   Fibre::Strand::Set::Tensor& fisher_info) {
                    throw Exception("should be implemented in derrived class.");
                }
                
                virtual double log_prob_and_fisher(const Fibre::Tractlet::Set& strands,
                                                   Fibre::Tractlet::Set& gradient,
                                                   Fibre::Tractlet::Set::Tensor& fisher_info) {
                    throw Exception("should be implemented in derrived class.");
                }
                
                template<typename T> double log_prob_tpl(const typename T::Set& fibres);

                template<typename T> double log_prob_tpl(const typename T::Set& fibres,
                                                         typename T::Set& gradient);

                template<typename T> double log_prob_tpl(const typename T::Set& fibres,
                                                         typename T::Set& gradient,
                                                         typename T::Set::Tensor& hessian);

                virtual double log_prob(double expected, double observed) = 0;

                virtual double b0_log_prob(double expected, double observed) = 0;

                virtual double log_prob(double expected, double observed, double& d_lprob) = 0;

                virtual double log_prob(double expected, double observed, double& d_lprob,
                                        double& d2_lprob2) = 0;

                void save_expected_image(const std::string& location)    //Used for debugging
                        {
                    exp_image->save(location);
                }
                
                const Image::Expected::Buffer& get_expected_image()    //Used for debugging
                {
                    return *exp_image;
                }
                
                const Image::Observed::Buffer& get_observed_image()    //Used for debugging
                {
                    return obs_image;
                }
                
                void set_enforce_bounds(bool flag);

                void set_assumed_snr(double snr, const std::string& ref_b0, double ref_signal);

            protected:
                
                //! Used to get the right used gradients for the templated type.
                template<typename T> typename Image::Container::Buffer<T>::Set& get_recycled_gradients();

                //! Used to get the right used hessians for the templated type.
                template<typename T> typename Image::Container::Buffer<typename T::Tensor>::Set& get_recycled_hessians();

                Likelihood() {
                }
                
                friend std::ostream& operator<<(std::ostream& stream, const Likelihood& likelihood);
                
        };
    
    }

}

#endif
