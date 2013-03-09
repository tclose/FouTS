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

#ifndef __bts_prob_likelihood_imagediff_gaussian_h__
#define __bts_prob_likelihood_imagediff_gaussian_h__

#include "bts/prob/prob.h"
#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood.h"

namespace BTS {
    
    namespace Prob {
        
        class Likelihood::Gaussian: public Likelihood {
                
            public:
                
                Gaussian(const Image::Observed::Buffer& obs_image,
                         Image::Expected::Buffer* const exp_image, double assumed_snr =
                                 Likelihood::SNR_DEFAULT,
                         const std::string& b0_include = Likelihood::B0_INCLUDE_DEFAULT,
                         double outside_scale = Likelihood::OUTSIDE_SCALE_DEFAULT,
                         const std::string& ref_b0 = Likelihood::REF_B0_DEFAULT, double ref_signal =
                                 NAN)
                        
                        : Likelihood(obs_image, exp_image, assumed_snr, b0_include, outside_scale,
                                ref_b0, ref_signal) {
                }
                
                Gaussian(const Gaussian& s)
                        : Likelihood(s) {
                }
                
                ~Gaussian() {
                }
                
                Gaussian& operator=(const Gaussian& s) {
                    this->Likelihood::operator=(s);
                    return *this;
                }
                
                using Likelihood::log_prob;

                double log_prob_and_fisher(
                        const Fibre::Strand::Set& strands, Fibre::Strand::Set& gradient,
                        Fibre::Strand::Set::Tensor& fisher_info,
                        std::vector<Fibre::Strand::Set::Tensor>& fisher_info_gradients) {
                    return log_prob_and_fisher_tpl<Fibre::Strand>(strands, gradient, fisher_info,
                            fisher_info_gradients);
                }
                
                double log_prob_and_fisher(
                        const Fibre::Tractlet::Set& tractlets, Fibre::Tractlet::Set& gradient,
                        Fibre::Tractlet::Set::Tensor& fisher_info,
                        std::vector<Fibre::Tractlet::Set::Tensor>& fisher_info_gradients) {
                    return log_prob_and_fisher_tpl<Fibre::Tractlet>(tractlets, gradient,
                            fisher_info, fisher_info_gradients);
                }
                
                template<typename T> double log_prob_and_fisher_tpl(
                        const typename T::Set& fibres, typename T::Set& gradient,
                        typename T::Set::Tensor& fisher_info,
                        std::vector<typename T::Set::Tensor>& fisher_info_gradients);

                double log_prob_and_fisher(const Fibre::Strand::Set& strands,
                                           Fibre::Strand::Set& gradient,
                                           Fibre::Strand::Set::Tensor& fisher_info) {
                    return log_prob_and_fisher_tpl<Fibre::Strand>(strands, gradient, fisher_info);
                }
                
                double log_prob_and_fisher(const Fibre::Tractlet::Set& strands,
                                           Fibre::Tractlet::Set& gradient,
                                           Fibre::Tractlet::Set::Tensor& fisher_info) {
                    return log_prob_and_fisher_tpl<Fibre::Tractlet>(strands, gradient, fisher_info);
                }
                
                template<typename T> double log_prob_and_fisher_tpl(
                        const typename T::Set& fibres, typename T::Set& gradient,
                        typename T::Set::Tensor& fisher_info);

                double log_prob(double expected, double observed) {
                    
                    double diff = expected - observed;
                    
                    double lprob = -MR::Math::pow2(diff) / (2.0 * this->sigma2);
                    
                    //Add normalising constant.
//          lprob +=       - MR::Math::log(MR::Math::sqrt(2.0 * M_PI * this->sigma2));
                    
                    return lprob;
                    
                }
                
                double b0_log_prob(double expected, double observed) {
                    
                    double diff = expected - observed;
                    
                    double lprob;
                    
                    if (b0_include == "full" || ((b0_include == "half") && diff > 0)) {
                        lprob = -MR::Math::pow2(diff) / (2.0 * this->sigma2);
                        //Add normalising constant.
//            lprob +=       - MR::Math::log(MR::Math::sqrt(2.0 * M_PI * this->sigma2));
                    } else
                        lprob = 0.0;
                    
                    return lprob;
                    
                }
                
                double log_prob(double expected, double observed, double& d_lprob) {
                    
                    double diff = expected - observed;
                    
                    d_lprob = -diff / this->sigma2;
                    
                    return log_prob(expected, observed);
                    
                }
                
                double log_prob(double expected, double observed, double& d_lprob,
                                double& d2_lprob2) {
                    
                    d2_lprob2 = -1.0 / this->sigma2;
                    
                    return log_prob(expected, observed, d_lprob);
                    
                }
                
        };
    
    }

}

#include "bts/prob/likelihood/gaussian.cpp.h"

#endif
