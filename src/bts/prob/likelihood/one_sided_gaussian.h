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

#ifndef __bts_prob_likelihood_imagediff_onesidedgaussian_h__
#define __bts_prob_likelihood_imagediff_onesidedgaussian_h__

#include "bts/prob/prob.h"
#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood.h"

namespace BTS {
    
    namespace Prob {
        
        class Likelihood::OneSidedGaussian: public Likelihood {
                
            public:
                
                OneSidedGaussian(const Image::Observed::Buffer& obs_image,
                                 Image::Expected::Buffer* const exp_image, double assumed_snr =
                                         Likelihood::SNR_DEFAULT,
                                 const std::string& b0_include = Likelihood::B0_INCLUDE_DEFAULT,
                                 double outside_scale = Likelihood::OUTSIDE_SCALE_DEFAULT,
                                 const std::string& ref_b0 = Likelihood::REF_B0_DEFAULT,
                                 double ref_signal = NAN,
                                 const Image::Double::Buffer& noise_map = Image::Double::Buffer())
                        
                        : Likelihood(obs_image, exp_image, assumed_snr, b0_include, outside_scale,
                                ref_b0, ref_signal, noise_map) {
                }
                
                OneSidedGaussian(const OneSidedGaussian& s)
                        : Likelihood(s) {
                }
                
                ~OneSidedGaussian() {
                }
                
                OneSidedGaussian& operator=(const OneSidedGaussian& s) {
                    this->Likelihood::operator=(s);
                    return *this;
                }
                
                using Likelihood::log_prob;

                double log_prob(double expected, double observed) {
                    
                    double diff = expected - observed;
                    
                    double lprob;
                    
                    if (diff > 0)
                        lprob = -MR::Math::pow2(diff) / (2.0 * this->sigma2);    // - MR::Math::log(MR::Math::sqrt(2.0 * M_PI * this->sigma2));
                    else
                        lprob = 0.0;
                    
                    return lprob;
                    
                }
                
                double b0_log_prob(double expected, double observed) {
                    
                    double diff = expected - observed;
                    
                    double lprob;
                    
                    if (this->b0_include == "full" || this->b0_include == "half") {
                        
                        if (diff > 0)
                            lprob = -MR::Math::pow2(diff) / (2.0 * this->sigma2);    // - MR::Math::log(MR::Math::sqrt(2.0 * M_PI * this->sigma2));
                        else
                            lprob = 0.0;
                        
                    } else
                        lprob = 0.0;
                    
                    return lprob;
                    
                }
                
                double log_prob(double expected, double observed, double& d_lprob) {
                    
                    double diff = expected - observed;
                    
                    if (diff < 0)
                        d_lprob = -diff / this->sigma2;
                    else
                        d_lprob = 0.0;
                    
                    return log_prob(expected, observed);
                    
                }
                
                double log_prob(double expected, double observed, double& d_lprob,
                                double& d2_lprob2) {
                    
                    double diff = expected - observed;
                    
                    if (diff < 0)
                        d2_lprob2 = -1.0 / this->sigma2;
                    else
                        d2_lprob2 = 0.0;
                    
                    return log_prob(expected, observed, d_lprob);
                    
                }
                
        };
    
    }

}

#endif
