/*
 Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

 Written by Thomas G Close, 5/05/09.

 This file is part of Fourier Tract Sampling (FouTS).

 FouTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FouTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FTS.  If not, see <http://www.gnu.org/licenses/>.

 */

#include "bts/prob/likelihood.h"
#include "bts/image/observed/buffer.h"
#include "bts/image/expected/trilinear/buffer.h"

#include "bts/prob/likelihood/gaussian.h"
#include "bts/prob/likelihood/one_sided_gaussian.h"
#include "bts/prob/likelihood/rician.h"

#include "bts/prob/likelihood/gaussian.cpp.h"

#include "bts/prob/inline_functions.h"
#include "bts/image/inline_functions.h"
#include "bts/prob/likelihood.cpp.h"
#include "bts/image/buffer.cpp.h"
#include "bts/diffusion/response.h"

namespace FTS {
    
    namespace Prob {
        
        const char* Likelihood::TYPE_DEFAULT = "gaussian";
        const double Likelihood::SNR_DEFAULT = 20.0;
        const char* Likelihood::B0_INCLUDE_DEFAULT = "none";    //Also "half" or "full"
        const double Likelihood::OUTSIDE_SCALE_DEFAULT = 1.0;
        const char* Likelihood::REF_B0_DEFAULT = "max";    // Also could be "average"
        
        const Image::Observed::Buffer Likelihood::dummy_observed_buffer;
        const Image::Expected::Trilinear::Buffer Likelihood::dummy_expected_buffer;
        
        Likelihood* Likelihood::factory(const std::string& type,
                                        const Image::Observed::Buffer& obs_image,
                                        Image::Expected::Buffer* exp_image, double assumed_snr,
                                        const std::string& b0_include, double outside_scale,
                                        const std::string& ref_b0, double ref_signal,
                                        const Image::Double::Buffer& noise_map) {
            
            //FIXME: Not sure why this happens, but before exp_image is cloned it must be cleared to produce consistent results.
            //Investigate and resolve.
            exp_image->clear();
            
            Likelihood* likelihood;
            
            if (type == "gaussian")
                
                likelihood = new Gaussian(obs_image, exp_image, assumed_snr, b0_include,
                        outside_scale, ref_b0, ref_signal, noise_map);
            
            else if (type == "one_sided_gaussian")
                
                likelihood = new OneSidedGaussian(obs_image, exp_image, assumed_snr, b0_include,
                        outside_scale, ref_b0, ref_signal, noise_map);
            
            else if (type == "rician") {
                
                const Diffusion::Model& diffusion_model = exp_image->get_diffusion_model();
                
                if (!diffusion_model.includes_isotropic())
                    throw Exception(
                            "Only diffusion models that include isotropic components can be used with Rician noise model.");
                
                if (!exp_image->bounds_are_enforced())
                    throw Exception(
                            "Enforce bounds must be set with Rician noise model to avoid -inf values.");
                
                likelihood = new Rician(obs_image, exp_image, assumed_snr, b0_include,
                        outside_scale, ref_b0, ref_signal, noise_map);
                
            } else
                throw Exception(
                        "Unrecognised likelihood type '" + type
                        + "' passed to option '-like_type'.");
            
            return likelihood;
            
        }
        
        Likelihood::Likelihood(const Image::Observed::Buffer& observed_image,
                               Image::Expected::Buffer* const expected_image, double assumed_snr,
                               const std::string& b0_include, double outside_scale,
                               const std::string& ref_b0, double ref_signal,
                               const Image::Double::Buffer& noise_map)
                : sigma2_map(noise_map), obs_image(observed_image), exp_image(expected_image->clone()),
                  b0_include(b0_include) {
            
            if (!expected_image->dims_match(obs_image))
                throw Exception(
                        "Expected image dimensions (" + str(expected_image->dims())
                        + ") do not match and observed image (" + str(obs_image.dims()) + ").");
            
            if (!sigma2_map.dim(X))
                set_assumed_snr(assumed_snr, ref_b0, ref_signal);
            else {
                sigma2_map *= noise_map; // Convert the noise map into a variance map
                // Set the maximum element as the default value for when the signal lies outside the
                // true image bounds.
                sigma2 = -INFINITY;
                for (Image::Double::Buffer::iterator it = sigma2_map.begin(); it != sigma2_map.end(); ++it)
                    if (it->second > sigma2)
                        sigma2 = it->second;
            }

            
            this->exp_image->zero();
            
        }
        
        Likelihood::Likelihood(const Likelihood& l)
                : sigma2(l.sigma2), sigma2_map(l.sigma2_map), obs_image(l.obs_image),
                  exp_image(l.exp_image->clone()), difference_mode(l.difference_mode),
                  b0_include(l.b0_include)
                   {
        }
        
        Likelihood& Likelihood::operator=(const Likelihood& l) {
            
            delete this->exp_image;
            
            sigma2 = l.sigma2;
            sigma2_map = l.sigma2_map;
            obs_image = l.obs_image;
            exp_image = l.exp_image->clone();
            b0_include = l.b0_include;
            difference_mode = l.difference_mode;
            
            return *this;
            
        }
        
        void Likelihood::set_assumed_snr(double assumed_snr, const std::string& ref_b0,
                                         double reference_signal) {
            
            if (!!sigma2_map)
                throw Exception("SNR cannot be set explicitly if noise map is already provided.");

            double ref_signal;

            if (!isnan(reference_signal))
                ref_signal = reference_signal;
            else if (obs_image.properties().count("noise_ref_signal"))
                ref_signal = to<double>(obs_image.properties()["noise_ref_signal"]);
            else if (ref_b0 == "average")
                ref_signal = obs_image.average_b0();
            else if (ref_b0 == "max")
                ref_signal = obs_image.max_b0();
            else
                throw Exception(
                        "Unrecognised value for '-like_ref_b0' ('" + ref_b0
                        + "'), can be either 'average', or 'max'.");

            if (!ref_signal)
                throw Exception("The " + ref_b0 + " b0 of observed image is 0.");

            sigma2 = MR::Math::pow2(ref_signal / assumed_snr);
            
        }
        
        void Likelihood::set_enforce_bounds(bool flag) {
            
            if (flag)
                exp_image->clear_and_enforce_bounds();
            else
                exp_image->relax_bounds();
            
        }
        
        //TODO: Collapse the workings of the different log_prob_tpl's into common subfunctions.
        double Likelihood::log_prob(Image::Expected::Buffer& image) {
            
            double lprob = 0.0;
            
            std::set<Image::Index> coords = image.non_empty();
            
            for (size_t x = 0; x < image.dim(X); ++x)
                for (size_t y = 0; y < image.dim(Y); ++y)
                    for (size_t z = 0; z < image.dim(Z); ++z)
                        coords.insert(Image::Index(x, y, z));
            
            for (std::set<Image::Index>::iterator index_it = coords.begin();
                    index_it != coords.end(); ++index_it) {
                
                const Image::Expected::Voxel& exp_voxel = image.operator()(*index_it);
                
                for (size_t encode_i = 0; encode_i < image.num_encodings(); encode_i++) {
                    
                    double observed;
                    
                    if (obs_image.in_bounds(*index_it))
                        observed = obs_image(*index_it)[encode_i];
                    else
                        observed = 0.0;
                    
                    if (image.encoding(encode_i).b_value())
                        lprob += log_prob(exp_voxel[encode_i], observed, *index_it);
                    else
                        lprob += b0_log_prob(exp_voxel[encode_i], observed);
                    
                }
                
            }
            
            return lprob;
            
        }
    
    }

}
