/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 14/07/2010.

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

#ifndef __bts_image_noise_h__
#define __bts_image_noise_h__

#define NOISE_PARAMETERS \
  Option ("noise_type", "Type of noise distribution. Can be either 'gaussian' or 'rician' (default 'gaussian').") \
   + Argument ("noise_type", "").type_text (Image::Noise::TYPE_DEFAULT), \
\
  Option ("noise_snr", "The signal-to-noise ratio of the noise.") \
   + Argument ("noise_snr", "").type_float (SMALL_FLOAT, Image::Noise::SNR_DEFAULT, LARGE_FLOAT), \
\
  Option ("noise_ref_signal", "The reference intensity that the noise variance for a given SNR will be calculated from.") \
   + Argument ("noise_ref_signal", "").type_float (SMALL_FLOAT, NAN, LARGE_FLOAT) \
\


//Loads the parameters into variables
#define SET_NOISE_PARAMETERS \
\
  std::string      noise_type           = Image::Noise::TYPE_DEFAULT; \
  double           noise_snr            = Image::Noise::SNR_DEFAULT; \
  double           noise_ref_signal     = NAN; \
\
  opt = get_options("noise_type"); \
  if (opt.size()) \
    noise_type = opt[0][0].c_str(); \
\
  opt = get_options("noise_snr"); \
  if (opt.size()) \
    noise_snr = opt[0][0]; \
\
  opt = get_options("noise_ref_signal"); \
  if (opt.size()) \
    noise_ref_signal = opt[0][0]; \


//Adds the  parameters to the properties to be saved with the data.
#define ADD_NOISE_PROPERTIES(properties) \
  properties["noise_type"]    = noise_type; \
  properties["noise_snr"]     = str(noise_snr); \
  if (!isnan(noise_ref_signal)) \
    properties["noise_ref_signal"] = str(noise_ref_signal);

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

#include "bts/image/observed/buffer.h"
#include "bts/image/expected/buffer.h"

namespace BTS {
    
    namespace Image {
        
        class Noise {
                
                //Public nested classes and typedefs
            public:
                
                class Gaussian;
                class Rician;

                //Public static constants
            public:
                
                const static char* TYPE_DEFAULT;
                const static double SNR_DEFAULT;

                //Protected member variables
            protected:
                
                gsl_rng* rand_gen;
                double snr;
                double ref_signal;

                //Public static member functions
            public:
                
                static Noise* factory(gsl_rng* rand_gen, const std::string& type, double snr,
                                      double ref_signal);

                //Public member functions
            public:
                
                Noise(gsl_rng* rand_gen, double snr, double ref_signal)
                        : rand_gen(rand_gen), snr(snr), ref_signal(ref_signal) {
                }
                
                virtual Image::Observed::Buffer& noisify(Image::Observed::Buffer& image) = 0;

                virtual Image::Expected::Buffer& noisify(Image::Expected::Buffer& image) = 0;

            protected:
                
                template<typename T> double noise_magnitude(T& image) const {
                    
                    double signal;
                    
                    if (!isnan(ref_signal))
                        signal = ref_signal;
                    else
                        signal = image.max_b0();
                    
                    return signal / snr;
                    
                }
                
        };
    
    }

}

#endif /* __bts_image_noise_h__ */
