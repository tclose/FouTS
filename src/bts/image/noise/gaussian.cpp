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

#include "math/math.h"

#include "bts/image/noise/gaussian.h"
#include "bts/image/observed/voxel.h"

#include "bts/image/inline_functions.h"

namespace FTS {
    
    namespace Image {
        
        Noise::Gaussian::Gaussian(gsl_rng* rand_gen, double snr, double ref_signal)
                : Noise(rand_gen, snr, ref_signal) {
        }
        
        template<typename T> void Noise::Gaussian::noisify_tpl(T& image) {
            
            double noise_mag = noise_magnitude(image);
            
            for (size_t x = 0; x < image.dim(X); x++)
                for (size_t y = 0; y < image.dim(Y); y++)
                    for (size_t z = 0; z < image.dim(Z); z++)
                        for (size_t orient_i = 0; orient_i < image.num_encodings(); orient_i++)
                            
                            image(x, y, z)[orient_i] += gsl_ran_gaussian(rand_gen, noise_mag);
            
        }
    
    }

}
