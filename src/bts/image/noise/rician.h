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

#ifndef __bts_image_noise_rician_h__
#define __bts_image_noise_rician_h__

#include "bts/image/noise.h"
#include "bts/image/observed/buffer.h"
#include "bts/image/expected/buffer.h"

namespace BTS {
    
    namespace Image {
        
        class Noise::Rician: public Noise {
                
                //Public static constants.
            public:
                
                //Protected member variables
            protected:
                
                //Public member functions
            public:
                
                Rician(gsl_rng* rand_gen, double snr, double ref_signal);

                Image::Observed::Buffer& noisify(Image::Observed::Buffer& image) {
                    noisify_tpl(image);
                    return image;
                }
                
                Image::Expected::Buffer& noisify(Image::Expected::Buffer& image) {
                    noisify_tpl(image);
                    return image;
                }
                
            protected:
                
                template<typename T> void noisify_tpl(T& image);
                
        };
    
    }

}

#endif /* __bts_image_noise_rician_h__ */
