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

#include "bts/image/noise.h"

#include "bts/image/noise/gaussian.h"
#include "bts/image/noise/rician.h"
#include "bts/image/buffer.cpp.h"

namespace BTS {
    
    namespace Image {
        
        const char* Noise::TYPE_DEFAULT = "gaussian";
        const double Noise::SNR_DEFAULT = 15;
        
        Noise* Noise::factory(gsl_rng* rand_gen, const std::string& type, double snr,
                              double ref_signal) {
            
            if (type == "gaussian")
                
                return new Noise::Gaussian(rand_gen, snr, ref_signal);
            
            else if (type == "rician")
                
                return new Noise::Rician(rand_gen, snr, ref_signal);
            
            else
                
                throw Exception("Unrecognised noise type '" + type + "'.");
            
        }
    
    }

}
