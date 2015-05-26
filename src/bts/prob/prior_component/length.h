/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

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

#ifndef __bts_prob_priorcomponent_length_h__
#define __bts_prob_priorcomponent_length_h__

#include "math/math.h"

#include "bts/fibre/strand/set.h"

namespace FTS {
    
    namespace Prob {
        
        namespace PriorComponent {
            
            class Length {
                    
                public:
                    
                    const static double SCALE_DEFAULT;
                    const static double MEAN_DEFAULT;
                    const static std::string NAME;

                    //Protected member variables.
                protected:
                    
                    double scale;
                    double mean;

                public:
                    
                    Length(double scale, double mean)
                            : scale(scale), mean(mean) {
                    }
                    
                    Length(const Length& m)
                            : scale(m.scale), mean(m.mean) {
                    }
                    
                    Length& operator=(const Length& m) {
                        scale = m.scale;
                        mean = m.mean;
                        return *this;
                    }
                    
                    ~Length() {
                    }
                    
                    Length* clone() const {
                        return new Length(*this);
                    }
                    
                    double log_prob(const Fibre::Strand strand);

                    double log_prob(const Fibre::Strand strand, Fibre::Strand gradient);

                    const std::string& get_name() {
                        return NAME;
                    }
            };
        
        }
    
    }

}

#endif
