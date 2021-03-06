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

#ifndef __bts_prob_priorcomponent_frequency_h__
#define __bts_prob_priorcomponent_frequency_h__

#include "math/math.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/strand/tensor.h"
#include "bts/fibre/strand/set/tensor.h"

namespace FTS {
    
    namespace Prob {
        
        namespace PriorComponent {
            
            class Frequency {
                    
                public:
                    
                    const static double SCALE_DEFAULT;
                    const static double AUX_SCALE_DEFAULT;
                    const static std::string NAME;

                    //Protected member variables.
                protected:
                    
                    double scale;
                    double aux_scale;

                public:
                    
                    Frequency(double scale, double aux_scale)
                            : scale(scale), aux_scale(aux_scale) {
                    }
                    
                    Frequency(const Frequency& m)
                            : scale(m.scale), aux_scale(m.aux_scale) {
                    }
                    
                    Frequency& operator=(const Frequency& m) {
                        scale = m.scale;
                        aux_scale = m.aux_scale;
                        return *this;
                    }
                    
                    ~Frequency() {
                    }
                    
                    Frequency* clone() const {
                        return new Frequency(*this);
                    }
                    
                    double log_prob(const Fibre::Strand strand, Fibre::Strand gradient);

                    double log_prob(const Fibre::Tractlet, Fibre::Tractlet gradient);

                    double log_prob(const Fibre::Strand strand, Fibre::Strand gradient,
                                    Fibre::Strand::Tensor hessian) {
                        throw Exception("Not implemented yet.");
                    }
                    
                    double log_prob_and_fisher(const Fibre::Strand strand, Fibre::Strand gradient,
                                               Fibre::Strand::Tensor fisher) {
                        throw Exception("Not implemented yet.");
                    }
                    
                    double log_prob_and_fisher(
                            const Fibre::Strand strand, Fibre::Strand gradient,
                            Fibre::Strand::Tensor fisher,
                            std::vector<Fibre::Strand::Tensor>& fisher_gradient) {
                        throw Exception("Not implemented yet.");
                    }
                    
                    const std::string& get_name() {
                        return NAME;
                    }
            };
        
        }
    
    }

}

#endif
