/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

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

#ifndef __bts_prob_priorcomponent_thinness_h__
#define __bts_prob_priorcomponent_thinness_h__

#include "math/math.h"

#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/tractlet/tensor.h"

namespace BTS {
    
    namespace Prob {
        
        namespace PriorComponent {
            
            class Thinness {
                    
                public:
                    
                    const static double SCALE_DEFAULT;
                    const static size_t POWER_DEFAULT;
                    const static std::string NAME;

                    //Protected member variables.
                protected:
                    
                    double scale;
                    size_t power;

                public:
                    
                    Thinness(double scale, size_t power)
                            : scale(scale), power(power) {
                    }
                    
                    Thinness(const Thinness& s)
                            : scale(s.scale), power(s.power) {
                    }
                    
                    Thinness& operator=(const Thinness& s) {
                        scale = s.scale;
                        power = s.power;
                        return *this;
                    }
                    
                    ~Thinness() {
                    }
                    
                    Thinness* clone() const {
                        return new Thinness(*this);
                    }
                    
                    double log_prob(const Fibre::Tractlet tract);

                    double log_prob(const Fibre::Tractlet tract, Fibre::Tractlet gradient);

                    double log_prob_and_fisher(const Fibre::Tractlet tractlet,
                                               Fibre::Tractlet gradient,
                                               Fibre::Tractlet::Tensor hessian) {
                        throw Exception("Not implemented");
                    }
                    
                    double log_prob_and_fisher(
                            const Fibre::Tractlet tractlet, Fibre::Tractlet gradient,
                            Fibre::Tractlet::Tensor fisher,
                            std::vector<Fibre::Tractlet::Tensor>& fisher_gradient) {
                        throw Exception("Not implemented");
                    }
                    
                    const std::string& get_name() {
                        return NAME;
                    }
            };
        
        }
    
    }

}

#endif
