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

#ifndef __bts_prob_prior_tractlet_density__
#define __bts_prob_prior_tractlet_density__

#include "bts/fibre/tractlet.h"
#include "bts/fibre/tractlet/tensor.h"

namespace BTS {
    
    namespace Prob {
        
        namespace PriorComponent {
            
            class Density {
                    
                public:
                    
                    const static double HIGH_SCALE_DEFAULT;
                    const static double LOW_SCALE_DEFAULT;
                    const static std::string NAME;
                    const static size_t NUM_POINTS_DEFAULT;

                protected:
                    
                    double high_scale;
                    double low_scale;
                    size_t num_points;

                public:
                    
                    Density() {
                    }
                    
                    Density(double high_scale, double low_scale, size_t num_points =
                            NUM_POINTS_DEFAULT)
                            : high_scale(high_scale), low_scale(low_scale), num_points(num_points) {
                    }
                    
                    Density(const Density& d)
                            : high_scale(d.high_scale), low_scale(d.low_scale), num_points(
                                      num_points) {
                    }
                    
                    Density& operator=(const Density& d) {
                        high_scale = d.high_scale;
                        low_scale = d.low_scale;
                        num_points = d.num_points;
                        return *this;
                    }
                    
                    Density* clone() const {
                        return new Density(*this);
                    }
                    
                    double log_prob(const Fibre::Tractlet, Fibre::Tractlet gradient);

                    double log_prob(const Fibre::Tractlet tractlet, Fibre::Tractlet gradient,
                                    Fibre::Tractlet::Tensor hessian) {
                        throw Exception("Not implemented");
                    }
                    
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
