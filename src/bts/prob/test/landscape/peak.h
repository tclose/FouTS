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

#ifndef __bts_prob_testpeaks_peak__
#define __bts_prob_testpeaks_peak__

#include "bts/prob/test/landscape.h"
#include "bts/mcmc/state/tensor.h"

namespace FTS {
    
    namespace Prob {
        
        namespace Test {
            
            class Landscape::Peak {
                    
                public:
                    
                    class Gaussian;
                    class Pyramid;
                    class Quartic;

                public:
                    
                    const static size_t GAUSSIAN = 0;
                    const static size_t PYRAMID = 1;

                public:
                    
                    class Reader;
                    class Writer;

                    friend class Reader;
                    friend class Writer;

                protected:
                    
                    size_t type;
                    size_t ndims;
                    double height;
                    MCMC::State centre;
                    MCMC::State::Tensor widths;

                public:
                    
                    Peak() {
                    }
                    
                    Peak(size_t ndims, const MCMC::State& centre, const MCMC::State::Tensor& widths,
                         double height, size_t type = GAUSSIAN);

                    double log_prob(const MCMC::State& point, MCMC::State& gradient);

                    double log_prob(const MCMC::State& point, MCMC::State& gradient,
                                    MCMC::State::Tensor& hessian);

                    double log_prob(const MCMC::State& point, MCMC::State& gradient,
                                    MCMC::State::Tensor& hessian,
                                    std::vector<MCMC::State::Tensor>& rank3_hessian);

                    double gaussian_log_prob(const MCMC::State& point, MCMC::State& gradient) {
                        MCMC::State::Tensor dummy;
                        return log_prob(point, gradient, dummy);
                    }
                    
                    double gaussian_log_prob(const MCMC::State& point, MCMC::State& gradient,
                                             MCMC::State::Tensor& hessian);

                    double gaussian_log_prob(const MCMC::State& point, MCMC::State& gradient,
                                             MCMC::State::Tensor& hessian,
                                             std::vector<MCMC::State::Tensor>& rank3_hessian);

                    double pyramid_log_prob(const MCMC::State& point, MCMC::State& gradient);

                    size_t num_dims() const {
                        return ndims;
                    }
                    
                    friend std::ostream& operator<<(std::ostream& stream, const Peak& peak);
                    
            };
        
        }
    
    }

}

#endif
