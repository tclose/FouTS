/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

 This file is part of Bayesian Tractlet Sampling (BTS).

 BTS is free software: you can reimageibute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 BTS is imageibuted in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with BTS.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef __bts_analysis_imagegradienttester_h__
#define __bts_analysis_imagegradienttester_h__

namespace BTS {
    
    namespace Analysis {
        
        template<typename ImageClass, typename State> class ImageGradientTester {
                
            public:
                
                //typedef the pointer to the function to test as 'Function'.
                typedef ImageClass& (ImageClass::*Function)(
                        const State&, typename Image::Container::Buffer<State>&);

            protected:
                
                ImageClass image;

            public:
                
                ImageGradientTester(const ImageClass& image)
                        : image(image) {
                    this->image.clear_and_enforce_bounds();
                }
                
                ~ImageGradientTester() {
                }
                
                void test(Function function, State& state, double step_size,
                          Image::Container::Buffer<State>& analytic_gradient,
                          Image::Container::Buffer<State>& numeric_gradient);
                
        };
    
    }

}

#include "progressbar.h"

#include "bts/common.h"
#include "bts/image/observed/buffer.h"

#include "bts/analysis/image_gradient_tester.h"

namespace BTS {
    
    namespace Analysis {
        
        template<typename ImageClass, typename State>
        void ImageGradientTester<ImageClass, State>::test(
                Function function, State& state, double step_size,
                Image::Container::Buffer<State>& analytic_gradient,
                Image::Container::Buffer<State>& numeric_gradient) {
            
            //Initialise gradient buffers.
            
            Image::Container::Buffer<State> dummy_gradient(image.dims(), image.num_encodings());
            
            analytic_gradient.reset(image.dims(), image.num_encodings());
            numeric_gradient.reset(image.dims(), image.num_encodings());
            
            //Calculate analytic hessian and unperturbed gradient
            
            image.zero();
            
            Image::Buffer image_intensities = (image.*function)(state, analytic_gradient);
            
            //Initialize state vector
            
            MR::Math::Vector<double>& state_vector = state;
            
            Image::Container::Buffer<MR::Math::Vector<double> > numeric_gradient_vectors(
                    image.dims(), image.num_encodings());
            
            for (size_t z = 0; z < image.dim(Z); z++)
                for (size_t y = 0; y < image.dim(Y); y++)
                    for (size_t x = 0; x < image.dim(X); x++)
                        for (size_t encode_i = 0; encode_i < image.num_encodings(); encode_i++)
                            numeric_gradient_vectors(x, y, z)[encode_i].resize(state_vector.size(),
                                    0.0);
            
            MR::ProgressBar progress_bar("Testing gradient calculations...", state_vector.size());
            
            for (size_t elem_i = 0; elem_i < state_vector.size(); ++elem_i) {
                
                state_vector[elem_i] += step_size;
                
                image.zero();
                
                ImageClass& perturbed_image_intensities = (image.*function)(state, dummy_gradient);
                
                for (Image::Buffer::iterator vox_it = image_intensities.begin();
                        vox_it != image_intensities.end(); ++vox_it)
                    for (size_t encode_i = 0; encode_i < image.num_encodings(); encode_i++)
                        if (!perturbed_image_intensities.is_empty(vox_it->first))
                            numeric_gradient_vectors(vox_it->first)[encode_i][elem_i] =
                                    (perturbed_image_intensities(vox_it->first)[encode_i] - vox_it->second[encode_i]) / step_size;
                
                state_vector[elem_i] -= step_size;
                
                ++progress_bar;
                
            }
            
            for (size_t z = 0; z < image.dim(Z); z++)
                for (size_t y = 0; y < image.dim(Y); y++)
                    for (size_t x = 0; x < image.dim(X); x++) {
                        
                        for (size_t encode_i = 0; encode_i < image.num_encodings(); encode_i++) {
                            numeric_gradient(x, y, z)[encode_i] = state;
                            numeric_gradient(x, y, z)[encode_i] =
                                    numeric_gradient_vectors(x, y, z)[encode_i];
                        }
                        
                        if (analytic_gradient.is_empty(x, y, z))
                            for (size_t encode_i = 0; encode_i < image.num_encodings(); encode_i++)
                                analytic_gradient(x, y, z)[encode_i] = State(state).zero();
                        
                    }
            
        }
    
    }

}

#endif
