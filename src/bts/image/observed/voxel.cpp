/*
 Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

 Written by Thomas G Close, 5/05/09.

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

#include "bts/image/voxel.h"
#include "bts/image/observed/buffer.h"
#include "bts/image/observed/voxel.h"

#include "bts/image/inline_functions.h"

#define LOOP(op) \
for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++) {\
  op \
}

namespace FTS {
    
    namespace Image {
        
        namespace Observed {
            
            Voxel::Voxel(Observed::Buffer& obs_image, const Index& coord)
                    : Image::Voxel<double>(obs_image.num_encodings()), obs_image(&obs_image), coordinate(
                              coord), centre_point(Observed::Buffer::voxel_centre(coordinate)) {
            }
            
            Voxel::Voxel(size_t num_encodings, const Index& coord)
                    : Image::Voxel<double>(num_encodings), obs_image(0), coordinate(coord), centre_point(
                              Observed::Buffer::voxel_centre(coordinate)) {
            }
            
            const Diffusion::Encoding& Voxel::encoding(size_t index) const {
                return obs_image->encoding(index);
            }
            
            // NB: This will not work for signal with isotropic component removed.
            double Voxel::sum_weak_difference(const Voxel& expected_image,
                                              Voxel& expected_image_gradient,
                                              bool include_b0s) const {
                
                double sum = 0.0;
                
                for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++) {
                    
                    // If direction is not a b=0 direction.
                    if (this->encoding(encode_i).b_value() || include_b0s) {
                        
                        double diff = this->operator[](encode_i) - expected_image[encode_i];
                        
                        if (diff > 0.0) {
                            
                            sum += diff * diff;
                            expected_image_gradient[encode_i] = -2.0 * diff;
                            
                        } else
                            expected_image_gradient[encode_i] = 0.0;
                        
                    } else
                        expected_image_gradient[encode_i] = 0.0;
                    
                }
                
                return sum;
            }
            
            double Voxel::sum_strong_difference(const Voxel& expected_image,
                                                Voxel& expected_image_gradient,
                                                bool include_b0s) const {
                
                double sum = 0;
                
                for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++) {
                    
                    // If direction is not a b=0 direction.
                    if (this->encoding(encode_i).b_value() || include_b0s) {
                        
                        double diff = this->operator[](encode_i) - expected_image[encode_i];
                        sum += diff * diff;
                        expected_image_gradient[encode_i] = -2.0 * diff;
                        
                    } else
                        expected_image_gradient[encode_i] = 0.0;
                    
                }
                
                return sum;
            }
            
            double Voxel::sum_square_signal(bool include_b0s) const {
                
                double sum = 0;
                
                for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++) {
                    
                    // If direction is not a b=0 direction.
                    if (this->encoding(encode_i).b_value() || include_b0s)
                        sum += MR::Math::pow2(this->operator[](encode_i));
                    
                }
                
                return sum;
                
            }
            
            double Voxel::b0() const {
                
                double sum_b0 = 0.0;
                size_t num_b0 = 0;
                
                for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++) {
                    
                    if (this->encoding(encode_i).b_value() == 0.0) {
                        sum_b0 += this->operator[](encode_i);
                        ++num_b0;
                    }
                }
                
                return sum_b0 / (double) num_b0;
            }
            
            double Voxel::isotropic() const {
                
                double sum = 0.0;
                
                size_t num_b0s = 0;
                
                for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++) {
                    
                    if (this->encoding(encode_i).b_value())
                        sum += this->operator[](encode_i);
                    else
                        ++num_b0s;
                    
                }
                
                return sum / (double) (num_encodings() - num_b0s);
                
            }
            
            double Voxel::rms() const {
                
                double sum = 0.0;
                
                for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++) {
                    
                    sum += MR::Math::pow2(this->operator[](encode_i));
                    
                }
                
                return MR::Math::sqrt(sum / (double) num_encodings());
                
            }
            
            void Voxel::remove_isotropic() {
                
                double iso = this->isotropic();
                
                for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++) {
                    
                    if (this->encoding(encode_i).b_value())
                        this->operator[](encode_i) -= iso;
                    
                }
                
            }
        
        }
    
    }

}

#undef LOOP
