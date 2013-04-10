/*
 Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

 Written by Thomas G Close, 5/05/09.

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

#ifndef __image_observed_voxel_h__
#define __image_observed_voxel_h__

#include <ostream>
#include <vector>

#include "bts/common.h"

#include "bts/image/observed/buffer.h"
#include "bts/diffusion/encoding.h"
#include "bts/image/voxel.h"
#include "bts/image/observed/buffer.h"
#include "bts/image/voxel.h"

namespace BTS {
    
    namespace Image {
        
        namespace Observed {
            
            class Voxel: public Image::Voxel<double> {
                    
                    //Public nested classes and typedefs
                public:
                    
                    typedef Observed::Buffer Buffer;

                protected:
                    
                    Buffer* obs_image;

                    Index coordinate;
                    Coord centre_point;

                public:
                    
                    Voxel() {
                    }
                    
                    Voxel(Buffer& buffer, const Index& coord);

                    Voxel(const Voxel& v)
                            : Image::Voxel<double>(v), obs_image(v.obs_image), coordinate(
                                      v.coordinate), centre_point(v.centre_point) {
                    }
                    
                    Voxel& operator=(const Voxel& v) {
                        
                        Image::Voxel<double>::operator=(v);
                        obs_image = v.obs_image;
                        coordinate = v.coordinate;
                        centre_point = v.centre_point;
                        
                        return *this;
                    }
                    
                    //Doesn't copy across the parent image or coordinate.
                    Voxel& copy_value(const Voxel& v) {
                        Image::Voxel<double>::operator=(v);
                        return *this;
                    }
                    
                    virtual ~Voxel() {
                    }
                    
                    virtual const Diffusion::Encoding& encoding(size_t index) const;

                    double b0() const;

                    const Coord& centre() const {
                        return centre_point;
                    }
                    
                    const Index& coord() const {
                        return coordinate;
                    }
                    
                    double sum_square_signal(bool include_b0s) const;

                    double mean_square_signal(bool include_b0s) const {
                        return sum_square_signal(include_b0s) / num_encodings();
                    }
                    
                    double sum_strong_difference(const Voxel& expected_image,
                                                 Voxel& expected_image_gradient,
                                                 bool include_b0s) const;

                    double sum_weak_difference(const Voxel& expected_image,
                                               Voxel& expected_image_gradient,
                                               bool include_b0s) const;

                    double isotropic() const;

                    void remove_isotropic();

                    double rms() const;

                    Voxel operator-() const {
                        Voxel answer(*this);
                        for (size_t encode_i = 0; encode_i < num_encodings(); encode_i++)
                            answer[encode_i] = -operator[](encode_i);
                        return answer;
                    }
                    
                protected:
                    
                    Voxel(size_t num_encodings, const Index& coord);

                    //      friend class Observed::Buffer;
                    
                    friend std::ostream& operator<<(std::ostream& stream,
                                                    const Observed::Voxel& voxel);
                    
            };
            
            inline std::ostream& operator<<(std::ostream& stream, const Observed::Voxel& voxel) {
                
                stream << "coord: " << voxel.coordinate << std::endl;
                stream << "centre point: " << voxel.centre_point << std::endl;
                stream << "intensities: " << std::endl;
                for (size_t encode_i = 0; encode_i < voxel.num_encodings(); encode_i++) {
                    stream << voxel[encode_i] << ", ";
                }
                stream << std::endl;
                
                return stream;
                
            }
        
        }
    
    }

}

#undef LOOP

#endif
