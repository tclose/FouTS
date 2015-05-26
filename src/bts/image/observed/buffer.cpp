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

#include "image/header.h"
#include "dataset/loop.h"

#include "bts/image/observed/buffer.h"

#include "bts/image/inline_functions.h"

namespace FTS {
    
    namespace Image {
        
        namespace Observed {
            
            std::ostream& operator<<(std::ostream& stream, const Buffer& buffer) {
                return buffer.to_stream(stream);
            }
            
            void Buffer::load(const std::string& location,
                              const Diffusion::Encoding::Set& default_encodings) {
                
                MR::Image::Header header(location);
                
                if (header.ndim() != 4)
                    throw Exception(
                            "dwi image should contain 4 dimensions, found " + str(header.ndim())
                            + ".");
                
                Diffusion::Encoding::Set encodings;
                if (header.get_DW_scheme().rows()) {
                    MR::Math::Matrix<double> encodings_matrix;
                    encodings.set(header.get_DW_scheme());
                } else if (default_encodings.size())
                    encodings = default_encodings;
                else
                    throw Exception(
                            "No encodings found in loaded image, and no defaults supplied.");
                
                if (header.dim(3) != (int) encodings.size())
                    throw Exception(
                            "Number of encoding directions in loaded image (" + str(header.dim(DW))
                            + ") do not match that in encoding file (" + str(encodings.size())
                            + ").");
                
                //Get image dimensions and voxel lengths.
                Triple<size_t> dimensions(header.dim(X), header.dim(Y), header.dim(Z));
                Triple<double> voxel_lengths(header.vox(X), header.vox(Y), header.vox(Z));
                
                //Get offset from transform.
                MR::Math::Matrix<double> transform = header.transform();
                Triple<double> mrtrix_offset(transform(0, 3), transform(1, 3), transform(2, 3));
                Triple<double> bts_offset = mrtrix_offset - voxel_lengths * 0.5;
                
                //Resize buffer to fit image.
                reset(dimensions, voxel_lengths, bts_offset, encodings);
                
                //Copy data from image to buffer.
                MR::Image::Voxel<double> voxel(header);
                MR::DataSet::Loop loop(0, 4);
                
                for (loop.start(voxel); loop.ok(); loop.next(voxel))
                    operator()(voxel[X], voxel[Y], voxel[Z])[voxel[DW]] = voxel.value();
                
                properties().insert(header.begin(), header.end());
                
            }
            
            Double::Buffer Buffer::isotropic(bool include_b0s) const {
                
                Double::Buffer mean_buffer(this->dims());
                
                for (const_iterator vox_it = this->begin(); vox_it != this->end(); ++vox_it) {
                    mean_buffer(vox_it->first) = vox_it->second.isotropic();
                }
                
                return mean_buffer;
                
            }
            
            void Buffer::remove_isotropic() {
                
                for (iterator vox_it = this->begin(); vox_it != this->end(); ++vox_it) {
                    vox_it->second.remove_isotropic();
                }
                
            }
            
            double Buffer::mean_square(bool include_b0s) const {
                
                double sum = 0.0;
                
                for (const_iterator vox_it = this->begin(); vox_it != this->end(); ++vox_it) {
                    sum += vox_it->second.mean_square_signal(include_b0s);
                }
                
                return sum / this->num_voxels_in_bounds();
                
            }
            
            Double::Buffer Buffer::b0() const {
                
                Double::Buffer b0_buffer(this->dims());
                
                for (const_iterator vox_it = this->begin(); vox_it != this->end(); ++vox_it) {
                    b0_buffer(vox_it->first) = vox_it->second.b0();
                }
                
                return b0_buffer;
                
            }
            
            double Buffer::average_b0() const {
                
                double b0_sum = 0;
                
                for (const_iterator vox_it = this->begin(); vox_it != this->end(); ++vox_it) {
                    b0_sum += vox_it->second.b0();
                }
                
                return b0_sum / this->num_voxels_in_bounds();
                
            }
        
        }
    
    }

}

#undef LOOP
#undef CONST_LOOP
