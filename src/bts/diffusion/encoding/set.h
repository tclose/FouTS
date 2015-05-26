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

#ifndef __bts_diffusion_encoding_settemplate_h__
#define __bts_diffusion_encoding_settemplate_h__

#include <vector>

#include "bts/diffusion/encoding.h"
#include "bts/diffusion/model.h"

namespace FTS {
    
    namespace Diffusion {
        
        class Encoding::Set {
                
            protected:
                
                std::vector<Encoding> encodings;

                //Public member functions
            public:
                
                Set() {
                }
                
                Set(const MR::Math::Matrix<double>& encodings_matrix) {
                    set(encodings_matrix);
                }
                
                Set(const std::string& location) {
                    load(location);
                }
                
                Set(const Set& s)
                        : encodings(s.encodings) {
                }
                
                Set(const Model& m);

                Set& operator=(const Set& s) {
                    
                    this->encodings = s.encodings;
                    
                    return *this;
                }
                
                void load(const std::string& location);

                void set(const MR::Math::Matrix<double>& encodings_matrix);

                Encoding& operator[](size_t index) {
                    return encodings[index];
                }
                const Encoding& operator[](size_t index) const {
                    return encodings[index];
                }
                
                size_t num_encodings() const {
                    return encodings.size();
                }
                size_t size() const {
                    return encodings.size();
                }
                
                const Encoding& encoding(size_t index) const {
                    return encodings[index];
                }
                
                friend std::ostream& operator<<(std::ostream& stream,
                                                const Encoding::Set& encodings);
        };
        
        std::ostream& operator<<(std::ostream& stream, const Encoding::Set& encodings);
    
    }
}

#endif
