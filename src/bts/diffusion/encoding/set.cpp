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

#include "bts/diffusion/encoding/set.h"

#include "bts/diffusion/inline_functions.h"

namespace BTS {
    
    namespace Diffusion {
        
        Encoding::Set::Set(const Model& m) {
            
            for (size_t encode_i = 0; encode_i < m.num_encodings(); encode_i++)
                encodings.push_back(m.encoding(encode_i));
            
        }
        
        void Encoding::Set::load(const std::string& location) {
            
            MR::Math::Matrix<double> encodings_matrix(location);
            set(encodings_matrix);
            
        }
        
        void Encoding::Set::set(const MR::Math::Matrix<double>& encodings_matrix) {
            
            encodings.clear();
            
            if (encodings_matrix.columns() != 4)
                throw Exception(
                        "Encoding matrix must have exactly 4 columns (" + str(
                                encodings_matrix.columns())
                        + " found).");
            
            for (size_t encode_i = 0; encode_i < encodings_matrix.rows(); encode_i++)
                encodings.push_back(Encoding(encodings_matrix.row(encode_i)));
            
        }
        
        std::ostream& operator<<(std::ostream& stream, const Encoding::Set& e_set) {
            
            for (std::vector<Encoding>::const_iterator encode_it = e_set.encodings.begin();
                    encode_it != e_set.encodings.end(); encode_it++)
                stream << *encode_it;
            
            return stream;
            
        }
    
    }

}
