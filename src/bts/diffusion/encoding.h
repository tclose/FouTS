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

#ifndef __bts_diffusion_encoding_h__
#define __bts_diffusion_encoding_h__



namespace BTS {

  namespace Diffusion {
  
    class Encoding;
  
  }
  
}

#include <vector>

#include "math/matrix.h"

#include "bts/triple.h"

#include "bts/common.h"

  


namespace BTS {

  namespace Diffusion {

    class Encoding {

      public:

        class Set;
    
      protected:

        double b;      

        Coord orient;
        
      public:
      
        Encoding(const MR::Math::Vector<double>::View& encoding) 
          : b(encoding[3]), orient(encoding[X], encoding[Y], encoding[Z]) {
          normalise_orient();
        }
          

        Encoding(const Triple<double>& orientation, double b_value)
          : b(b_value), orient(orientation) { normalise_orient(); }


        Encoding(const Triple<double>& direction)
          : b(direction.norm()), orient(direction) { normalise_orient(); }
                                        

        Encoding(const Encoding& e) 
          : b(e.b), orient(e.orient) {}
          
          
        Encoding&                       operator= (const Encoding& e) 
          { b = e.b; orient = e.orient; return *this; }
          
          
        double&                         operator[] (size_t dim_i)
          { return orient[dim_i]; }

        const double&                   operator[] (size_t dim_i) const
          { return orient[dim_i]; }

        Coord                           orientation() const
          { return orient; }
          
          
        double                          b_value() const
          { return b; }
        

        static std::vector<Encoding>    parse_encodings(const MR::Math::Matrix<double>::View& encodings_matrix);
        
        
        static MR::Math::Matrix<double> encodings_matrix(const std::vector<Encoding>& encodings);
    
      protected:

        /*! Normalises the orientation, checking to see if the norm of the orientation vector isn't zero first */
        void normalise_orient() {

          double norm = orient.norm();

          //Check to see if the orientation vector isn't all zero before normalising.
          if (norm) {
            orient[X] /= norm;
            orient[Y] /= norm;
            orient[Z] /= norm;
          }

        }


      friend std::ostream&              operator<< (std::ostream& stream, const Encoding& e);

    };
  
  
    inline std::vector<Encoding>        parse_encodings(const MR::Math::Matrix<double>::View& encodings_matrix) {
    
      std::vector<Encoding> encodings;
      
      for (size_t encode_i = 0; encode_i < encodings_matrix.rows(); ++encode_i)
        encodings.push_back(Encoding(encodings_matrix.row(encode_i)));
        
      return encodings;  
      
    }


    inline MR::Math::Matrix<double>          encodings_matrix(const std::vector<Encoding>& encodings) {
    
      MR::Math::Matrix<double> encodings_matrix(encodings.size(), 4);
      
      for (size_t encode_i = 0; encode_i < encodings.size(); ++encode_i) {
        encodings_matrix(encode_i, X) = encodings[encode_i].orientation()[X];
        encodings_matrix(encode_i, Y) = encodings[encode_i].orientation()[Y];
        encodings_matrix(encode_i, Z) = encodings[encode_i].orientation()[Z];        
        encodings_matrix(encode_i, 3) = encodings[encode_i].b_value();
      }
      
      return encodings_matrix;  
      
    }


    inline std::ostream&                operator<< (std::ostream& stream, const Encoding& e) {

      stream << "b=" << e.b << ", orient=" << e.orient << std::endl;

      return stream;

    }
  
  }
  
}


#endif
