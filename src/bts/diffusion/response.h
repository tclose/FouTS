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


#ifndef __bts_diffusion_response_h__
#define __bts_diffusion_response_h__

namespace BTS {

  namespace Diffusion {
  
    class Response;
    
  }
  
}

#include "bts/common.h"

#include "bts/triple.h"
#include "bts/coord.h"
#include "bts/diffusion/encoding.h"


namespace BTS {

  namespace Diffusion {
  
    class Response : public Encoding {

      protected:
        
        MR::Math::Vector<double>   coeffs;
        size_t index;

        
#ifdef OPTIMISED
//        // Used to store the weighting and the gradient of the weighting for the current Fibre::Strand::Section
//        // to save it being recalculated at every voxel it affects.
//        double precalc_weighting;
//
//        Triple<double> precalc_gradient;
//        Coord::Tensor precalc_hessian;
#endif

        //NB: b value is incorporated into the response coefficients.
        
      public:
        
        Response (const MR::Math::Vector<double>::View& encoding_row, MR::Math::Vector<double>& coeffs, size_t index)
          : Encoding(encoding_row), coeffs(coeffs), index(index) {}

        Response (const Triple<double>& orientation, double b_value, MR::Math::Vector<double>& coeffs, size_t index)
          : Encoding(orientation, b_value), coeffs(coeffs), index(index) {}
        
        Response (const Triple<double>& direction, MR::Math::Vector<double>& coeffs, size_t index)
          : Encoding(direction), coeffs(coeffs), index(index) {}

        Response (const Response& r)
          : Encoding(r), coeffs(r.coeffs), index(r.index) {}
        
        Response&   operator= (const Response& r)
          { Encoding::operator= (r); coeffs = r.coeffs; index = r.index; return *this; }

        
        double       weighting (const Coord& tangent) const;

        double       weighting (const Coord& tangent, Coord& tangent_gradient) const;

        double       weighting (const Coord& tangent, Coord& tangent_gradient, Coord::Tensor& tangent_hessian) const;

//--------------------------------------------------------------------------------------------------//
//These functions are only used to get the signature right for the GradientTester::Function classes.
//--------------------------------------------------------------------------------------------------//

        double       weighting (const Coord& tangent)
          { const Response& const_this = *this; return const_this.weighting(tangent); }

        double       weighting (const Coord& tangent, Coord& tangent_gradient)
          { const Response& const_this = *this; return const_this.weighting(tangent, tangent_gradient); }

        double       weighting (const Coord& tangent, Coord& tangent_gradient, Coord::Tensor& tangent_hessian)
          { const Response& const_this = *this; return const_this.weighting(tangent, tangent_gradient, tangent_hessian); }
//--------------------------//


        size_t         get_index() const
          { return index; }

#ifdef OPTIMISED
//        void        precalculate_weighting(const Triple<double>& tangent)
//          { this->precalc_weighting = weighting(tangent); }
//
//        void        precalculate_weighting_and_gradient(const Triple<double>& tangent)
//          { this->precalc_weighting = weighting(tangent, this->precalc_gradient); }
//
//        void        precalculate_weighting_gradient_and_hessian(const Triple<double>& tangent)
//          { this->precalc_weighting = weighting(tangent, this->precalc_gradient, this->precalc_hessian); }
//
//        double       precalculated_weighting() const
//          { return precalc_weighting; }
//
//
//        const Triple<double>&         precalculated_gradient() const
//          { return precalc_gradient; }
//
//        const Coord::Tensor& precalculated_hessian() const
//          { return precalc_hessian; }

#endif

        void        scale_coeffs(double scalar)
          { this->coeffs *= scalar; }


      
        
      friend std::ostream& operator<< (std::ostream& stream, const Diffusion::Response& response);

    };
    
    
//    inline
    std::ostream& operator<< (std::ostream& stream, const Diffusion::Response& response);


  }
  
  
}


#endif
