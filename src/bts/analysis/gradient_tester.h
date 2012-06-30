/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia
 
 Created by Tom Close on 13/03/09.
 
 This file is part of Bayesian Tractlet Sampling (BTS).
 
 BTS is free software: you can reobjectibute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 BTS is objectibuted in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with BTS.  If not, see <http://www.gnu.org/licenses/>.
 
 */


#ifndef __bts_analysis_gradienttester_h__
#define __bts_analysis_gradienttester_h__


namespace BTS {

  namespace Analysis {
  
  
    template <typename Object, typename State> class GradientTester {
    
    
      public:
        
        //typedef the pointer to the function to test as 'Function'.
        typedef double (Object::*Function)(const State, State);
      
      protected:
      
        Object        *object;
        
      public:

        GradientTester(Object& object)
          : object(&object) {}
          
        ~GradientTester() {}
        
        void                          test(Function function, State state, double step_size, State analytic_gradient, State numeric_gradient);

    
    };
    
    
  }

}



#include "progressbar.h"
#include "math/vector.h"


namespace BTS {

  namespace Analysis {


    template <typename Object, typename State> void     GradientTester<Object, State>::test(Function function, State state, double step_size, State analytic_gradient, State numeric_gradient) {

      analytic_gradient = state;
      analytic_gradient.invalidate();

      double px = (*object.*function)(state, analytic_gradient);

      numeric_gradient = state;
      numeric_gradient.invalidate();

      MR::Math::Vector<double>& state_vector = state;
      MR::Math::Vector<double>& gradient_vector = numeric_gradient;

      State dummy_gradient(state);

      MR::ProgressBar progress_bar ("Testing gradient calculations...", state_vector.size());

      for (size_t elem_i = 0; elem_i < state_vector.size(); ++elem_i) {

        state_vector[elem_i] += step_size;

        gradient_vector[elem_i] = ((*object.*function)(state, dummy_gradient) - px) / step_size;

        state_vector[elem_i] -= step_size;

        ++progress_bar;

      }


    }


  }

}



#endif
