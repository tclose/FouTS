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


#ifndef __bts_analysis_hessiantester_h__
#define __bts_analysis_hessiantester_h__


#include "bts/mcmc/state/tensor.h"
#include "bts/triple.h"



namespace BTS {

  namespace Analysis {


    template <typename Object, typename State> class HessianTester {
    
    
      public:
        
        //typedef the pointer to the function to test as 'Function'.
        typedef double (Object::*Function)(const State&, State&, typename State::Tensor&);
      
      protected:
      
        Object                  *object;
        
      public:

        HessianTester(Object& object)
          : object(&object) {}

          
        ~HessianTester() {}
        
        void                          test(Function function, State& state, double step_size, typename State::Tensor& analytic_hessian, typename State::Tensor& numeric_hessian);

    };
    
    
  }

}



#include "progressbar.h"
#include "math/vector.h"
#include "math/matrix.h"


namespace BTS {

  namespace Analysis {

    template <typename Object, typename State> void     HessianTester<Object, State>::test(Function function, State& state, double step_size, typename State::Tensor& analytic_hessian, typename State::Tensor& numeric_hessian) {

      State perturbed_state(state), gradient(state), perturbed_gradient(state);

      gradient.invalidate();

      typename State::Tensor dummy_hessian(state);

      (analytic_hessian = typename State::Tensor(state)).invalidate();

      (*object.*function)(state, gradient, analytic_hessian);

      MR::Math::Vector<double>& state_vector = state;

      numeric_hessian.invalidate();

      MR::ProgressBar progress_bar ("Testing hessian calculations...", state_vector.size());

      for (size_t elem_i = 0; elem_i < state_vector.size(); ++elem_i) {

        perturbed_gradient.invalidate();

        state_vector[elem_i] += step_size;

        (*object.*function)(state, perturbed_gradient, dummy_hessian);

        state_vector[elem_i] -= step_size;

        perturbed_gradient -= gradient;
        perturbed_gradient /= step_size;

        numeric_hessian.row(elem_i) = perturbed_gradient;

        ++progress_bar;

      }

    }


  }

}



#endif
