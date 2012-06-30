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


#ifndef __bts_analysis_fishergradientmatrixtester_h__
#define __bts_analysis_fishergradientmatrixtester_h__


#include "bts/mcmc/state/tensor.h"
#include "bts/triple.h"



namespace BTS {

  namespace Analysis {


    template <typename Object, typename State> class FisherGradientMatrixTester {
    
    
      public:
        
        //typedef the pointer to the function to test as 'Function'.
        typedef double (Object::*Function)(const State&, State&, typename State::Tensor&, std::vector<typename State::Tensor>&);
      
      protected:
      
        Object                  *object;
        
      public:

        FisherGradientMatrixTester(Object& object)
          : object(&object) {}

          
        ~FisherGradientMatrixTester() {}
        
        void                          test(Function function, State& state, double step_size, std::vector<typename State::Tensor>& analytic_fisher_gradient, std::vector<typename State::Tensor>& numeric_fisher_gradient);

    };
    
    
  }

}



#include "progressbar.h"
#include "math/vector.h"
#include "math/matrix.h"


namespace BTS {

  namespace Analysis {

    template <typename Object, typename State> void     FisherGradientMatrixTester<Object, State>::test(Function function, State& state, double step_size, std::vector<typename State::Tensor>& analytic_fisher_gradient, std::vector<typename State::Tensor>& numeric_fisher_gradient) {

      State dummy_gradient(state);
      typename State::Tensor fisher(state), perturbed_fisher(state);
      std::vector<typename State::Tensor> dummy_fisher_gradient;

      fisher.invalidate();

      numeric_fisher_gradient.clear();

      (*object.*function)(state, dummy_gradient, fisher, analytic_fisher_gradient);

      MR::Math::Vector<double>& state_vector = state;

      MR::ProgressBar progress_bar ("Testing fisher_gradient calculations...", state_vector.size());

      for (size_t elem_i = 0; elem_i < state_vector.size(); ++elem_i) {

        perturbed_fisher.invalidate();

        state_vector[elem_i] += step_size;

        (*object.*function)(state, dummy_gradient, perturbed_fisher, dummy_fisher_gradient);

        state_vector[elem_i] -= step_size;

        numeric_fisher_gradient.push_back(perturbed_fisher);
        numeric_fisher_gradient.back() -= fisher;
        numeric_fisher_gradient.back() /= step_size;

        ++progress_bar;

      }

    }


  }

}



#endif
