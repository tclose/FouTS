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


#ifndef __bts_analysis_rank3hessiantester_h__
#define __bts_analysis_rank3hessiantester_h__

#include <vector>


#include "bts/mcmc/state/tensor.h"
#include "bts/triple.h"




namespace BTS {

  namespace Analysis {


    template <typename Object, typename State> class Rank3HessianTester {
    
    
      public:
        
        //typedef the pointer to the function to test as 'Function'.
        typedef double (Object::*Function)(const State&, State&, typename State::Tensor&, std::vector<typename State::Tensor>&);
      
      protected:
      
        Object                  *object;
        
      public:

        Rank3HessianTester(Object& object)
          : object(&object) {}

          
        ~Rank3HessianTester() {}
        
        void                          test(Function function, State& state, double step_size, std::vector<typename State::Tensor>& analytic_rank3hessian, std::vector<typename State::Tensor>& numeric_rank3hessian);

    };
    
    
  }

}



#include "progressbar.h"
#include "math/vector.h"
#include "math/matrix.h"


namespace BTS {

  namespace Analysis {

    template <typename Object, typename State> void     Rank3HessianTester<Object, State>::test(Function function, State& state, double step_size, std::vector<typename State::Tensor>& analytic_rank3hessian, std::vector<typename State::Tensor>& numeric_rank3hessian) {

      State perturbed_state(state), dummy_gradient(state);

      typename State::Tensor hessian(state), perturbed_hessian(state);

      hessian.invalidate();
      perturbed_hessian.invalidate();

      std::vector<typename State::Tensor> dummy_rank3hessian;

      (*object.*function)(state, dummy_gradient, hessian, analytic_rank3hessian);

      MR::Math::Vector<double>& state_vector = state;

      MR::Math::Matrix<double> template_rank3hessian_matrix(state_vector.size(), state_vector.size());
      template_rank3hessian_matrix = NAN;


      MR::ProgressBar progress_bar ("Testing rank3hessian calculations...", state_vector.size());

      for (size_t elem_i = 0; elem_i < state_vector.size(); ++elem_i) {

        perturbed_hessian.invalidate();

        state_vector[elem_i] += step_size;

        perturbed_state.from_vector(state_vector);

        (*object.*function)(perturbed_state, dummy_gradient, perturbed_hessian, dummy_rank3hessian);

        state_vector[elem_i] -= step_size;

        MR::Math::Matrix<double>& perturbed_hessian_matrix = perturbed_hessian;

        perturbed_hessian_matrix -= hessian;
        perturbed_hessian_matrix /= step_size;

        numeric_rank3hessian.push_back(perturbed_hessian_matrix);

        ++progress_bar;

      }

    }


  }

}



#endif
