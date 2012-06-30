/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 01/10/2010.

    This file is part of MRtrix.

    MRtrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRtrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef __bts_prob_test_bayeslogregression_h__
#define __bts_prob_test_bayeslogregression_h__

#include <cmath>
#include <vector>

#include "math/vector.h"
#include "math/matrix.h"
#include "math/math.h"

#include "bts/mcmc/state.h"

namespace BTS {

	namespace Prob {

	  namespace Test {

	    class BayesLogRegression {

	      public:

	        const static char* DATA_LOCATION_DEFAULT;
	        const static size_t POLY_ORDER_DEFAULT;
	        const static float PRIOR_VARIANCE_DEFAULT;

//	      protected:

	        MR::Math::Matrix<double> X;
	        MR::Math::Matrix<double> XX;
	        MR::Math::Vector<double> t;

	        MR::Math::Vector<double> f;
	        MR::Math::Matrix<double> working;
	        MR::Math::Vector<double> v;
	        MR::Math::Vector<double> p;
	        MR::Math::Matrix<double> V;
	        MR::Math::Matrix<double> Z;
	        MR::Math::Matrix<double> Z1;
	        MR::Math::Matrix<double> Z2;

	        size_t poly_order;
	        double alpha;
	        size_t D;
	        size_t N;


	      public:

	        BayesLogRegression(const std::string& location, size_t poly_order = 1, double alpha = 100.0);


	         double log_prob (const MR::Math::Vector<double>& w);


	        double log_prob (const MR::Math::Vector<double>& w, MR::Math::Vector<double>& d_w);


          double log_prob_and_fisher (const MR::Math::Vector<double>& w, MR::Math::Vector<double>& d_w, MR::Math::Matrix<double>& G);


          double log_prob_and_fisher (const MR::Math::Vector<double>& w, MR::Math::Vector<double>& d_w, MR::Math::Matrix<double>& G, std::vector<MCMC::State::Tensor>& d_G);


	    };

	  }

	}

}

#endif /* __bts_prob_test_bayeslogregression_h__ */
