/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, Sep 27, 2010.

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

#include "math/matrix.h"
#include "math/cholesky.h"
#include "bts/mcmc/naninf_exception.h"

#include "bts/math/common.h"

#include "bts/mcmc/proposal/momentum/weighted/non_separable.h"

namespace BTS {

	namespace MCMC {

	  namespace Proposal {


//      double        Momentum::Weighted::NonSeparable::log_kinetic_energy(const MR::Math::Matrix<double>& weights_chol) const {
//
//        MR::Math::Vector<double> weights_momen(size());
//
//        MR::Math::Cholesky::solve(weights_momen, weights_chol, momen);
//        double energy = MR::Math::dot(weights_momen, momen) / 2.0;
//
//        energy += log_normalization(MR::Math::Cholesky::det_from_decomp(weights_chol));
//
//        return energy;
//
//      }


      double             Momentum::Weighted::NonSeparable::predicted_change(const MR::Math::Vector<double>& gradient, const MR::Math::Matrix<double>& fisher_chol, double time_direction) const {

        double pred_change = 0.0;

        MR::Math::Vector<double> working(gradient.size());

        MR::Math::Cholesky::solve(working, fisher_chol, momen);

        for (size_t elem_i = 0; elem_i < size(); elem_i++)
          pred_change += working[elem_i] * gradient[elem_i] * step[elem_i];

        return pred_change * time_direction;

      }

	  }

	}

}
