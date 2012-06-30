/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 15/09/2010.

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

#include "bts/mcmc/proposal/momentum/weighted.h"
#include "math/cholesky.h"

namespace BTS {

	namespace MCMC {

	  namespace Proposal {


	    Momentum::Weighted::Weighted(Distribution* const proposal_distribution,
                                   const MR::Math::Vector<double>& relative_step_sizes)
        : Momentum(proposal_distribution, relative_step_sizes),
          momen_step(relative_step_sizes.size()) {}




      //!The kinetic energy is the weighted inner product of the momentum, normalized by 1/2 * log(2*pi)^D |weights|.
      double        Momentum::Weighted::log_kinetic_energy(const MR::Math::Matrix<double>& weights_chol) const {

        MR::Math::Vector<double> weights_momen(size());

        MR::Math::Cholesky::solve(weights_momen, weights_chol, momen);
        double energy = MR::Math::dot(weights_momen, momen) / 2.0;

        energy += log_normalization(MR::Math::Cholesky::det_from_decomp(weights_chol));

        return energy;

      }



      //!The kinetic energy is the weighted inner product of the momentum, normalized by 1/2 * log(2*pi)^D |weights|.
      double             Momentum::Weighted::log_normalization(double inv_weights_determinant) const {

        return 0.5 * (LOG_2 + LOG_PI * size() +  MR::Math::log(inv_weights_determinant));

      }

      //!At the start of each leapfrog chain the momentum is randomized to have the a normal distribution with covariance matrix equal to the weighting matrix.
      void               Momentum::Weighted::randomize(const MR::Math::Matrix<double>& cholesky) {

        MR::Math::Vector<double> random_vec(size());

        for (size_t elem_i = 0; elem_i < size(); ++elem_i)
          random_vec[elem_i] = prop_distr->sample(0.0, 1.0);


//        MR::Math::Vector<double> random_vec("/home/tclose/data/mcmc/state/random_vec.txt");

        MR::Math::Matrix<double> lower_cholesky = cholesky;

        for (size_t elem_i = 0; elem_i < size(); ++elem_i)
          for (size_t elem_i2 = elem_i+1; elem_i2 < size(); ++elem_i2)
            lower_cholesky(elem_i, elem_i2) = 0.0;


//        MR::Math::mult(momen, CblasNoTrans, CblasLower, cholesky, random_vec);

        MR::Math::mult(momen, lower_cholesky, random_vec);



      }


      void               Momentum::Weighted::update_state(MR::Math::Vector<double>& state, const MR::Math::Matrix<double>& weights_chol, double time_direction) const {


        MR::Math::Vector<double> working(state.size());
        MR::Math::Cholesky::solve(working, weights_chol, momen);

        for (size_t elem_i = 0; elem_i < size(); elem_i++)
          state[elem_i] += time_direction * working[elem_i] * step[elem_i];


      }

      double             Momentum::Weighted::predicted_change(const MR::Math::Vector<double>& gradient, const MR::Math::Matrix<double>& weights_chol, double time_direction) const {

        double pred_change = 0.0;

        MR::Math::Vector<double> working(gradient.size());
        MR::Math::Cholesky::solve(working, weights_chol, momen);

        for (size_t elem_i = 0; elem_i < size(); elem_i++)
          pred_change += time_direction * working[elem_i] * gradient[elem_i] * step[elem_i];

        return pred_change;

      }


	  }

	}

}
