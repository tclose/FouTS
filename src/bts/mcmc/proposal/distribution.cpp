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

#include "bts/mcmc/proposal/distribution.h"
#include "bts/mcmc/proposal/distribution/gaussian.h"
#include "bts/common.h"

#include "bts/mcmc/inline_functions.h"

namespace BTS {

  namespace MCMC {

    namespace Proposal {

      const char*              Distribution::TYPE_DEFAULT                      = "gaussian";


      Distribution*             Distribution::factory(const std::string& type, gsl_rng* rand_gen) {

        Proposal::Distribution* distribution;

        if (type == "gaussian")

          distribution = new Proposal::Distribution::Gaussian (rand_gen);

        else
          throw Exception ("Unrecognised proposal type '" + type + "' passed to option '-prop_distr_type'.");

        return distribution;

      }

//      const double               PROPOSAL_DEGREE_ARRAY_DEFAULT[] = {1.0, 1.5, 0.2, 0.01, 0.002, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0};
//      const std::vector<double>  Distribution::DEGREE_DEFAULTS (PROPOSAL_DEGREE_ARRAY_DEFAULT, &(PROPOSAL_DEGREE_ARRAY_DEFAULT[12]));
//
//
//      const double               PROPOSAL_SECONDARY_DEGREE_ARRAY_DEFAULT[] = {0.2, 0.1, 0.1, 0.01, 0.002, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0};
//      const std::vector<double>  Distribution::SECONDARY_DEGREE_DEFAULTS (PROPOSAL_SECONDARY_DEGREE_ARRAY_DEFAULT, &(PROPOSAL_SECONDARY_DEGREE_ARRAY_DEFAULT[12]));


    }

  }

}
