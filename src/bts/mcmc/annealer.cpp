/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 24/01/2011.

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

#include "bts/mcmc/annealer.h"

namespace BTS {

	namespace MCMC {


    const double Annealer::START_FRACTION_DEFAULT = 0.01;


    Annealer::Annealer(size_t num_steps, double start_fraction) {

      t = MR::Math::log(start_fraction);
      t_inc = - t / (double)num_steps;

    }

    double Annealer::factor() const {

      return MR::Math::exp(t);

    }

    void   Annealer::increment() {

      t += t_inc;

    }


	}

}
