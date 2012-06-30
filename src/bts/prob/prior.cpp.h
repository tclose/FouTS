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

#include "bts/prob/prior.h"

#include "bts/prob/inline_functions.h"


namespace BTS {

  namespace Prob {

    const char* Prior::TYPE_DEFAULT = "default";
    const double Prior::SCALE_DEFAULT = 1;


    template <typename T> std::map<std::string, std::string>     get_component_values(const T& fibre) {

      std::map<std::string, std::string> component_values;

      component_values["magnitude"] = magnitude.log_prob(fibre);
      component_values["hook"] = hook.log_prob(fibre);
      component_values["density"] = density.log_prob(fibre);

      return component_values;

    }


  }
}
