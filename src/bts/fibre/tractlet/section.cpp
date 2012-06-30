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

#include "bts/fibre/tractlet.h"

#include "bts/fibre/tractlet/section.h"



namespace BTS {

  namespace Fibre {

    const std::string            Tractlet::Section::FILE_EXTENSION = "tsc";
    const Coord                  Tractlet::Section::FILE_SEPARATOR = Coord (-INFINITY, NAN, NAN);


    std::ostream&         operator<< (std::ostream& stream , const BTS::Fibre::Tractlet::Section& section) {

      stream << (const Fibre::Strand::BasicSection&)section;

      stream << "width1:         " << section.width1() << std::endl;
      stream << "witdh2:         " << section.width2() << std::endl;

      return (stream);

    }



  }
}
