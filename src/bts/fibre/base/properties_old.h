/*
    Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

    Written by Thomas G Close, 04/06/2010.

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

#ifndef __bts_fibre_properties_h__
#define __bts_fibre_properties_h__

#include <string>
#include <map>


namespace BTS {

	namespace Fibre {

		
    class Properties : public std::map<std::string, std::string> {

      public:

        class Extended;
        class Header;
        template <typename T> class Intrinsic_tpl;
        class Intrinsic;

      public:

        Properties() {}

        Properties(const std::map<std::string, std::string>& properties)
          : std::map<std::string, std::string>(properties) {}


        Properties(const MR::DWI::Tractography::Properties& mr_properties)
          : std::map<std::string, std::string>(mr_properties) {}


        Properties&   operator=(const std::map<std::string, std::string>& properties)
          { std::map<std::string, std::string>::operator=(properties); return *this; }


        ~Properties() {}


    };


	}

}

#endif /* __bts_fibre_properties_h__ */
