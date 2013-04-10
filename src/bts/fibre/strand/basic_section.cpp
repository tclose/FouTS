/*
 *     Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia
 *
 *     Written by Thomas G Close, 5/05/09.
 *
 *     This file is part of Bayesian Tractlet Sampling (BTS).
 *
 *     BTS is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     BTS is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with BTS.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "bts/fibre/strand.h"
#include "bts/fibre/strand/basic_section.h"

namespace BTS {
    
    namespace Fibre {
        
        const char* Strand::BasicSection::INTENSITY_PROP = "intensity";
        const char* Strand::BasicSection::PROPS_LIST[] = { Strand::BasicSection::INTENSITY_PROP,
                                                           PROPS_LIST_END };
        const std::vector<const char*> Strand::BasicSection::PROPS(&(PROPS_LIST[0]),
                &(PROPS_LIST[1]));
        
        std::ostream& operator<<(std::ostream& stream,
                                 const BTS::Fibre::Strand::BasicSection& seg) {
            
            stream << "position:      " << seg.position() << std::endl;
            stream << "tangent:       " << seg.tangent() << std::endl;
            stream << "intensity:     " << seg.intensity() << std::endl;
            
            return (stream);
            
        }
    
    }

}
