/*
 *     Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia
 *
 *     Written by Thomas G Close, 5/05/09.
 *
 *     This file is part of Fourier Tract Sampling (FouTS).
 *
 *     FouTS is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     FouTS is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with FTS.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "bts/fibre/strand.h"
#include "bts/fibre/strand/section.h"

namespace FTS {
    
    namespace Fibre {
        
        const std::string Strand::Section::FILE_EXTENSION = "ssc";
        const Coord Strand::Section::FILE_SEPARATOR = Coord(-INFINITY, NAN, NAN);
    
    }

}
