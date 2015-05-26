/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Written by Thomas G. Close, 04/03/2009.

 This file is part of Fourier Tract Sampling (FouTS).

 FouTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FouTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FTS.  If not, see <http://www.gnu.org/licenses/>.

 */

#include "bts/cmd.h"

#include "bts/common.h"

#include "bts/fibre/strand/set.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"
#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/tractlet/section.h"
#include "bts/fibre/base/set_reader.h"
#include "bts/fibre/base/set_writer.h"

#include "bts/inline_functions.h"

using namespace FTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Print strands to screen",
    "",
    NULL
};

ARGUMENTS= {
    Argument()
};

OPTIONS= {

    Option()};

using MR::Math::Vector;

EXECUTE {
    
        Vector<double> v(10);
        for (size_t i = 0; i < 10; ++i)
            v[i] = (double) i;
        
        bool sw = false;
        
        Vector<double> v2(sw ? v.view() : v);
        
        std::cout << v2 << std::endl;
        
    }
