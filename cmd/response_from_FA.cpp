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

#include "bts/diffusion/model.h"

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

    DIFFUSION_PARAMETERS,

    Option("b_value", "The b value used to generate the response vector")
    + Argument().type_float(SMALL_FLOAT, 3000, LARGE_FLOAT),

    Option()};

using MR::Math::Vector;

EXECUTE {
    
    SET_DIFFUSION_PARAMETERS;

    double b_value = 3000.0;
    Options opt = get_options("b_value");
    if (opt.size())
        b_value = opt[0][0];

    MR::Math::Vector<double> sh = Diffusion::Model::tensor_m0_SH(diff_adc, diff_fa, b_value);
    std::cout << sh << std::endl;
}
