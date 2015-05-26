/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

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

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

#include "bts/cmd.h"
#include "progressbar.h"

#include "bts/common.h"
#include "bts/image/expected/buffer.h"
#include "bts/image/observed/buffer.h"

#include "bts/image/noise.h"
#include "bts/image/noise/gaussian.h"

#include "bts/inline_functions.h"

using namespace FTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Finds the maximum b0 intensity of an image",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("image", "The image.").type_file(),

    Argument()
};

OPTIONS= {


    Option("percentile", "The percentile of the ""maximum"" b=0 value.")
    + Argument("").type_float(-1.01, -1.0, 100.0),

    Option("threshold", "The cut off value above which the percentile is calculated above.")
        + Argument("").type_float(0.0, 0.0, LARGE_FLOAT),

    Option()

};

EXECUTE {
    
    double percentile = -1.0;
    Options opt = get_options("percentile");
    if (opt.size())
        percentile = opt[0][0];

    double threshold = -1.0;
    opt = get_options("threshold");
    if (opt.size())
        threshold = opt[0][0];

    std::string input_location = argument[0];

    Image::Observed::Buffer image(input_location);

    std::cout << image.max_b0(percentile, threshold);

}
    
