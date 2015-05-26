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

#include "bts/cmd.h"
#include "progressbar.h"

#include "bts/common.h"

#include "bts/fibre/track/set.h"

#include "phantom/interface.h"

#include "phantom/shared/strand_collection.h"
#include "phantom/resample/resample.h"

#include "bts/inline_functions.h"

using namespace FTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Resamples control points along tracks.",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("input", "The tracks to be resampled.").type_file (),
    Argument ("output", "The returned resampled tracks.").optional().type_file(),
    Argument()
};

OPTIONS= {

    Option ("step_size", "Defines the length of the new inter control-point intervals")
    + Argument ("step_size", "").type_float (0.00001, 0.02, 10000.0),

    Option ("num_points", "The number of points that will be generated along the strand path if the loaded tracks are fourier descriptors")
    + Argument ("num_points", "").type_integer (1, 100, 2000),

    Option()};

EXECUTE {
    
        std::string input_path = argument[0];
        std::string output_path;

        if (argument.size() > 1)
            output_path = argument[1].c_str();
        else
            output_path = input_path;
        
        double step_size = 0.02;
        size_t num_points = 0;
        
        Options opt = get_options("step_size");
        if (opt.size())
            step_size = opt[0][0];
        
        opt = get_options("num_points");
        if (opt.size())
            num_points = opt[0][0];
        
        FTS::Fibre::Track::Set tcks(input_path, num_points);
        FTS::Fibre::Track::Set resampled_tcks(tcks.get_extend_props());
        
        for (size_t tck_i = 0; tck_i < tcks.size(); ++tck_i)
            resampled_tcks.push_back(tcks[tck_i].resample(step_size));

        resampled_tcks.save(output_path);
        
    }
    
