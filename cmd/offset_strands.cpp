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

#include "bts/common.h"

#include "bts/cmd.h"
#include "progressbar.h"

#include "bts/fibre/track.h"
#include "bts/fibre/track/set.h"

#include "bts/inline_functions.h"

using namespace FTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Trims paths to a sphere of specified radius.",
    "The sections of paths that lie outside the radius of the sphere are removed and sections that reenter are split into new paths.",
    NULL
};

ARGUMENTS= {
    Argument ("input", "The tracks to be offset.").type_file (),
    Argument ("output", "The returned offset tracks.").optional().type_file(),
    Argument()
};

OPTIONS= {
    Option ("num_points", "The number of points that will be generated along the strand location")
    + Argument ("num_points", "The number of points that will be generated along the strand location").type_integer (1, 100, 2000),

    Option ("offset", "Offset that will be applied to the tracks")
    + Argument ("offset", "Offset that will be applied to the tracks"),

    Option()};

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location;
        
        if (argument.size() > 1)
            output_location = argument[1].c_str();
        else
            output_location = input_location;
        
        size_t num_points = 0;
        Triple<double> offset = Triple<double>::Zeros;
        
        Options opt = get_options("num_points");
        if (opt.size())
            num_points = opt[0][0];
        
        opt = get_options("offset");
        if (opt.size())
            offset = parse_triple<double>(std::string(opt[0][0]));
        
        FTS::Fibre::Track::Set tcks(input_location, num_points);
        
        for (size_t tck_i = 0; tck_i < tcks.size(); tck_i++)
            for (size_t point_i = 0; point_i < tcks[tck_i].size(); point_i++)
                tcks[tck_i][point_i] += offset;
        
        tcks.set_extend_prop("origin offset", str(-offset));
        
        tcks.save(output_location);
        
        //MR::ProgressBar::done();
    }
    
