/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Sep 27, 2010.

 This file is part of MRtrix.

 MRtrix is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 MRtrix is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

 */

#include "bts/cmd.h"

#include "bts/common.h"

#include "bts/inline_functions.h"
#include "bts/fibre/track/set.h"

using namespace FTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

const size_t NUM_LENGTH_SECTIONS_DEFAULT = 4;

DESCRIPTION = {
    "remove_short_tracks",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("input_location","").type_file (),
    Argument ("output_location", "").optional().type_file (),
    Argument()
};

OPTIONS= {

    Option("num_points","The number of points below which a track will be removed.")
    + Argument("num_points","").type_integer(0,NUM_LENGTH_SECTIONS_DEFAULT,LARGE_INT),

    Option()};

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location;
        
        if (argument.size() == 2)
            output_location = argument[1].c_str();
        else
            output_location = input_location;
        
        size_t num_points = NUM_LENGTH_SECTIONS_DEFAULT;
        
        Options opt;
        
        opt = get_options("num_points");
        if (opt.size())
            num_points = opt[0][0];
        
        Fibre::Track::Set in_tracks(input_location);
        Fibre::Track::Set out_tracks(in_tracks.get_extend_props());
        out_tracks.add_extend_elem_props(in_tracks);
        
        for (size_t tck_i = 0; tck_i < in_tracks.size(); ++tck_i)
            if (in_tracks[tck_i].size() >= num_points)
                out_tracks.push_back(in_tracks[tck_i], in_tracks.get_extend_elem_prop_row(tck_i));
        
        out_tracks.save(output_location);
        
    }
