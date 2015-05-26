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

#include "bts/math/odr.h"

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
//  Argument ("input", "The tracks to be offset.").type_file (),
//  Argument ("output", "The returned offset tracks.").optional().type_file(),
    Argument()
};

OPTIONS= {

    Option ("degree", "The degree of the merged strand.")
    + Argument ("degree", "").type_integer (1, 4, LARGE_INT),

    Option()};

EXECUTE {
    
//  std::string input_location = argument[0];
//  std::string output_location;
//
//  if (argument.size() > 1)
//    output_location = argument[1];
//  else
//    output_location = input_location;
    
        size_t degree = 4;
        
        Options opt = get_options("degree");
        if (opt.size())
            degree = opt[0][0];
        
        Fibre::Strand::Set strands("/home/tclose/thesis/MCMC/data/merge/both.str");
        Fibre::Strand::Set combined;
        
        MR::ProgressBar pb("Merging...");
        
        Fibre::Strand merged = strands[0].merge(strands[1], 100, degree);
        
        combined.push_back(merged);
        
        combined.save("/home/tclose/thesis/MCMC/data/merge/merged.str");
        
        //Generate ODR line
        Fibre::Track::Set tcks = strands.to_tracks(10);
        
        std::vector<Coord> points;
        
        for (size_t tck_i = 0; tck_i < 2; tck_i++)
            for (size_t point_i = 0; point_i < tcks[tck_i].size(); ++point_i)
                points.push_back(tcks[tck_i][point_i]);
        
        Math::Line odr_line = Math::odr(points);
        
        std::ofstream fout("/home/tclose/thesis/MCMC/data/merge/odr.txt");
        
        fout << odr_line.centre[0] << " " << odr_line.centre[1] << " " << odr_line.centre[2] << " "
             << std::endl;
        fout << odr_line.direction[0] << " " << odr_line.direction[1] << " "
             << odr_line.direction[2] << " " << std::endl;
        
    }
    
