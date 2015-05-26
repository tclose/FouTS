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

#include "bts/fibre/strand/set.h"

#include "bts/inline_functions.h"

using namespace FTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Projects file locations to a sphere of supplied radius (defaults to 1.0).",
    "Projects file locations to a sphere of supplied radius (defaults to 1.0).",
    NULL
};

ARGUMENTS= {
    Argument ("input", "The file locations to be projected.").type_file (),
    Argument ("output", "The returned projected file locations.").optional().type_file(),
    Argument()
};

OPTIONS= {
    Option ("sphere_radius", "The sphere radius the file locations will be projected onto")
    + Argument ("sphere_radius", "The sphere radius the file locations will be projected onto").type_float (SMALL_FLOAT, 1.0, LARGE_FLOAT),

    Option ("num_points", "The number of points that will be generated along the strand location")
    + Argument ("points", "The number of points that will be generated along the strand location").type_integer (1, 100, 2000),

    Option ("degree", "The degree of the Strand coefficients used to describe the strands")
    + Argument ("degree", "The degree of the outputted Strand coefficients").type_integer (1, 3, 200),

    Option()};

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location;
        
        if (argument.size() > 1)
            output_location = argument[1].c_str();
        else
            output_location = input_location;
        
        double sphere_radius = 1.0;
        size_t degree = 0;
        size_t num_points = 0;
        
        Options opt = get_options("sphere_radius");
        if (opt.size())
            sphere_radius = opt[0][0];
        
        opt = get_options("num_points");
        if (opt.size())
            num_points = opt[0][0];
        
        opt = get_options("degree");
        if (opt.size())
            degree = opt[0][0];
        
//----------------------//
//  Load original file locations
//----------------------//
        
        Fibre::Strand::Set strands(input_location, degree);
        
        MR::ProgressBar progress_bar(
                "Projecting strands onto a sphere of radius " + str(sphere_radius) + " ...");
        
        strands.project_to_sphere(sphere_radius);
        
        strands.save(output_location, num_points);
        
        //MR::ProgressBar::done();
        
    }
    
