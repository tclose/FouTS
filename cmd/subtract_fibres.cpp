/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

 This file is part of Bayesian Tractlet Sampling (BTS).

 BTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 BTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with BTS.  If not, see <http://www.gnu.org/licenses/>.

 */

#include "bts/common.h"

#include "bts/cmd.h"
#include "progressbar.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/track/set.h"

#include "bts/inline_functions.h"

using namespace BTS;
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
    Argument ("input", "The fibres to be subtracted from.").type_file (),
    Argument ("subtract", "The fibres to subtract.").type_file (),
    Argument ("output", "The resultant fibres.").optional().type_file(),
    Argument()
};

OPTIONS= {

    Option ("no_acs", "Don't acs in the subtraction."),

    Option()};

EXECUTE {
    
        std::string input_location = argument[0];
        std::string subtract_location = argument[1];
        std::string output_location;
        
        if (argument.size() > 2)
            output_location = argument[2].c_str();
        else
            output_location = input_location;
        
        if (File::extension(input_location) != File::extension(subtract_location))
            throw Exception(
                    "Extensions of input files do not match ('" + File::extension(input_location)
                    + "' and '" + File::extension(input_location) + "').");
        
        if (File::has_extension<Fibre::Strand>(input_location)) {
            
            MR::ProgressBar progress_bar("Subtracting strands...");
            
            BTS::Fibre::Strand::Set input(input_location);
            BTS::Fibre::Strand::Set subtract(subtract_location);
            
            if (input.size() != subtract.size())
                throw Exception("Number of points don't match.");
            
            double base_intensity = input.base_intensity();
            
            for (size_t tck_i = 0; tck_i < input.size(); tck_i++) {
                
                if (input[tck_i].degree() != subtract[tck_i].degree())
                    throw Exception("Degrees don't match.");
                
                double acs = input[tck_i].acs();
                
                input[tck_i] -= subtract[tck_i];
                
                input[tck_i].set_acs(acs);
                
            }
            
            input.set_base_intensity(base_intensity);
            
            input.save(output_location);
            
        } else if (File::has_extension<Fibre::Tractlet>(input_location)) {
            
            MR::ProgressBar progress_bar("Subtracting tractlets...");
            
            BTS::Fibre::Tractlet::Set input(input_location);
            BTS::Fibre::Tractlet::Set subtract(subtract_location);
            
            if (input.size() != subtract.size())
                throw Exception("Number of points don't match.");
            
            double base_intensity = input.base_intensity();
            
            for (size_t tck_i = 0; tck_i < input.size(); tck_i++) {
                
                if (input[tck_i].degree() != subtract[tck_i].degree())
                    throw Exception("Degrees don't match.");
                
                double acs = input[tck_i].acs();
                
                input[tck_i] -= subtract[tck_i];
                
                input[tck_i].set_acs(acs);
                
            }
            
            input.set_base_intensity(base_intensity);
            
            input.save(output_location);
            
        } else if (File::has_extension<Fibre::Track>(input_location)) {
            
            MR::ProgressBar progress_bar("Subtracting tracks...");
            
            BTS::Fibre::Track::Set input(input_location);
            BTS::Fibre::Track::Set subtract(subtract_location);
            
            double base_intensity = input.base_intensity();
            
            if (input.size() != subtract.size())
                throw Exception("Number of points don't match.");
            
            for (size_t tck_i = 0; tck_i < input.size(); tck_i++) {
                
                if (input[tck_i].num_points() != subtract[tck_i].num_points())
                    throw Exception("Number of points don't match.");
                
                input[tck_i] -= subtract[tck_i];
                
            }
            
            input.set_base_intensity(base_intensity);
            
            input.save(output_location);
            
        } else
            throw Exception("Unrecognised extension ('" + File::extension(input_location));
        
    }
