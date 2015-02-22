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

#include "bts/cmd.h"

#include "progressbar.h"

#include "bts/common.h"

#include "bts/fibre/track/set.h"
#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"

#include "bts/inline_functions.h"

using namespace BTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Peels fibres.",
    "Removes all fibres that do not pass withinness an inner 'inclusion radius'.",
    NULL
};

ARGUMENTS= {
    Argument ("input", "The paths to be peeled.").type_file (),
    Argument ("output", "The returned peeled paths.").optional().type_file(),
    Argument()
};

OPTIONS= {
    Option ("distance", "The inner radius the paths must pass through to be included")
    + Argument ("distance", "The inner radius the paths must pass through to be included").type_float (SMALL_FLOAT, 0.8, LARGE_FLOAT),

    Option ("num_length_sections", "The number of points that will be generated along the strand location if the loaded paths are strand descriptors")
    + Argument ("num_length_sections", "").type_integer (1, 100, 2000),

    Option ("num_width_sections", "Number of width sections, only relevant for tractlet class.")
    + Argument ("num_width_sections", "").type_integer (1, 3, 200),

    Option("cube", "Peel to cube instead of sphere."),

    Option("by_bundles", "Peel to cube instead of sphere."),

    Option()};

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location;
        
        if (argument.size() > 1)
            output_location = argument[1].c_str();
        else
            output_location = input_location;
        
        double distance = 0.0;
        bool cube = false;
        bool by_bundles = false;
        size_t num_length_sections = 100;
        size_t num_width_sections = 6;
        
        Options opt = get_options("distance");
        if (opt.size())
            distance = opt[0][0];
        
        opt = get_options("cube");
        if (opt.size())
            cube = true;
        
        opt = get_options("by_bundles");
        if (opt.size())
            by_bundles = true;

        opt = get_options("num_length_sections");
        if (opt.size())
            num_length_sections = opt[0][0];
        
        opt = get_options("num_width_sections");
        if (opt.size())
            num_width_sections = opt[0][0];
        
        if (File::extension(input_location) != File::extension(output_location))
            throw Exception(
                    "Input ('" + input_location + "') and output_location ('" + output_location
                    + "') extensions need to match.");
        
//----------------------//
//  Load original paths
//----------------------//
        
        if (File::has_extension<Fibre::Strand>(input_location) || File::has_extension<Fibre::Track>(
                    input_location)) {
            
            Fibre::Track::Set all_tcks(input_location, num_length_sections);
            
            MR::ProgressBar progress_bar(
                    "Removing tracks/strands that do not pass through centre of sphere...");
            
            Fibre::Track::Set include_tcks;
            
            include_tcks = all_tcks.peel(distance, cube, by_bundles);
            
            include_tcks.save(output_location);
            
            //MR::ProgressBar::done();
            
        } else if (File::has_extension<Fibre::Tractlet>(input_location)) {
            
            Fibre::Tractlet::Set all_tractlets(input_location);
            
            MR::ProgressBar progress_bar(
                    "Removing tractlets that do not pass through centre of sphere...");
            
            Fibre::Tractlet::Set include_tractlets;
            
            include_tractlets = all_tractlets.peel(distance, num_length_sections,
                    num_width_sections, cube);
            
            include_tractlets.save(output_location);
            
            //MR::ProgressBar::done();
            
        } else if (File::has_extension<Fibre::Strand::Set>(input_location)) {
            
            Fibre::Strand::Set::Reader reader(input_location);
            Fibre::Strand::Set::Writer writer(output_location, reader, reader.get_extend_props());
            
            Fibre::Strand::Set all_strands;
            std::map<std::string, std::string> properties_row;
            
            MR::ProgressBar progress_bar(
                    "Removing strands that do not pass through centre of sphere...");
            
            while (reader.next(all_strands)) {
                writer.append(all_strands.peel(distance, num_length_sections, cube));
                progress_bar++;
            }
            
            reader.close();
            writer.close();
            
        } else if (File::has_extension<Fibre::Tractlet::Set>(input_location)) {
            
            Fibre::Tractlet::Set::Reader reader(input_location);
            Fibre::Tractlet::Set::Writer writer(output_location, reader, reader.get_extend_props());
            
            Fibre::Tractlet::Set all_tractlets;
            std::map<std::string, std::string> properties_row;
            
            MR::ProgressBar progress_bar(
                    "Removing tractlets that do not pass through centre of sphere...");
            
            while (reader.next(all_tractlets)) {
                writer.append(
                        all_tractlets.peel(distance, num_length_sections, num_width_sections,
                                cube));
                progress_bar++;
            }
            
            reader.close();
            writer.close();
            
        } else if (File::has_extension<Fibre::Track::Set>(input_location)) {
            
            Fibre::Track::Set::Reader reader(input_location);
            Fibre::Track::Set::Writer writer(output_location, reader, reader.get_extend_props());
            
            Fibre::Track::Set all_tcks;
            std::map<std::string, std::string> properties_row;
            
            MR::ProgressBar progress_bar(
                    "Removing tcks that do not pass through centre of sphere...");
            
            while (reader.next(all_tcks)) {
                writer.append(all_tcks.peel(distance, cube));
                progress_bar++;
            }
            
            reader.close();
            writer.close();
            
        } else
            throw Exception("Unrecognised extension '" + File::extension(input_location) + "'.");
        
    }
    
