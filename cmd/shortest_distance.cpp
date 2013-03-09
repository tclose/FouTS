/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Mar 14, 2011.

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

#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"

#include "bts/inline_functions.h"

using namespace BTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Map colours to reference set by switching bundle indices to closest match.",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("reference","Reference fibres the match fibres will be mapped to.").type_file (),
    Argument ("samples","Fibre samples that stats will collected from.").type_file (),
    Argument ("output","Output filename, if not provided will be same as samples with 'stat' extension.").optional().type_file (),
    Argument()
};

OPTIONS= {

    Option()};

EXECUTE {
    
        std::string reference_location = argument[0];
        std::string samples_location = argument[1];
        std::string output_location = File::strip_extension(samples_location) + ".shrt."
                                      + File::extension(samples_location);
        
        if (argument.size() == 3)
            output_location = argument[2].c_str();
        
        if (File::has_extension<Fibre::Strand>(reference_location)) {
            
            if (!File::has_extension<Fibre::Strand::Set>(samples_location))
                throw Exception(
                        "Samples location ('" + samples_location
                        + "') does not match reference location ('" + reference_location + "').");
            
            Fibre::Strand::Set reference(reference_location);
            
            Fibre::Strand::Set strands;
            
            Fibre::Strand::Set::Reader reader(samples_location);
            Fibre::Strand::Set::Writer writer(output_location, reader, reader.get_extend_props());
            
            MR::ProgressBar progress_bar("Calculating shortest distances...",
                    to<size_t>(reader.get_extend_props()["count"]));
            
            size_t count = 0;
            while (reader.next(strands)) {
                
                if (strands.size() != reference.size())
                    throw Exception(
                            "Size of tractlet set " + str(count) + " (" + str(strands.size())
                            + ") does not match that of 'reference' set (" + str(reference.size())
                            + ")");
                
                Fibre::Strand::Set shortest_dist = strands.smallest_distance_set(reference);
                writer.append(shortest_dist);
                
                ++count;
                ++progress_bar;
                
            }
            
        } else if (File::has_extension<Fibre::Tractlet>(reference_location)) {
            
            if (!File::has_extension<Fibre::Tractlet::Set>(samples_location))
                throw Exception(
                        "Samples location ('" + samples_location
                        + "') does not match reference location ('" + reference_location + "').");
            
            Fibre::Tractlet::Set reference(reference_location);
            
            Fibre::Tractlet::Set tractlets;
            
            Fibre::Tractlet::Set::Reader reader(samples_location);
            Fibre::Tractlet::Set::Writer writer(output_location, reader, reader.get_extend_props());
            
            MR::ProgressBar progress_bar("Calculating shortest distances...",
                    to<size_t>(reader.get_extend_props()["count"]));
            
            size_t count = 0;
            while (reader.next(tractlets)) {
                
                if (tractlets.size() != reference.size())
                    throw Exception(
                            "Size of tractlet set " + str(count) + " (" + str(tractlets.size())
                            + ") does not match that of 'reference' set (" + str(reference.size())
                            + ")");
                
                Fibre::Tractlet::Set shortest_dist = tractlets.smallest_distance_set(reference);
                writer.append(shortest_dist);
                
                ++count;
                ++progress_bar;
                
            }
            
        } else
            throw Exception(
                    "Usupported output file type ('" + output_location
                    + "') for match of type tractlet.");
        
    }
