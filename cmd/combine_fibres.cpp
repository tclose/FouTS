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

#include "progressbar.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"
#include "bts/fibre/tractlet/set.h"

#include "bts/inline_functions.h"

using namespace FTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Merge multiple fibres files into one file.",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("input", "The input strands to be merged into the output file.").allow_multiple().type_file (),
    Argument ("output", "The merged strands file in FTS format").type_file(),
    Argument()
};

OPTIONS= {

    Option()};

template<typename T> void combine_fibres(const std::vector<std::string>& input_locations,
                                         const std::string& output_location,
                                         MR::ProgressBar& progress_bar);

template<typename T> void add_fibres(typename T::Writer& writer, const std::string& location,
                                     MR::ProgressBar& progress_bar);

EXECUTE {
    
        std::vector<std::string> input_locations;
        
        if (!argument.size())
            throw Exception("No input_locations provided.");
        
        for (size_t arg_i = 0; arg_i < argument.size() - 1; arg_i++)
            input_locations.push_back(argument[arg_i]);
        
        std::string output_location = argument.back();
        
        MR::ProgressBar progress_bar("Combining fibres ...");
        
        if (File::has_extension<Fibre::Strand>(output_location)) {
            
            combine_fibres<Fibre::Strand>(input_locations, output_location, progress_bar);
            
        } else if (File::has_extension<Fibre::Tractlet>(output_location)) {
            
            combine_fibres<Fibre::Tractlet>(input_locations, output_location, progress_bar);
            
        } else if (File::has_extension<Fibre::Track>(output_location)) {
            
            combine_fibres<Fibre::Track>(input_locations, output_location, progress_bar);
            
        } else if (File::has_extension<Fibre::Strand::Set>(output_location)) {
            
            combine_fibres<Fibre::Strand::Set>(input_locations, output_location, progress_bar);
            
        } else if (File::has_extension<Fibre::Tractlet::Set>(output_location)) {
            
            combine_fibres<Fibre::Tractlet::Set>(input_locations, output_location, progress_bar);
            
        } else if (File::has_extension<Fibre::Track::Set>(output_location)) {
            
            combine_fibres<Fibre::Track::Set>(input_locations, output_location, progress_bar);
            
        } else
            throw Exception("Unknown file extension in output file " + output_location + ".");
        
        //MR::ProgressBar::done();
        
    }
    
    template<typename T> void combine_fibres(const std::vector<std::string>& input_locations,
                                             const std::string& output_location,
                                             MR::ProgressBar& progress_bar) {
        
        std::map<std::string, std::string> properties;
        std::vector<std::string> set_header, elem_header;
        
        typename T::Reader sneak_peek_reader(input_locations[0]);
        
        sneak_peek_reader.close();
        
        const Fibre::Properties& sneak_props = sneak_peek_reader.get_extend_props();
//        const std::string step_size_key = str("step_size");
        Fibre::Properties props;
        if (sneak_props.count("step_size"))
            props["step_size"] = "0.05";

        typename T::Writer writer(output_location, sneak_peek_reader, props);
        
        for (typename std::vector<std::string>::const_iterator input_it = input_locations.begin();
                input_it != input_locations.end(); ++input_it)
            add_fibres<T>(writer, *input_it, progress_bar);
        
    }
    
    template<typename T> void add_fibres(typename T::Writer& writer, const std::string& location,
                                         MR::ProgressBar& progress_bar) {
        
        //If location is directory recursively enter and add all
        if (File::is_dir(location)) {
            
            std::vector<std::string> filenames = File::list_filenames(location, "",
                    T::FILE_EXTENSION);
            
            for (size_t filename_i = 0; filename_i < filenames.size(); filename_i++)
                add_fibres<T>(writer, File::join(location, filenames[filename_i]), progress_bar);
            
        } else {
            
            if (!File::has_extension<T>(location))
                throw Exception("Incompatible file extension found at " + location + ".");
            
            std::map<std::string, std::string> properties;
            std::vector<std::string> set_header, elem_header;
            
            typename T::Reader reader(location);
            
            T fibre;
            std::map<std::string, std::string> row;
            
            while (reader.next(fibre, row))
                writer.append(fibre, row);
            
            progress_bar++;
            
        }
        
    }
