/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Written by Thomas G. Close, 04/03/2009.

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

#include "progressbar.h"

#include "bts/cmd.h"

#include "bts/triple.h"

#include "bts/common.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"

const double RADII_DEFAULT = 0.03;
const size_t NUM_LENGTH_SECTIONS_DEFAULT = 100;

using namespace BTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Export MRtrix tracks file to NFG format.",
    "The output ascii files will consist of 3xN matrices, corresponding to the [ X Y Z ] coordinates of the N points making up the path. All paths will be written to the output folder.",
    NULL
};

ARGUMENTS= {
    Argument ("input", "The strands to be exported.").type_file (),
    Argument ("output", "The output directory where the ASCII files will be written to.").type_file(),
    Argument()
};

OPTIONS= {

    Option ("radii","The radii of the exported tracks.")
    + Argument ("radii","").type_float (SMALL_FLOAT, RADII_DEFAULT, LARGE_FLOAT),

    Option ("num_length_sections","Number of length sections to convert strands into tracks before exporting.")
    + Argument ("num_length_sections","").type_integer (0,NUM_LENGTH_SECTIONS_DEFAULT,LARGE_INT),

    Option()};

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location = argument[1];
        
        double radii = RADII_DEFAULT;
        size_t num_length_sections = NUM_LENGTH_SECTIONS_DEFAULT;
        
        bool use_loaded_radii = true;
        
        Options opt = get_options("radii");
        if (opt.size()) {
            radii = opt[0][0];
            use_loaded_radii = false;
        }
        
        opt = get_options("num_length_sections");
        if (opt.size())
            num_length_sections = opt[0][0];
        
        Fibre::Track::Set tcks(input_location, num_length_sections);
        
        MR::ProgressBar progress_bar("Exporting MRtrix tracks file to NFG format...", tcks.size());
        
        File::mkdir(output_location);
        
        size_t count = 0;
        
        size_t max_bundle_index = 0;
        
        if (tcks.has_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP)) {
            
            for (size_t tck_i = 0; tck_i < tcks.size(); ++tck_i) {
                
                size_t bundle_index = to<size_t>(
                        tcks.get_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, tck_i));
                
                if (bundle_index > max_bundle_index)
                    max_bundle_index = bundle_index;
                
            }
            
        } else
            max_bundle_index = tcks.size();
        
        size_t num_tck_i_dec_places = num_dec_places(tcks.size());
        size_t num_bundle_dec_places = num_dec_places(max_bundle_index);
        
        for (size_t tck_i = 0; tck_i < tcks.size(); ++tck_i) {
            
            size_t bundle_index;
            
            if (tcks.has_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP))
                bundle_index = to<size_t>(
                        tcks.get_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, tck_i));
            else
                bundle_index = tck_i;
            
            double track_radius;
            
            if (tcks.has_extend_elem_prop(Fibre::Track::RADIUS_PROP) && use_loaded_radii)
                track_radius = to<double>(
                        tcks.get_extend_elem_prop(Fibre::Track::RADIUS_PROP, tck_i));
            else
                track_radius = radii;
            
            std::stringstream filename_stream;
            
            filename_stream << output_location << std::setfill('0') << "/strand_"
                            << std::setw(num_tck_i_dec_places) << tck_i << "-"
                            << std::setw(num_bundle_dec_places) << bundle_index << "-r"
                            << track_radius << ".txt";
            std::string filename = filename_stream.str();
            
            std::ofstream out(filename.c_str());
            
            if (!out) {
                throw Exception(
                        "error opening ascii file \"" + filename + "\": " + strerror(errno));
            }
            
            for (size_t point_i = 0; point_i < tcks[tck_i].size(); ++point_i)
                out << tcks[tck_i][point_i][X] << " " << tcks[tck_i][point_i][Y] << " "
                    << tcks[tck_i][point_i][Z] << "\n";
            
            out.close();
            count++;
            ++progress_bar;
            
        }
        
    }
