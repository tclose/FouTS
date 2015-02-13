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

#include "phantom/interface.h"

#include "phantom/trim/trim.h"
#include "bts/fibre/track.h"
#include "bts/fibre/track/set.h"
#include "bts/fibre/tractlet.h"
#include "bts/fibre/tractlet/set.h"

#include "bts/inline_functions.h"

const size_t NUM_LENGTH_SECTIONS_DEFAULT = 50;
const size_t NUM_WIDTH_SECTIONS_DEFAULT = 40;

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
    Argument ("input", "The fibres to be trimmed.").type_file (),
    Argument ("output", "The returned trimmed fibres.").optional().type_file(),
    Argument()
};

OPTIONS= {

    Option ("degree", "Limit the degree of the Strand coefficients used to describe the strands to the first - degrees")
    + Argument ("degree", "Limit the degree of the Strand coefficients used to describe the strands to the first - degrees").type_integer (1, 200, 3),

    Option ("distance", "The radius of the sphere the strands will be trimmed to")
    + Argument ("new_sphere_r", "").type_float (SMALL_FLOAT, 1.0, LARGE_FLOAT),

    Option ("length_reject_threshold", "The length below which new strands will be rejected")
    + Argument ("length_reject_threshold", "").type_float (SMALL_FLOAT, 1.0, LARGE_FLOAT),

    Option ("num_length_sections", "The number of length sections used in conversion between tractlets or strands and tracks.")
    + Argument ("num_length_sections", "").type_integer (1, NUM_LENGTH_SECTIONS_DEFAULT, LARGE_INT),

    Option ("num_width_sections", "The number of width sections used in conversion between tractlets and tracks.")
    + Argument ("num_width_sections", "").type_integer (1, NUM_WIDTH_SECTIONS_DEFAULT, LARGE_INT),

    Option ("offsets", "an offset by which the fibres should be shifted before they are trimmed to the sphere/cube")
    + Argument("offsets", "").type_text ("[0.0,0.0,0.0]"), \

    Option ("cube", "Trim the set to a cube instead of the default sphere."),

    Option()};

BTS::Fibre::Track::Set trim_tcks(BTS::Fibre::Track::Set& input_tcks, double sphere_r,
                                 double length_reject_threshold, bool trim_to_cube,
                                 const Triple<double>& offsets);

BTS::Fibre::Tractlet::Set trim_tractlets(BTS::Fibre::Tractlet::Set& input_tractlets,
                                         double distance, double length_reject_threshold,
                                         bool trim_to_cube, size_t num_length_sections,
                                         size_t num_width_sections,
                                         const Triple<double>& offsets);

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location;
        
        if (argument.size() > 1)
            output_location = argument[1].c_str();
        else
            output_location = input_location;
        
        size_t num_length_sections = NUM_LENGTH_SECTIONS_DEFAULT;
        size_t num_width_sections = NUM_WIDTH_SECTIONS_DEFAULT;
        size_t degree = 3;
        double distance = 1.0;
        double length_reject_threshold = 0.05;
        bool cube = false;
        Triple<double> offsets = Triple<double>::Zeros;
        
        Options opt = get_options("num_length_sections");
        if (opt.size())
            num_length_sections = opt[0][0];
        
        opt = get_options("num_width_sections");
        if (opt.size())
            num_width_sections = opt[0][0];
        
        opt = get_options("degree");
        if (opt.size())
            degree = opt[0][0];
        
        opt = get_options("distance");
        if (opt.size())
            distance = opt[0][0];
        
        opt = get_options("length_reject_threshold");
        if (opt.size())
            length_reject_threshold = opt[0][0];
        
        opt = get_options("cube");
        if (opt.size())
            cube = true;
        
        opt = get_options("offsets");
        if (opt.size())
          offsets = parse_triple<double>(std::string(opt[0][0]));

        if (File::has_or_txt_extension<BTS::Fibre::Strand>(input_location) || File::has_or_txt_extension<
                    BTS::Fibre::Track>(input_location)) {
            
            BTS::Fibre::Track::Set input(input_location, num_length_sections);
            
            MR::ProgressBar progress_bar("Trimming strands/tracks...");
            
            BTS::Fibre::Track::Set output = trim_tcks(input, distance, length_reject_threshold,
                    cube, offsets);
            
            output.save(output_location, degree);
            
        } else if (File::has_or_txt_extension<BTS::Fibre::Tractlet>(input_location)) {
            
            MR::ProgressBar progress_bar("Trimming tractlets...");
            
            BTS::Fibre::Tractlet::Set input(input_location);
            
            BTS::Fibre::Tractlet::Set output = trim_tractlets(input, distance,
                    length_reject_threshold, cube, num_length_sections,
                    num_width_sections, offsets);
            
            output.save(output_location);
            
        } else if (File::has_or_txt_extension<BTS::Fibre::Strand::Set>(input_location)) {
            
            MR::ProgressBar progress_bar("Trimming strand sets...");
            
            BTS::Fibre::Strand::Set::Reader reader(input_location);
            
            std::map<std::string, std::string> props = reader.get_extend_props();
            
            props["method"] = "trim";
            props["trim_dist"] = str(distance);
            props["cube"] = str(cube);
            
            std::vector<std::string> elem_header;
            
            elem_header.push_back(BTS::Fibre::Track::BUNDLE_INDEX_EPROP);
            
            BTS::Fibre::Strand::Set::Writer writer(output_location, reader,
                    reader.extend_prop_keys(), elem_header, props);
            
            BTS::Fibre::Strand::Set strands;
            
            while (reader.next(strands)) {
                
                BTS::Fibre::Track::Set tcks = strands.to_tracks(num_length_sections);
                
                std::map<std::string, std::string> props = strands.get_extend_props();
                
                BTS::Fibre::Track::Set trimmed_tcks = trim_tcks(tcks, distance,
                        length_reject_threshold, cube, offsets);
                
                //FIXME: Remove explicit as it should be handled by the "loaded_degree" property
                strands = trimmed_tcks.to_strands(strands[0].degree());
                
                strands.set_extend_props(props);
                
                writer.append(strands);
                
                ++progress_bar;
                
            }
            
        } else if (File::has_or_txt_extension<BTS::Fibre::Track::Set>(input_location)) {
            
            MR::ProgressBar progress_bar("Trimming track sets...");
            
            BTS::Fibre::Track::Set::Reader reader(input_location);
            BTS::Fibre::Track::Set::Writer writer(output_location, reader,
                    reader.extend_prop_keys(), reader.extend_elem_prop_keys(),
                    reader.get_extend_props());
            
            BTS::Fibre::Track::Set tracks;
            
            while (reader.next(tracks)) {
                writer.append(trim_tcks(tracks, distance, length_reject_threshold, cube, offsets));
                
                ++progress_bar;
            }
            
        } else if (File::has_or_txt_extension<BTS::Fibre::Tractlet::Set>(input_location)) {
            
            MR::ProgressBar progress_bar("Trimming tractlet sets...");
            
            BTS::Fibre::Tractlet::Set::Reader reader(input_location);
            
            std::map<std::string, std::string> props = reader.get_extend_props();
            
            props["method"] = "trim";
            props["trim_dist"] = str(distance);
            props["cube"] = str(cube);
            
            std::vector<std::string> elem_header;
            
            elem_header.push_back(BTS::Fibre::Track::BUNDLE_INDEX_EPROP);
            
            BTS::Fibre::Tractlet::Set::Writer writer(output_location, reader,
                    reader.extend_prop_keys(), elem_header, props);
            
            BTS::Fibre::Tractlet::Set tractlets;
            
            while (reader.next(tractlets)) {
                writer.append(
                        trim_tractlets(tractlets, distance, length_reject_threshold, cube,
                                num_length_sections, num_width_sections, offsets));
                ++progress_bar;
            }
            
        } else
            throw Exception(
                    "Unrecognised extension for input file '" + File::extension(input_location)
                    + "'.");
        
    }
    
    BTS::Fibre::Track::Set trim_tcks(BTS::Fibre::Track::Set& input_tcks, double distance,
                                     double length_reject_threshold, bool trim_to_cube,
                                     const Triple<double>& offsets) {
        
        BTS::Fibre::Track::Set trimmed_tcks(input_tcks.get_extend_props());
        
        std::vector<Triple<double> > pre_points;
        std::vector<Triple<double> > post_points;
        
        for (size_t tck_i = 0; tck_i < input_tcks.size(); ++tck_i)
            input_tcks[tck_i] -= offsets;

        Strand_collection c, trimmed_c;
        
        input_tcks.remove_short_tracks(3);
        
        generate_pre_points(input_tcks, pre_points);
        
        generate_post_points(input_tcks, post_points);
        
        convert_mr_to_nfg(&c, input_tcks, pre_points, post_points);
        
        trim(&trimmed_c, &c, distance, length_reject_threshold, 0, trim_to_cube);
        
        convert_nfg_to_mr(trimmed_tcks, pre_points, post_points, &trimmed_c);
        
        return trimmed_tcks;
        
    }
    
    BTS::Fibre::Tractlet::Set trim_tractlets(BTS::Fibre::Tractlet::Set& input_tractlets,
                                             double distance, double length_reject_threshold,
                                             bool trim_to_cube, size_t num_length_sections,
                                             size_t num_width_sections,
                                             const Triple<double>& offsets) {
        
        //FIXME: Remove this as it should be handled by the "loaded_degree" property
        size_t degree = input_tractlets[0].degree();
        
        BTS::Fibre::Track::Set input_tcks = input_tractlets.to_tracks(num_length_sections,
                num_width_sections);
        
        BTS::Fibre::Track::Set trimmed_tcks = trim_tcks(input_tcks, distance,
                length_reject_threshold, trim_to_cube, offsets);
        
        BTS::Fibre::Tractlet::Set trimmed_tractlets = trimmed_tcks.to_tractlets(degree);
        
        trimmed_tractlets.set_extend_props(input_tractlets.get_extend_props());
        
        return trimmed_tractlets;
        
    }
    
