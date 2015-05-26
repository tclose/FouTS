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

#include "progressbar.h"

#include "bts/common.h"
#include "bts/file.h"

#include "bts/fibre/track/set.h"
#include "bts/fibre/strand/set.h"

#include "bts/inline_functions.h"

using namespace FTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Import ascii path files into FTS format.",
    "The input ascii files should consist of 3xN matrices, corresponding to the [ X Y Z ] coordinates of the N points making up the path. All the input paths will be included into the same output FTS path file.",
    NULL
};

ARGUMENTS= {
    Argument ("input", "The input strand text files to be imported into FTSTrix format.").allow_multiple().type_file (),
    Argument ("output", "The output strands file in FTS format").type_file(),
    Argument()
};

OPTIONS= {

    Option ("prefix", "Match files in the input directory with the following prefix")
    + Argument ("prefix", ""),

    Option ("nfg_format", "Import the NFG format and save the extended track properties"),

    Option()};

void import(Fibre::Strand::Set& s_set, const std::string& location, const std::string& prefix,
            MR::ProgressBar& progress_bar);
void import(Fibre::Track::Set& t_set, const std::string& location, const std::string& prefix,
            bool nfg_format, MR::ProgressBar& progress_bar);

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location = argument.back();
        
        bool nfg_format = false;
        std::string prefix = "";
        
        Options opt = get_options("prefix");
        if (opt.size())
            prefix = opt[0][0].c_str();
        
        opt = get_options("nfg_format");
        
        if (opt.size()) {
            
            nfg_format = true;
            
            if (File::has_extension<Fibre::Strand>(output_location))
                throw Exception("'-nfg_format' can only be used with the '.tck' file format.");
            
            if (prefix.size())
                throw Exception(
                        "'-prefix' option cannot be used in conjunction with '-nfg_format'.");
            
            prefix = "strand_";
            
        }
        
        if (!File::exists(input_location))
            throw Exception("Input location, " + input_location + ", does not exist");
        
        MR::ProgressBar progress_bar("adding ASCII strand files");
        
        if (File::has_extension<Fibre::Strand>(output_location)) {
            
            Fibre::Strand::Set s_set;
            
            std::string input_locations = "";
            for (size_t n = 0; n < argument.size() - 1; n++) {
                import(s_set, argument[n], prefix, progress_bar);
                input_locations += argument[n].c_str();
            }
            
            std::map<std::string, std::string> properties;
            properties["import_dir"] = input_locations;
            
            s_set.save(output_location);
            
        } else if (File::has_extension<Fibre::Track>(output_location)) {
            
            Fibre::Track::Set t_set;
            
            if (nfg_format) {
                t_set.add_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, "0");
                
                t_set.add_extend_elem_prop(Fibre::Track::RADIUS_PROP, "0.0");
            }
            
            std::string input_locations = "";
            for (size_t n = 0; n < argument.size() - 1; n++) {
                
                input_locations += argument[n].c_str();
                import(t_set, argument[n], prefix, nfg_format, progress_bar);
                
            }
            
            t_set.set_extend_prop("import_dir", input_locations);
            
            t_set.save(output_location);
            
        } else
            throw Exception("Unknown file extension in output file " + output_location + ".");
        
        //MR::ProgressBar::done();
    }
    
    void import(Fibre::Strand::Set& s_set, const std::string& location, const std::string& prefix,
                MR::ProgressBar& progress_bar) {
        
        if (location.length() >= 21) {
            if (!location.substr(location.length() - 21).compare("isotropic_regions.txt"))
                throw Exception(
                        "Found 'isotropic_regions.txt' file, which is assumed to be from Numerical Fibre Generator (NFG) output, please use '-nfg_format' option to import these tracks, and use the '.tck' for the imported tracks file extension.");
        }
        
        if (File::is_dir(location)) {
            
            std::vector<std::string> filenames = File::list_filenames(location, prefix, "txt");
            
            for (size_t filename_i = 0; filename_i < filenames.size(); filename_i++)
                import(s_set, File::join(location, filenames[filename_i]), prefix, progress_bar);
            
        } else {
            
            MR::Math::Matrix<double> M;
            
            M.load(location);
            if (M.columns() != 3)
                throw Exception(
                        "WARNING: file \"" + location + "\" does not contain 3 columns - ignored");
            
            Fibre::Strand strand;
            
            for (size_t i = 0; i < M.rows(); i++) {
                Triple<double> d(M(i, 0), M(i, 1), M(i, 2));
                strand.push_back(d);
            }
            s_set.push_back(strand);
            
            progress_bar++;
            
        }
        
    }
    
    void import(Fibre::Track::Set& t_set, const std::string& location, const std::string& prefix,
                bool nfg_format, MR::ProgressBar& progress_bar) {
        
        if (location.length() >= 21) {
            if (!location.substr(location.length() - 21).compare("isotropic_regions.txt") && !nfg_format)
                throw Exception(
                        "Found 'isotropic_regions.txt' file, which is assumed to be from Numerical Fibre Generator (NFG) output, please use '-nfg_format' option to import these tracks.");
        }
        
        if (File::is_dir(location)) {
            
            std::vector<std::string> filenames = File::list_filenames(location, prefix, "txt");
            
            for (size_t filename_i = 0; filename_i < filenames.size(); filename_i++)
                import(t_set, File::join(location, filenames[filename_i]), prefix, nfg_format,
                        progress_bar);
            
        } else {
            
            MR::Math::Matrix<double> M;
            
            M.load(location);
            if (M.columns() != 3)
                throw Exception(
                        "WARNING: file \"" + location + "\" does not contain 3 columns - ignored");
            
            Fibre::Track tck;
            
            for (size_t i = 0; i < M.rows(); i++) {
                tck.push_back(Coord(M(i, 0), M(i, 1), M(i, 2)));
            }
            
            t_set.push_back(tck);
            
            if (nfg_format) {
                
                int bundle_index;
                double track_radius;
                
                sscanf(File::basename(location).c_str(), "strand_%*d-%d-r%lf.txt", &bundle_index,
                        &track_radius);
                
                t_set.set_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, str(bundle_index),
                        t_set.size() - 1);
                
                t_set.set_extend_elem_prop(Fibre::Track::RADIUS_PROP, str(track_radius),
                        t_set.size() - 1);
                
            }
            
            progress_bar++;
            
        }
        
    }
