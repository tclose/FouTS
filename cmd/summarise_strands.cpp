/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Written by Thomas G. Close, 04/03/2009.

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
//

#include "progressbar.h"

#include "bts/common.h"

#include "bts/fibre/strand.h"
#include "bts/fibre/strand/set.h"

#include "bts/inline_functions.h"

using namespace BTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Extractlets salient features from a set of Strand descriptors",
    "Features are extractleted from the strand descriptor representation of a collection of strands.  These features can then be compared with priors on these features.",
    NULL
};

ARGUMENTS= {
    Argument ("input", "The strand coefficients to be analyzed.").allow_multiple().type_file (),
    Argument ("output", "output directory where the feature files will be written").type_file(),
    Argument()
};

OPTIONS= {

    Option ("all", "Saves all the features. Overrides all other flags"),

    Option ("v0_norm", "Save the norm of the v0 vector"),

    Option ("v1_norm", "Save the norm of the v1 vector"),

    Option ("v2_norm", "Save the norm of the v2 vector"),

    Option ("dot_v1_v2", "Save the dot product of the v1 and v2 vectors"),

    Option ("angle_v1_v2", "Save the angle between the v1 and v2 vectors"),

    Option ("dot_v0_v1", "Save the dot product of the v0 and v1 vectors"),

    Option ("degree", "The degree of the Strand coefficients used to describe the strands")
    + Argument ("degree", "").type_integer (1, 3, 200),

    Option()};

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_dir = argument[1];
        
        bool save_all = get_options("save_all").size();
        bool save_v0_norm = get_options("save_v0_norm").size() || save_all;
        bool save_v1_norm = get_options("save_v1_norm").size() || save_all;
        bool save_v2_norm = get_options("save_v2_norm").size() || save_all;
        bool save_dot_v1_v2 = get_options("save_dot_v1_v2").size() || save_all;
        bool save_angle_v1_v2 = get_options("save_angle_v1_v2").size() || save_all;
        bool save_dot_v0_v1 = get_options("save_dot_v0_v1").size() || save_all;
        size_t degree = 0;
        
        Options opt = get_options("degree");
        if (opt.size())
            degree = opt[0][0];
        
        Fibre::Strand::Set strands(input_location, degree);
        
        int count = 0;
        
        if (!(save_all || save_v0_norm || save_v1_norm || save_v2_norm || save_dot_v0_v1
              || save_dot_v1_v2 || save_angle_v1_v2))
            throw Exception("No options supplied, please select a feature to save");
        
        std::ofstream v0_norm_out, v1_norm_out, v2_norm_out, dot_v1_v2_out, angle_v1_v2_out,
                dot_v0_v1_out;
        
        if (save_v0_norm)
            v0_norm_out.open((output_dir + "/v0_norm.txt").c_str(), std::ios::out);
        
        if (save_v1_norm)
            v1_norm_out.open((output_dir + "/v1_norm.txt").c_str(), std::ios::out);
        
        if (save_v2_norm)
            v2_norm_out.open((output_dir + "/v2_norm.txt").c_str(), std::ios::out);
        
        if (save_dot_v1_v2)
            dot_v1_v2_out.open((output_dir + "/dot_v1_v2.txt").c_str(), std::ios::out);
        
        if (save_angle_v1_v2)
            angle_v1_v2_out.open((output_dir + "/angle_v1_v2.txt").c_str(), std::ios::out);
        
        if (save_dot_v0_v1)
            dot_v0_v1_out.open((output_dir + "/dot_v0_v1.txt").c_str(), std::ios::out);
        
        MR::ProgressBar progress_bar("extractleting features from strand descriptors...",
                strands.size());
        
        for (size_t strand_i = 0; strand_i < strands.size(); ++strand_i) {
            
            Fibre::Strand strand = strands[strand_i];
            
            Coord v0 = strand[0];
            Coord v1 = strand[1];
            Coord v2 = strand[2];
            
            double v0_norm = v0.norm();
            double v1_norm = v1.norm();
            double v2_norm = v2.norm();
            double dot_v1_v2 = v1.dot(v2) / (v1_norm * v2_norm);
            double dot_v0_v1 = v0.dot(v1) / (v1_norm * v0_norm);
            
            if (save_v0_norm)
                v0_norm_out << v0_norm << std::endl;
            
            if (save_v1_norm)
                v1_norm_out << v1_norm << std::endl;
            
            if (save_v2_norm)
                v2_norm_out << v2_norm << std::endl;
            
            if (save_dot_v1_v2)
                dot_v1_v2_out << dot_v1_v2 << std::endl;
            
            if (save_angle_v1_v2)
                angle_v1_v2_out << acos(dot_v1_v2) * 180 / M_PI << std::endl;
            
            if (save_dot_v0_v1)
                dot_v0_v1_out << dot_v0_v1 << std::endl;
            
            count++;
            progress_bar++;
            
        }
        
        if (save_v0_norm)
            v0_norm_out.close();
        
        if (save_v1_norm)
            v1_norm_out.close();
        
        if (save_v2_norm)
            v2_norm_out.close();
        
        if (save_dot_v1_v2)
            dot_v1_v2_out.close();
        
        if (save_angle_v1_v2)
            angle_v1_v2_out.close();
        
        if (save_dot_v0_v1)
            dot_v0_v1_out.close();
        
    }
