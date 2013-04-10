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

#include "phantom/interface.h"

#include "phantom/shared/strand_collection.h"
#include "phantom/resample/resample.h"

#include "bts/inline_functions.h"

using namespace BTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Resamples control points along tracks.",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("input", "The tracks to be resampled.").type_file (),
    Argument ("output", "The returned resampled tracks.").optional().type_file(),
    Argument()
};

OPTIONS= {

    Option ("resample_length", "Defines the length of the new inter control-point intervals")
    + Argument ("resample_length", "").type_float (0.00001, 0.02, 10000.0),

    Option ("double_back_angle_threshold", "The angle between successive segments above which the second segment is considered to have doubled back on itself")
    + Argument ("double_back_angle_threshold", "").type_float (90.0, 150.0, 180.0),

    Option ("forward_angle_threshold", "The angle after a double back has been detected above which the subsequent segment is considered to have travelled back in a forward direction again.")
    + Argument ("forward_angle_threshold","").type_float (5.0, 80.0, 180.0),

    Option ("num_points", "The number of points that will be generated along the strand path if the loaded tracks are fourier descriptors")
    + Argument ("num_points", "").type_integer (1, 100, 2000),

    Option ("degree", "The degree of the Fourier coefficients used to describe the ")
    + Argument ("degree", "").type_integer (1, 8, 200),

    Option()};

EXECUTE {
    
        std::string input_path = argument[0];
        std::string output_path;
        
        if (argument.size() > 1)
            output_path = argument[1].c_str();
        else
            output_path = input_path;
        
        double resample_length = 0.02;
        double double_back_angle_threshold = 150.0;
        double forward_angle_threshold = 80.0;
        size_t num_points = 0;
        size_t degree = 0;
        
        Options opt = get_options("resample_length");
        if (opt.size())
            resample_length = opt[0][0];
        
        opt = get_options("double_back_angle_threshold");
        if (opt.size())
            double_back_angle_threshold = opt[0][0];
        
        opt = get_options("forward_angle_threshold");
        if (opt.size())
            forward_angle_threshold = opt[0][0];
        
        opt = get_options("num_points");
        if (opt.size())
            num_points = opt[0][0];
        
        opt = get_options("degree");
        if (opt.size())
            degree = opt[0][0];
        
        BTS::Fibre::Track::Set tcks(input_path, num_points);
        BTS::Fibre::Track::Set resampled_tcks(tcks.get_extend_props());
        
        std::vector<Triple<double> > pre_points;
        std::vector<Triple<double> > post_points;
        
        Strand_collection c, resampled_c;
        
        generate_pre_points(tcks, pre_points);
        
        generate_post_points(tcks, post_points);
        
        convert_mr_to_nfg(&c, tcks, pre_points, post_points);
        
        MR::ProgressBar progress_bar("Resampling tracks...");
        resample_collection(&c, &resampled_c, resample_length, double_back_angle_threshold,
                forward_angle_threshold);
        //MR::ProgressBar::done();
        
        convert_nfg_to_mr(resampled_tcks, pre_points, post_points, &resampled_c);
        
        resampled_tcks.save(output_path, degree);
        
    }
    
