/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 10/09/2010.

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
#include "bts/file.h"
#include "image/header.h"
#include "bts/prob/test/landscape.h"

#include "bts/inline_functions.h"

using namespace FTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "test_scan",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("input_peaks","Location of the input peaks file").type_file (),
    Argument ("output_image","Location of the output image").optional().type_image_out (),
    Argument()
};

double NUM_STEPS_DEFAULT = 100;
double AXIS_SCALE_DEFAULT = 15.0;
OPTIONS= {

    Option ("num_steps1", "Number of taken along axis.")
    + Argument ("num_steps", "").type_integer (1, NUM_STEPS_DEFAULT, LARGE_INT),

    Option ("num_steps2", "Number of taken along second axis.")
    + Argument ("num_steps", "").type_integer (1, NUM_STEPS_DEFAULT, LARGE_INT),

    Option ("axis1_scale", "Scale applied to axis1, used to quickly scale binary masks.")
    + Argument ("axis1_scale", "").type_float (SMALL_FLOAT, AXIS_SCALE_DEFAULT, LARGE_FLOAT),

    Option ("axis2_scale", "Scale applied to axis2, used to quickly scale binary masks.")
    + Argument ("axis2_scale", "").type_float (SMALL_FLOAT, AXIS_SCALE_DEFAULT, LARGE_FLOAT),

    Option()};

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location;
        
        if (argument.size() == 2)
            output_location = argument[0].c_str();
        else
            output_location = File::strip_extension(input_location) + ".mif";
        
        size_t num_steps1 = NUM_STEPS_DEFAULT;
        ;
        size_t num_steps2;
        
        double axis1_scale = AXIS_SCALE_DEFAULT;
        double axis2_scale;
        
        Options opt = get_options("num_steps1");
        if (opt.size())
            num_steps1 = opt[0][0];
        
        opt = get_options("num_steps2");
        if (opt.size())
            num_steps2 = opt[0][0];
        else
            num_steps2 = num_steps1;
        
        opt = get_options("axis1_scale");
        if (opt.size())
            axis1_scale = opt[0][0];
        
        opt = get_options("axis2_scale");
        if (opt.size())
            axis2_scale = opt[0][0];
        else
            axis2_scale = axis1_scale;
        
        Prob::Test::Landscape peaks(input_location);
        
        if (peaks.num_dims() != 2)
            throw Exception(
                    "Can only produce scan image for peaks of dimension 2 (" + str(peaks.num_dims())
                    + " supplied).");
        
        MCMC::State point(2);
        
        MR::Image::Header header;
        
        header.set_ndim(2);
        
        header.set_dim(0, num_steps1);
        header.set_dim(1, num_steps2);
        
        header.set_description(0, File::basename("Axis 1"));
        header.set_description(1, File::basename("Axis 2"));
        
        header.set_vox(0, 2.0 * axis1_scale / (double) num_steps1);
        header.set_vox(1, 2.0 * axis2_scale / (double) num_steps2);
        
        header["num_steps1"] = str(num_steps1);
        header["num_steps2"] = str(num_steps2);
        
        header["axis1_scale"] = str(axis1_scale);
        header["axis2_scale"] = str(axis2_scale);
        
        File::clear_path(output_location);
        
        header.create(output_location);
        
        double inc1 = 2.0 / (double) (num_steps1 - 1);
        double inc2 = 2.0 / (double) (num_steps2 - 1);
        
        MR::Image::Voxel<double> pixel(header);
        
        //Nicety, makes sure the loop goes out of scope, and prints its close message before the MATLAB plot command is issued.
        if (1) {
            
            MR::DataSet::Loop loop("Scanning over 2 dimensions...", 0, 2);
            loop.start(pixel);
            
            for (double frac1 = -1.0; frac1 <= 1.0 + inc1 / 2.0; frac1 += inc1) {    // The +inc/2.0 safeguards against rounding errors.
                    
                for (double frac2 = -1.0; frac2 <= 1.0 + inc2 / 2.0; frac2 += inc2) {
                    
                    point[0] = frac1 * axis1_scale;
                    point[1] = frac2 * axis2_scale;
                    
                    pixel.value() = peaks.log_prob(point);
                    
                    loop.next(pixel);
                    
                }
                
            }
        }
        
        std::cout << std::endl << "\nMATLAB plot command:\n\nplot_test_scan " << output_location
                  << "\n" << std::endl;
        
    }
