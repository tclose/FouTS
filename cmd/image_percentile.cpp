/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

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
#include <numeric>

#include "bts/common.h"
#include "bts/cmd.h"
#include "math/matrix.h"
#include "image/voxel.h"
#include "dwi/gradient.h"
#include "dataset/loop.h"

#include "bts/inline_functions.h"

using namespace FTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Finds the maximum b0 intensity of an image",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("dwi", "the input diffusion-weighted image.").type_image_in(),
    Argument ("mask", "the input mask.").type_image_in(),
    Argument()
};

OPTIONS= {

    Option("percentile", "Used to select the method used to get the b0, one of 'average', 'max', or 'median' (default).")
    + Argument().type_float(0.0, 50.0, 100.0),

    Option()

};

EXECUTE {
    
        MR::Image::Header header(argument[0]);
        
        if (header.ndim() != 3)
            throw Exception("image should contain 3 dimensions (found " + str(header.ndim()) + ")");
        
        MR::Math::Matrix<float> grad;
        double percentile = 50.0;
        
        Options opt = get_options("percentile");
        if (opt.size())
            percentile = opt[0][0];
        if (percentile < 0.0 || percentile > 100.0)
            throw Exception("Percentile out of range (" + str(percentile) + ")");
        
        MR::Image::Header mask_header(argument[1]);
        if (mask_header.dim(0) != header.dim(0) || mask_header.dim(1) != header.dim(1)
            || mask_header.dim(2) != header.dim(2))
            throw Exception("mask & DWI image dimensions do not match");
        // Loop through all voxels in the image and push them to a vector
        MR::Image::Voxel<float> vox(header);
        MR::Image::Voxel<bool> mask(mask_header);
        std::vector<double> values;
        values.reserve(header.dim(X) * header.dim(Y) * header.dim(Z));
        MR::DataSet::Loop loop;
        for (loop.start(mask, vox); loop.ok(); loop.next(mask, vox)) {
            if (mask.value()) {
                double value = vox.value();
                values.push_back(value);
            }
        }
        // Calculate the desired percentile
        std::vector<double>::iterator first = values.begin();
        std::vector<double>::iterator last = values.end();
        size_t index = (size_t)(round((double)(last - first - 1) * percentile / 100.0));
        std::vector<double>::iterator target = first + index;
        std::nth_element(first, target, last);
        // Print percentile to standard output
        std::cout << *target << std::endl;
    }
    
