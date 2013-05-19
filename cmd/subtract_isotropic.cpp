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
#include <numeric>

#include "bts/common.h"
#include "bts/cmd.h"
#include "math/matrix.h"
#include "image/voxel.h"
#include "dwi/gradient.h"
#include "dataset/loop.h"
#include "dataset/subset.h"
#include "math/eigen.h"

#include "bts/math/svd.h"

#include "bts/image/expected/buffer.h"
#include "bts/image/expected/trilinear/buffer.h"
#include "bts/image/expected/gaussian/buffer.h"
#include "bts/image/expected/quartic/buffer.h"
#include "bts/diffusion/model.h"

#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood/one_sided_gaussian.h"
#include "bts/prob/likelihood/gaussian.h"

#include "bts/inline_functions.h"

using namespace BTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Finds the estimate noise map from the standard deviation between b=0 encoding images",
    "",
    NULL
};

ARGUMENTS= {
    Argument("dwi", "the input diffusion-weighted image.").type_image_in(),
    Argument("new", "the noise levels estimated at each voxel.").type_image_out(),
    Argument()
};

OPTIONS= {

    Option ("grad", "specify the diffusion-weighted gradient scheme used in the acquisition. "
            "The program will normally attempt to use the encoding stored in the image "
            "header. This should be supplied as a 4xN text file with each line is in "
            "the format [ X Y Z b ], where [ X Y Z ] describe the direction of the "
            "applied gradient, and b gives the b-value in units (1000 s/mm^2).")
    + Argument ("encoding").type_file(),

    Option()

};

EXECUTE {
        //------------------------------------------------------------------------------------------
        //Load the DWI image header
        MR::Image::Header dwi_header(argument[0]);
        if (dwi_header.ndim() != 4)
            throw Exception("dwi image should contain 4 dimensions");

        //------------------------------------------------------------------------------------------
        // Create an image to contain the estimated intensities (not actually used except for sanity
        // check).
        MR::Image::Header noiso_header = dwi_header;
        noiso_header.create (argument[1]);

        //------------------------------------------------------------------------------------------
        // Supply the DW gradient scheme if required
        MR::Math::Matrix<float> grad;
        Options opt = get_options("grad");
        if (opt.size())
            grad.load(opt[0][0]);
        else {
            if (!dwi_header.DW_scheme().is_set())
                throw Exception("no diffusion encoding found in image \"" + dwi_header.name() + "\"");
            grad = dwi_header.DW_scheme();
        }
        if (grad.rows() < 7 || grad.columns() != 4)
            throw Exception("unexpected diffusion encoding matrix dimensions");
        if (dwi_header.dim(3) != (int)grad.rows())
            throw Exception("number of studies in base image does not match that in encoding file");

        //------------------------------------------------------------------------------------------
        // Get the matrix of non b=0 encoding directions
        MR::DWI::normalise_grad(grad);
        std::vector<int> bzeros, dwis;
        MR::DWI::guess_DW_directions(dwis, bzeros, grad);

        //------------------------------------------------------------------------------------------
        // Loop through all voxels and calculate the base intensities that would produce the
        // least difference from the voxel generated by the reference tract
        MR::Image::Voxel<float> dwi_vox(dwi_header), noiso_vox(noiso_header);
        MR::DataSet::Loop vox_loop("Subtracting isotropic components...", 0, 3);
        MR::DataSet::Loop encode_loop(3);
        for (vox_loop.start (dwi_vox, noiso_vox); vox_loop.ok();
                vox_loop.next (dwi_vox, noiso_vox)) {

            float average = 0.0;
            //--------------------------------------------------------------------------------------
            for (std::vector<int>::iterator it = dwis.begin(); it != dwis.end(); ++it) {
                dwi_vox[3] = *it;
                average += dwi_vox.value();
            }

            average /= (float)dwis.size();

            for (std::vector<int>::iterator it = dwis.begin(); it != dwis.end(); ++it) {
                dwi_vox[3] = *it;
                noiso_vox[3] = *it;
                noiso_vox.value() = dwi_vox.value() - average;
            }

        }

    }
