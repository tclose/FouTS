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

const double MIDDLE_DEFAULT = 0.5;
const double CURVATURE_DEFAULT = 0.0;
const size_t NUM_LENGTH_SECTIONS_DEFAULT = 10;
const size_t NUM_WIDTH_SECTIONS_DEFAULT = 4;

DESCRIPTION = {
    "Finds the maximum b0 intensity of an image",
    "",
    NULL
};

ARGUMENTS= {
    Argument("dwi", "the input diffusion-weighted image.").type_image_in(),
    Argument("orientations", "The orientations used to align the reference tracts to.").type_image_in(),
    Argument("mask", "only perform computation within the specified binary brain mask image.").type_image_in (),
    Argument("intensities", "a sanity check of the intensities estimated at each voxel.").optional().type_image_out(),
    Argument()
};

OPTIONS= {

    Option("curvature", "The relative curvature that is applied to the reference tract")
    + Argument("curvature").type_float(SMALL_FLOAT, CURVATURE_DEFAULT, LARGE_FLOAT),

    Option ("grad", "specify the diffusion-weighted gradient scheme used in the acquisition. "
            "The program will normally attempt to use the encoding stored in the image "
            "header. This should be supplied as a 4xN text file with each line is in "
            "the format [ X Y Z b ], where [ X Y Z ] describe the direction of the "
            "applied gradient, and b gives the b-value in units (1000 s/mm^2).")
    + Argument ("encoding").type_file(),


    DIFFUSION_PARAMETERS,

    EXPECTED_IMAGE_PARAMETERS,

    LIKELIHOOD_PARAMETERS,

    Option()

};

EXECUTE {
        //------------------------------------------------------------------------------------------
        //Load the DWI image header
        MR::Image::Header dwi_header(argument[0]);
        if (dwi_header.ndim() != 4)
            throw Exception("dwi image should contain 4 dimensions");

        //------------------------------------------------------------------------------------------
        //Load the orientations image header
        MR::Image::Header orient_header(argument[1]);
        if (orient_header.ndim() != 4)
            throw Exception("orientations image should contain 4 dimensions");

        //------------------------------------------------------------------------------------------
        //Load the orientations image header
        MR::Image::Header mask_header(argument[2]);
        if (mask_header.ndim() != 3)
            throw Exception("mask image should contain 3 dimensions");

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
        // The relative curvature of the reference tract used to calculate the base intensity from
        Fibre::Tractlet::Set reference_tract;
        opt = get_options("reference_tract");
        if (opt.size())
            reference_tract.load(opt[0][0].c_str());

        //------------------------------------------------------------------------------------------
        // The relative curvature of the reference tract used to calculate the base intensity from
        double curvature = CURVATURE_DEFAULT;
        opt = get_options("curvature");
        if (opt.size())
            curvature = opt[0][0];

        //------------------------------------------------------------------------------------------
        // Loads parameters to construct Diffusion::Model ('diff_' prefix)
        SET_DIFFUSION_PARAMETERS;
        // Loads extra parameters to construct Image::Expected::*::Buffer ('exp_' prefix)
        SET_EXPECTED_IMAGE_PARAMETERS
        // Overide the default number of width and length sections if they haven't been explicitly
        // provided
        opt = get_options("exp_num_length_sections");
        if (!opt.size())
            exp_num_width_sections = NUM_LENGTH_SECTIONS_DEFAULT;
        opt = get_options("exp_num_width_sections");
        if (!opt.size())
            exp_num_width_sections = NUM_WIDTH_SECTIONS_DEFAULT;
        // Set the diffusion model to isotropic (overiding the default)

        //------------------------------------------------------------------------------------------
        // Get the matrix of non b=0 encoding directions
        MR::DWI::normalise_grad(grad);
        std::vector<int> bzeros, dwis;
        MR::DWI::guess_DW_directions(dwis, bzeros, grad);
        if (!bzeros.size())
            throw Exception("No b=0 encodings found in gradient encoding scheme");
        MR::Math::Matrix<float> all_encodings = dwi_header.get_DW_scheme();
        double estimated_intensity;

        MR::Math::Matrix<double> nonb0_encodings(dwis.size(), 4);
        size_t num_nonb0_encodings = 0;
        for (size_t encode_i = 0; encode_i < dwi_header.get_DW_scheme().rows(); ++encode_i)
            if (std::find(bzeros.begin(), bzeros.end(), encode_i) == bzeros.end())
                nonb0_encodings.row(num_nonb0_encodings++) = all_encodings.row(encode_i);
        // Reference the b_values for readibility
        const MR::Math::Vector<float>& b_values = nonb0_encodings.column(DW);

        //------------------------------------------------------------------------------------------
        // Create the expected image voxel
        Diffusion::Model diffusion_model = Diffusion::Model::factory(nonb0_encodings,
                diff_response_SH, diff_adc, diff_fa, diff_isotropic, diff_warn_b_mismatch);

        //------------------------------------------------------------------------------------------
        // Create an image to contain the estimated intensities (not actually used except for sanity
        // check).
        MR::Image::Header intens_header = dwi_header;
        intens_header.reset_scaling();
        intens_header.set_ndim (3);
        for (size_t dim_i = 0; dim_i < 3; ++dim_i)
            intens_header.set_dim (dim_i, dwi_header.dim(dim_i));
        intens_header.set_datatype (MR::DataType::Float32);
        intens_header.create (argument[3]);

        Triple<double> vox_lengths(dwi_header.vox(X), dwi_header.vox(Y), dwi_header.vox(Z));
        Triple<size_t> dims(1.0, 1.0, 1.0);
        Triple<double> offsets(0.0, 0.0, 0.0);
        Image::Expected::Buffer* exp_image = Image::Expected::Buffer::factory(exp_type, dims,
                vox_lengths, diffusion_model, exp_num_length_sections, exp_num_width_sections,
                exp_interp_extent, offsets, exp_enforce_bounds, exp_half_width);

        //------------------------------------------------------------------------------------------
        // Loop through all voxels and calculate the base intensities that would produce the
        // least difference from the voxel generated by the reference tract
//        std::vector<double> intensities;
        MR::Image::Voxel<float> dwi_vox(dwi_header), orient_vox(orient_header), mask_vox(mask_header),
                intens_vox(intens_header);
        MR::DataSet::Loop vox_loop("Estimating intensities...", 0, 3);
        MR::DataSet::Loop vox_loop2;
        MR::DataSet::Loop encode_loop(3);
        vox_loop2.start(intens_vox); // Can't do this in vox_loop because it only takes up to three
        for (vox_loop.start (dwi_vox, orient_vox, mask_vox); vox_loop.ok();
                vox_loop.next (dwi_vox, orient_vox, mask_vox)) {

            if (mask_vox.value() > 0.5) { // If current voxel is within the mask

                Triple<double> orient;
                size_t dim_i = 0;
                for (encode_loop.start(orient_vox); encode_loop.ok(); encode_loop.next(orient_vox))
                    orient[dim_i++] = orient_vox.value();
                orient.normalise();

                // Here get perpendicular vectors for use creating tract
                Triple<double> perp1 = orient.cross(orient.min_axis());
                Triple<double> perp2 = orient.cross(perp1);

                //------------------------------------------------------------------------------
                // Create a tract that spans the space that voxel will draw signal from and is
                // aligned with the principle eigenvector of the estimated diffusion tensor,
                // with a ratio betwen the m=1 and m=2 coefficient vectors of the primary
                // axis provided by the '--curvature' option.
                Coord tract_extent = exp_interp_extent * vox_lengths * 2.0;
                Fibre::Tractlet::Set tcts(1, 3);
                tcts.zero();
                // Centre the tract in the middle of the voxel
                tcts[0](0, 0) = vox_lengths / 2.0;
                // Set the main orientation of the tract to be in the direction of the principal
                // eigenvector
                tcts[0](0, 1) = orient;
                tcts[0](0, 1) *= tract_extent;
                // Set the width of the tract so that it extents past the extent of the voxel
                tcts[0](1, 0) = perp1;
                tcts[0](1, 0) *= tract_extent;
                tcts[0](2, 0) = perp2;
                tcts[0](2, 0) *= tract_extent;
                // Set the curvature of the tract along the second diffusion eigenvector
                // with relative magnitude to the main orientation of '--curvature'
                tcts[0](0, 2) = perp1;
                tcts[0](0, 2) *= tract_extent * curvature;
                // Normalize the density of the tract and set the base_intensity of the set to
                // 1.0, to calculate the required base intensity value to match that of the
                // reference.
                tcts.normalise_densities();
                tcts.set_base_intensity(1.0);
                exp_image->expected_image(tcts);
                // The estimated base intensity is then the average scaling required to get the
                // expected intensity to match the observed intensity.
                double intensity = 0.0;
                size_t encode_i = 0;
                for (std::vector<int>::iterator it = dwis.begin(); it != dwis.end(); ++it) {
                    dwi_vox[3] = *it;
                    float observed = dwi_vox.value();
                    float expected = (*exp_image)(0, 0, 0)[encode_i++];
                    intensity += observed / expected;
                }
                intensity /= (double)num_nonb0_encodings;
//
//                intensities.push_back(intensity);

                // Set the intensity in the intensities image
                intens_vox.value() = intensity;
            }
            vox_loop2.next(intens_vox);
        }
//        //------------------------------------------------------------------------------------------
//        // Get the median intensity and print it to the command line.
//        std::sort(intensities.begin(), intensities.end());
//        estimated_intensity = intensities[intensities.size() / 2];
//        // Print the estimated intensity to the terminal for use with other commands
//        std::cout << estimated_intensity << std::endl;
    }
