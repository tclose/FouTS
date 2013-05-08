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
const double FA_THRESHOLD_DEFAULT = 0.7;
const double DENSITY_PERCENTILE_DEFAULT = 95.0;
const size_t NUM_AREA_SAMPLES_DEFAULT = 100;
const size_t NUM_LENGTH_SECTIONS_DEFAULT = 10;
const size_t NUM_WIDTH_SECTIONS_DEFAULT = 4;
const double RANGE_STEP_DEFAULT = 0.001;
const double RANGE_LOW_DEFAULT = 0.75;
const double RANGE_HIGH_DEFAULT = 1.25;


DESCRIPTION = {
    "Finds the maximum b0 intensity of an image",
    "",
    NULL
};

ARGUMENTS= {
    Argument("dwi", "the input diffusion-weighted image.").type_image_in(),
    Argument("intensities", "a sanity check of the intensities estimated at each voxel.").type_image_out(),
    Argument()
};

OPTIONS= {

    Option("middle", "the fraction of voxels which should be included in the estimation "
            "taken from the centre of the range of estimated intensities above the FA threshold")
    + Argument("middle").type_float(SMALL_FLOAT, MIDDLE_DEFAULT, LARGE_FLOAT),

    Option("curvature", "The relative curvature that is applied to the reference tract")
    + Argument("curvature").type_float(SMALL_FLOAT, CURVATURE_DEFAULT, LARGE_FLOAT),

    Option("fa_threshold", "the threshold above which the current voxel is counted as a white "
            "matter voxel")
    + Argument("fa_threshold").type_float(SMALL_FLOAT, FA_THRESHOLD_DEFAULT, LARGE_FLOAT),

    Option ("grad", "specify the diffusion-weighted gradient scheme used in the acquisition. "
            "The program will normally attempt to use the encoding stored in the image "
            "header. This should be supplied as a 4xN text file with each line is in "
            "the format [ X Y Z b ], where [ X Y Z ] describe the direction of the "
            "applied gradient, and b gives the b-value in units (1000 s/mm^2).")
    + Argument ("encoding").type_file(),

    Option ("reference_tract", "Instead of using straight tracks to fit the base intensity, a "
            "set of reference tracts can be provided instead.")
    + Argument ("reference").type_file(),

    Option("num_area_samples", "The number of samples taken along the reference tract when "
            "calculating the estimated area.")
    + Argument("num_area_samples").type_integer(1, NUM_AREA_SAMPLES_DEFAULT, LARGE_INT),

    Option("density_percentile", "The percentile used to normalise the density to (only applicable "
            "with 'ref_tracts' provided (fraction between 0 and 1).")
    + Argument("fa_threshold").type_float(0.0, DENSITY_PERCENTILE_DEFAULT, 100.0),

    Option("range", "The percentile used to normalise the density to (only applicable "
            "with 'ref_tracts' provided (fraction between 0 and 1).")
    + Argument("low").type_float(0.0, RANGE_LOW_DEFAULT, 100.0)
    + Argument("step").type_float(0.0, RANGE_STEP_DEFAULT, 100.0)
    + Argument("high").type_float(0.0, RANGE_HIGH_DEFAULT, 100.0),

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
        // Create an image to contain the estimated intensities (not actually used except for sanity
        // check).
        MR::Image::Header intens_header = dwi_header;
        intens_header.reset_scaling();
        intens_header.set_ndim (3);
        for (size_t dim_i = 0; dim_i < 3; ++dim_i)
            intens_header.set_dim (dim_i, dwi_header.dim(dim_i));
        intens_header.set_datatype (MR::DataType::Float32);
        intens_header.create (argument[1]);

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
        // Set the fraction of voxels to include in the final estimate (the middle '--fraction' is used)
        double middle = MIDDLE_DEFAULT;
        opt = get_options("middle");
        if (opt.size())
            middle = opt[0][0];
        if ((middle <= 0.0) || (middle > 1.0))
            throw Exception("'--middle' must be within (0-1] as it is the fraction of voxels to include "
                    "in the estimation (found " + str(middle) + ")");

        //------------------------------------------------------------------------------------------
        // The relative curvature of the reference tract used to calculate the base intensity from
        double curvature = CURVATURE_DEFAULT;
        opt = get_options("curvature");
        if (opt.size())
            curvature = opt[0][0];

        //------------------------------------------------------------------------------------------
        // The relative curvature of the reference tract used to calculate the base intensity from
        double fa_threshold = FA_THRESHOLD_DEFAULT;
        opt = get_options("fa_threshold");
        if (opt.size())
            fa_threshold = opt[0][0];

        //------------------------------------------------------------------------------------------
        // The relative curvature of the reference tract used to calculate the base intensity from
        double density_percentile = DENSITY_PERCENTILE_DEFAULT;
        opt = get_options("density_percentile");
        if (opt.size())
          density_percentile = opt[0][0];

        //------------------------------------------------------------------------------------------
        // The relative curvature of the reference tract used to calculate the base intensity from
        size_t num_area_samples = NUM_AREA_SAMPLES_DEFAULT;
        opt = get_options("num_length_sections");
        if (opt.size())
            num_area_samples = opt[0][0];

        //------------------------------------------------------------------------------------------
        // The relative curvature of the reference tract used to calculate the base intensity from
        double range_low = RANGE_LOW_DEFAULT;
        double range_step = RANGE_STEP_DEFAULT;
        double range_high = RANGE_HIGH_DEFAULT;
        opt = get_options("range");
        if (opt.size()) {
          range_low = opt[0][0];
          range_step = opt[0][1];
          range_high = opt[0][2];
        }

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
        diff_isotropic = true;

        //------------------------------------------------------------------------------------------
        // Get the matrix of non b=0 encoding directions
        MR::DWI::normalise_grad(grad);
        std::vector<int> bzeros, dwis;
        MR::DWI::guess_DW_directions(dwis, bzeros, grad);
        if (!bzeros.size())
            throw Exception("No b=0 encodings found in gradient encoding scheme");
        MR::Math::Matrix<float> all_encodings = dwi_header.get_DW_scheme();
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

        double estimated_intensity;

        if (reference_tract.size()) {

            if (reference_tract.size() > 1)
                throw Exception("The reference tracts is only meant to contain a single tract "
                        "(found " + str(reference_tract.size()) + ").");

            Image::Observed::Buffer obs_image(argument[0]);

            Image::Expected::Buffer* exp_image = Image::Expected::Buffer::factory(exp_type, obs_image,
                            diffusion_model, exp_num_length_sections, exp_num_width_sections,
                            exp_interp_extent, exp_enforce_bounds, exp_half_width);

            // Normalise reference tract density to the "density_percentile" and set the base
            // intensity to 1.0.
            reference_tract[0].normalise_density(num_area_samples, density_percentile);

            double intens_low = range_low * reference_tract.base_intensity();
            double intens_step = range_step * reference_tract.base_intensity();
            double intens_high = range_high * reference_tract.base_intensity();

            // Loop through the range of intensities and determine the optimum one.
            double max_likelihood = -INFINITY;
            for (double intens = intens_low; intens < intens_high; intens += intens_step) {
                //Generate image
                reference_tract.set_base_intensity(intens);
                exp_image->expected_image(reference_tract);

                double likelihood = 0.0;
                // Get the average ratio between the observed image and the expected image and set that
                // to be the estimated intensity.
                for (size_t z = 0; z < obs_image.dim(Z); ++z) {
                    for (size_t y = 0; y < obs_image.dim(Y); ++y) {
                        for (size_t x = 0; x < obs_image.dim(X); ++x) {
                            for (size_t encode_i = 0; encode_i < dwis.size(); ++encode_i) {
                                likelihood -= MR::Math::pow2(obs_image(x,y,z)[dwis[encode_i]] -
                                        (*exp_image)(x,y,z)[dwis[encode_i]]);
                            }
                        }
                    }
                }

                if (likelihood > max_likelihood) {
                    max_likelihood = likelihood;
                    estimated_intensity = intens;
                }
            }

        // If no reference tracts are provided estimate the intensity from the fit of straight tracts
        // to each single fibre voxel.
        } else {


            Triple<double> vox_lengths(dwi_header.vox(X), dwi_header.vox(Y), dwi_header.vox(Z));
            Triple<size_t> dims(1.0, 1.0, 1.0);
            Triple<double> offsets(0.0, 0.0, 0.0);
            Image::Expected::Buffer* exp_image = Image::Expected::Buffer::factory(exp_type, dims,
                    vox_lengths, diffusion_model, exp_num_length_sections, exp_num_width_sections,
                    exp_interp_extent, offsets, exp_enforce_bounds, exp_half_width);

            //------------------------------------------------------------------------------------------
            // Create the matrix and its singular-value decomposition used to calculate the diffusion
            // tensor at each voxel (for the alignment of the reference Fourier tract). The order of the
            // rows is d_xx, d_yy, d_zz, d_xy, d_xz, dyz.
            MR::Math::Matrix<double> tensor_encodings(num_nonb0_encodings, 6);
            for (size_t dim_i = 0; dim_i < 3; ++dim_i) {
                // Do the d_xx, d_yy and d_zz columns
                tensor_encodings.column(dim_i) = nonb0_encodings.column(dim_i);
                tensor_encodings.column(dim_i) *= nonb0_encodings.column(dim_i);
                // Make the d_xy and d_yz columns
                if (dim_i < 2) {
                    tensor_encodings.column(dim_i * 2 + 3) = nonb0_encodings.column(dim_i);
                    tensor_encodings.column(dim_i * 2 + 3) *= nonb0_encodings.column(dim_i + 1);
                }
            }
            // Finish of with the d_xz column
            tensor_encodings.column(4) = nonb0_encodings.column(X);
            tensor_encodings.column(4) *= nonb0_encodings.column(Z);
            // Caluculate the singular value decomposition on the tensor encodings matrix
            Math::USV usv = Math::svd(tensor_encodings);

            //------------------------------------------------------------------------------------------
            // Loop through all voxels and calculate the base intensities that would produce the
            // least difference from the voxel generated by the reference tract
            std::vector<double> intensities;
            MR::Image::Voxel<float> dwi_vox(dwi_header), intens_vox(intens_header);
            MR::DataSet::Loop vox_loop("Estimating intensities...", 0, 3);
            MR::DataSet::Loop encode_loop(3);
            for (vox_loop.start (dwi_vox, intens_vox); vox_loop.ok();
                    vox_loop.next (dwi_vox, intens_vox)) {

                //--------------------------------------------------------------------------------------
                // Initialise the intensity to 0 (also the default if it can't be estimated from the
                // current voxel)
                double intensity = 0.0;

                //--------------------------------------------------------------------------------------
                // Read non-(b==0) encodings into MR::Math::Vector and get the average b==0 value
                MR::Math::Vector<double> observed(num_nonb0_encodings);
                size_t all_encode_i = 0;
                size_t nonb0_encode_i = 0;
                double b0 = 0.0;
                double max_value = 0.0;
                for (encode_loop.start (dwi_vox); encode_loop.ok(); encode_loop.next (dwi_vox)) {
                    double value = (double)dwi_vox.value();
                    if (value > max_value)
                        max_value = value;
                    if (value < 0)
                        throw Exception("Negative value found in observed image (" + str(value) + "). "
                                        "NB: The diffusion tensor can not be found for images with no "
                                        "istropic components");
                    if (std::find(bzeros.begin(), bzeros.end(), all_encode_i) == bzeros.end())
                        observed[nonb0_encode_i++] = value;
                    else
                        b0 += value;
                    ++all_encode_i;
                }
                b0 /= bzeros.size();

                if (max_value > 0.0) {
                    //----------------------------------------------------------------------------------
                    // Calculate the diffusion tensor from the observed intens_header
                    MR::Math::Vector<double> observed_log_ratio(num_nonb0_encodings), tensor_vec(6);
                    for (size_t encode_i = 0; encode_i < num_nonb0_encodings; ++encode_i)
                        observed_log_ratio[encode_i] = -MR::Math::log(observed[encode_i] / b0)
                                / b_values[encode_i];
                    Math::solve_psuedo_inverse(observed_log_ratio, tensor_vec, usv);
                    MR::Math::Matrix<double> tensor(3, 3);
                    for (size_t dim_i = 0; dim_i < 3; ++dim_i) {
                        tensor(dim_i, dim_i) = tensor_vec[dim_i];
                        // Make the d_xy and d_yz columns
                        if (dim_i < 2)
                            tensor(dim_i, dim_i + 1) = tensor(dim_i + 1, dim_i) = tensor_vec[dim_i * 2 + 3];
                    }
                    tensor(X, Z) = tensor(Z, X) = tensor_vec[4];
                    // Get the eigenvectors of the tensor matrix
                    MR::Math::Matrix<double> evec(3, 3);
                    MR::Math::Vector<double> eval(3);
                    MR::Math::Eigen::SymmV<double> symmv(3);
                    symmv(eval, tensor, evec);
                    double trace = (eval[0] + eval[1] + eval[2]) / 3.0;
                    double fa =
                            MR::Math::sqrt(
                                    3.0 / 2.0 * (MR::Math::pow2(eval[0] - trace)
                                            + MR::Math::pow2(eval[1] - trace)
                                            + MR::Math::pow2(eval[2] - trace))
                                    / (MR::Math::pow2(eval[0]) + MR::Math::pow2(eval[1])
                                       + MR::Math::pow2(eval[2])));

                    if (fa > fa_threshold) {
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
                        tcts[0](0, 1) = evec.column(0);
                        tcts[0](0, 1) *= tract_extent;
                        // Set the width of the tract so that it extents past the extent of the voxel
                        tcts[0](1, 0) = evec.column(1);
                        tcts[0](1, 0) *= tract_extent;
                        tcts[0](2, 0) = evec.column(2);
                        tcts[0](2, 0) *= tract_extent;
                        // Set the curvature of the tract along the second diffusion eigenvector
                        // with relative magnitude to the main orientation of '--curvature'
                        tcts[0](0, 2) = evec.column(1);
                        tcts[0](0, 2) *= tract_extent * curvature;
                        // Normalize the density of the tract and set the base_intensity of the set to
                        // 1.0, to calculate the required base intensity value to match that of the
                        // reference.
                        tcts.normalise_densities();
                        tcts.set_base_intensity(1.0);
                        exp_image->expected_image(tcts);
                        // The estimated base intensity is then the average scaling required to get the
                        // expected intensity to match the observed intensity.
                        for (size_t encode_i = 0; encode_i < num_nonb0_encodings; ++encode_i)
                            intensity += observed[encode_i] / (*exp_image)(0, 0, 0)[encode_i];
                        intensity /= (double)num_nonb0_encodings;

                        intensities.push_back(intensity);
                    }
                }
                // Set the intensity in the intensities image
                intens_vox.value() = intensity;
            }

            //------------------------------------------------------------------------------------------
            // Get the median intensity and print it to the command line.
            if (!intensities.size())
                throw Exception("No voxels were above the FA threshold (" + str(fa_threshold) + ")");
            std::sort(intensities.begin(), intensities.end());
            estimated_intensity = intensities[intensities.size() / 2];
        }
        // Print the estimated intensity to the terminal for use with other commands
        std::cout << estimated_intensity << std::endl;
    }
