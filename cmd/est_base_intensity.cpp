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

const double FA_THRESHOLD_DEFAULT = 0.6;

DESCRIPTION = {
    "Finds the maximum b0 intensity of an image",
    "",
    NULL
};

ARGUMENTS= {
    Argument("dwi", "the input diffusion-weighted image.").type_image_in(),
    Argument()
};

OPTIONS= {
    Option("index", "the voxel indices in the image where the base intensity "
            "is calculated from. Should be in regions with known "
            "single fibre populations.").allow_multiple()
    + Argument("X", "The index along the x dimension").type_integer()
    + Argument("Y", "The index along the y dimension").type_integer()
    + Argument("Z", "The index along the z dimension").type_integer(),

    Option("seed", "the voxel indices in the image where the base intensity "
            "is calculated from. Should be in regions with known "
            "single fibre populations.").allow_multiple()
    + Argument("X", "The x coordinate").type_float()
    + Argument("Y", "The y coordinate").type_float()
    + Argument("Z", "The z coordinate").type_float(),

    Option ("grad", "specify the diffusion-weighted gradient scheme used in the acquisition. "
            "The program will normally attempt to use the encoding stored in the image "
            "header. This should be supplied as a 4xN text file with each line is in "
            "the format [ X Y Z b ], where [ X Y Z ] describe the direction of the "
            "applied gradient, and b gives the b-value in units (1000 s/mm^2).")
    + Argument ("encoding").type_file(),

    Option ("fa_threshold", "The minimum threshold for the calculated Fractional Anisotropy (FA) "
            "below which the provided will be omitted (with warning).")
    + Argument ().type_float(SMALL_FLOAT, FA_THRESHOLD_DEFAULT, LARGE_FLOAT),

    DIFFUSION_PARAMETERS,

    EXPECTED_IMAGE_PARAMETERS,

    LIKELIHOOD_PARAMETERS,

    Option()

};

EXECUTE {

        MR::Image::Header dwi(argument[0]);
        if (dwi.ndim() != 4)
            throw Exception("dwi image should contain 4 dimensions");
        // Set the script defaults
        MR::Math::Matrix<float> grad;
        double fa_threshold = FA_THRESHOLD_DEFAULT;
        std::vector<Triple<size_t> > indices;
        // Loads parameters to construct Diffusion::Model ('diff_' prefix)
        SET_DIFFUSION_PARAMETERS;
        // Loads extra parameters to construct Image::Expected::*::Buffer ('exp_' prefix)
        SET_EXPECTED_IMAGE_PARAMETERS
        // Set the diffusion model to isotropic (overiding the default) because diffusion tensors
        // only work with isotropic images
        diff_isotropic = true;
        // Supply the DW gradient scheme if required
        Options opt = get_options("grad");
        if (opt.size())
            grad.load(opt[0][0]);
        else {
            if (!dwi.DW_scheme().is_set())
                throw Exception("no diffusion encoding found in image \"" + dwi.name() + "\"");
            grad = dwi.DW_scheme();
        }
        if (grad.rows() < 7 || grad.columns() != 4)
            throw Exception("unexpected diffusion encoding matrix dimensions");
        if (dwi.dim(3) != (int)grad.rows())
            throw Exception("number of studies in base image does not match that in encoding file");
        // Set the FA anisotropy threshold
        opt = get_options("fa_threshold");
        if (opt.size())
            fa_threshold = opt[0][0];
        // Get the indices of the b=0 encoding directions
        MR::DWI::normalise_grad(grad);
        std::vector<int> bzeros, dwis;
        MR::DWI::guess_DW_directions(dwis, bzeros, grad);
        if (!bzeros.size())
            throw Exception("No b=0 encodings found in gradient encoding scheme");
        // Load the reference indices from the supplied '-index' and '-seed' options
        opt = get_options("index");
        for (size_t i = 0; i < opt.size(); ++i) {
            size_t x_ind = opt[i][X];
            size_t y_ind = opt[i][Y];
            size_t z_ind = opt[i][Z];
            indices.push_back(Triple<size_t>(x_ind, y_ind, z_ind));
        }
        opt = get_options("seed");
        for (size_t i = 0; i < opt.size(); ++i) {
            double x = opt[i][X];
            double y = opt[i][Y];
            double z = opt[i][Z];
            size_t x_ind = (size_t)((x - dwi.transform()(X, 3)) / dwi.vox(X));
            size_t y_ind = (size_t)((y - dwi.transform()(Y, 3)) / dwi.vox(Y));
            size_t z_ind = (size_t)((z - dwi.transform()(Z, 3)) / dwi.vox(Z));
            indices.push_back(Triple<size_t>(x_ind, y_ind, z_ind));
        }
        if (!indices.size())
            throw Exception("At least one reference voxel needs to be specified, using either the "
                            "'-index' or '-seed' options");
        // Get the matrix of non b=0 encoding directions
        Triple<double> vox_lengths(dwi.vox(X), dwi.vox(Y), dwi.vox(Z));
        Triple<size_t> dims(1.0, 1.0, 1.0);
        Triple<double> offsets(0.0, 0.0, 0.0);
        MR::Math::Matrix<float> all_encodings = dwi.get_DW_scheme();
        MR::Math::Matrix<double> nonb0_encodings(dwis.size(), 4);
        size_t num_nonb0_encodings = 0;
        for (size_t encode_i = 0; encode_i < dwi.get_DW_scheme().rows(); ++encode_i)
            if (std::find(bzeros.begin(), bzeros.end(), encode_i) == bzeros.end())
                nonb0_encodings.row(num_nonb0_encodings++) = all_encodings.row(encode_i);
        // Create the matrix used to calculate the diffusion tensor
        // The order of the rows is d_xx, d_yy, d_zz, d_xy, d_xz, dyz
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
        // Reference the b_values for readibility
        const MR::Math::Vector<float>& b_values = nonb0_encodings.column(DW);
        // Create the expected image voxel
        Diffusion::Model diffusion_model = Diffusion::Model::factory(nonb0_encodings,
                diff_response_SH, diff_adc, diff_fa, diff_isotropic, diff_warn_b_mismatch);
        Image::Expected::Buffer* exp_image = Image::Expected::Buffer::factory(exp_type, dims,
                vox_lengths, diffusion_model, exp_num_length_sections, exp_num_width_sections,
                exp_interp_extent, offsets, exp_enforce_bounds, exp_half_width);
        double est_base_intensity = 0.0;
        // Loop through all the provided single fibre indices and calculate the optimum b0
        for (std::vector<Triple<size_t> >::iterator index_it = indices.begin(); index_it != indices.end(); ++index_it) {
            // The following is quite a clunky way to get an observed buffer consisting of a single
            // voxel from the index given as a triple. NB the first argument to this script is the
            // image location hence the offset by one.
            MR::Math::Vector<double> observed(num_nonb0_encodings);
            size_t index[4] = {0, 0, 0, 0};
            (*index_it).copy_to(index);
            MR::Image::Voxel<float> vox(dwi);
            size_t loop_dims[4] = {1, 1, 1, all_encodings.rows()};
            MR::DataSet::Subset<MR::Image::Voxel<float> > vox_subset(vox, index, loop_dims);
            MR::DataSet::Loop loop(0, 4);
            size_t all_encode_i = 0;
            size_t nonb0_encode_i = 0;
            double b0 = 0.0;
            double max_value = 0.0;
            for (loop.start(vox_subset); loop.ok(); loop.next(vox_subset)) {
                double value = (double)vox_subset.value();
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
            if (max_value == 0.0)
                throw Exception("Selected voxel '" + str(*index_it) + "' contains no signal.");
            // Calculate the diffusion tensor from the observed intensities
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
                // Create a tract that spans the space that voxel will draw signal from and is aligned
                // with the principle eigenvector of the estimated diffusion tensor.
                Coord tract_extent = exp_interp_extent * vox_lengths * 2.0;
                Fibre::Tractlet::Set tcts(1, 2);
                tcts.zero();
                // Centre the tract in the middle of the voxel
                tcts[0](0, 0) = vox_lengths / 2.0;
                // Set the orientation of the tract to be in the direction of the principal eigenvector
                tcts[0](0, 1) = evec.column(0);
                tcts[0](0, 1) *= tract_extent;
                tcts[0](1, 0) = evec.column(1);
                tcts[0](1, 0) *= tract_extent;
                tcts[0](2, 0) = evec.column(2);
                tcts[0](2, 0) *= tract_extent;
                // Normalize the density of the tract and set the base_intensity of the set to 1.0,
                // to calculate the required base intensity value to match that of the reference.
                tcts.normalise_densities();
                tcts.set_base_intensity(1.0);
                exp_image->expected_image(tcts);
                // The estimated base intensity is then the average scaling required to get the expected
                // intensity to match the observed intensity.
                for (size_t encode_i = 0; encode_i < num_nonb0_encodings; ++encode_i)
                    est_base_intensity += observed[encode_i] / (*exp_image)(0, 0, 0)[encode_i];
            }
        }
        // Divide the sum of the scaling factor required to scale each
        est_base_intensity /= (double)(num_nonb0_encodings * indices.size());
        // Print the estimated intensity for use with other commands
        std::cout << est_base_intensity;
    }
