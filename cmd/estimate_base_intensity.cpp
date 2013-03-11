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

const size_t ONE_VOX[3] = {1, 1, 1}
const double FA_THRESHOLD_DEFAULT = 0.6;

DESCRIPTION = {
    "Finds the maximum b0 intensity of an image",
    "",
    NULL
};

ARGUMENTS= {
    Argument("dwi", "the input diffusion-weighted image.").type_image_in(),
    Argument("indices", "the voxel indices in the image where the base intensity "
            "is calculated from. Should be in regions with known "
            "single fibre populations.").type_text().allow_multiple(),
    Argument()
};

OPTIONS= {

    Option ("grad", "specify the diffusion-weighted gradient scheme used in the acquisition. "
            "The program will normally attempt to use the encoding stored in the image "
            "header. This should be supplied as a 4xN text file with each line is in "
            "the format [ X Y Z b ], where [ X Y Z ] describe the direction of the "
            "applied gradient, and b gives the b-value in units (1000 s/mm^2).")
    + Argument ("encoding").type_file(),

    Option ("fa_threshold", "The minimum threshold for the calculated Fractional Anisotropy (FA) "
            "below which the provided will be omitted (with warning).")

    DIFFUSION_PARAMETERS,

    EXPECTED_IMAGE_PARAMETERS,

    LIKELIHOOD_PARAMETERS,

    Option()

};

EXECUTE {

        if (argument.size() < 2)
            throw Exception(
                    "Found less than two arguments (" + str(argument.size()) + "). The reference"
                    " image and at least one " "single fibre" " index needs to be provided");

        MR::Image::Header dwi(argument[0]);

        if (dwi.ndim() != 4)
            throw Exception("dwi image should contain 4 dimensions");

        MR::Math::Matrix<float> grad;
        double fa_threshold = FA_THRESHOLD_DEFAULT;

        // Loads parameters to construct Diffusion::Model ('diff_' prefix)
        SET_DIFFUSION_PARAMETERS;

        // Loads extra parameters to construct Image::Expected::*::Buffer ('exp_' prefix)
        SET_EXPECTED_IMAGE_PARAMETERS

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
        MR::DWI::normalise_grad(grad);
        std::vector<int> bzeros, dwis;
        MR::DWI::guess_DW_directions(dwis, bzeros, grad);
        // Set the FA anisotropy threshold
        opt = get_options("fa_threshold");
        if (opt.size())
            fa_threshold = opt[0][0];
        // Create a 1-voxel buffer to hold the single fibre voxels the base intensity will be estimated
        // against
        Triple<double> vox_lengths(dwi.vox(X), dwi.vox(Y), dwi.vox(Z));
        Triple<size_t> dims(1.0, 1.0, 1.0);
        Triple<double> offsets(0.0, 0.0, 0.0);
        MR::Math::Matrix<float> encodings = dwi.get_DW_scheme();
        size_t num_encodings = encodings.rows();
        // Create the matrix used to calculate the diffusion tensor
        // The order of the rows is d_xx, d_yy, d_zz, d_xy, d_xz, dyz
        MR::Math::Matrix<double> tensor_encodings(num_encodings, 6);
        for (size_t dim_i = 0; dim_i < 3; ++dim_i) {
            // Do the d_xx, d_yy and d_zz columns
            tensor_encodings.column(dim_i) = encodings.column(dim_i);
            tensor_encodings.column(dim_i) *= encodings.column(dim_i);
            // Make the d_xy and d_yz columns
            if (dim_i < 2) {
                tensor_encodings.column(dim_i * 2 + 3) = encodings.column(dim_i);
                tensor_encodings.column(dim_i * 2 + 3) *= encodings.column(dim_i + 1);
            }
        }
        // Finish of with the d_xz column
        tensor_encodings.column(4) = encodings.column(X);
        tensor_encodings.column(4) *= encodings.column(Z);
        // Caluculate the singular value decomposition on the tensor encodings matrix
        Math::USV usv = Math::svd(tensor_encodings);
        // Reference the b_values for readibility
        MR::Math::Vector<float>& b_values = encodings.column(DW);
        // Create the expected image voxel
        Diffusion::Model diffusion_model = Diffusion::Model::factory(encodings, diff_response_SH,
                diff_adc, diff_fa, diff_isotropic, diff_warn_b_mismatch);
        Image::Expected::Buffer* exp_image = Image::Expected::Buffer::factory(exp_type, dims,
                vox_lengths, diffusion_model, exp_num_length_sections, exp_num_width_sections,
                exp_interp_extent, offsets, exp_enforce_bounds, exp_half_width);
        double est_base_intensity = 0.0;
        // Loop through all the provided single fibre indices and calculate the optimum b0
        size_t num_samples = argument.size() - 1;
        for (size_t sample_i = 0; sample_i < (argument.size() - 1); ++sample_i) {
            // The following is quite a clunky way to get an observed buffer consisting of a single
            // voxel from the index given as a triple. NB the first argument to this script is the
            // image location hence the offset by one.
            MR::Math::Vector<double> observed(num_encodings);
            Triple<size_t> ref_index = parse_triple<size_t>(std::string(argument[sample_i + 1]));
            size_t index[3];
            ref_index.copy_to(index);
            MR::Image::Voxel<float> vox(dwi);
            MR::DataSet::Subset<MR::Image::Voxel<float> > vox_subset(vox, index, ONE_VOX);
            MR::DataSet::Loop loop(0, 3);
            size_t encode_i = 0;
            for (loop.start(vox_subset); loop.ok(); loop.next(vox_subset))
                observed[encode_i] = (double)vox_subset.value();
            // Calculate the diffusion tensor from the observed intensities
            MR::Math::Vector<double> observed_scaled_logged(num_encodings), tensor_vec(6);
            for (size_t encode_i = 0; encode_i < num_encodings; ++encode_i)
                observed_scaled_logged[encode_i] = -MR::Math::log(observed[encode_i])
                        / b_values[encode_i];
            Math::solve_psuedo_inverse(observed_scaled_logged, tensor_vec, usv);
            MR::Math::Matrix<double> tensor(3, 3);
            for (size_t dim_i = 0; dim_i < 3; ++dim_i) {
                tensor(dim_i, dim_i) = tensor_vec[dim_i];
                // Make the d_xy and d_yz columns
                if (dim_i < 2)
                    tensor(dim_i, dim_i + 1) = tensor(dim_i + 1, dim_i) = tensor_vec[dim_i * 2 + 3];
            }
            tensor(X, Z) = tensor(Z, X) = tensor_vec[4];
            // Get the eigenvectors of the tensor matrix
            MR::Math::Matrix<double> evec(3);
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
            if (fa < fa_threshold)
                throw Exception("Sample index '" + str(ref_index) + "' has a FA of less than "
                           + str(fa_threshold) + " (" + str(fa) + "). It will be omitted from the "
                           "estimation.")
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
            for (size_t encode_i = 0; encode_i < num_encodings; ++encode_i)
                est_base_intensity += observed[encode_i] / exp_image(0, 0, 0)[encode_i];
        }

        est_base_intensity /= num_encodings * (argument.size() - 1);

        std::cout << est_base_intensity;
    }
