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

#include "bts/image/expected/buffer.h"
#include "bts/image/expected/trilinear/buffer.h"
#include "bts/image/expected/gaussian/buffer.h"
#include "bts/image/expected/quartic/buffer.h"
#include "bts/diffusion/model.h"

#include "bts/inline_functions.h"

using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

const size_t ONE_VOX[3] = {1, 1, 1}

DESCRIPTION = {
  "Finds the maximum b0 intensity of an image",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("dwi", "the input diffusion-weighted image.").type_image_in(),
  Argument ("indices", "the voxel indices in the image where the base intensity "
                       "is calculated from. Should be in regions with known "
                       "single fibre populations.").type_text().allow_multiple(),
  Argument()
};


OPTIONS = {

    Option ("grad",
    "specify the diffusion-weighted gradient scheme used in the acquisition. "
    "The program will normally attempt to use the encoding stored in the image "
    "header. This should be supplied as a 4xN text file with each line is in "
    "the format [ X Y Z b ], where [ X Y Z ] describe the direction of the "
    "applied gradient, and b gives the b-value in units (1000 s/mm^2).")
    + Argument ("encoding").type_file(),

    DIFFUSION_PARAMETERS,

    EXPECTED_IMAGE_PARAMETERS,

    Option ("img_dims", "The number of voxels along each dimension")
     + Argument ("img_dims", "").type_text (str(DIMS_DEFAULT).c_str()),

    Option()

};



EXECUTE {

  MR::Image::Header dwi (argument[0]);

  if (dwi.ndim() != 4)
    throw Exception ("dwi image should contain 4 dimensions");

  MR::Math::Matrix<float> grad;

  // Loads parameters to construct Diffusion::Model ('diff_' prefix)
  SET_DIFFUSION_PARAMETERS;

  // Loads extra parameters to construct Image::Expected::*::Buffer ('exp_' prefix)
  SET_EXPECTED_IMAGE_PARAMETERS;

  Triple<size_t> dims = DIMS_DEFAULT;

  Options img_opt = get_options("img_dims");
  if (img_opt.size())
    dims = parse_triple<size_t>(std::string(img_opt[0][0]));

  Options opt = get_options ("grad");
  if (opt.size())
    grad.load (opt[0][0]);
  else {
    if (!dwi.DW_scheme().is_set())
      throw Exception ("no diffusion encoding found in image \"" + dwi.name() + "\"");
    grad = dwi.DW_scheme();
  }

  if (grad.rows() < 7 || grad.columns() != 4)
    throw Exception ("unexpected diffusion encoding matrix dimensions");

  if (dwi.dim(3) != (int) grad.rows())
    throw Exception ("number of studies in base image does not match that in encoding file");

  MR::DWI::normalise_grad (grad);

  std::vector<int> bzeros, dwis;
  MR::DWI::guess_DW_directions (dwis, bzeros, grad);

  MR::Image::Voxel<float> vox (dwi);


  for (size_t ref_i = 1; ref_i < argument.size(); ++ref_i) {
    Triple<size_t> ref_index = parse_triple<size_t>(std::string(argument[ref_i]));
    size_t index[3];
    ref_index.copy_to(index);
    MR::DataSet::Subset<MR::Image::Voxel<float> > vox_subset(vox, index, ONE_VOX);
    MR::DataSet::Loop loop(0,3);
    BTS::Image::Observed::Buffer buff(Triple<size_t>::ONES, )
    for (loop.start(vox); loop.ok(); loop.next(vox)) {
    }

    Coord interp_length = this->interp_extent * this->vox_lengths();
    // Create tract spans the interpolation length of the interpolation kernel which is centred on the
    // bottom left voxel.
    Fibre::Tractlet::Set tcts (1,2);
    tcts.zero();
    tcts[0](0,0) = this->corner_offsets + this->vox_lengths() / 2.0;
    tcts[0](0,1) = Coord(interp_length[0] / M_SQRT2, 0.0, 0.0);
    tcts[0](1,0) = Coord(0.0, interp_length[1] * M_SQRT2, 0.0);
    tcts[0](2,0) = Coord(0.0, 0.0, interp_length[2] * M_SQRT2);
    // Normalize the density of the tract and set the base_intensity of the set to 1.0, to calculate the required
    // base intensity value to match that of the reference.
    tcts.normalise_densities();
    tcts.set_base_intensity(1.0);
    expected_image<Fibre::Tractlet>(tcts);
    // Divide the reference b0 by the value in the test voxel in the bottom left corner.
    base_intensity = ref_b0 / this->operator()(0,0,0).b0();
    clear();





  }

  double estimated_base_intensity = 0.0;


  std::cout << estimated_base_intensity << std::endl;

}


