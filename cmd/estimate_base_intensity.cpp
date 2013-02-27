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

#include "bts/inline_functions.h"

using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

const char* METHOD_DEFAULT = "median";

DESCRIPTION = {
  "Finds the maximum b0 intensity of an image",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("dwi", "the input diffusion-weighted image.").type_image_in(),
  Argument ("samples", "the locations in the image where the base intensity "
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

    Option("method", "Used to select the method used to get the b0, one of 'average', 'max', or 'median' (default).")
    + Argument().type_text(METHOD_DEFAULT),

    Option()

};



EXECUTE {

  MR::Image::Header dwi_header (argument[0]);

  if (dwi_header.ndim() != 4)
    throw Exception ("dwi image should contain 4 dimensions");

  MR::Math::Matrix<float> grad;
  std::string method = METHOD_DEFAULT;

  Options opt = get_options ("grad");
  if (opt.size())
    grad.load (opt[0][0]);
  else {
    if (!dwi_header.DW_scheme().is_set())
      throw Exception ("no diffusion encoding found in image \"" + dwi_header.name() + "\"");
    grad = dwi_header.DW_scheme();
  }

  opt = get_options("method");
  if (opt.size())
    method = opt[0][0].c_str();

  if (grad.rows() < 7 || grad.columns() != 4)
    throw Exception ("unexpected diffusion encoding matrix dimensions");

  if (dwi_header.dim(3) != (int) grad.rows())
    throw Exception ("number of studies in base image does not match that in encoding file");

  MR::DWI::normalise_grad (grad);

  std::vector<int> bzeros, dwis;
  MR::DWI::guess_DW_directions (dwis, bzeros, grad);

  MR::Image::Header mask_header (argument[1]);
  if (mask_header.dim(0) != dwi_header.dim(0) ||
      mask_header.dim(1) != dwi_header.dim(1) ||
      mask_header.dim(2) != dwi_header.dim(2))
    throw Exception ("mask & DWI image dimensions do not match");

  MR::Image::Voxel<float> dwi (dwi_header);
  MR::Image::Voxel<float> mask (mask_header);

  std::vector<double> all_bzeros;
  // Just to make it a little quicker reserve some places for the vector.
  all_bzeros.reserve(bzeros.size() * 1e4);
  MR::DataSet::Loop loop(0,3);
  for (loop.start(mask, dwi); loop.ok(); loop.next(mask, dwi)) {

  }
  double estimated_base_intensity;


  std::cout << estimated_base_intensity << std::endl;

}


