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


extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

#include "bts/common.h"
#include "bts/cmd.h"
#include "progressbar.h"

#include "bts/common.h"
#include "bts/image/expected/buffer.h"
#include "bts/image/observed/buffer.h"

#include "bts/image/noise.h"
#include "bts/image/noise/gaussian.h"



#include "bts/inline_functions.h"

using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Selects a rectangular region of interest from a larger image.",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("initial", "The initial image the second image will be subtracted from.").type_image_in(),

  Argument ("output_image", "The resulting image").type_file (),

  Argument()
};


OPTIONS = {

  Option ("dims", "The number of voxels along each dimension")
   + Argument ("dims", "").type_text ("[3,3,3]"),

  Option ("offsets", "The offset of the centre of the image from the origin (0,0,0).")
   + Argument ("offsets", "").type_text ("auto"),


Option()

};



EXECUTE {


  Triple<size_t>    dims (3,3,3);
  Triple<size_t>    offsets (0,0,0);

  Options opt = get_options("dims");
  if (opt.size())
    dims = parse_triple<size_t>(std::string(opt[0][0]));

  opt = get_options("offsets");
  if (opt.size())
    offsets = parse_triple<size_t>(std::string(opt[0][0]));


  MR::Image::Header in (argument[0]);
  std::string output_location   = argument[1];

  MR::Math::Matrix<float> dw_scheme = in.get_DW_scheme();

  Triple<double> spatial_offset(in.vox(X) * offsets[X], in.vox(Y) * offsets[Y], in.vox(Z) * offsets[Z]);

  Image::Observed::Buffer out(dims, Triple<double>(in.vox(X), in.vox(Y), in.vox(Z)), spatial_offset, dw_scheme);

  Image::Voxel<float> dwi (in);

  for (size_t x = 0; x < dims[X]; ++x)
    for (size_t y = 0; y < dims[Y]; ++y)
      for (size_t z = 0; z < dims[Z]; ++z)
        out(x,y,z) = in(x + offsets[X], y + offsets[Y], z + offsets[Z]);


  out.properties()["original_image"] = argument[0].c_str();
  out.properties()["offsets"] = str(offsets);


//-------------//
//  Save Image //
//-------------//

  out.save(output_location);


}


