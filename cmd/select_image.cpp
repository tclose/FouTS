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

#include "image/header.h"
#include "image/voxel.h"
#include "dataset/subset.h"
#include "dataset/loop.h"

#include "bts/inline_functions.h"

using namespace BTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
  "Selects a rectangular region of interest from a larger image.",
  "",
  NULL
};

ARGUMENTS= {
  Argument ("initial", "The initial image the second image will be subtracted from.").type_image_in(),

  Argument ("output_image", "The resulting image").type_file (),

  Argument()
};

OPTIONS= {

  Option ("dims", "The number of voxels along each dimension")
  + Argument ("dims", "").type_text ("[3,3,3]"),

  Option ("offsets", "The offset of the centre of the image from the origin (0,0,0).")
  + Argument ("offsets", "").type_text ("auto"),

  Option("remove_isotropic", "Removes the isotropic component of the signal (assumes that the signal is acquired at a \
      constant b value"),

  Option()

};

EXECUTE {

    Triple<size_t> dims(3, 3, 3);
    Triple<size_t> offsets(0, 0, 0);
    bool remove_isotropic = false;

    Options opt = get_options("dims");
    if (opt.size())
      dims = parse_triple<size_t>(std::string(opt[0][0]));

    opt = get_options("offsets");
    if (opt.size())
      offsets = parse_triple<size_t>(std::string(opt[0][0]));

    opt = get_options("remove_isotropic");
    if (opt.size())
      remove_isotropic = true;

    MR::Image::Header in(argument[0]);
    std::string output_location = argument[1];

    Diffusion::Encoding::Set dw_scheme(in.get_DW_scheme());

    Triple<double> spatial_offset(in.vox(X) * offsets[X] + in.transform()(X,3),
        in.vox(Y) * offsets[Y] + in.transform()(Y,3),
        in.vox(Z) * offsets[Z] + in.transform()(Z,3));

    Image::Observed::Buffer out(dims,
        Triple<double>(in.vox(X), in.vox(Y), in.vox(Z)), spatial_offset,
        dw_scheme);

    size_t dim_array[] =
        { dims[X], dims[Y], dims[Z], dw_scheme.num_encodings() };
    size_t offset_array[] = { offsets[X], offsets[Y], offsets[Z], 0 };

    MR::Image::Voxel<float> vox(in);
    MR::DataSet::Subset<MR::Image::Voxel<float> > vox_subset(vox, offset_array,
        dim_array);

    MR::ProgressBar progress_bar("Selecting subset of image...",
        dw_scheme.num_encodings());
    MR::DataSet::Loop loop(0, 4);
    size_t x, y, z, encoding_i;
    x = y = z = encoding_i = 0;
    for (loop.start(vox_subset); loop.ok(); loop.next(vox_subset)) {
      if (x >= dims[X]) {
        x = 0;
        ++y;
        if (y >= dims[Y]) {
          y = 0;
          ++z;
          if (z >= dims[Z]) {
            z = 0;
            ++encoding_i;
            ++progress_bar;
            if (encoding_i >= dw_scheme.num_encodings())
              throw Exception("Incremented past the end of the number of "
                              "encodings dimension, something has gone wrong.");
          }
        }
      }
      out(x, y, z)[encoding_i] = (float) vox_subset.value();
      ++x;
    }

    out.properties()["original_image"] = argument[0].c_str();
    out.properties()["offset_from_original"] = str(offsets);

    if (remove_isotropic)
      out.remove_isotropic();

//-------------//
//  Save Image //
//-------------//
    out.save(output_location);

  }

