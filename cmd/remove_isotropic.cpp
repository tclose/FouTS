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
  "Subtracts the isotropic components from the image.",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("initial", "The initial image the second image will be subtracted from.").type_file(),

  Argument ("output_image", "The resulting image").optional().type_file (),

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

  std::string input_location    = argument[0];
  std::string output_location   = input_location;

  if (argument.size() > 1)
    output_location   = argument[1].c_str();

  Image::Observed::Buffer image(input_location);

  MR::ProgressBar progress_bar ("Removing isotropic components from image...");

  image.remove_isotropic();

  image.properties()["istropic_removed"] = "1";


//-------------//
//  Save Image //
//-------------//

  image.save(output_location);


}


