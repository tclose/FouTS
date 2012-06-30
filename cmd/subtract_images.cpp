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




#include "bts/cmd.h"
#include "progressbar.h"

#include "bts/common.h"
#include "bts/image/expected/buffer.h"
#include "bts/image/observed/buffer.h"




#include "bts/inline_functions.h"


using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Subtractlets one image from another",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("initial", "The initial image the second image will be subtracted from.").type_image_in(),

  Argument ("subtract", "The image that will be subtracted from the initial image.").type_image_in(),

  Argument ("output_image", "The resulting image").type_image_out (),

  Argument()
};


OPTIONS = {

Option()

};



EXECUTE {


//------------------------------//
//  Load Input Image Buffer //
//------------------------------//


  std::string image_location = argument[0];
  std::string subtract_image_location = argument[1];
  std::string output_location = argument[2];

  Image::Observed::Buffer image(image_location);
  Image::Observed::Buffer subtract_image (subtract_image_location);

  image -= subtract_image;

  image.save(output_location);


}






