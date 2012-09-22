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
  "Selects a rectangular ROI from a larger image",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("initial", "The initial image the second image will be subtracted from.").type_file(),

  Argument ("output_image", "The resulting image").type_file (),

  Argument()
};


OPTIONS = {

  NOISE_PARAMETERS,

  Option ("seed", "Random seed used to initialise the random generator.")
   + Argument ("", "").type_float (1e-9, INFINITY, 0.1),

  Option ("dim", "dimensions of the ")
    + Argument ("", "").type_float (1e-9, INFINITY, 0.1),

Option()

};



EXECUTE {


//------------------------------//
//  Load Input Image Buffer //
//------------------------------//



  size_t seed = time(NULL);

  Options opt;

  SET_NOISE_PARAMETERS;


  opt = get_options("seed");
  if (opt.size())
    seed = opt[0][0];


//--------------------------//
//  Load Input Image Buffer //
//--------------------------//

  std::string input_location    = argument[0];
  std::string output_location   = argument[1];

  Image::Observed::Buffer image(input_location);

//------------//
//  Add noise //
//------------//


  gsl_rng* rand_gen = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rand_gen, seed);

  Image::Noise* noise = Image::Noise::factory(rand_gen, noise_type, noise_snr, noise_ref_signal);

  noise->noisify(image);

//---------------------------//
//  Add Properties to header //
//---------------------------//

  ADD_NOISE_PROPERTIES(image.properties());

  image.properties()["seed"] = str(seed);

//-------------//
//  Save Image //
//-------------//

  image.save(output_location);


}


