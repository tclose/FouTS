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

#include "bts/math/common.h"

#include "image/header.h"
#include "image/header.h"
#include "dwi/gradient.h"
#include "progressbar.h"

#include "bts/common.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"

#include "bts/image/observed/buffer.h"
#include "bts/image/expected/buffer.h"
#include "bts/image/expected/trilinear/buffer.h"
#include "bts/image/expected/gaussian/buffer.h"
#include "bts/image/expected/quartic/buffer.h"
#include "bts/diffusion/model.h"

#include "bts/image/noise/gaussian.h"
#include "bts/image/noise/rician.h"



#include "bts/inline_functions.h"


using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Generates a set of DW-BTS images from a given strand or tractlet configuration.",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("input", "The strands file the image set will be expected_image from.").type_file (),
  Argument ("output", "The generated DW-MR image set").optional().type_image_out (),
  Argument()
};


OPTIONS = {
  
  Option ("clean", "No noise added to the image."),

  Option ("seed", "The random seed that is passed to the random generator")
   + Argument ("seed", ""),

  Option ("degree", "Degree of strands when converted from tracks")
   + Argument ("degree", "").type_integer (0, Fibre::Strand::DEFAULT_DEGREE, LARGE_INT),

  DIFFUSION_PARAMETERS,

  IMAGE_PARAMETERS,

  EXPECTED_IMAGE_PARAMETERS,
  
  NOISE_PARAMETERS,

  COMMON_PARAMETERS,  
  

Option() };



EXECUTE {


//--------------//
//  Get options //
//--------------//

  std::string input_location = argument[0];
  std::string output_location;
  
  if (argument.size() > 1)
    output_location = argument[1].c_str();
  else
    output_location = File::strip_extension (input_location) + ".mif";

  bool clean = false;
  size_t seed = time(0);
  size_t degree = Fibre::Strand::DEFAULT_DEGREE;


  Options opt = get_options("clean");
  if (opt.size())
    clean = true;
  
  opt = get_options("seed");
   if (opt.size()) {
     std::string seed_str = opt[0][0];
     seed = to<size_t>(seed_str);
   }

  opt = get_options("degree");
  if (opt.size())
    degree = opt[0][0];

  // Loads parameters to construct Diffusion::Model ('diff_' prefix)
  SET_DIFFUSION_PARAMETERS;

  // Loads parameters to construct Image::Expected::*::Buffer that are inherited from Image::Observed::Buffer ('img_' prefix)
  SET_IMAGE_PARAMETERS;

  // Loads extra parameters to construct Image::Expected::*::Buffer ('img_' prefix)
  SET_EXPECTED_IMAGE_PARAMETERS;

  // Loads parameters to construct Image::Noise::* ('noise_' prefix)
  SET_NOISE_PARAMETERS;

  // Loads parameters that are common to all commands.
  SET_COMMON_PARAMETERS;


//----------------------------------------------------//
//  Auto generate parameters if not set or irrelevant //
//----------------------------------------------------//

  if (!img_offsets.valid())
    img_offsets = Image::Observed::Buffer::default_corner_offset(img_dims, img_vox_lengths);


  exp_enforce_bounds = true;

  if (!diff_isotropic && (noise_type == "rician"))
    throw Exception ("Rician noise cannot be added to images with isotropic components removed.");

//-------------------------//
//  Set up diffusion model //
//-------------------------//

  Diffusion::Model diffusion_model = Diffusion::Model::factory (diff_encodings,
                                                                diff_response_SH,
                                                                diff_adc,
                                                                diff_fa,
                                                                diff_isotropic);


//----------------------------------//
//  Initialize expected_image image //
//----------------------------------//



  Image::Expected::Buffer* image = Image::Expected::Buffer::factory(exp_type,
                                                        img_dims,
                                                        img_vox_lengths,
                                                        diffusion_model,
                                                        exp_num_length_sections,
                                                        exp_num_width_sections,
                                                        exp_interp_extent,
                                                        img_offsets,
                                                        exp_enforce_bounds,
                                                        exp_half_width);



//-----------------//
// Fabricate image //
//-----------------//


  MR::ProgressBar progress_bar ("Generating image from provided fibres...");

  if (File::has_or_txt_extension<Fibre::Strand>(input_location) || File::has_or_txt_extension<Fibre::Track>(input_location)) {

    Fibre::Strand::Set strands (input_location, degree);

    //exp_base_intensity is recorded in header so needs to be synchronised with base_intensity of strands
    if (exp_base_intensity < 0)
      exp_base_intensity = strands.base_intensity();
    else
      strands.set_base_intensity(exp_base_intensity);

    image->expected_image(strands);

  } else if (File::has_or_txt_extension<Fibre::Tractlet>(input_location)) {

    Fibre::Tractlet::Set tractlets (input_location);

    //exp_base_intensity is recorded in header so needs to be synchronised with base_intensity of tractlets
    if (exp_base_intensity < 0)
      exp_base_intensity = tractlets.base_intensity();
    else
      tractlets.set_base_intensity(exp_base_intensity);

    image->expected_image(tractlets);

  } else
    throw Exception ("Unrecognised extension '" + input_location + "'.");

  image->properties()["type"]                  = "synthetic";
  image->properties()["method"]                = "generate_image";
  image->properties()["state_location"]        = input_location;

  image->properties()["software version"]  = version_number_string();
  image->properties()["datetime"]          = current_datetime();


  ADD_DIFFUSION_PROPERTIES(image->properties());

  ADD_IMAGE_PROPERTIES(image->properties());

  ADD_EXPECTED_IMAGE_PROPERTIES(image->properties());

  ADD_COMMON_PROPERTIES(image->properties());


  //If '-clean' option is not selected noise is be added to the image.
  if (!clean) {

    gsl_rng* rand_gen = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rand_gen, seed);

    if (isnan(noise_ref_signal))
      noise_ref_signal = image->max_b0();

    Image::Noise* noise = Image::Noise::factory(rand_gen, noise_type, noise_snr, noise_ref_signal);

    noise->noisify(*image);

    ADD_NOISE_PROPERTIES(image->properties());
  }


  //------------//
  // Save image //
  //------------//


  image->save(output_location);



}




