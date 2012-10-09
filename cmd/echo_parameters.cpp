/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, Jul 6, 2010.

    This file is part of MRtrix.

    MRtrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRtrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

*/
/*

#include "bts/cmd.h"

#include "bts/common.h"
#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/strand/set.h"
#include "bts/image/observed/buffer.h"

#include "bts/mcmc/riemannian.h"
#include "bts/mcmc/hamiltonian.h"
#include "bts/mcmc/metropolis.h"
#include "bts/analysis/scan.h"
#include "bts/analysis/gradient_tester.h"
#include "bts/image/expected/buffer.h"

#include "bts/prob/likelihood.h"

#include "bts/prob/prior.h"





#include "bts/mcmc/proposal/momentum.h"
#include "bts/mcmc/proposal/distribution.h"
#include "bts/mcmc/proposal/distribution/gaussian.h"

#include "bts/fibre/tractlet/set/momentum.h"
#include "bts/fibre/strand/set/momentum.h"
#include "bts/fibre/strand/set/walker.h"

#include "bts/fibre/tractlet/set/walker.h"


#include "image/header.h"
#include "bts/image/observed/buffer.h"


#include "bts/inline_functions.h"


using namespace BTS;


void                              check_parameter(std::map<std::string, std::string>& echo_options,
                                                  const std::string& key,
                                                  const std::string& default_value);

void                              check_parameters_group ( std::map<std::string, std::string>& echo_options,
                                                           const std::string& name,
                                                           std::vector< std::pair<std::string, std::string> > key_default_value_pairs,
                                                           bool default_active);


template <typename T> void        load_and_check_parameters(std::map<std::string, std::string>& echo_options,
                                                            const std::string& key,
                                                            const std::string& default_value_location = "");

void                              load_and_check_state_parameters ( std::map<std::string, std::string>& echo_options,
                                                                    const std::string& key,
                                                                    const std::vector<std::string>& default_value_locations = std::vector<std::string>(3));

bool                               is_riemannian(const std::string& method)
  { return method == "riemannian"; }

bool                               is_hamiltonian(const std::string& method)
  { return method == "hamiltonian"; }

bool                               is_metropolis(const std::string& method)
  { return method == "metropolis"; }

bool                               is_scan(const std::string& method)
  { return method == "scan" || method == "scan_batch"; }

bool                               is_gradient_test(const std::string& method)
  { return method == "test_gradient"; }

bool                               is_mcmc(const std::string& method)
  { return is_metropolis(method) || is_hamiltonian(method) || is_riemannian(method); }

bool                               is_analysis(const std::string& method)
  { return is_scan(method) || is_gradient_test(method); }



SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Echo Parameters",
  "Determine which parameters in a previous run need to be set to exactly copy the run",
  NULL
};

ARGUMENTS = {
  Argument ("template_location", "Template location", "Location of template file, from which to copy the parameters.").type_file (),

  Argument()
};


OPTIONS = {

  Option ("method", "The output method", "The output method to print if not the original method.")
   + Argument ("method", ""),

  Option ("output_location", "Output location", "Override the location of the output.")
   + Argument ("output_location", ""),

  Option ("state_location", "New state location", "Override the location of the state.")
   + Argument ("state_location", ""),

  Option ("seed_copy", "Copy seed", "Also copy the seed."),

  Option ("last_place_output", "Put the output location last", "Put the output location last for easy access to change the name."),

  Option ("origin_location", "New origin location", "Override the location of the origin.")
   + Argument ("origin_location", ""),

  Option ("axis1_location", "Axis 1 location", "Pass the location for axis1 for use in mcmc -> scan echo's.")
   + Argument ("axis1_location", ""),

  Option ("axis2_location", "Axis 2 location", "Pass the location for axis2 for use in mcmc -> scan echo's.")
   + Argument ("axis2_location", ""),

  Option ("axis3_location", "Axis 3 location", "Pass the location for axis3 for use in mcmc -> scan echo's.")
   + Argument ("axis3_location", ""),
Option() };


EXECUTE {


  std::vector<std::string> echo_arguments;
  std::map<std::string, std::string> echo_options;






//------------------------------//
//  Load the template echo_options.  //
//------------------------------//

  std::string template_location = argument[0];


  if (File::has_extension<Fibre::Strand::Set>(template_location)) {

    Fibre::Strand::Set::Reader reader(template_location);

    echo_options = reader.get_extend_props();

  } else if (File::has_extension<Fibre::Tractlet::Set>(template_location)) {

    Fibre::Tractlet::Set::Reader reader(template_location);

    echo_options = reader.get_extend_props();

  } else if (File::has_extension<Image::Expected::Buffer>(template_location)) {

    MR::Image::Header header;

    header.open(template_location);

    echo_options = header;

  } else
    throw Exception ("Unrecognised extension '" + File::extension(template_location) + "'.");




  //-----------------------//
  //  Legacy conversions.  //
  //-----------------------//

  if (echo_options.count("origin_state")) {
    echo_options["origin"] = echo_options["origin_state"];
    echo_options.erase("origin_state");
  }

  if (echo_options.count("momen_rand_scales"))
    echo_options.erase("momen_rand_scales");

  if (echo_options.count("momen_step_sizes")) {
    echo_options["momen_step"] = echo_options["momen_step_sizes"];
    echo_options.erase("momen_step_sizes");
  }


  if (echo_options.count("leapfrog")) {
    echo_options["momen_step"] = echo_options["leapfrog"];
    echo_options.erase("leapfrog");
  }

  echo_options.erase("momen_over_rand_scale");
  echo_options.erase("momen_rand_scales_location");
  echo_options.erase("momen_over_randomize");
  echo_options.erase("momen_randomize_location");
  echo_options.erase("momen_randomize");

  if (echo_options.count("momen_over_step_size")) {
    echo_options["momen_step_scale"] = echo_options["momen_over_step_size"];
    echo_options.erase("momen_over_step_size");
  }


  if (echo_options.count("momen_step_sizes_location")) {
    echo_options["momen_step_location"] = echo_options["momen_step_sizes_location"];
    echo_options.erase("momen_step_sizes_location");
  }


  if (echo_options.count("leapfrog_scale")) {
    echo_options["momen_step_scale"] = echo_options["leapfrog_scale"];
    echo_options.erase("leapfrog_scale");
  }


  if (echo_options.count("leapfrog_location")) {
    echo_options["momen_step_location"] = echo_options["leapfrog_location"];
    echo_options.erase("leapfrog_location");
  }


  if (echo_options.count("reference_image")) {
    echo_options["obs_image"] = echo_options["reference_image"];
    echo_options.erase("reference_image");
  }

  if (echo_options.count("ref_image")) {
    echo_options["obs_image"] = echo_options["ref_image"];
    echo_options.erase("ref_image");
  }

//---------------------------------------//
//  Set the template and output method.  //
//---------------------------------------//

  std::string loaded_method;
  loaded_method = echo_options["method"];


  //Check for legacy names for methods.
  if (loaded_method == "Hamiltonian")
    loaded_method = "hamiltonian";

  else if (loaded_method == "Metropolis Hastings")
    loaded_method = "metropolis";

  else if (loaded_method == "Scan")
    loaded_method = "scan";

  else if (loaded_method == "Container Test")
    loaded_method = "test_gradient";


  std::string method = loaded_method;

  //Get the method.
  Options opt = get_options("method");
  if (opt.size())
    method = opt[0][0];


  echo_options.erase("method");



//----------------------------//
//  Get the Output location.  //
//----------------------------//

  std::string output_location;

  //If it is the same method, default to overide previous file.
  if (loaded_method == method)
    output_location = template_location;
  else
    output_location = "__output_location_place_holder__";

  opt = get_options("output_location");
  if (opt.size())
    output_location = opt[0][0];



//------------------------------------//
//  Get the Observed Image location.  //
//------------------------------------//

  std::string obs_image_location;
  std::string origin_location;

  if (echo_options.count("obs_image")) {
    obs_image_location = echo_options["obs_image"];
    echo_options.erase("obs_image");

  } else
    obs_image_location = "__obs_image_place_holder__";


  //Load the observed image to check the auto generated parameters.
  Image::Observed::Buffer obs_image;

  if (obs_image_location != "__obs_image_place_holder__") {

    try {

      obs_image.load(obs_image_location);
      origin_location = obs_image.properties()["state_location"];

    } catch(...) {

      std::cout << "Warning!! Could not load image from observed image location." << std::endl;

    }

  }


//---------------------------//
//  Get the state location.  //
//---------------------------//

  std::string state_location, state_key;
  bool overide_state_location = false;

  if (is_mcmc(loaded_method))
    state_key = "initial_state";
  else if (is_scan(loaded_method))
    state_key = "origin";

  opt = get_options("overide_state_location");
  if (opt.size()) {
    state_location = opt[0][0];
    overide_state_location = true;
    echo_options.erase(state_key);
    echo_options.erase(state_key + "_location");
  } else {

    state_location = echo_options[state_key + "_location"];

    if (File::has_extension<Fibre::Strand>(state_location))
      load_and_check_parameters<Fibre::Strand::Set>(echo_options, state_key);
    else if (File::has_extension<Fibre::Tractlet>(state_location))
      load_and_check_parameters<Fibre::Tractlet::Set>(echo_options, state_key);

  }


//------------------------------------//
//  Set the random seed if required.  //
//------------------------------------//

  bool seed_copy = false;
  opt = get_options("seed_copy");
  if (opt.size())
    seed_copy = true;


  if (echo_options.count("random_seed")) {
    echo_options["seed"] = echo_options["random_seed"];
    echo_options.erase("random_seed");
  }

  if (seed_copy) {

    if (!echo_options.count("seed"))
      echo_options["seed"] = "__seed_place_holder__";

  } else
    echo_options.erase("seed");


//----------------------------------------------//
//  Get the remaining options for this program  //
//----------------------------------------------//

  bool set_output_last = 0;
  std::string axis1_location;
  std::string axis2_location;
  std::string axis3_location;

  opt = get_options("last_place_output");
  if (opt.size())
    set_output_last = 1;

  opt = get_options("origin_location");
  if (opt.size())
    origin_location = opt[0][0];

  opt = get_options("axis1_location");
  if (opt.size())
    axis1_location = opt[0][0];

  opt = get_options("axis2_location");
  if (opt.size())
    axis2_location = opt[0][0];

  opt = get_options("axis3_location");
  if (opt.size())
    axis3_location = opt[0][0];
//-----------------------------------//
//  Erase non-parameter properties.  //
//-----------------------------------//

  echo_options.erase("datetime");
  echo_options.erase("count");
  echo_options.erase("total_count");
  echo_options.erase("method");
  echo_options.erase("software version");


//-----------------------------------//
//  Check base density parameter.  //
//-----------------------------------//

  if (obs_image_location != "__obs_image_place_holder__" && state_location != "") {

    MR::Math::Matrix<double> diff_encodings;
    Image::Expected::Buffer::Buffer* exp_image;
    MR::Math::Vector<double> diff_response_SH;
    bool can_check_base_intensity = true;


    bool load_encodings = false;

    if (echo_options.count("diff_encodings"))
      load_encodings = echo_options["diff_encodings"] == "from observed image" || echo_options["diff_encodings"] == "from reference image";
    else if (echo_options.count("diff_encodings_location"))
      load_encodings = echo_options["diff_encodings_location"] == "from observed image" || echo_options["diff_encodings_location"] == "from reference image";

    if (load_encodings) {

      try {

        MR::Image::Header header;

        header.open(obs_image_location);

        //If gradient scheme is included in reference image header, use that instead of default (NB: Will override any gradients passed to '-diff_encodings' option).
        diff_encodings = header.get_DW_scheme();

      } catch(...) {

        std::cout << "Warning!! Could not load image from observed image location." << std::endl;

      }

    } else if (echo_options.count("diff_encodings_location")) {

      diff_encodings.load(echo_options["diff_encodings_location"]);

    } else
      can_check_base_intensity = false;


    Diffusion::Model diffusion_model;
    bool diff_isotropic = false;

    if (echo_options.count("diff_isotropic"))
      diff_isotropic = to<bool>(echo_options["diff_isotropic"]);

    if (echo_options.count("diff_response_SH_location")) {

      try {

        diff_response_SH.load(echo_options["diff_response_SH_location"]);

        diffusion_model = Diffusion::Model (diff_encodings, diff_response_SH, diff_isotropic);

      } catch (Exception e) {

        std::cout << "WARNING!! Could not load diffusion responses from location '" << echo_options["diff_response_SH_location"] << "'." << std::endl;
        can_check_base_intensity = false;
      }

    } else if (echo_options.count("diff_adc") && echo_options.count("diff_fa") && echo_options.count("diff_b_value")) {

      diffusion_model = Diffusion::Model (diff_encodings,
                                          to<double>(echo_options["diff_adc"]),
                                          to<double>(echo_options["diff_fa"]),
                                          diff_isotropic);

    } else
      can_check_base_intensity = false;


    if (can_check_base_intensity && state_location.size()) {


      exp_image = Image::Expected::Buffer::Buffer::factory( echo_options["exp_type"],
                                            obs_image,
                                            diffusion_model,
                                            to<size_t>(echo_options["exp_num_length_sections"]),
                                            to<size_t>(echo_options["exp_num_width_sections"]),
                                            to<double>(echo_options["exp_interp_extent"]),
                                            to<size_t>(echo_options["exp_enforce_bounds"]),
                                            to<double>(echo_options["exp_half_width"]));

      double calc_exp_base_intensity = exp_image->base_intensity_default(obs_image, state_location);

      check_parameter(echo_options, "exp_base_intensity", str(calc_exp_base_intensity));

    }

  }

//-----------------------------//
//  Check general parameters.  //
//-----------------------------//

  //General image parameters.
  check_parameter(echo_options, "img_dims", str(Image::Expected::Buffer::DIMS_DEFAULT));
  check_parameter(echo_options, "img_vox_lengths", str(Image::Expected::Buffer::VOX_LENGTHS_DEFAULT));
  check_parameter(echo_options, "img_offsets", str(- obs_image.vox_lengths() * obs_image.dims() * 0.5));

  //Diffusion response parameters.
  check_parameter(echo_options, "diff_encodings_location", "from observed image");
  check_parameter(echo_options, "diff_encodings_location", "from reference image");
  check_parameter(echo_options, "diff_isotropic", str(Diffusion::Model::ISOTROPIC_DEFAULT));


  if (echo_options.count("diff_response_SH_location"))
    load_and_check_parameters< MR::Math::Vector<double> >(echo_options, "diff_response_SH", "");
  else {
    check_parameter(echo_options, "diff_adc", str(Diffusion::Model::ADC_DEFAULT));
    check_parameter(echo_options, "diff_fa", str(Diffusion::Model::FA_DEFAULT));
  }

  //Expected image parameters.
  check_parameter(echo_options, "exp_num_length_sections", str(Image::Expected::Buffer::NUM_LENGTH_SECTIONS_DEFAULT));
  check_parameter(echo_options, "exp_num_width_sections", str(Image::Expected::Buffer::NUM_WIDTH_SECTIONS_DEFAULT));
  check_parameter(echo_options, "exp_type", str(Image::Expected::Buffer::TYPE_DEFAULT));
  check_parameter(echo_options, "exp_interp_extent", str(Image::Expected::Buffer::INTERP_EXTENT_DEFAULT));
  check_parameter(echo_options, "exp_enforce_bounds", str(Image::Expected::Buffer::ENFORCE_BOUNDS_DEFAULT));
  check_parameter(echo_options, "exp_half_width", str(Image::Expected::Buffer::HALF_WIDTH_DEFAULT));



  //Legacy conversion.
  if (echo_options["like_type"] == "strong")
    echo_options["like_type"] = "gaussian";

  //Likelihood parameters.
  check_parameter(echo_options, "like_snr", str(Prob::Likelihood::SNR_DEFAULT));
  check_parameter(echo_options, "like_type", str(Prob::Likelihood::TYPE_DEFAULT));
  check_parameter(echo_options, "like_b0_include", str(Prob::Likelihood::B0_INCLUDE_DEFAULT));
  check_parameter(echo_options, "like_outside_scale", str(Prob::Likelihood::OUTSIDE_SCALE_DEFAULT));

  //Proposal parameters.
  if (File::has_extension<Fibre::Strand>(state_location)) {
    check_parameter(echo_options, "walk_type", str(Fibre::Strand::Set::Walker::TYPE_DEFAULT));
    check_parameter(echo_options, "walk_scale", str(Fibre::Strand::Set::Walker::STEP_SCALE_DEFAULT));
    check_parameter(echo_options, "walk_intens_scale", str(Fibre::Strand::Walker::Constrained::INTENSITY_SCALE_DEFAULT));
    check_parameter(echo_options, "walk_base_intens_scale", str(Fibre::Strand::Set::Walker::Constrained::BASE_INTENSITY_SCALE_DEFAULT));
  } else if (File::has_extension<Fibre::Tractlet>(state_location)) {
    check_parameter(echo_options, "walk_type", str(Fibre::Tractlet::Set::Walker::TYPE_DEFAULT));
    check_parameter(echo_options, "walk_scale", str(Fibre::Tractlet::Set::Walker::STEP_SCALE_DEFAULT));
    check_parameter(echo_options, "walk_intens_scale", str(Fibre::Tractlet::Walker::Constrained::INTENSITY_SCALE_DEFAULT));
    check_parameter(echo_options, "walk_base_intens_scale", str(Fibre::Tractlet::Set::Walker::Constrained::BASE_INTENSITY_SCALE_DEFAULT));
  }

  check_parameter(echo_options, "walk_midpoint_scale", str(Fibre::Strand::Walker::Constrained::MIDPOINT_SCALE_DEFAULT));
  check_parameter(echo_options, "walk_endpoint_scale", str(Fibre::Strand::Walker::Constrained::ENDPOINT_SCALE_DEFAULT));
  check_parameter(echo_options, "walk_sheer_scale", str(Fibre::Tractlet::Walker::Constrained::SHEER_SCALE_DEFAULT));
  check_parameter(echo_options, "walk_basewidth_scale", str(Fibre::Tractlet::Walker::Constrained::BASEWIDTH_SCALE_DEFAULT));
  check_parameter(echo_options, "prop_distr_type", str(MCMC::Proposal::Distribution::TYPE_DEFAULT));
  check_parameter(echo_options, "walk_sphere_radius", str(Fibre::Strand::Walker::Constrained::radius_default(obs_image)));


  std::vector<std::string> walk_step_sizes_defaults;
  walk_step_sizes_defaults.push_back(str(Fibre::Strand::Set::Walker::STEP_LOCATION_DEFAULT));
  walk_step_sizes_defaults.push_back(str(Fibre::Tractlet::Set::Walker::STEP_LOCATION_DEFAULT));
  walk_step_sizes_defaults.push_back("");

  load_and_check_state_parameters(echo_options, "walk_step", walk_step_sizes_defaults);

  std::vector<std::string> momen_step_sizes_defaults;
  momen_step_sizes_defaults.push_back(str(Fibre::Strand::Set::Momentum::STEP_LOCATION_DEFAULT));
  momen_step_sizes_defaults.push_back(str(Fibre::Tractlet::Set::Momentum::STEP_LOCATION_DEFAULT));
  momen_step_sizes_defaults.push_back("");

  load_and_check_state_parameters(echo_options, "momen_step", momen_step_sizes_defaults);

  //Prior parameters.
  check_parameter(echo_options, "prior_scale", str(Prob::Prior::SCALE_DEFAULT));
  check_parameter(echo_options, "prior_only", str(0));


  //-----------------------------------//
  //  Check general group parameters.  //
  //-----------------------------------//

  if (echo_options.count("prior_curv_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_curv_scale", str(Prob::Prior::Frequency::SCALE_DEFAULT)));

    check_parameters_group(echo_options, "prior_curv", key_values, Prob::Prior::Frequency::ACTIVE_DEFAULT);


  } else if (Prob::Prior::Frequency::ACTIVE_DEFAULT == Active)
    echo_options["prior_curv"] = "0";


  if (echo_options.count("prior_sheer_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_sheer_scale", str(Prob::Prior::Sheer::SCALE_DEFAULT)));

    check_parameters_group(echo_options, "prior_sheer", key_values, Prob::Prior::Sheer::ACTIVE_DEFAULT);


  } else if (Prob::Prior::Sheer::ACTIVE_DEFAULT == Active)
    echo_options["prior_sheer"] = "0 1";


  if (echo_options.count("prior_end_in_sphere_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_end_in_sphere_scale", str(Prob::Prior::EndInSphere::SCALE_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_end_in_sphere_radius", str(Prob::Prior::EndInSphere::radius_default(obs_image))));


    check_parameters_group(echo_options, "prior_end_in_sphere", key_values, Prob::Prior::EndInSphere::ACTIVE_DEFAULT);


  } else if (Prob::Prior::EndInSphere::ACTIVE_DEFAULT == Active)
    echo_options["prior_end_in_sphere"] = "0 1";


  if (echo_options.count("prior_end_out_sphere_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_end_out_sphere_scale", str(Prob::Prior::EndOutSphere::SCALE_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_end_out_sphere_sd", str(Prob::Prior::EndOutSphere::STD_DEV_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_end_out_sphere_radius", str(Prob::Prior::EndOutSphere::radius_default(obs_image))));


    check_parameters_group(echo_options, "prior_end_out_sphere", key_values, Prob::Prior::EndOutSphere::ACTIVE_DEFAULT);


  } else if (Prob::Prior::EndOutSphere::ACTIVE_DEFAULT == Active)
    echo_options["prior_end_out_sphere"] = "0 1";

  if (echo_options.count("prior_end_on_sphere_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_end_on_sphere_scale", str(Prob::Prior::EndOnSphere::SCALE_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_end_on_sphere_radius", str(Prob::Prior::EndOnSphere::radius_default(obs_image))));


    check_parameters_group(echo_options, "prior_end_on_sphere", key_values, Prob::Prior::EndOnSphere::ACTIVE_DEFAULT);


  } else if (Prob::Prior::EndOnSphere::ACTIVE_DEFAULT == Active)
    echo_options["prior_end_on_sphere"] = "0 1";


  if (echo_options.count("prior_mid_sphere_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_mid_sphere_scale", str(Prob::Prior::MidpointInSphere::SCALE_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_mid_sphere_sd", str(Prob::Prior::MidpointInSphere::STD_DEV_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_mid_sphere_radius", str(Prob::Prior::MidpointInSphere::radius_default(obs_image))));


    check_parameters_group(echo_options, "prior_mid_sphere", key_values, Prob::Prior::MidpointInSphere::ACTIVE_DEFAULT);


  } else if (Prob::Prior::MidpointInSphere::ACTIVE_DEFAULT == Active)
    echo_options["prior_mid_sphere"] = "0 1";

  if (echo_options.count("prior_hook_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_hook_scale", str(Prob::Prior::Hook::SCALE_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_mid_cube_sd", str(Prob::Prior::Hook::STD_DEV_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_hook_num_points", str(Prob::Prior::Hook::side_length_default(obs_image))));


    check_parameters_group(echo_options, "prior_mid_cube", key_values, Prob::Prior::Hook::ACTIVE_DEFAULT);


  } else if (Prob::Prior::Hook::ACTIVE_DEFAULT == Active)
    echo_options["prior_mid_cube"] = "0 1";


  if (echo_options.count("prior_all_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_all_scale", str(Prob::Prior::AllInCube::SCALE_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_all_sd", str(Prob::Prior::AllInCube::STD_DEV_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_all_num_length_sections", str(Prob::Prior::AllInCube::NUM_LENGTH_SECTIONS_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_all_num_width_sections", str(Prob::Prior::AllInCube::NUM_WIDTH_SECTIONS_DEFAULT)));

    check_parameters_group(echo_options, "prior_all", key_values, Prob::Prior::Hook::ACTIVE_DEFAULT);


  } else if (Prob::Prior::Hook::ACTIVE_DEFAULT == Active)
    echo_options["prior_all"] = "0 1";


  if (echo_options.count("prior_tract_sheer_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_tract_sheer_scale", str(Prob::Prior::Sheer::SCALE_DEFAULT)));

    check_parameters_group(echo_options, "prior_tract_sheer", key_values, Prob::Prior::Sheer::ACTIVE_DEFAULT);


  } else if (Prob::Prior::Sheer::ACTIVE_DEFAULT == Active)
    echo_options["prior_tract_sheer"] = "0 1";


  if (echo_options.count("prior_prior_tract_freq_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_prior_tract_freq_scale", str(Prob::Prior::Frequency::SCALE_DEFAULT)));

    check_parameters_group(echo_options, "prior_prior_tract_mag", key_values, Prob::Prior::Frequency::ACTIVE_DEFAULT);


  } else if (Prob::Prior::Frequency::ACTIVE_DEFAULT == Active)
    echo_options["prior_prior_tract_mag"] = "0 1";


  if (echo_options.count("prior_tract_avg_width_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_tract_avg_width_scale", str(Prob::Prior::BaseWidth::SCALE_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_tract_avg_width_sd", str(Prob::Prior::BaseWidth::STD_DEV_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_tract_avg_width_mean", str(Prob::Prior::BaseWidth::MEAN_DEFAULT)));

    check_parameters_group(echo_options, "prior_tract_basewidth", key_values, Prob::Prior::BaseWidth::ACTIVE_DEFAULT);


  } else if (Prob::Prior::BaseWidth::ACTIVE_DEFAULT == Active)
    echo_options["prior_tract_basewidth"] = "0 1";


  if (echo_options.count("prior_intens_gauss_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_intens_gauss_scale", str(Prob::Prior::GaussDensity::SCALE_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_intens_gauss_mean", str(Prob::Prior::GaussDensity::MEAN_DEFAULT)));

    check_parameters_group(echo_options, "prior_intens_gauss", key_values, Prob::Prior::GaussDensity::ACTIVE_DEFAULT);


  } else if (Prob::Prior::GaussDensity::ACTIVE_DEFAULT == Active)
    echo_options["prior_intens_gauss"] = "0 1";


  if (echo_options.count("prior_b_intens_gauss_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_b_intens_gauss_scale", str(Prob::Prior::BaseIntensity::SCALE_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_b_intens_gauss_sd_frac", str(Prob::Prior::BaseIntensity::STD_DEV_FRACTION_DEFAULT)));

    check_parameters_group(echo_options, "prior_b_intens_gauss", key_values, Prob::Prior::BaseIntensity::ACTIVE_DEFAULT);


  } else if (Prob::Prior::BaseIntensity::ACTIVE_DEFAULT == Active)
    echo_options["prior_b_intens_gauss"] = "0 1";


  if (echo_options.count("prior_intens_flattop_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_b_intens_flattop_scale", str(Prob::Prior::BaseIntensity::FlatTopQuadratic::SCALE_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_b_intens_flattop_sd", str(Prob::Prior::BaseIntensity::FlatTopQuadratic::STD_DEV_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_b_intens_flattop_low", str(Prob::Prior::BaseIntensity::FlatTopQuadratic::LOW_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_b_intens_flattop_high", str(Prob::Prior::BaseIntensity::FlatTopQuadratic::HIGH_DEFAULT)));

    check_parameters_group(echo_options, "prior_intens_flattop", key_values, Prob::Prior::BaseIntensity::FlatTopQuadratic::ACTIVE_DEFAULT);


  } else if (Prob::Prior::BaseIntensity::FlatTopQuadratic::ACTIVE_DEFAULT == Active)
    echo_options["prior_intens_flattop"] = "0 1";


  if (echo_options.count("prior_intens_log_scale")) {

    std::vector< std::pair<std::string, std::string> > key_values;

    key_values.push_back(std::pair<std::string, std::string> ("prior_b_intens_log_scale", str(Prob::Prior::BaseIntensity::Log::SCALE_DEFAULT)));
    key_values.push_back(std::pair<std::string, std::string> ("prior_b_intens_log_in_scale", str(Prob::Prior::BaseIntensity::Log::INPUT_SCALE_DEFAULT)));

    check_parameters_group(echo_options, "prior_intens_log", key_values, Prob::Prior::BaseIntensity::Log::ACTIVE_DEFAULT);


  } else if (Prob::Prior::BaseIntensity::Log::ACTIVE_DEFAULT == Active)
    echo_options["prior_intens_log"] = "0 1";

  //------------------------------//
  //  Method Specific Paramters.  //
  //------------------------------//


  if (is_mcmc(method)) {

    echo_arguments.push_back(obs_image_location);

    if (overide_state_location) {
      echo_arguments.push_back(state_location);
      echo_options.erase("initial_state");
      echo_options.erase("initial_state_location");

    } else if (is_mcmc(loaded_method)) {
      load_and_check_state_parameters(echo_options, "initial_state");
      echo_arguments.push_back(echo_options["initial_state_location"]);
      echo_options.erase("initial_state_location");

    } else
      echo_arguments.push_back("__initial_state_location_place_holder__");

  } else if (is_mcmc(loaded_method)) {

    echo_options.erase("initial_state");
    echo_options.erase("initial_state_location");
    echo_options.erase("sample_period");
    echo_options.erase("burn_sample_period");
    echo_options.erase("burn_snr");

  }


  if (is_analysis(method)) {

    if (obs_image_location != "__obs_image_place_holder__")
      echo_options["obs_image"] = obs_image_location;

  }

  //Hamiltonian method
  if (is_riemannian(method)) {

    check_parameter(echo_options,"num_samples", str(MCMC::Riemannian::NUM_SAMPLES_DEFAULT));
    check_parameter(echo_options,"num_leapfrog_steps", str(MCMC::Hamiltonian::NUM_LEAPFROG_STEPS_DEFAULT));
    check_parameter(echo_options,"num_newton_steps", str(MCMC::Riemannian::NUM_NEWTON_STEPS_DEFAULT));
    check_parameter(echo_options,"burn_num_samples", str(MCMC::Riemannian::BURN_NUM_SAMPLES_DEFAULT));
    check_parameter(echo_options,"burn_num_leapfrog_steps", str(MCMC::Riemannian::BURN_NUM_LEAPFROG_STEPS_DEFAULT));
    check_parameter(echo_options,"burn_num_newton_steps", str(MCMC::Riemannian::BURN_NUM_NEWTON_STEPS_DEFAULT));
    check_parameter(echo_options,"burn_snr", str(MCMC::Riemannian::BURN_SNR_DEFAULT));
    check_parameter(echo_options,"precondition", str(MCMC::Riemannian::PRECONDITION_DEFAULT));

    std::vector<std::string> leapfrog_defaults;

//    leapfrog_defaults.push_back(str(Fibre::Strand::Set::Momentum::STEP_LOCATION_DEFAULT));
//    leapfrog_defaults.push_back(str(Fibre::Tractlet::Set::Momentum::STEP_LOCATION_DEFAULT));
//    leapfrog_defaults.push_back("");
//
//    load_and_check_state_parameters (echo_options, "momen_step", leapfrog_defaults);


  } else if (is_hamiltonian(loaded_method)) {

    echo_options.erase("num_samples");
    echo_options.erase("burn_num_samples");

    echo_options.erase("leapfrog");
    echo_options.erase("leapfrog_location");
    echo_options.erase("leapfrog_scale");

  }


  //Hamiltonian method
  if (is_hamiltonian(method)) {

    check_parameter(echo_options,"num_samples", str(MCMC::Hamiltonian::NUM_SAMPLES_DEFAULT));
    check_parameter(echo_options,"num_leapfrog_steps", str(MCMC::Hamiltonian::NUM_LEAPFROG_STEPS_DEFAULT));
    check_parameter(echo_options,"burn_num_samples", str(MCMC::Hamiltonian::BURN_NUM_SAMPLES_DEFAULT));
    check_parameter(echo_options,"burn_num_leapfrog_steps", str(MCMC::Hamiltonian::BURN_NUM_LEAPFROG_STEPS_DEFAULT));
    check_parameter(echo_options,"burn_snr", str(MCMC::Hamiltonian::BURN_SNR_DEFAULT));

    std::vector<std::string> leapfrog_defaults;

//    leapfrog_defaults.push_back(str(Fibre::Strand::Set::Momentum::STEP_LOCATION_DEFAULT));
//    leapfrog_defaults.push_back(str(Fibre::Tractlet::Set::Momentum::STEP_LOCATION_DEFAULT));
//    leapfrog_defaults.push_back("");
//
//    load_and_check_state_parameters (echo_options, "leapfrog", leapfrog_defaults);


  } else if (is_hamiltonian(loaded_method)) {

    echo_options.erase("num_samples");
    echo_options.erase("burn_num_samples");

    echo_options.erase("leapfrog");
    echo_options.erase("leapfrog_location");
    echo_options.erase("leapfrog_scale");

  }


  //Metropolis method
  if (is_metropolis(method)) {

    check_parameter(echo_options,"num_iterations", str(MCMC::Metropolis::NUM_ITERATIONS_DEFAULT));
    check_parameter(echo_options,"sample_period", str(MCMC::Metropolis::SAMPLE_PERIOD_DEFAULT));
    check_parameter(echo_options,"burn_num_iterations", str(MCMC::Metropolis::BURN_NUM_ITERATIONS_DEFAULT));
    check_parameter(echo_options,"burn_sample_period", str(MCMC::Metropolis::BURN_SAMPLE_PERIOD_DEFAULT));
    check_parameter(echo_options,"burn_snr", str(MCMC::Metropolis::BURN_SNR_DEFAULT));

  } else if (is_metropolis(loaded_method)) {

    echo_options.erase("num_iterations");
    echo_options.erase("burn_num_iterations");

  }



  //Scan method
  if (is_scan(method)) {

    check_parameter(echo_options, "only_first", "0");

    //Append the origin location to the argument list.
    if (origin_location.size())
      echo_arguments.push_back(origin_location);

    else if (is_scan(loaded_method)) {

      load_and_check_state_parameters(echo_options, "origin");
      echo_arguments.push_back(echo_options["origin_location"]);

    } else
      echo_arguments.push_back("__origin_location_place_holder__");


    //Append the axes locations to the argument list.
    if (axis1_location.size()) {

      echo_arguments.push_back(axis1_location);

      if (axis2_location.size()) {

        echo_arguments.push_back(axis2_location);

        if (axis3_location.size())
          echo_arguments.push_back(axis3_location);

      }

    } else {

      if (echo_options.count("axis1")) {
        load_and_check_state_parameters(echo_options, "axis1");
        echo_arguments.push_back(echo_options["axis1_location"]);
      } else
        echo_arguments.push_back("__axes_place_holder__");

      if (echo_options.count("axis2")) {
        load_and_check_state_parameters(echo_options, "axis2");
        echo_arguments.push_back(echo_options["axis2_location"]);
      }

      if (echo_options.count("axis3")) {
        load_and_check_state_parameters(echo_options, "axis3");
        echo_arguments.push_back(echo_options["axis3_location"]);
      }

    }


  } else if (is_scan(loaded_method)) {

    echo_options.erase("only_first");

  }

  echo_options.erase("origin");
  echo_options.erase("origin_location");
  echo_options.erase("axis1");
  echo_options.erase("axis1_scale");
  echo_options.erase("axis1_location");
  echo_options.erase("axis2");
  echo_options.erase("axis2_scale");
  echo_options.erase("axis2_location");
  echo_options.erase("axis3");
  echo_options.erase("axis3_scale");
  echo_options.erase("axis3_location");




  //Test gradient method
  if (is_gradient_test(method)) {

    if (state_location.size())
      echo_arguments.push_back(state_location);
    else if (is_gradient_test(loaded_method))
      echo_arguments.push_back(echo_options["state_location"]);
    else
      echo_arguments.push_back("__state_location_place_holder__");

    echo_options["obs_image"] = obs_image_location;

  } else if (is_gradient_test(loaded_method)) {

    echo_options.erase("state_location");

  }




  //All methods
  echo_arguments.push_back(output_location);



//--------------------------------//
//  Echo the required parameters  //
//--------------------------------//


  //Print preceding newlines.
  std::cout << std::endl << std::endl;

  //print method name.
  std::cout << method;

  //print echo_arguments.
  for (size_t arg_i = 0; arg_i < echo_arguments.size() - set_output_last; arg_i++)
    std::cout << " " << echo_arguments[arg_i];

  //print echo_options.
  for (std::map<std::string, std::string>::iterator opt_it = echo_options.begin(); opt_it != echo_options.end(); ++opt_it) {

    //If option actuall has value (in the checking procedure empty values may be added).
    if (opt_it->second.size())
      std::cout << " -" << opt_it->first << " " << opt_it->second;

  }

  //If output argument is to be printed last, then print it now.
  if (set_output_last)
    std::cout << " " << echo_arguments.back();

  //print following newlines.
  std::cout << std::endl << std::endl;

}



void                              check_parameter(std::map<std::string, std::string>& echo_options,
                                                  const std::string& key,
                                                  const std::string& default_value) {

  if (echo_options.count(key)) {

    if (echo_options[key] == default_value)
      echo_options.erase(key);

  }

}


void                            check_parameters_group (std::map<std::string, std::string>& echo_options,
                                                        const std::string& name,
                                                        std::vector< std::pair<std::string, std::string> > key_default_value_pairs,
                                                        bool default_active) {

  int num_append = 0;

  for (size_t pair_i = 0; pair_i < key_default_value_pairs.size(); pair_i++) {

    if (echo_options.count(key_default_value_pairs[pair_i].first)) {

      if (echo_options[key_default_value_pairs[pair_i].first] != key_default_value_pairs[pair_i].second) {
        key_default_value_pairs[pair_i].second = echo_options[key_default_value_pairs[pair_i].first];
        num_append = pair_i+1;
      }

      echo_options.erase(key_default_value_pairs[pair_i].first);

    }

  }

  if (num_append || !default_active) {

    echo_options[name] = "";

    for (size_t append_i = 0; append_i < key_default_value_pairs.size(); append_i++) //  for (size_t append_i = 0; append_i < num_append; append_i++)
      echo_options[name] += " " + key_default_value_pairs[append_i].second;

  }

}


template <typename T> void        load_and_check_parameters ( std::map<std::string, std::string>& echo_options,
                                                              const std::string& key,
                                                              const std::string& default_value_location = "") {


  T default_value, scaled_default_value, loaded_value, scaled_loaded_value;
  std::string value_str = echo_options[key];
  std::string place_holder = "\"__" + key + "(" + value_str + ")_place_holder__\"";

  if (value_str.size()) {

    double scale = 1.0;
    bool matches_default = false;

    if (echo_options.count(key + "_scale"))
      scale = to<double>(echo_options[key + "_scale"]);

    std::string loaded_value_location = echo_options[key + "_location"];

    echo_options.erase(key);
    echo_options.erase(key+"_scale");
    echo_options.erase(key+"_location");


    //First, check to see if saved value matches default value or a scaled version of the defaul value
    if (default_value_location.size()) {

      try {

        default_value.load (default_value_location);
        scaled_default_value = default_value;
        scaled_default_value *= scale;

        //Loaded value still matches saved value.
        if (default_value.matlab_str() == value_str)

          //'loaded value' == 'Default_value' no paramters need to be set.
          matches_default = true;

        else if (scaled_default_value.matlab_str() == value_str) {

          //'loaded value' == scaled version of 'Default_value' only scale parameter needs to be set.
          echo_options[key+"_scale"] = str(scale);
          matches_default = true;

        }

      } catch (Exception e) {

        std::cout << "WARNING!! Could not load default parameter from location '" + default_value_location + "'." << std::endl << std::endl;

      }

    }

    //If saved value does not equal default, check to see if it still matches the saved location
    if (!matches_default) {

      if (!loaded_value_location.size()) {
        echo_options[key] = place_holder;

      } else {

        //If it doesn't check the loaded value to see if that is still the same.
        try {

          loaded_value.load(loaded_value_location);
          scaled_loaded_value = loaded_value;
          scaled_loaded_value *= scale;

          //Loaded value still matches saved value.
          if (loaded_value.matlab_str() == value_str)

            echo_options[key+"_location"] = loaded_value_location;

          else if (scaled_loaded_value.matlab_str() == value_str) {

            //Scale is added back into the echo_options
            echo_options[key+"_location"] = loaded_value_location;
            echo_options[key+"_scale"] = str(scale);


          } else

            echo_options[key] = place_holder;


        } catch (Exception e) {

          std::cout << "WARNING!! Could not load parameter from location '" + loaded_value_location + "', skipping." << std::endl << std::endl;

          echo_options[key] = place_holder;

        }

      }

    }

  }

}


void                              load_and_check_state_parameters ( std::map<std::string, std::string>& echo_options,
                                                                    const std::string& key,
                                                                    const std::vector<std::string>& default_value_locations) {


  if (echo_options.count(key)) {

    std::string state_location = echo_options[key + "_location"];

    if (File::has_extension<Fibre::Strand>(state_location))
      load_and_check_parameters<Fibre::Strand::Set>(echo_options, key, default_value_locations[0]);

    else if (File::has_extension<Fibre::Tractlet>(state_location))
      load_and_check_parameters<Fibre::Tractlet::Set>(echo_options, key, default_value_locations[1]);

    else if (File::has_extension<MCMC::State>(state_location))
      load_and_check_parameters<MCMC::State>(echo_options, key, default_value_locations[2]);

    else {
      echo_options[key + "_location"] =  "__" + key + "_place_holder_(" + echo_options[key] + ")__";
      echo_options.erase(key);
      echo_options.erase(key+"_scale");
    }

  }

}
*/

int main() {

  return 0;

}
