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




#include <ctime>

#include "bts/cmd.h"

#include "math/matrix.h"

#include "progressbar.h"


#include "bts/common.h"

#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/strand/set.h"


#include "bts/diffusion/model.h"
#include "bts/image/expected/trilinear/buffer.h"
#include "bts/image/expected/gaussian/buffer.h"
#include "bts/image/expected/quartic/buffer.h"

#include "bts/image/observed/buffer.h"

#include "bts/prob/uniform.h"
#include "bts/prob/prior.h"





#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood/one_sided_gaussian.h"
#include "bts/prob/likelihood/gaussian.h"
#include "bts/analysis/scan.h"

#include "bts/file.h"

#include "bts/math/common.h"


#include "bts/inline_functions.h"


extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

using namespace BTS;

double NUM_STEPS_DEFAULT = 30;

double AXIS_SCALE_DEFAULT = 1.0;


SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Runs Metropolis-Hastings sampling on a given image from a given starting configuration of tractlets or strands.",
  "",
  NULL
};

ARGUMENTS = {

  Argument ("start_fibres", "The starting point of the scan.").type_file (),

  Argument ("end_fibres", "The ending point of the scan.").optional().type_file (),

  Argument ("second_end_fibres", "The ending point of the scan.").optional().type_file (),

  Argument ("third_end_fibres", "The ending point of the scan.").optional().type_file (),

  Argument ("output", "The output file containing the function scan. Will be either a text file if it is a linear scan or an image file if it is multi-dimension scan.").optional().type_file (),

  Argument()
};


OPTIONS = {


  Option ("num_steps", "Number of taken along axis.")
   + Argument ("num_steps", "").type_integer (1, NUM_STEPS_DEFAULT, LARGE_INT),


  Option ("num_steps2", "Number of taken along second axis.")
   + Argument ("num_steps", "").type_integer (1, NUM_STEPS_DEFAULT, LARGE_INT),
  

  Option ("num_steps3", "Number of taken along third axis.")
   + Argument ("num_steps", "").type_integer (1, NUM_STEPS_DEFAULT, LARGE_INT),
    

  Option ("centred", "For 1D case (multi-D cases are automatically centred) first set is interpreted as origin of image second set is interpreted as a perturbation axes."),


  Option ("axis1_scale", "Scale applied to axis1, used to quickly scale binary masks.")
     + Argument ("axis1_scale", "").type_float (SMALL_FLOAT, AXIS_SCALE_DEFAULT, LARGE_FLOAT),


  Option ("axis2_scale", "Scale applied to axis2, used to quickly scale binary masks.")
   + Argument ("axis2_scale", "").type_float (SMALL_FLOAT, AXIS_SCALE_DEFAULT, LARGE_FLOAT),


  Option ("axis3_scale", "Scale applied to axis3, used to quickly scale binary masks.")
   + Argument ("axis3_scale", "").type_float (SMALL_FLOAT, AXIS_SCALE_DEFAULT, LARGE_FLOAT),


  Option ("only_prob", "Only use the following probability (sets the rest to Prob::Uniform).")
   + Argument("only_prob", "").type_text(),


  Option ("no_prob", "Don't caculate probability, used to create a range of states."),

  Option ("save_gradient", "Save the gradient as well as the log probability for single axis scan."),

  Option ("only_first", "Use the first fibre in the input files instead of the whole set."),

  Option ("no_map", "Don't map the set to its smallest permutation before the scale."),

  DIFFUSION_PARAMETERS,

  IMAGE_PARAMETERS,

  EXPECTED_IMAGE_PARAMETERS,

  LIKELIHOOD_PARAMETERS,

  PRIOR_PARAMETERS,

  COMMON_PARAMETERS,

//  Option ("num_ref_sections", "Number of sections used to generate the reference images")
//   + Argument ("num_ref_sections", "").type_integer (1, LARGE_INT, 1000),

  Option ("obs_image", "The location of the reference image that is to be set")
   + Argument ("obs_image", "").type_image_in(),

Option() };






EXECUTE {

//-----------------//
//  Load Arguments //
//-----------------//

  std::string origin_location  = argument[0];

  std::string output_location, output_gradient_location, output_hessian_location;
  std::string axis1_location, axis2_location, axis3_location;


  int num_axes;

  if (argument.size() == 1)
    num_axes = 0;

  else if (argument.size() == 2) {
    num_axes = -1;
    output_location = argument[1].c_str();

    if (File::extension(origin_location) != File::extension(output_location))
      throw Exception ("Output location extension must match origin location extension (and must be a 'set' extension) when only two arguments are provided.");

  } else
    num_axes = argument.size() - 2;

  if (num_axes > 3)
    throw Exception ("Max number of axes is 3, " + str(num_axes) + " supplied.");

  if (num_axes >= 1) {

    axis1_location = argument[1].c_str();

    output_location   = argument[num_axes+1].c_str();
    output_gradient_location = File::strip_extension(output_location) + ".gradient." + File::extension(output_location);

    if (File::extension(origin_location) != File::extension(axis1_location))
      throw Exception ("Extension of axis1 location (2nd argument), '" + File::extension(axis1_location) + "', does not match extension of origin location (1st argument), '" + File::extension(origin_location) + "'.");
  }

  if (num_axes >= 2) {
    axis2_location = argument[2].c_str();

    if (File::extension(origin_location) != File::extension(axis2_location))
      throw Exception ("Extension of second axis location (3rd argument), '" + File::extension(axis2_location) + "', does not match extension of origin location (1st argument), '" + File::extension(origin_location) + "'.");
  }

  if (num_axes == 3) {
    axis3_location = argument[3].c_str();

    if (File::extension(origin_location) != File::extension(axis3_location))
      throw Exception ("Extension of third axis location (4th argument), '" + File::extension(axis3_location) + "', does not match extension of origin location (1st argument), '" + File::extension(origin_location) + "'.");

  }

  if ( (num_axes > 1) && !File::has_extension<Image::Expected::Buffer>(output_location))
    throw Exception ("Output file must be a mrtrix image ('mif' extension) for more scans with multiple axes.");


//----------------------------------//
//  Get and Set Optional Parameters //
//----------------------------------//


  size_t num_steps    = NUM_STEPS_DEFAULT;
  size_t num_steps2;
  size_t num_steps3;
  bool centred = false;
  double axis1_scale = AXIS_SCALE_DEFAULT;
  double axis2_scale;
  double axis3_scale;
  std::string only_prob = "";
  bool no_prob = false;
  bool save_gradient = false;
  bool only_first = false;
  std::string prior_b_intens_gauss_mean_type;
  std::string obs_image_name;
  bool no_map = false;




  Options opt = get_options("num_steps");
  if (opt.size())
    num_steps = opt[0][0];
    
  opt = get_options("num_steps2");
  if (opt.size())
    num_steps2 = opt[0][0];
  else
    num_steps2 = num_steps;
    
  opt = get_options("num_steps3");
  if (opt.size())
    num_steps3 = opt[0][0];
  else
    num_steps3 = num_steps;    

  opt = get_options("centred");
  if (opt.size())
    centred = true;

  opt = get_options("axis1_scale");
  if (opt.size())
    axis1_scale = opt[0][0];

  opt = get_options("axis2_scale");
  if (opt.size())
    axis2_scale = opt[0][0];
  else
    axis2_scale = axis1_scale;

  opt = get_options("axis3_scale");
  if (opt.size())
    axis3_scale = opt[0][0];
  else
    axis3_scale = axis1_scale;

  opt = get_options("only_prob");
  if (opt.size())
    only_prob = opt[0][0].c_str();

  opt = get_options("no_prob");
  if (opt.size())
    no_prob = true;

  opt = get_options("save_gradient");
  if (opt.size())
    save_gradient = true;

  opt = get_options("only_first");
  if (opt.size())
    only_first = true;

  opt = get_options("no_map");
  if (opt.size())
    no_map = true;

  // Loads parameters to construct Diffusion::Model ('diff_' prefix)
  SET_DIFFUSION_PARAMETERS;

  // Loads parameters to construct Image::Expected::*::Buffer that are inherited from Image::Observed::Buffer ('img_' prefix)
  SET_IMAGE_PARAMETERS;

  // Loads parameters to construct Image::Expected::*::Buffer ('img_' prefix)
  SET_EXPECTED_IMAGE_PARAMETERS;

  // Loads parameters to construct Prob::Likelihood ('like_' prefix)
  SET_LIKELIHOOD_PARAMETERS;

  // Loads parameters to construct Prob::Prior ('prior_' prefix)
  SET_PRIOR_PARAMETERS;

  // Loads parameters that are common to all commands.
  SET_COMMON_PARAMETERS;


//  size_t num_ref_sections = exp_num_length_sections;

//  opt = get_options("num_ref_sections");
//  if (opt.size())
//    num_ref_sections = opt[0][0];



  if (num_axes > 1) {

    if (centred)
      std::cout << "Warning: option '-centred' is not required with used with multi-dimension scans as non-centred scans are not available.";

    if (no_prob)
      throw Exception ("The '-no_prob' option is meant to be used to generate a sequence of states can can't be used with multi-dimension scans.");

  }


//----------------------//
//  Load observed image //
//----------------------//

  Image::Observed::Buffer obs_image;

  opt = get_options("obs_image");
  if (opt.size()) {
    obs_image_name = opt[0][0].c_str();
    obs_image.load(obs_image_name);

    if (obs_image.properties().count("diff_response_SH")) {

      if (Math::matlab_str(diff_response_SH) != obs_image.properties()["diff_response_SH"])
        std::cout << std::endl << "Warning! Diffusion response function harmonics (" << Math::matlab_str(diff_response_SH)  << ") do not match reference image (" << obs_image.properties()["diff_response_SH"] + ")!" << std::endl;

    }

  }


 //----------------------------//
 //  Initialize Expected Image //
 //----------------------------//


  Diffusion::Model diffusion_model = Diffusion::Model::factory (diff_encodings,
                                                                diff_response_SH,
                                                                diff_adc,
                                                                diff_fa,
                                                                diff_isotropic);

  Image::Expected::Buffer* exp_image = Image::Expected::Buffer::factory(exp_type,
                                                         img_dims,
                                                         img_vox_lengths,
                                                         diffusion_model,
                                                         exp_num_length_sections,
                                                         exp_num_width_sections,
                                                         exp_interp_extent,
                                                         img_offsets,
                                                         exp_enforce_bounds,
                                                         exp_half_width);


//----------------------------------------------------------------//
// Auto-generate base intensity initial value/b_intens_gauss_mean //
//----------------------------------------------------------------//

  if ((exp_base_intensity < 0) && !no_prob && (!only_prob.size() || only_prob.substr(0,16) == "Prob::Likelihood"))
     exp_base_intensity = exp_image->base_intensity_default(obs_image, origin_location);


//-----------------------//
// Initialize Likelihood //
//-----------------------//

  if (no_prob || (only_prob.size() && only_prob.substr(0,16) != "Prob::Likelihood"))
    like_snr = 0.0;

  Prob::Likelihood* likelihood = Prob::Likelihood::factory(like_type, obs_image, exp_image, like_snr, like_b0_include, like_outside_scale, like_ref_b0, like_ref_signal);



//------------------//
// Initialize Prior //
//------------------//

  if (no_prob || only_prob == "Prob::Likelihood")
    prior_scale = 0.0;

  Prob::Prior prior (prior_scale,
                    prior_freq_scale,
                    prior_freq_aux_scale,
                    prior_hook_scale,
                    prior_hook_num_points,
                    prior_hook_num_width_sections,
                    prior_density_high_scale,
                    prior_density_low_scale,
                    prior_density_num_points,
                    prior_acs_scale,
                    prior_acs_mean,
                    prior_length_scale,
                    prior_length_mean,
                    prior_thinness_scale,
                    prior_thinness_power);

//----------------//
// Set Properties //
//----------------//
  
  std::map<std::string, std::string> properties;

  
  properties["method"]              = "scan";
  properties["obs_image"]           = obs_image_name;
  if (num_axes != -1)
    properties["origin"]              = Fibre::Base::Object::load_matlab_str(origin_location);
  properties["origin_location"]     = origin_location;

  if (num_axes >= 1) {
    properties["axis1"]               = Fibre::Base::Object::load_matlab_str(axis1_location, axis1_scale);
    properties["axis1_scale"]         = str(axis1_scale);
    properties["axis1_location"]      = axis1_location;
    properties["num_steps"]           = str(num_steps);
  }

  if (num_axes >= 2) {
    properties["axis2"]               = Fibre::Base::Object::load_matlab_str(axis2_location, axis2_scale);
    properties["axis2_scale"]         = str(axis2_scale);
    properties["axis2_location"]      = axis2_location;
    properties["num_steps"]           = str(num_steps);
  }

  if (num_axes == 3) {
    properties["axis3"]               = Fibre::Base::Object::load_matlab_str(axis3_location, axis3_scale);
    properties["axis3_scale"]         = str(axis3_scale);
    properties["axis3_location"]      = axis3_location;
    properties["num_steps"]           = str(num_steps);
  }

  properties["centred"]             = str(centred);
  properties["only_first"]          = str(only_first);


  if (no_prob)
    properties["no_prob"]             = str(no_prob);

  else {

    properties["only_prob"]           = only_prob;

    ADD_DIFFUSION_PROPERTIES(properties);

    if (!obs_image.dim(X)) {
      ADD_IMAGE_PROPERTIES(properties);
    }

    ADD_EXPECTED_IMAGE_PROPERTIES(properties);

    ADD_LIKELIHOOD_PROPERTIES(properties);

    ADD_PRIOR_PROPERTIES(properties);

    ADD_COMMON_PROPERTIES(properties);

  }



//----------------------------//
// Actually Perform the scans //
//----------------------------//


  if (only_first) {

    if (File::has_extension<Fibre::Strand>(origin_location)) {

      Fibre::Strand::Set origin_set(origin_location);
      Fibre::Strand origin(origin_set[0]);

      if (num_axes == 0) {

        Fibre::Strand prior_gradient;
        prior_gradient = origin;

        double total_prior  = prior.log_prob(origin, prior_gradient);

        std::cout << "Total Prior: " << total_prior << std::endl;
        std::cout << "Total Prior Gradient:" << prior_gradient << std::endl << std::endl;

        std::map<std::string, double> component_values = prior.get_component_values(origin);

        for (std::map<std::string, double>::const_iterator comp_it = component_values.begin(); comp_it != component_values.end(); ++comp_it)
          std::cout << "Prior " << comp_it->first << ": " << comp_it->second << std::endl;

        std::cout << std::endl;



      } else if (num_axes == 1) {

        Fibre::Strand::Set axis1_set(axis1_location);
        Fibre::Strand axis1(axis1_set[0]);

        //Convert the 'start' and 'finish' sets into an origin and perturbation around this origin.
        if (!centred) {

          Fibre::Strand::Set start_set(origin_location), end_set(axis1_location);
          Fibre::Strand start(start_set[0]), end(end_set[0]);

          if (!no_map)
            end = end.smallest_distance_set(start);

          axis1 = (end - start) / 2.0;
          origin = start + axis1;
        }

        axis1 *= axis1_scale;

        Analysis::scan<Fibre::Strand>(*likelihood, prior, origin, axis1, num_steps, output_location, properties, save_gradient);


      } else if (num_axes == 2) {

        Fibre::Strand::Set axis1_set(axis1_location), axis2_set(axis2_location);
        Fibre::Strand axis1(axis1_set[0]), axis2(axis2_set[0]);

        axis1 *= axis1_scale;
        axis2 *= axis2_scale;


        Analysis::scan<Fibre::Strand>(*likelihood, prior, origin, axis1, axis2, num_steps, num_steps2, output_location, properties);


      } else if (num_axes == 3) {


        Fibre::Strand::Set axis1_set(axis1_location), axis2_set(axis2_location), axis3_set(axis3_location);
        Fibre::Strand axis1(axis1_set[0]),axis2(axis2_set[0]), axis3(axis3_set[0]);

        axis1 *= axis1_scale;
        axis2 *= axis2_scale;
        axis3 *= axis3_scale;

        Analysis::scan<Fibre::Strand>(*likelihood, prior, origin, axis1, axis2, axis3, num_steps, num_steps2, num_steps3, output_location, properties);



      } else if (num_axes == -1)
        throw Exception ("Can't do pre-generated set (implied by 2 input arguments) for individual set.");


    } else if (File::has_extension<Fibre::Tractlet>(origin_location)) {

      Fibre::Tractlet::Set origin_set(origin_location);

      Fibre::Tractlet origin(origin_set[0]);

      if (num_axes == 0) {


        Fibre::Tractlet prior_gradient(origin);

        double total_prior  = prior.log_prob(origin, prior_gradient);


        std::cout << "Total Prior: " << total_prior << std::endl;
        std::cout << "Total Prior Gradient:" << prior_gradient << std::endl << std::endl;

        std::map<std::string, double> component_values = prior.get_component_values(origin);

        for (std::map<std::string, double>::const_iterator comp_it = component_values.begin(); comp_it != component_values.end(); ++comp_it)
          std::cout << "Prior " << comp_it->first << ": " << comp_it->second << std::endl;

        std::cout << std::endl;



      } else if (num_axes == 1) {

        Fibre::Tractlet::Set axis1_set(axis1_location);
        Fibre::Tractlet axis1(axis1_set[0]);

        //Convert the 'start' and 'finish' sets into an origin and perturbation around this origin.
        if (!centred) {

          Fibre::Tractlet::Set start_set(origin_location), end_set(axis1_location);
          Fibre::Tractlet start(start_set[0]), end(end_set[0]);

          if (!no_map)
            end = end.smallest_distance_set(start);

          axis1 = (end - start) / 2.0;
          origin = start + axis1;
        }

        axis1 *= axis1_scale;

        Analysis::scan<Fibre::Tractlet>(*likelihood, prior, origin, axis1, num_steps, output_location, properties, save_gradient);


      } else if (num_axes == 2) {

        Fibre::Tractlet::Set axis1_set(axis1_location), axis2_set(axis2_location);
        Fibre::Tractlet axis1(axis1_set[0]), axis2(axis2_set[0]);

        axis1 *= axis1_scale;
        axis2 *= axis2_scale;

        Analysis::scan<Fibre::Tractlet>(*likelihood, prior, origin, axis1, axis2, num_steps, num_steps2, output_location, properties);


      } else if (num_axes == 3) {


        Fibre::Tractlet::Set axis1_set(axis1_location), axis2_set(axis2_location), axis3_set(axis3_location);
        Fibre::Tractlet axis1(axis1_set[0]), axis2(axis2_set[0]), axis3(axis3_set[0]);

        axis1 *= axis1_scale;
        axis2 *= axis2_scale;
        axis3 *= axis3_scale;

        Analysis::scan<Fibre::Tractlet>(*likelihood, prior, origin, axis1, axis2, axis3, num_steps, num_steps2, num_steps3, output_location, properties);



      } else if (num_axes == -1)
        throw Exception ("Can't do pre-generated set (implied by 2 input arguments) for individual set.");

    } else
      throw Exception ("Unrecognised extension '" + File::extension(origin_location));
    


//********************************************
// End only first
//*********************************************



  } else {

    if (File::has_extension<Fibre::Strand>(origin_location)) {

      Fibre::Strand::Set origin(origin_location);

      if (num_axes == 0) {


        Fibre::Strand::Set prior_gradient(origin), likelihood_gradient(origin), gradient(origin);
        Fibre::Strand::Set::Tensor likelihood_fisher(origin); //, prior_hessian(origin);
        std::vector< MR::Math::Matrix<double> > likelihood_fisher_gradient;

        double total_prior;
        double likelihood_px;

        if (save_gradient) {
          total_prior  = prior.log_prob(origin, prior_gradient);//, prior_hessian);
          likelihood_px   = likelihood->log_prob_and_fisher(origin, likelihood_gradient, likelihood_fisher);

          gradient = prior_gradient;
          gradient += likelihood_gradient;

          gradient.save(File::strip_extension(origin_location) + ".gradient." + File::extension(origin_location));
          prior_gradient.save(File::strip_extension(origin_location) + ".prior_gradient." + File::extension(origin_location));
          likelihood_gradient.save(File::strip_extension(origin_location) + ".likelihood_gradient." + File::extension(origin_location));
          likelihood_fisher.save(File::strip_extension(origin_location) + ".likelihood_hessian." + File::extension(origin_location) + ".tnr");

        } else {
          total_prior  = prior.log_prob(origin);//, prior_hessian);
          likelihood_px   = likelihood->log_prob(origin);
        }

        std::cout << "Total Prior: " << total_prior << std::endl;
        std::cout << "Likelihood: " << likelihood_px << std::endl;

        std::map<std::string, double> component_values = prior.get_component_values(origin);

        for (std::map<std::string, double>::const_iterator comp_it = component_values.begin(); comp_it != component_values.end(); ++comp_it)
          std::cout << "Prior " << comp_it->first << ": " << comp_it->second << std::endl;

        std::cout << std::endl;


        std::cout << "Total Prior Gradient:" << prior_gradient << std::endl << std::endl;

      } else if (num_axes == 1) {


        Fibre::Strand::Set axis1(axis1_location);

        //Convert the 'start' and 'finish' sets into an origin and perturbation around this origin.
        if (!centred) {

          Fibre::Strand::Set start(origin_location), end(axis1_location);

          if (!no_map)
            end = end.smallest_distance_set(start);

          axis1 = (end - start) / 2.0;

          origin = start + axis1;
        }

        axis1 *= axis1_scale;

        Analysis::scan<Fibre::Strand::Set>(*likelihood, prior, origin, axis1, num_steps, output_location, properties, save_gradient);


      } else if (num_axes == 2) {

        Fibre::Strand::Set axis1(axis1_location), axis2(axis2_location);

        axis1 *= axis1_scale;
        axis2 *= axis2_scale;

        Analysis::scan<Fibre::Strand::Set>(*likelihood, prior, origin, axis1, axis2, num_steps, num_steps2, output_location, properties);

      } else if (num_axes == 3) {

        Fibre::Strand::Set axis1(axis1_location), axis2(axis2_location), axis3(axis3_location);

        axis1 *= axis1_scale;
        axis2 *= axis2_scale;
        axis3 *= axis3_scale;

        Analysis::scan<Fibre::Strand::Set>(*likelihood, prior, origin, axis1, axis2, axis3, num_steps, num_steps2, num_steps3, output_location, properties);


      } else if (num_axes == -1)
        throw Exception ("Can't do pre-generated set (implied by 2 input arguments) for individual set.");


    } else if (File::has_extension<Fibre::Tractlet>(origin_location)) {

      Fibre::Tractlet::Set origin(origin_location);

      if (num_axes == 0) {


        Fibre::Tractlet::Set prior_gradient(origin), likelihood_gradient(origin), gradient(origin);
        Fibre::Tractlet::Set::Tensor likelihood_fisher(origin); //prior_hessian(origin);

        double total_prior;
        double likelihood_px;

        if (save_gradient) {
          total_prior  = prior.log_prob(origin, prior_gradient);//, prior_hessian);
          likelihood_px   = likelihood->log_prob_and_fisher(origin, likelihood_gradient, likelihood_fisher);

          gradient = prior_gradient;
          gradient += likelihood_gradient;

          gradient.save(File::strip_extension(origin_location) + ".gradient." + File::extension(origin_location));
          prior_gradient.save(File::strip_extension(origin_location) + ".prior_gradient." + File::extension(origin_location));
          likelihood_gradient.save(File::strip_extension(origin_location) + ".likelihood_gradient." + File::extension(origin_location));
          likelihood_fisher.save(File::strip_extension(origin_location) + ".likelihood_hessian." + File::extension(origin_location) + ".tnr");

        } else {
          total_prior  = prior.log_prob(origin);//, prior_hessian);
          likelihood_px   = likelihood->log_prob(origin);
        }

        std::cout << "Total Prior: " << total_prior << std::endl;
        std::cout << "Likelihood: " << likelihood_px << std::endl;

        std::map<std::string, double> component_values = prior.get_component_values(origin);

        for (std::map<std::string, double>::const_iterator comp_it = component_values.begin(); comp_it != component_values.end(); ++comp_it)
          std::cout << "Prior " << comp_it->first << ": " << comp_it->second << std::endl;

        std::cout << std::endl;

      } else if (num_axes == 1) {

        Fibre::Tractlet::Set axis1(axis1_location);

        //Convert the 'start' and 'finish' sets into an origin and perturbation around this origin.
        if (!centred) {

          Fibre::Tractlet::Set start(origin_location), end(axis1_location);

          if (!no_map)
            end = end.smallest_distance_set(start);

          axis1 = (end - start) / 2.0;
          origin = start + axis1;
        }

        axis1 *= axis1_scale;

        Analysis::scan<Fibre::Tractlet::Set>(*likelihood, prior, origin, axis1, num_steps, output_location, properties, save_gradient);

      
      } else if (num_axes == 2) {
      
        Fibre::Tractlet::Set axis1(axis1_location), axis2(axis2_location);

        axis1 *= axis1_scale;
        axis2 *= axis2_scale;

        Analysis::scan<Fibre::Tractlet::Set>(*likelihood, prior, origin, axis1, axis2, num_steps, num_steps2, output_location, properties);

      } else if (num_axes == 3) {
      
        Fibre::Tractlet::Set axis1(axis1_location), axis2(axis2_location), axis3(axis3_location);

        axis1 *= axis1_scale;
        axis2 *= axis2_scale;
        axis3 *= axis3_scale;

        Analysis::scan<Fibre::Tractlet::Set>(*likelihood, prior, origin, axis1, axis2, axis3, num_steps, num_steps2, num_steps3, output_location, properties);

      } else if (num_axes == -1)
        throw Exception ("Can't do pre-generated set (implied by 2 input arguments) for individual set.");

    } else if (File::has_extension<Fibre::Strand::Set>(origin_location)) {

      //Special case where the origin is actually a precalculated set of fibres
      if (num_axes != -1)
        throw Exception ("Set origin only allowed when exactly two arguments are provided.");


//      Fibre::Strand::Set true_strands ("/home/tclose/data/fibre/strand/single/z.str");

      std::vector<std::string> writer_header;
      writer_header.push_back("px");

      Fibre::Strand::Set::Reader reader(origin_location);
      Fibre::Strand::Set::Writer writer(output_location, reader, writer_header);

      Fibre::Strand::Set strands;

      MR::ProgressBar progress_bar ("Scanning through strand sets...", to<size_t>(reader.get_extend_props()["count"]));

      while (reader.next(strands)) {

        strands.set_extend_prop("px", str(likelihood->log_prob(strands) + prior.log_prob(strands)));
        writer.append(strands);
        ++progress_bar;
      }

    } else if (File::has_extension<Fibre::Tractlet::Set>(origin_location)) {

      //Special case where the origin is actually a precalculated set of fibres
      if (num_axes != -1)
        throw Exception ("Set origin only allowed when exactly two arguments are provided.");

      std::vector<std::string> writer_header;
      writer_header.push_back("px");

      Fibre::Tractlet::Set::Reader reader(origin_location);
      Fibre::Tractlet::Set::Writer writer(output_location, reader, writer_header);

      Fibre::Tractlet::Set tractlets;

      MR::ProgressBar progress_bar ("Scanning through tractlet sets...", to<size_t>(reader.get_extend_props()["count"]));

      while (reader.next(tractlets)) {

        tractlets.set_extend_prop("px", str(likelihood->log_prob(tractlets) + prior.log_prob(tractlets)));
        writer.append(tractlets);

        ++progress_bar;
      }


    } else
      throw Exception ("Unrecognised extension '" + File::extension(origin_location));


  }

  delete likelihood;
  delete exp_image;


}

