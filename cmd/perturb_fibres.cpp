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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "bts/cmd.h"
#include "progressbar.h"

#include "bts/common.h"


#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"
#include "bts/mcmc/proposal/distribution/gaussian.h"




#include "bts/inline_functions.h"


using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Perturbs fibres, useful for debugging.",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("input", "The paths to be perturbed.").type_file (),
  Argument ("output", "The output perturbed paths.").optional().type_file(),
  Argument()
};


OPTIONS = {

  Option ("stdev", "The base standard deviation the parameters will be perturbed by.")
    + Argument ("stdev", "").type_float (SMALL_FLOAT, 0.05, LARGE_FLOAT),

  Option ("scales_location", "The relative scales of the perturbations.")
   + Argument ("scales_location", "").type_file (),

  Option ("seed", "The random seed that is passed to the random generator")
   + Argument ("seed", ""),

  Option("flat_position", "Use a flat distribution for a position parameters."),

Option() };


template <typename T> void      perturb_fibres(const std::string& input_location,
                                                const std::string& output_location,
                                                double stdev,
                                                const std::string& scale_location,
                                                gsl_rng* rand_gen,
                                                const MR::Math::Vector<double>& flat_mask);

EXECUTE {


  std::string input_location = argument[0];
  std::string output_location;
  
  if (argument.size() > 1)
    output_location = argument[1].c_str();
  else
    output_location = input_location;

  double stdev = 0.05;
  std::string scales_location;
  size_t seed = time(NULL);
  bool flat_position = false;

  Options opt = get_options("stdev");
  if (opt.size())
    stdev = opt[0][0];

  opt = get_options("scales_location");
  if (opt.size())
    scales_location = opt[0][0].c_str();

  opt = get_options("seed");
  if (opt.size()) {
    std::string seed_string = opt[0][0];
    seed = to<size_t>(seed_string);
  } else
   std::cout << "No random seed supplied. Using timestamp: " << seed << std::endl;


  opt = get_options("flat_position");
  if (opt.size())
    flat_position = true;


//-------------------------------//
// Setup Random number Generator //
//-------------------------------//


  gsl_rng *rand_gen;
  rand_gen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rand_gen, seed);


//----------------------------------------------//
//  Loop through paths and perturb parameters  //
//----------------------------------------------//


	if (File::has_extension<Fibre::Strand>(input_location)) {

	  Fibre::Strand::Set mask (input_location);

	  mask.zero();

	  if (flat_position)
	    for (size_t strand_i = 0; strand_i < mask.size(); ++strand_i)
	      mask[strand_i][0].set(1,1,1);

	  MR::Math::Vector<double>& vector_mask = mask;

	  perturb_fibres<Fibre::Strand>(input_location, output_location, stdev, scales_location, rand_gen, vector_mask);

	} else if (File::has_extension<Fibre::Tractlet>(input_location)) {

    Fibre::Tractlet::Set mask (input_location);

    mask.zero();

    if (flat_position)
      for (size_t tract_i = 0; tract_i < mask.size(); ++tract_i)
        mask[tract_i][0][0].set(1,1,1);

    MR::Math::Vector<double>& vector_mask = mask;

	  perturb_fibres<Fibre::Tractlet>(input_location, output_location, stdev, scales_location, rand_gen, vector_mask);

	} else
	  throw Exception ("Unrecognised extension '" + File::extension(input_location) + "'.");

}


template <typename T> void      perturb_fibres(const std::string& input_location,
                                                const std::string& output_location,
                                                double stdev,
                                                const std::string& scales_location,
                                                gsl_rng* rand_gen,
                                                const MR::Math::Vector<double>& flat_mask) {

  typename T::Set fibres(input_location);
  typename T::Set scales;

  if (scales_location.size()) {

    scales.load(scales_location);

    if (scales.size() == 1) {
      T tmpl_scale = scales[0];
      scales = fibres;

      for (size_t fibre_i = 0; fibre_i < fibres.size(); ++fibre_i)
        scales[fibre_i] = tmpl_scale;

    } else if (scales.size() != fibres.size())
      throw Exception ("scales' size do not match fibres.");

  } else {

    scales = fibres;
    scales.set(1.0);

  }

  scales *= stdev;

  MR::Math::Vector<double>& fibres_vector = fibres;
  MR::Math::Vector<double>& scales_vector = scales;

  if (fibres_vector.size() != scales_vector.size())
    throw Exception ("Scales vector does not match fibre vector.");

  MR::ProgressBar progress_bar  ("Perturbing strands ...", fibres_vector.size());

  for (size_t elem_i = 0; elem_i < fibres_vector.size(); elem_i++) {
    if (flat_mask[elem_i])
      fibres_vector[elem_i] += (double)gsl_ran_flat(rand_gen, -scales_vector[elem_i], scales_vector[elem_i]);
    else
      fibres_vector[elem_i] += (double)gsl_ran_gaussian(rand_gen, scales_vector[elem_i]);
    progress_bar++;
  }


  fibres.save(output_location);
  

}



