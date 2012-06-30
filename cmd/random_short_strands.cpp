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

#include "bts/fibre/strand/set.h"




#include "bts/inline_functions.h"


using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Trims paths to a sphere of specified radius.",
  "The sections of paths that lie outside the radius of the sphere are removed and sections that reenter are split into new paths.",
  NULL
};

ARGUMENTS = {
  Argument ("output", "The returned trimmed paths.").type_file(),
  Argument()
};


OPTIONS = {

  Option ("degree", "Limit the degree of the Strand coefficients used to describe the strands to the first - degrees")
   + Argument ("degree", "Limit the degree of the Strand coefficients used to describe the strands to the first - degrees").type_integer (1, 200, 3),

  Option ("sphere_r", "The radius of the sphere the strands will be trimmed to")
   + Argument ("sphere_r", "").type_float (SMALL_FLOAT, 1.0, LARGE_FLOAT),

  Option ("num_strands", "The number of strands to be created")
   + Argument ("num_strands", "").type_integer(1, 10, LARGE_INT),

  Option ("strand_length", "The length of the short strands")
   + Argument ("strand_length", "").type_float (SMALL_FLOAT, 0.01, LARGE_FLOAT),

  Option ("num_points", "The number of points that will be generated along the strand location")
   + Argument ("points", "The number of points that will be generated along the strand location").type_integer (1, 100, 2000),

Option() };


EXECUTE {


  std::string output_location = argument[0];;



  size_t degree = 3;
  double sphere_radius = 1.0;
  size_t num_strands = 0.05;
  double strand_length = 0.01;
  size_t num_points = 0;
  

  Options opt = get_options(0);
  if (opt.size())
    degree = opt[0][0];

  opt = get_options("sphere_radius");
  if (opt.size())
    sphere_radius = opt[0][0];

  opt = get_options("num_strands");
  if (opt.size())
    num_strands = opt[0][0];

  opt = get_options("strand_length");
  if (opt.size())
    strand_length = opt[0][0];
    
  opt = get_options("num_points");
  if (opt.size())
    num_points = opt[0][0];    

  if (degree < 2)
    throw Exception( "DegreeVector (" + str(degree) + ") must be greater than 2.");

  gsl_rng* rand_gen = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rand_gen, time(NULL));

  Fibre::Strand::Set strands;

  MR::ProgressBar progress_bar ("Generating " + str(num_strands) + " short strands on sphere_surface...");
  
  for (size_t strand_i = 0; strand_i < num_strands; strand_i++) {
    Fibre::Strand strand(degree, 0.0);
    
    strand[0] = Triple<double>((double)gsl_ran_gaussian(rand_gen, 1.0), (double)gsl_ran_gaussian(rand_gen, 1.0), (double)gsl_ran_gaussian(rand_gen, 1.0)).normalise() * sphere_radius;   

    strand[1] = strand[0].cross(strand[0].min_axis()).normalise() * strand_length;
  
    strands.push_back(strand);
  }

  strands.save(output_location, num_points);

  //MR::ProgressBar::done();
}


