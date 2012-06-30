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


#include "phantom/interface.h"

#include "phantom/rand_init/rand_init.h"
#include "bts/fibre/track.h"
#include "bts/fibre/track/set.h"




#include "bts/inline_functions.h"

const size_t DEGREE_DEFAULT = 3;
const double SPHERE_RADIUS_DEFAULT = 0.225;
const double BASE_WIDTH_DEFAULT = 0.1;
const double STRAND_RADIUS_LBOUND_DEFAULT =  0.12;
const double STRAND_RADIUS_UBOUND_DEFAULT = 0.3;
const double CONTROL_POINT_FREQUENCY_DEFAULT = 50.0;
const double STRAND_RADIUS_BUFFER_RATIO_DEFAULT = 1.12;


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
  Argument ("output", "Random strands").type_file(),
  Argument()
};


OPTIONS = {

  Option ("degree", "Limit the degree of the Strand coefficients used to describe the strands to the first - degrees")
   + Argument ("degree", "").type_integer (1, DEGREE_DEFAULT, 200),

  Option ("sphere_r", "The radius of the sphere the strands will be trimmed to")
   + Argument ("sphere_r", "").type_float (SMALL_FLOAT, SPHERE_RADIUS_DEFAULT, LARGE_FLOAT),

  Option ("control_point_freq", "The frequency of the control-points generated between the strand's end-points")
   + Argument ("control_point_freq", "").type_float (SMALL_FLOAT, CONTROL_POINT_FREQUENCY_DEFAULT, LARGE_FLOAT),

  Option ("strand_r_lbound", "The minimum strand radius")
   + Argument ("strand_r_lbound", "The minimum strand radius").type_float (SMALL_FLOAT, STRAND_RADIUS_LBOUND_DEFAULT, LARGE_FLOAT),

  Option ("strand_r_ubound", "The maximum strand radius")
   + Argument ("strand_r_ubound", "The maximum strand radius").type_float (SMALL_FLOAT, STRAND_RADIUS_UBOUND_DEFAULT, LARGE_FLOAT),

  Option ("strand_r_buffer_ratio", "Scales the 'buffer radius' around the endpoints of each strand that other strands cannot overlap.")
   + Argument ("strand_r_buffer_ratio", "").type_float (SMALL_FLOAT, STRAND_RADIUS_BUFFER_RATIO_DEFAULT, LARGE_FLOAT),

  Option ("max_attempts", "Max. attempts to add another strand")
   + Argument ("max_attempts", "").type_integer (1, 1000000, LARGE_INT),

  Option ("seed", "Seed for the random generation")
   + Argument ("seed", ""),


Option() };


EXECUTE {


  std::string output_location = argument[0];

  size_t degree = DEGREE_DEFAULT;
  double sphere_r = SPHERE_RADIUS_DEFAULT;
  double control_point_freq = CONTROL_POINT_FREQUENCY_DEFAULT;
  double strand_r_lbound = STRAND_RADIUS_LBOUND_DEFAULT;
  double strand_r_ubound = STRAND_RADIUS_UBOUND_DEFAULT;
  double strand_r_buffer_ratio = STRAND_RADIUS_BUFFER_RATIO_DEFAULT;
  size_t max_attempts = 1000000;
  size_t seed = time(NULL);

  Options opt = get_options("degree");
  if (opt.size())
    degree = opt[0][0];

  opt = get_options("sphere_r");
  if (opt.size())
    sphere_r = opt[0][0];

  opt = get_options("control_point_freq");
  if (opt.size())
    control_point_freq = opt[0][0];

  opt = get_options("strand_r_lbound");
  if (opt.size())
    strand_r_lbound = opt[0][0];

  opt = get_options("strand_r_ubound");
  if (opt.size())
    strand_r_ubound = opt[0][0];

  opt = get_options("strand_r_buffer_ratio");
  if (opt.size())
    strand_r_buffer_ratio = opt[0][0];

  opt = get_options("max_attempts");
  if (opt.size())
    max_attempts = opt[0][0];

  opt = get_options("seed");
  if (opt.size()) {
    std::string seed_str = opt[0][0];
    seed = to<size_t>(seed_str);
  }

  gsl_rng* rand_gen = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rand_gen, seed);

  MR::ProgressBar progress_bar ("Generating random strands...");

  Strand_collection c;

  rand_init_collection(&c, sphere_r, control_point_freq, strand_r_lbound, strand_r_ubound, strand_r_buffer_ratio, rand_gen, max_attempts, 0, NULL);

  BTS::Fibre::Track::Set tcks;
  std::vector<Triple<double> > pre_points;
  std::vector<Triple<double> > post_points;

  convert_nfg_to_mr(tcks, pre_points, post_points, &c);

  tcks.save(output_location, degree);

  //MR::ProgressBar::done();
  
  std::cout << "Generated " + str(tcks.size()) + " strands." << std::endl;
}


