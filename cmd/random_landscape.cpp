/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 10/09/2010.

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


#include "bts/cmd.h"

#include "bts/common.h"
#include "bts/prob/test/landscape.h"
#include "bts/inline_functions.h"

using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "random_peaks",
  "Generates a random probability landscape",
  NULL
};

ARGUMENTS = {
  Argument ("output_location", "Location of the output file.").type_file (),
  Argument()
};


OPTIONS = {

    Option ("location", "")
      + Argument ("location", ""),

    Option ("num_dims", "")
      + Argument ("num_dims", "").type_integer (1, Prob::Test::Landscape::NUM_DIMS_DEFAULT, LARGE_INT),

    Option ("num_peaks", "")
      + Argument ("num_peaks", "").type_integer (1, Prob::Test::Landscape::NUM_PEAKS_DEFAULT, LARGE_INT),

    Option ("width_mu", "")
     + Argument ("width_mu", "").type_float (SMALL_FLOAT, Prob::Test::Landscape::WIDTH_MU_DEFAULT, LARGE_FLOAT),

    Option ("width_sigma", "")
     + Argument ("width_sigma", "").type_float (SMALL_FLOAT, Prob::Test::Landscape::WIDTH_SIGMA_DEFAULT, LARGE_FLOAT),

    Option ("width_min", "")
      + Argument ("width_min", "").type_float (SMALL_FLOAT, Prob::Test::Landscape::WIDTH_MIN_DEFAULT, LARGE_FLOAT),

    Option ("height_var", "")
      + Argument ("height_var", "").type_float (SMALL_FLOAT, Prob::Test::Landscape::HEIGHT_VAR_DEFAULT, LARGE_FLOAT),

    Option ("roi_radius", "")
      + Argument ("roi_radius", "").type_float (SMALL_FLOAT, Prob::Test::Landscape::ROI_RADIUS_DEFAULT, LARGE_FLOAT),

    Option ("fraction_pyramid", "")
      + Argument ("fraction_pyramid", "").type_float (SMALL_FLOAT, Prob::Test::Landscape::FRACTION_PYRAMID_DEFAULT, LARGE_FLOAT),

    Option ("barrier_rate", "")
      + Argument ("barrier_rate", "").type_float (SMALL_FLOAT, Prob::Test::Landscape::BARRIER_RATE_DEFAULT, LARGE_FLOAT),

    Option ("seed", "The random seed that is passed to the random generator")
     + Argument ("seed", ""),

Option() };


EXECUTE {

  std::string output_location = argument[0];

  size_t seed = time(NULL);

  Options opt;

  size_t num_dims             = Prob::Test::Landscape::NUM_DIMS_DEFAULT;
  size_t  num_peaks           = Prob::Test::Landscape::NUM_PEAKS_DEFAULT;
  double width_mu           = Prob::Test::Landscape::WIDTH_MU_DEFAULT;
  double width_sigma        = Prob::Test::Landscape::WIDTH_SIGMA_DEFAULT;
  double width_min          = Prob::Test::Landscape::WIDTH_MIN_DEFAULT;
  double height_var         = Prob::Test::Landscape::HEIGHT_VAR_DEFAULT;
  double roi_radius         = Prob::Test::Landscape::ROI_RADIUS_DEFAULT;
  double fraction_pyramid   = Prob::Test::Landscape::FRACTION_PYRAMID_DEFAULT;
  double barrier_rate       = Prob::Test::Landscape::BARRIER_RATE_DEFAULT;

  opt = get_options("num_dims");
  if (opt.size())
    num_dims = opt[0][0];

  opt = get_options("num_peaks");
  if (opt.size())
    num_peaks = opt[0][0];

  opt = get_options("width_mu");
  if (opt.size())
    width_mu = opt[0][0];

  opt = get_options("width_sigma");
  if (opt.size())
    width_sigma = opt[0][0];

  opt = get_options("width_min");
  if (opt.size())
    width_min = opt[0][0];

  opt = get_options("height_var");
  if (opt.size())
    height_var = opt[0][0];

  opt = get_options("roi_radius");
  if (opt.size())
    roi_radius = opt[0][0];

  opt = get_options("fraction_pyramid");
  if (opt.size())
    fraction_pyramid = opt[0][0];

  opt = get_options("barrier_rate");
  if (opt.size())
    barrier_rate = opt[0][0];

  opt = get_options("seed");
  if (opt.size()) {
    std::string seed_string = opt[0][0];
    seed = to<size_t>(seed_string);
  } else
    std::cout << "No random seed supplied. Using timestamp: " << seed << std::endl;


  std::map<std::string, std::string> properties;

  properties["num_dims"]          = str(num_dims);
  properties["barrier_rate"]      = str(barrier_rate);
  properties["roi_radius"]        = str(roi_radius);
  properties["num_peaks"]         = str(num_peaks);
  properties["width_mu"]          = str(width_mu);
  properties["width_sigma"]       = str(width_sigma);
  properties["width_min"]         = str(width_min);
  properties["height_var"]        = str(height_var);
  properties["fraction_pyramid"]  = str(fraction_pyramid);



  Prob::Test::Landscape peaks = Prob::Test::Landscape::randomly_generate( num_dims,
                                                                  num_peaks,
                                                                  width_mu,
                                                                  width_sigma,
                                                                  width_min,
                                                                  height_var,
                                                                  roi_radius,
                                                                  fraction_pyramid,
                                                                  barrier_rate, seed);

  peaks.save(output_location);

  std::cout << "Generated " << num_peaks << " peaks." << std::endl;

}
