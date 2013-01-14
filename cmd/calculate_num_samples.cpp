/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Written by Thomas G. Close, 04/03/2009.

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

#include "bts/common.h"
#include "bts/math/common.h"

#include "bts/fibre/strand/set.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"
#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/tractlet/section.h"
#include "bts/fibre/base/set_reader.h"
#include "bts/fibre/base/set_writer.h"

#include "bts/inline_functions.h"

using namespace BTS;

const float SAMPLES_PER_ACS_DEFAULT = 1000;
const float SAMPLES_PER_LENGTH_DEFAULT = 1000;
const size_t MIN_WIDTH_SAMPLES_DEFAULT = 4;
const size_t MIN_LENGTH_SAMPLES_DEFAULT = 4;
const char* STRATEGY_DEFAULT = "max";
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
  "Print strands to screen",
  "",
  NULL
};

ARGUMENTS= {
  Argument ("input", "The strands to be displayed.").type_file (),
  Argument ("output", "The prefix for the outputted files (the two output "
      "files will be appended with 'length.txt' and "
      "'width.txt' respectively.").type_file (),
  Argument()
};

OPTIONS= {

  Option ("samples_per_acs","")
  + Argument ("samples_per_acs","").type_float (0,SAMPLES_PER_ACS_DEFAULT,
      LARGE_FLOAT),

  Option ("samples_per_length","")
  + Argument ("samples_per_length")
  .type_float (0,SAMPLES_PER_LENGTH_DEFAULT,LARGE_FLOAT),

  Option ("min_length_samples","")
  + Argument ("min_length_samples")
  .type_integer (1,MIN_LENGTH_SAMPLES_DEFAULT,LARGE_INT),

  Option ("min_width_samples","")
  + Argument ("min_width_samples")
  .type_integer (1,MIN_WIDTH_SAMPLES_DEFAULT,LARGE_INT),

  Option ("strategy","The strategy used to handle multiple Fourier tracts."
      "Can be either 'max', 'average' or 'median' (defaults to "
      "'max').")
  + Argument ("strategy").type_text(STRATEGY_DEFAULT),
  Option()};

EXECUTE {
  
  std::string input_location = argument[0];
  std::string output_prefix = argument[1];

  double samples_per_acs = SAMPLES_PER_ACS_DEFAULT;
  double samples_per_length = SAMPLES_PER_LENGTH_DEFAULT;
  size_t min_width_samples = MIN_WIDTH_SAMPLES_DEFAULT;
  size_t min_length_samples = MIN_LENGTH_SAMPLES_DEFAULT;
  std::string strategy = STRATEGY_DEFAULT;

  Options opt = get_options("samples_per_acs");
  if (opt.size())
    samples_per_acs = opt[0][0];

  opt = get_options("samples_per_length");
  if (opt.size())
    samples_per_length = opt[0][0];

  opt = get_options("min_width_samples");
  if (opt.size())
    min_width_samples = opt[0][0];

  opt = get_options("min_length_samples");
  if (opt.size())
    min_length_samples = opt[0][0];

  opt = get_options("strategy");
  if (opt.size())
    strategy = opt[0][0].c_str();

  Fibre::Tractlet::Set tcts(input_location);
  size_t num_tracts = tcts.size();

  std::vector<double> acs_vector;
  std::vector<double> length_vector;

  for (size_t tract_i = 0; tract_i < tcts.size(); ++tract_i) {
    acs_vector.push_back(tcts[tract_i].acs());
    length_vector.push_back(tcts[tract_i].path_length());
  }

  double acs, length;

  if (strategy == "max") {
    acs = *std::max(acs_vector.begin(), acs_vector.end());
    length = *std::max(length_vector.begin(), length_vector.end());
  } else if (strategy == "average") {
    length = acs = 0.0;
    for (size_t tract_i = 0; tract_i < tcts.size(); ++tract_i) {
      acs += acs_vector[tract_i];
      length += length_vector[tract_i];
    }
    acs /= (double)num_tracts;
    length /= (double)num_tracts;
  } else if (strategy == "median") {
    acs = Math::median(acs_vector);
    length = Math::median(length_vector);
  } else
    throw Exception(
        "Unrecognised argument to 'strategy' option '" + strategy
            + "' (can be either 'max', 'average', or 'median').");

  size_t num_width_samples = round(acs * samples_per_acs);
  size_t num_length_samples = round(length * samples_per_length);

  if (num_width_samples < min_width_samples)
    num_width_samples = min_width_samples;
  if (num_length_samples < min_length_samples)
    num_length_samples = min_length_samples;

  std::ofstream out;

  // Write width value to file
  out.open((output_prefix + "width.txt").c_str(), std::ios::out);
  out << num_width_samples << std::endl;
  out.close();

  // Write length value to file
  out.open((output_prefix + "length.txt").c_str(), std::ios::out);
  out << num_length_samples << std::endl;
  out.close();


}
