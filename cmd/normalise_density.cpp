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

#include "bts/common.h"

#include "bts/cmd.h"
#include "progressbar.h"

#include "bts/fibre/tractlet/set.h"

#include "bts/inline_functions.h"


using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Normalises the density of the given Fourier tracts.",
  "Normalises the density of the given Fourier tracts.",
  NULL
};

ARGUMENTS = {
  Argument ("input", "The Fourier tracts to be normalised.").type_file (),
  Argument ("output", "Output location of the normalised tracts (if omitted the origin location will be used)").optional().type_file(),
  Argument()
};


OPTIONS = {

    Option ("width_epsilon", "The scaling of the width epsilon parameter used to provide a lower bound on the ACS that is tied to the size of the 0th degree auxiliary vectors")
     + Argument ("width_epsilon", "The scaling of the width epsilon parameter used to provide a lower bound on the ACS that is tied to the size of the 0th degree auxiliary vectors").type_float (SMALL_FLOAT, Fibre::Tractlet::WIDTH_EPSILON_DEFAULT, LARGE_FLOAT),

    Option ("length_epsilon", "The scaling of the width epsilon parameter used to provide a lower bound on the ACS that is tied to the size of the 1st degree primary vector")
    + Argument ("length_epsilon", "The scaling of the width epsilon parameter used to provide a lower bound on the ACS that is tied to the size of the 1st degree primar vector").type_float (SMALL_FLOAT, Fibre::Tractlet::LENGTH_EPSILON_DEFAULT, LARGE_FLOAT),

    Option ("num_points", "The number of points that will be generated along the strand location")
     + Argument ("num_points", "The number of points that will be generated along the strand location").type_integer (1, 100, 2000),

Option() };


EXECUTE {


  std::string input_location = argument[0];
  std::string output_location;

  if (argument.size() > 1)
    output_location = argument[1].c_str();
  else
    output_location = input_location;

  double width_epsilon = Fibre::Tractlet::WIDTH_EPSILON_DEFAULT;
  double length_epsilon = Fibre::Tractlet::LENGTH_EPSILON_DEFAULT;
  size_t num_points = 100;

  Options opt = get_options("width_epsilon");
  if (opt.size())
    width_epsilon = opt[0][0];

  opt = get_options("length_epsilon");
  if (opt.size())
    length_epsilon = opt[0][0];

  opt = get_options("num_points");
  if (opt.size())
    num_points = opt[0][0];

  Fibre::Tractlet::Set tcts(input_location);
  
  if (!tcts.has_elem_prop(Fibre::Tractlet::ALPHA_PROP))
      tcts.add_elem_prop(Fibre::Tractlet::ALPHA_PROP);

  MR::ProgressBar progress_bar ("Normalising densities of Fourier tracts...");

  tcts.normalise_densities(width_epsilon, length_epsilon, num_points);

  tcts.set_extend_prop("width_epsilon", str(width_epsilon));
  tcts.set_extend_prop("length_epsilon", str(length_epsilon));
  tcts.save(output_location);

}








