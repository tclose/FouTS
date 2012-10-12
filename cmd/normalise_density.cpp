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

  size_t num_points = 100;

  Options opt = get_options("num_points");
  if (opt.size())
    num_points = opt[0][0];

  Fibre::Tractlet::Set tcts(input_location);
  
  if (!tcts.has_elem_prop(Fibre::Tractlet::ACS_SQRT_PROP))
      tcts.add_elem_prop(Fibre::Tractlet::ACS_SQRT_PROP);

  MR::ProgressBar progress_bar ("Normalising densities of Fourier tracts...", tcts.size());

  for (size_t tct_i = 0; tct_i < tcts.size(); tct_i++) {
    std::vector<double> areas = tcts[tct_i].cross_sectional_areas(num_points);
    double avg_area = 0.0;
    for (size_t area_i = 0; area_i < num_points; ++area_i)
      avg_area += areas[area_i];
    avg_area /= (double)num_points;
    tcts[tct_i].set_acs(avg_area);
  }
 
  tcts.save(output_location);

}









