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

#include "bts/fibre/track/set.h"
#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"



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
  Argument ("input", "The input fibres.").type_file(),
  Argument ("output", "The output fibres.").optional().type_file(),
  Argument()
};


OPTIONS = {

  Option ("degree", "Limit the degree of the Strand coefficients used to describe the strands to the first - degrees")
   + Argument ("degree", "Limit the degree of the Strand coefficients used to describe the strands to the first - degrees").type_integer (1, 200, 3),

  Option ("min_length", "The length of the short strands")
   + Argument ("min_length", "").type_float (-1, -1, LARGE_FLOAT),

  Option ("min_acs", "The minimum tract volume.")
   + Argument ("min_acs", "").type_float (-1, -1, LARGE_FLOAT),

  Option ("num_length_sections", "The number of points that will be generated along the strand location")
   + Argument ("points", "").type_integer (1, 100, 2000),

  Option ("num_width_sections", "The number of points that will be generated along the strand location")
   + Argument ("points", "").type_integer (1, 100, 2000),

Option() };


EXECUTE {

  std::string input_location = argument[0];
  std::string output_location = input_location;

  if (argument.size() > 1)
    output_location = argument[1].c_str();


//  size_t degree = 3;
  double min_length = -1;
  double min_acs = -1;
//  size_t num_length_sections = 100;
//  size_t num_width_sections = 6;

  Options opt = get_options("min_length");
  if (opt.size())
    min_length = opt[0][0];

  opt = get_options("min_acs");
  if (opt.size())
    min_acs = opt[0][0];

//  opt = get_options("degree");
//  if (opt.size())
//    degree = opt[0][0];
//  opt = get_options("num_length_sections");
//  if (opt.size())
//    num_length_sections = opt[0][0];
//  opt = get_options("num_width_sections");
//  if (opt.size())
//    num_width_sections = opt[0][0];


  if ((File::has_extension<Fibre::Strand>(input_location) || File::has_extension<Fibre::Track>(input_location)) && (min_length > 0) ) {
    Fibre::Track::Set in(input_location);
    Fibre::Track::Set out(in.get_extend_props());
    out.add_extend_elem_props(in);

    MR::ProgressBar progress_bar("Removing short strands/tracks....");

    for (size_t tck_i = 0; tck_i < in.size(); ++tck_i) {

      Fibre::Track tck = in[tck_i];

      double length = 0.0;

      for (size_t point_i = 1; point_i < tck.size(); ++point_i) {
       length += (tck[point_i] - tck[point_i - 1]).norm();
      }

      if (length > min_length)
        out.push_back(tck, in.get_extend_elem_prop_row(tck_i));

      ++progress_bar;

    }

    in = out;

    out.save(output_location);

  }


  if (File::has_extension<Fibre::Tractlet>(input_location) && (min_acs > 0)) {

    Fibre::Tractlet::Set in(input_location);
    Fibre::Tractlet::Set out(in.get_extend_props());
    out.add_extend_elem_props(in);


    MR::ProgressBar progress_bar("Removing thin tractlets....");

    for (size_t tract_i = 0; tract_i < in.size(); ++tract_i) {

      if (in[tract_i].acs() > min_acs)
        out.push_back(in[tract_i], in.get_extend_elem_prop_row(tract_i));

      ++progress_bar;

    }

    in = out;

    out.save(output_location);

  }



}


