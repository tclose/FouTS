/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 18/06/2010.

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

#include "progressbar.h"
#include "bts/cmd.h"

#include "bts/common.h"

#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"


#include "bts/inline_functions.h"

using namespace BTS;

const double ACS_DEFAULT = 1.0;
const double BASE_INTENSITY_DEFAULT = 0.0;
const size_t  DEGREE_DEFAULT = 3;
const size_t  NUM_LENGTH_SECTIONS_DEFAULT = 10;
const size_t  SET_SIZE_DEFAULT = 1;
const size_t  NUM_SETS_DEFAULT = 1;
const double INITIAL_VALUE_DEFAULT = 0.0;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Generates new fibres [sets] of arbitrary size initialised to an arbitrary value",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("output","Output Location"),
  Argument()
};


OPTIONS = {

	Option("acs","Apparrent Connection Strength (ACS) of the generated fibres.")
	 + Argument("acs","").type_float(SMALL_FLOAT,ACS_DEFAULT,LARGE_FLOAT),

  Option("base_intensity","Base intensity of the generated fibres.")
   + Argument("base_intensity","").type_float(SMALL_FLOAT,ACS_DEFAULT,LARGE_FLOAT),

	Option("degree","The degree of the generated fibres.")
	 + Argument("degree","").type_integer(0,DEGREE_DEFAULT,LARGE_INT),

  Option("num_points","The degree of the generated fibres.")
   + Argument("degree","").type_integer(0,NUM_LENGTH_SECTIONS_DEFAULT,LARGE_INT),

  Option("set_size","Size of the generated fibre sets.")
   + Argument("set_size","").type_integer(0,SET_SIZE_DEFAULT,LARGE_INT),

	Option("num_sets","Number of generated fibre sets.")
	 + Argument("num_sets","").type_integer(0,NUM_SETS_DEFAULT,LARGE_INT),

  Option("initial_value","Value all the elements will be initialized to.")
   + Argument("initial_value","").type_float(-INFINITY,ACS_DEFAULT,INFINITY),


Option() };


EXECUTE {


  std::string output_location = argument[0];


  double acs      = 1.0;
  double base_intensity = 1.0;
  size_t degree         = DEGREE_DEFAULT;
  size_t num_points     = NUM_LENGTH_SECTIONS_DEFAULT;
  size_t set_size       = SET_SIZE_DEFAULT;
  size_t num_sets       = NUM_SETS_DEFAULT;
  double initial_value = INITIAL_VALUE_DEFAULT;


  Options opt;

  opt = get_options("acs");
  if (opt.size())
    acs = opt[0][0];

  opt = get_options("base_intensity");
  if (opt.size())
    base_intensity = opt[0][0];

  opt = get_options("degree");
  if (opt.size())
    degree = opt[0][0];

  opt = get_options("num_points");
  if (opt.size())
    num_points = opt[0][0];

  opt = get_options("set_size");
  if (opt.size())
    set_size = opt[0][0];

  opt = get_options("num_sets");
  if (opt.size())
    num_sets = opt[0][0];

  opt = get_options("initial_value");
  if (opt.size())
    initial_value = opt[0][0];


  MR::ProgressBar progress_bar ("Generating new fibres...");

  std::vector<const char*> props;
  if (base_intensity >= 0.0)
    props.push_back(Fibre::Strand::Set::BASE_INTENSITY_PROP);

  std::vector<const char*> elem_props;
  if (acs >= 0.0)
    elem_props.push_back(Fibre::Strand::ALPHA_PROP);


  if (File::has_or_txt_extension<Fibre::Strand>(output_location) || File::has_extension<Fibre::Strand::Set>(output_location) || File::has_txt_extension<Fibre::Strand::Set>(output_location)) {

    Fibre::Strand::Set strands(set_size, degree, props, elem_props);

    strands = initial_value;

    if (base_intensity >= 0.0)
      strands.set_base_intensity(base_intensity);

    if (acs >= 0.0)
      for (size_t strand_i = 0; strand_i < strands.size(); ++strand_i)
        strands[strand_i].set_acs(acs);

    if (File::has_or_txt_extension<Fibre::Strand>(output_location)) {

      strands.save(output_location);

    } else if (File::has_extension<Fibre::Strand::Set>(output_location)) {

      Fibre::Strand::Set::Writer writer(output_location, strands);

      for (size_t set_i = 0; set_i < num_sets; set_i++)
        writer.append(strands);

    } else if (File::has_txt_extension<Fibre::Strand::Set>(output_location)) {

      Fibre::Strand::Set::TextWriter writer(output_location, strands);

      for (size_t set_i = 0; set_i < num_sets; set_i++)
        writer.append(strands);

    } else
      throw Exception ("should not be able to get here.");

  } else if (File::has_or_txt_extension<Fibre::Tractlet>(output_location) || File::has_extension<Fibre::Tractlet::Set>(output_location) || File::has_txt_extension<Fibre::Tractlet::Set>(output_location)) {

    Fibre::Tractlet::Set tractlets(set_size, degree, props, elem_props);

    tractlets = initial_value;

    if (base_intensity >= 0.0)
      tractlets.set_base_intensity(base_intensity);

    if (acs >= 0.0)
      for (size_t tractlet_i = 0; tractlet_i < tractlets.size(); ++tractlet_i)
        tractlets[tractlet_i].set_acs(acs);

    if (File::has_or_txt_extension<Fibre::Tractlet>(output_location)) {

      tractlets.save(output_location);

    } else if (File::has_extension<Fibre::Tractlet::Set>(output_location)) {

      Fibre::Tractlet::Set::Writer writer(output_location, tractlets);

      for (size_t set_i = 0; set_i < num_sets; set_i++)
        writer.append(tractlets);

    } else if (File::has_txt_extension<Fibre::Tractlet::Set>(output_location)) {

      Fibre::Tractlet::Set::TextWriter writer(output_location, tractlets);

      for (size_t set_i = 0; set_i < num_sets; set_i++)
        writer.append(tractlets);

    } else
      throw Exception ("should not be able to get here.");

  } else if (File::has_or_txt_extension<Fibre::Track>(output_location) || File::has_extension<Fibre::Track::Set>(output_location) || File::has_txt_extension<Fibre::Track::Set>(output_location)) {

    Fibre::Track::Set tcks(set_size, num_points, props, elem_props);

    tcks = initial_value;

    if (base_intensity >= 0.0)
      tcks.set_base_intensity(base_intensity);

    if (acs >= 0.0)
      for (size_t track_i = 0; track_i < tcks.size(); ++track_i)
        tcks[track_i].set_acs(acs);

    if (File::has_or_txt_extension<Fibre::Track>(output_location)) {

      tcks.save(output_location);

    } else if (File::has_extension<Fibre::Track::Set>(output_location)) {

      Fibre::Track::Set::Writer writer(output_location, tcks);

      for (size_t set_i = 0; set_i < num_sets; set_i++)
        writer.append(tcks);

    } else if (File::has_txt_extension<Fibre::Track::Set>(output_location)) {

      Fibre::Track::Set::TextWriter writer(output_location, tcks);

      for (size_t set_i = 0; set_i < num_sets; set_i++)
        writer.append(tcks);

    } else
      throw Exception ("should not be able to get here.");

  } else
    throw Exception ("Unrecognised file extension '" + File::extension(output_location) + "'.");


  //MR::ProgressBar::done();

}
