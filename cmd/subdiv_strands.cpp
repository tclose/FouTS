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

#include "bts/cmd.h"

#include "progressbar.h"

#include "bts/common.h"

#include "bts/fibre/track.h"
#include "bts/fibre/strand.h"
#include "bts/fibre/track/set.h"





#include "bts/inline_functions.h"

#include "phantom/subdiv/subdiv.h"


#include "phantom/interface.h"



using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Subdivides thicker strands into more thinnessner strands",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("input", "The paths to be subdivided.").type_file (),
  Argument ("output", "The subdivided paths.").optional().type_file(),
  Argument()
};


OPTIONS = {
  Option ("num_points", "The number of points that will be generated along the strand location")
   + Argument ("num_points", "").type_integer (1, 100, 2000),

  Option ("old_strand_r", "The old radius of the strands")
   + Argument ("old_strand_r", "").type_float (SMALL_FLOAT, 1.0, LARGE_FLOAT),

  Option ("strands_per_acs", "Instead of specifying the radius of the subdivided strands, the strands per area can be specified instead.")
   + Argument ("strands_per_acs", "").type_float (-1.0, -1.0, LARGE_FLOAT),

  Option ("new_strand_r", "The new radius of the strands")
   + Argument ("new_strand_r", "").type_float (SMALL_FLOAT, 1.0, LARGE_FLOAT),

  Option ("degree", "DegreeVector of the output strands if converting from track to strand")
   + Argument ("degree", "").type_integer (1, 100, 2000),


Option() };


EXECUTE {

  std::string input_location = argument[0];
  std::string output_location;
  
  if (argument.size() > 1)
    output_location = argument[1].c_str();
  else
    output_location = input_location;


  size_t num_points = 0;
  size_t degree = 0;
  double old_strand_r = 1.0;
  double new_strand_r = 1.0;
  double strands_per_acs = -1.0;

  Options opt = get_options("num_points");
  if (opt.size())
    num_points = opt[0][0];

  opt = get_options("old_strand_r");
  if (opt.size())
    old_strand_r = opt[0][0];

  opt = get_options("new_strand_r");
  if (opt.size())
    new_strand_r = opt[0][0];

  opt = get_options("strands_per_acs");
  if (opt.size()) {
    if (get_options("new_strand_r").size())
      throw Exception ("'-strands_per_acs' option cannot be used in conjunction with '-new_strand_r' option.");

    strands_per_acs = opt[0][0];

  }

  opt = get_options("degree");
  if (opt.size())
    degree = opt[0][0];


  if (strands_per_acs > 0)
    new_strand_r = old_strand_r / MR::Math::sqrt(strands_per_acs);


  BTS::Fibre::Track::Set tcks(input_location, num_points);
  BTS::Fibre::Track::Set subdivided_tcks(tcks.get_extend_props());

  std::vector<Triple<double> > pre_points;
  std::vector<Triple<double> > post_points;

  Strand_collection c, subdivided_c;
  
  generate_pre_points(tcks, pre_points);
  
  generate_post_points(tcks, post_points);

  MR::ProgressBar progress_bar ("Subdividing strands...", 0);

  tcks.add_extend_elem_prop(BTS::Fibre::Track::RADIUS_PROP, str(old_strand_r));

  convert_mr_to_nfg(&c, tcks, pre_points, post_points);
  progress_bar++;
  
  subdivide_collection(&subdivided_c, &c, new_strand_r);

  progress_bar++;
  convert_nfg_to_mr(subdivided_tcks, pre_points, post_points, &subdivided_c);

  if (tcks.has_extend_elem_prop(BTS::Fibre::Track::Set::ORIGINAL_DEGREE_PROP)) {

    subdivided_tcks.add_extend_elem_prop(BTS::Fibre::Track::Set::ORIGINAL_DEGREE_PROP, str(BTS::Fibre::Strand::DEFAULT_DEGREE));

    for (size_t tck_i = 0; tck_i < tcks.size(); ++tck_i) {

      size_t bundle_index = to<size_t>(subdivided_tcks.get_extend_elem_prop(BTS::Fibre::Track::BUNDLE_INDEX_EPROP, tck_i));
      std::string original_degree = tcks.get_extend_elem_prop(BTS::Fibre::Track::Set::ORIGINAL_DEGREE_PROP, bundle_index);
      subdivided_tcks.set_extend_elem_prop(BTS::Fibre::Track::Set::ORIGINAL_DEGREE_PROP, original_degree, tck_i);

    }
  }

  progress_bar++;
  subdivided_tcks.save(output_location, degree);

  //MR::ProgressBar::done();
}
